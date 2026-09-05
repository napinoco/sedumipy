"""Port of sedumi.m: the top-level predictor-corrector interior-point
driver that ties together every other Phase 3-c/3-d piece into one
solve of MINIMIZE c'*x SUCH THAT A*x=b, x in K.

SCOPE (v1): LP, second-order-cone (Lorentz K.q/K.r -- pretransfo.py
already folds rotated K.r cones into standard K.q ones), and PSD (K.s)
problems. sedumi.m's own main loop branches on `sum(K.s)==0` for how it
rebuilds ADA every iteration: the `==0` branch uses getada.py directly;
the nonzero branch uses getada_psd.py's build_aord()/getada_psd()
(getada1->getada2->getada3 orchestration), with its own one-time
pre-loop Aord/getsymbada/symbchol setup -- original sedumi.m runs that setup
unconditionally for both branches, but this port's `K.s==0` path keeps
its existing (already-tested) simpler post-`sdinit` ordering instead of
being rewritten to match; the two branches genuinely differ here.

Dense-column preconditioning IS wired in: `getdense()` flags a small
proportion of LP/Lorentz columns as dense right after `pretransfo`,
those rows are zeroed out of `A2` exactly once (before the pre-loop
`Aord`/`symbchol`/`symbcholden` setup, the main loop, and the final-tasks
diagnostics below -- all of which reuse this same permanently-zeroed
`A2`, matching upstream sedumi.m's own single `A(dense.cols,:) = 0.0`),
and `symbcholden()` builds the `symLden` structure `deninfac()` needs
alongside `symbchol()`'s ordinary `Lsym`.

Also not ported (cosmetic/diagnostic, no effect on the returned
(x,y,info)): the console progress printout (my_fprintf/pars.fid),
pars.vplot's v-plot, pars.stopat's interactive debug break, the
optional pre-solve rank/infeasibility diagnostic (a warning heuristic),
and the origcoeff DIMACS error-measures block (info.err).
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from . import _native
from .amul import amul
from .checkpars import checkpars
from .cone import asmDxq, frameit, maxeigK, psdscale
from .deninfac import deninfac
from .getada import getada
from .getada_psd import build_aord, build_q_pattern, getada_psd
from .getdatm import getDAtm
from .getdense import getdense
from .getsymbada import getsymbada
from .optstep import optstep
from .posttransfo import posttransfo
from .pretransfo import pretransfo
from .sdfactor import sdfactor
from .sdinit import sdinit
from .symbchol import symbchol
from .symbcholden import symbcholden
from .updtransfo import updtransfo
from .wregion import wregion


def sedumi(A, b, c, K: dict, pars: dict | None = None, **pars_kwargs):
    """x, y, info = sedumi(A, b, c, K, pars=None, **pars_kwargs)

    `pars_kwargs` lets individual options be overridden as keyword
    arguments (e.g. `sedumi(A, b, c, K, eps=1e-9)`) instead of building a
    `pars` dict by hand; entries in `pars_kwargs` win over the same key
    in `pars` when both are given. See checkpars.py for the full list of
    recognized `pars` fields and their defaults."""
    if pars_kwargs:
        pars = {**(pars or {}), **pars_kwargs}
    A2, b2, c2, K2, prep, _origcoeff = pretransfo(A, b, c, K, pars or {})
    b2 = np.asarray(b2, dtype=np.float64).ravel()
    c2 = np.asarray(c2, dtype=np.float64).ravel()

    has_psd = len(K2.get("s", [])) > 0

    lponly = int(K2["l"]) == len(c2)
    pars = checkpars(pars)

    # ---- Remove dense columns (if any) -- sedumi.m lines ~352-364. This
    # zeroing of A2 happens exactly once, here, and the same A2 is reused
    # for the rest of this function (pre-loop setup, main loop, and the
    # final-tasks diagnostics), matching upstream's single A(dense.cols,:)
    # = 0.0 assignment.
    Ablkjc = _native.partitA(A2, K2["mainblks"])
    dense, DAtdenq = getdense(A2, Ablkjc, K2, pars)
    if dense["cols"].size:
        dense["A"] = A2.tocsc()[dense["cols"] - 1, :].T.tocsc()
        A2 = A2.tolil()
        A2[dense["cols"] - 1, :] = 0.0
        A2 = A2.tocsc()
        Ablkjc = _native.partitA(A2, K2["mainblks"])
    else:
        dense["A"] = sp.csc_matrix((b2.size, 0))

    if has_psd:
        is_dense = False  # getada_psd() ignores DAt["q"]'s representation
        Ablkjc, Aord, ADA = build_aord(A2, K2, dense)
        Lsym = symbchol(ADA)
        symLden = symbcholden(Lsym, dense, {"denq": DAtdenq})
        d, v, vfrm, y, y0, R = sdinit(A2, b2, c2, dense, K2, pars)
    else:
        # ---- one-time symbolic Cholesky pattern for ADA ----
        # Original sedumi.m builds this via getsymbada.m *before* sdinit.m
        # ever runs, unconditionally for both the K.s==0 and K.s!=0
        # cases (this module's SCOPE docstring already notes the
        # has_psd branch above replicates that unconditional setup via
        # build_aord()/getsymbada() -- a purely structural, 0/1
        # computation from A's own sparsity pattern). This branch used
        # to instead build the pattern numerically, from getada() on a
        # placeholder d with d["q2"] forced to all-ones (needed because
        # sdinit.py always starts the real d["q2"] at exactly 0 --
        # sdinit.m's own d.q2 = zeros(...) -- which would otherwise miss
        # the rank-1 cross-constraint coupling a Lorentz block's d["q2"]
        # term contributes once later iterations make it nonzero).
        # That placeholder fix was confirmed against the real Octave
        # build on vendor/sedumi-upstream/examples/nb.mat (394 Lorentz
        # blocks), but a numeric matrix product can ALSO under-cover the
        # true pattern through ordinary floating-point cancellation --
        # confirmed on DIMACS's nb_L2.mat (839 Lorentz blocks feeding
        # just 123 constraints, so many blocks routinely land in the
        # same ADA[i,j] cell): the placeholder's particular q1/l/det
        # values happened to cancel some of those sums to exactly 0 at
        # iteration 1, silently dropping positions later iterations
        # really need and degrading PCG's preconditioner until it
        # stopped converging within the iteration cap (numerr=2), while
        # the real Octave/MEX build solves the same file cleanly
        # (numerr=0). getsymbada()'s spones-based approach unions
        # structural patterns, never numeric values, so it can't suffer
        # that kind of cancellation by construction.
        #
        # Built via build_q_pattern()+getsymbada() directly rather than
        # the has_psd branch's build_aord(): this branch only needs the
        # resulting structural ADA pattern, not Aord's other fields
        # (lqperm/sperm/dz), and build_aord() computes those unconditionally
        # via incorder() -- a Python-level loop over every one of A's m
        # columns regardless of how small the PSD range is, which
        # dominates runtime for no benefit at nql180/qssp180 scale
        # (m ~ 1.3e5, no PSD blocks at all).
        Aord = None
        q_pattern = build_q_pattern(A2, Ablkjc, K2, dense)
        ADA_symbolic = getsymbada(A2, Ablkjc, q_pattern, K2["sblkstart"])
        # getDAtm()/getada() rebuild ADA every iteration of the main loop
        # below, so which in-memory form (sparse vs dense) they use for
        # it matters for runtime, not just memory: reuse getsymbada()'s
        # own density verdict (same >0.9 threshold it already applies to
        # its own fallback) rather than adding a second, inconsistent
        # size/density cutoff. Computed once, from a matrix already in
        # hand -- O(1) beyond the getsymbada() call above regardless of
        # m -- and threaded through every getDAtm() call for this solve
        # (getada() itself just follows whatever representation DAt["q"]
        # comes in as). See getada.py's SCOPE docstring for why this
        # matters: sparse @ sparse avoids OOM on large near-sparse ADA
        # (DIMACS's nql180/qssp180, m ~ 1.3e5) but is pure overhead
        # relative to a dense BLAS matmul when ADA comes out dense/
        # near-dense anyway (e.g. nb.mat, m=123).
        m2 = A2.shape[1]
        is_dense = (ADA_symbolic.nnz / (m2 * m2) if m2 else 0.0) > 0.9
        d, v, vfrm, y, y0, R = sdinit(A2, b2, c2, dense, K2, pars)
        DAt = getDAtm(A2, Ablkjc, dense, DAtdenq, d, K2, is_dense=is_dense)
        DAtdenq = DAt["denq"]

        ADA, _absd0 = getada(A2, K2, d, DAt)
        Lsym = symbchol(ADA_symbolic)
        symLden = symbcholden(Lsym, dense, DAt)

    n = vfrm["lab"].size
    Kl = int(K2["l"])

    merit = (float(np.sum(R["w"])) + max(R["sd"], 0.0)) ** 2 * y0 / R["b0"]

    STOP = 0
    iter_ = 0
    wr = {"delta": 0.0, "desc": 1}
    feasratio = 0.0
    xsol = ysol = None
    Lnum = None
    Lsd = {"kcg": 0}
    err = {"kcg": 0}

    while STOP == 0:
        iter_ += 1

        if pars["stepdif"] == 2 and (
            iter_ > 20
            or (iter_ > 1 and (err["kcg"] + Lsd["kcg"] > 3))
            or (iter_ > 5 and abs(1 - feasratio) < 0.05)
        ):
            pars["stepdif"] = 1

        # ---- ADA update + factorization ----
        DAt = getDAtm(A2, Ablkjc, dense, DAtdenq, d, K2, is_dense=is_dense)
        DAtdenq = DAt["denq"]
        if has_psd:
            ADA, absd = getada_psd(ADA, A2, Ablkjc, Aord, DAt, d, K2)
        else:
            ADA, absd = getada(A2, K2, d, DAt)

        fact = _native.numeric_cholesky(Lsym, ADA, pars["chol"], absd)
        Lnum = {
            "L": fact["L"], "d": fact["d"], "skip": fact["skip"],
            "perm": Lsym["perm"], "xsuper": Lsym["xsuper"], "tmpsiz": Lsym["tmpsiz"],
        }
        Lden, Lnum["d"] = deninfac(symLden, Lnum, dense, DAt, d, absd, K2.get("qblkstart"), pars["chol"])

        Lsd = sdfactor(Lnum, Lden, dense, DAt, d, v, y, A2, c2, K2, R, y0, pars)

        y0Old = y0
        xscl, yNxt, zscl, y0Nxt, w, relt, dxmdz, err, wr = wregion(
            Lnum, Lden, Lsd, d, v, vfrm, A2, DAt, dense, R, K2, y, y0, b2, pars, wr
        )

        if y0Nxt > 0:
            R["b"] = R["b"] + err["b"] / y0Nxt
            R["sd"] = R["sd"] + err["g"] / y0Nxt
            R["b0"] = R["b0"] + err["db0"] / y0Nxt
            y0 = y0Nxt
        else:
            R["b"] = (y0Nxt * R["b"] + err["b"]) / y0Old
            R["sd"] = (y0Nxt * R["sd"] + err["g"]) / y0Old
            R["b0"] = (y0Nxt * R["b0"] + err["db0"]) / y0Old
            R["w"][1] = abs(y0Nxt / y0Old) * R["w"][1]
            R["c"] = (y0Nxt / y0Old) * R["c"]
            R["maxRc"] = float(np.linalg.norm(R["c"], np.inf)) if R["c"].size else 0.0
            y0 = y0Old

        R["maxRb"] = float(np.linalg.norm(R["b"], np.inf)) if R["b"].size else 0.0
        R["w"][0] = 2 * pars["w"][0] * R["maxRb"] / (1 + R["maxb"])
        meritOld = merit
        merit = (float(np.sum(R["w"])) + max(R["sd"], 0.0)) ** 2 * y0 / R["b0"]
        rate = merit / meritOld

        if rate >= 0.9999 and wr["desc"] == 1:
            STOP = -1
            iter_ -= 1
            y0 = y0Old
            break

        feasratio = float(dxmdz[0] / v[0])

        y = yNxt
        by = float(np.sum(b2 * y))
        d, vfrm = updtransfo(xscl, zscl, w, d, K2)
        v = frameit(vfrm["lab"], vfrm["q"], vfrm["s"], K2)
        x0 = float(np.sqrt(d["l"][0]) * v[0])

        if lponly and rate < 0.05:
            xsol_try, ysol_try = optstep(
                A2, b2, c2, y0, y, d, v, dxmdz, K2, Lnum, symLden, dense, feasratio, R, pars
            )
            if xsol_try is not None:
                STOP = 2
                feasratio = 1 - 2 * (xsol_try[0] == 0)
                xsol, ysol = xsol_try, ysol_try
                break
        elif by > 0 and abs(1 + feasratio) < 0.05 and R["b0"] * y0 < 0.5:
            if maxeigK(amul(A2, dense, y, transp=True), K2) <= pars["eps"] * by:
                STOP = 3
                break

        r0 = float(np.sum(R["w"]))
        cx = by + y0 * R["sd"] - x0 / d["l"][0]
        rgap = max(cx - by, 0.0) / max(abs(cx), abs(by), 1e-3 * x0)
        precision1 = y0 * r0 / (1 + x0)
        precision2 = (y0 * r0 + rgap) / x0
        if precision1 < pars["eps"]:
            if precision2 < pars["eps"]:
                STOP = 1
                break
            elif y0 * R["maxRb"] + x0 * R["maxb"] < -pars["eps"] * cx:
                STOP = 1
                break
            elif y0 * R["maxRc"] + x0 * R["maxc"] < pars["eps"] * by:
                STOP = 1
                break

        if iter_ >= pars["maxiter"]:
            STOP = -1

    # ************************************************************
    # FINAL TASKS
    # ************************************************************
    info = {"iter": iter_, "feasratio": feasratio, "pinf": 0, "dinf": 0, "numerr": 0, "r0": np.inf}

    if STOP == 2:
        x = xsol
        y = ysol
    elif STOP == 3:
        x = np.zeros(len(c2))
    else:
        x = np.concatenate([np.sqrt(d["l"]) * v[:Kl], asmDxq(d, v, K2), psdscale(d, v, K2, transp=True)])

    x0 = float(x[0])
    cx = float(np.sum(c2 * x))
    abscx = float(np.sum(np.abs(c2) * np.abs(x)))
    by = float(np.sum(b2 * y))
    Ax = amul(A2, dense, x, transp=False)
    Ay = amul(A2, dense, y, transp=True)
    normy = float(np.linalg.norm(y))
    normx = float(np.linalg.norm(x[1:]))

    pinf = float(np.linalg.norm(x0 * b2 - Ax))
    dinf = float(maxeigK(Ay - x0 * c2, K2))
    if x0 > 0:
        relinf = max(pinf / (1 + R["maxb"]), dinf / (1 + R["maxc"])) / x0
        if relinf > pars["eps"]:
            pdirinf = float(np.linalg.norm(Ax))
            ddirinf = float(maxeigK(Ay, K2))
            reldirinf = pdirinf / (-cx) if cx < 0.0 else np.inf
            if by > 0.0:
                reldirinf = min(reldirinf, ddirinf / by)
            if reldirinf < pars["eps"] or relinf > max(pars["bigeps"], reldirinf):
                x0 = 0.0
                pinf = pdirinf
                dinf = ddirinf

    if x0 > 0:
        x = x / x0
        y = y / x0
        pinf = pinf / x0
        dinf = dinf / x0
        cx = cx / x0
        by = by / x0
        normx = normx / x0
        normy = normy / x0
        if cx <= by:
            r0 = 0.0
        elif cx == 0.0:
            r0 = -by / (R["maxb"] * normy + 1e-10 * x0)
        elif by == 0.0:
            r0 = cx / (R["maxc"] * normx + 1e-10 * x0)
        else:
            r0 = (cx - by) / (abs(by) + 1e-5 * (x0 + abscx))

        denom = np.array([1.0, 1 + R["maxb"] + 1e-3 * R["maxRb"], 1 + R["maxc"] + 1e-3 * R["maxRc"]])
        info["r0"] = float(np.max(np.array([r0, pinf, dinf]) / denom))
        if STOP == -1:
            if info["r0"] > pars["bigeps"]:
                info["numerr"] = 2
            elif info["r0"] > pars["eps"]:
                info["numerr"] = 1
            else:
                info["numerr"] = 0
        else:
            info["r0"] = min(info["r0"], pars["eps"])
    else:
        if pinf < -pars["bigeps"] * cx:
            info["r0"] = abs(pinf / cx)
            info["dinf"] = 1
            abscx = -cx
            pinf = pinf / abscx
            normx = normx / abscx
            x = x / abscx
        if dinf < pars["bigeps"] * by:
            info["r0"] = abs(dinf / by)
            info["pinf"] = 1
            dinf = dinf / by
            normy = normy / by
            y = y / by
        if info["pinf"] + info["dinf"] == 0:
            info["numerr"] = 2
        elif STOP == -1:
            if pinf > -pars["eps"] * cx and dinf > pars["eps"] * by:
                info["numerr"] = 1
            else:
                info["numerr"] = 0

    x, y, _K_out = posttransfo(x, y, prep, K2)

    return x, y, info
