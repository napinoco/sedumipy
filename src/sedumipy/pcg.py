"""Port of loopPcg.m and wrapPcg.m: the preconditioned-conjugate-gradient
solver for A*P(d)*A' * y = rhs that sits at the heart of every interior-
point iteration (sddir.m/sdfactor.m both delegate to wrapPcg), plus
sparfwslv.m/sparbwslv.m (thin wrappers SeDuMi itself uses around
fwblkslv/bwblkslv -- i.e. this port's own _native.fwsolve/bwsolve).

The dense-column preconditioning correction (dense.q/DAt.denq, from
deninfac.py's product-form Lden) IS implemented: loopPcg's residual
update includes the `DAt.denq*DApq(dense.q)` term (loopPcg.m line
~113-114) alongside the always-applied Lorentz-cone `DAt.q'*DApq` term.
wrapPcg needs no separate dense-column term of its own -- its residual
updates go entirely through amul() (already dense-column-aware) and
loopPcg().
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from . import _native
from .amul import amul
from .cone import PopK, asmDxq, psdscale


def sparfwslv(L: dict, b, ysymb=None):
    """y = sparfwslv(L,b[,ysymb]): forward-solve y := L\\b(L.perm) (L a
    blkchol.m-style dict: L["L"] the unit-lower-triangular CSC factor,
    L["xsuper"] the supernode boundaries, L["perm"] the elimination-order
    permutation -- 0-indexed in this port). Matches real fwblkslv.c's
    mexFunction exactly: `b` is gathered by `L["perm"]` (real/original
    index space -> the storage order L["L"]/L["d"] are expressed in)
    before forward-substituting; the result is left in storage order (no
    scatter on the way out -- that's sparbwslv's job on the way back).

    _native.fwsolve() mutates its `y` argument IN PLACE (by design, it
    wraps fwblkslv.c's own in-place solve) -- but MATLAB's sparfwslv/
    fwblkslv have ordinary value semantics (the caller's own `b` is
    never touched). Copying `b` here keeps that contract for every
    caller in this port (loopPcg/wrapPcg both reuse their `r` after
    calling sparfwslv) instead of requiring each call site to remember
    to copy -- forgetting this once (silently corrupting `r` under
    `Lr = fwdpr1(Lden, sparfwslv(L, r))`, then using the now-mutated `r`
    again for `rnew = r - alpha*tmp2`) is exactly what caused this
    port's first loopPcg/wrapPcg oracle mismatch.

    `b` may also be a 2-D dense array or a scipy-sparse matrix (multiple
    right-hand-side columns at once) -- the shape deninfac.py's
    `sparfwslv(L, Ad, symLden.LAD)` call needs to forward-solve the
    dense-column matrix `Ad` through `L`. Real fwblkslv.c's own MEX
    wrapper handles a sparse `b` with a dedicated sparsity-aware
    selfwsolve() kernel that needs a 3rd argument `y` (here: `ysymb`)
    giving the *exact sparsity pattern* of the output -- but restricted
    to that pattern, selfwsolve()'s VALUES are identical to an ordinary
    dense forward-solve (it's a performance optimization, skipping
    supernodes that can't affect the requested output positions, not a
    different computation -- confirmed by reading fwblkslv.c's
    mexFunction in full). So a sparse `b` is solved here by densifying,
    gathering rows by `L["perm"]`, and forward-solving column-by-column
    via the same raw fwsolve() kernel, then restricting the result to
    `ysymb`'s exact (row, col) support.
    """
    perm = np.asarray(L["perm"], dtype=np.int64).ravel()
    if sp.issparse(b):
        if ysymb is None:
            raise ValueError("sparfwslv: ysymb (3rd argument) is required for sparse b")
        B = np.asarray(b.todense(), dtype=np.float64)[perm, :]
        m, n = B.shape
        Y = np.empty((m, n), dtype=np.float64)
        for j in range(n):
            Y[:, j] = _native.fwsolve(L["L"], L["xsuper"], np.array(B[:, j], copy=True))
        ysymb_csc = ysymb.tocsc()
        rows = ysymb_csc.indices
        cols = np.repeat(np.arange(n), np.diff(ysymb_csc.indptr))
        data = Y[rows, cols]
        return sp.csc_matrix((data, rows.copy(), ysymb_csc.indptr.copy()), shape=(m, n))

    b_arr = np.asarray(b, dtype=np.float64)
    if b_arr.ndim == 1:
        return _native.fwsolve(L["L"], L["xsuper"], b_arr[perm].copy())
    m, n = b_arr.shape
    Bp = b_arr[perm, :]
    Y = np.empty((m, n), dtype=np.float64)
    for j in range(n):
        Y[:, j] = _native.fwsolve(L["L"], L["xsuper"], np.array(Bp[:, j], copy=True))
    return Y


def sparbwslv(L: dict, b):
    """y = sparbwslv(L,b): backward-solve y(L.perm) := L'\\b, i.e. `b` is
    already in storage order (typically straight from sparfwslv/fwdpr1/
    bwdpr1) and the result is scattered back to real/original index space
    by `L["perm"]` -- matching real bwblkslv.c's mexFunction exactly. See
    sparfwslv()'s docstring for why `b` is copied before the in-place
    native call."""
    perm = np.asarray(L["perm"], dtype=np.int64).ravel()
    z = _native.bwsolve(L["L"], L["xsuper"], np.array(b, dtype=np.float64, copy=True))
    y = np.empty_like(z)
    y[perm] = z
    return y


def _DAy_from_y(At, dense, d, K, y):
    """[sqrt(d.l).*Ap(1:K.l); asmDxq(d,Ap,K); psdscale(d,Ap,K)] for Ap =
    vecsym(Amul(At,dense,y,1),K) -- the common "DA'y" tail shared by
    loopPcg's k>1 path (once or twice, for a hi/lo quadadd pair)."""
    Kl = int(K["l"])
    Ap = _native.vecsym(amul(At, dense, y, transp=True), K)
    return np.concatenate([np.sqrt(d["l"]) * Ap[:Kl], asmDxq(d, Ap, K), psdscale(d, Ap, K)])


def loopPcg(L, Lden, At, dense, d, DAt, K, b, p, ssqrNew, cgpars, restol):
    """[y,k,DAy] = loopPcg(...): solve A*P(d)*A' * y = b via PCG with L
    as the (block sparse Cholesky) preconditioner. p=None starts PCG
    from scratch; ssqrNew is only used when p is not None (continuing a
    previous step's search direction)."""
    b = np.asarray(b, dtype=np.float64).ravel()
    Kl = int(K["l"])
    lorN = len(K.get("q", []))
    dense_q = np.asarray(dense.get("q", np.zeros(0))).ravel().astype(np.int64)

    k = 0
    r = b.copy()
    finew = 0.0
    y = None
    y_lo = None  # None means y (if set) is a plain array, not a hi/lo quad-precision pair
    normrmin = float(np.max(np.abs(r))) if r.size else 0.0
    ymin = ymin_lo = None
    alpha = None
    STOP = 0
    Ap = DApq = DAps = None

    while STOP == 0:
        Lr = _native.fwdpr1(Lden, sparfwslv(L, r))
        tmp = Lr / L["d"]
        if p is None:
            ssqrNew = float(Lr @ tmp)
            p = sparbwslv(L, _native.bwdpr1(Lden, tmp))
        else:
            ssqrOld = ssqrNew
            ssqrNew = float(Lr @ tmp)
            p = (ssqrNew / ssqrOld) * p
            p = p + sparbwslv(L, _native.bwdpr1(Lden, tmp))

        Ap = _native.vecsym(amul(At, dense, p, transp=True), K)
        DDAp, DApq, DAps, ssqrDAp = PopK(d, Ap, K)
        if ssqrDAp > 0.0:
            k += 1
            alpha = ssqrNew / ssqrDAp
            if y is not None:
                if y_lo is not None:
                    y, y_lo = _native.quadadd(y, y_lo, alpha * p)
                else:
                    y = y + alpha * p
            elif cgpars.get("qprec", 1) > 0:
                y = alpha * p
                y_lo = np.zeros_like(p)
            else:
                y = alpha * p

            tmp2 = amul(At, dense, DDAp, transp=False)
            if lorN:
                tmp2 = tmp2 + DAt["q"].T @ DApq
                if dense_q.size:
                    tmp2 = tmp2 + DAt["denq"] @ DApq[dense_q - 1]
            r = r - alpha * tmp2

            fiprev = finew
            if y_lo is not None:
                finew = float((b + r) @ y + (b + r) @ y_lo)
            else:
                finew = float((b + r) @ y)
            normr = float(np.max(np.abs(r))) if r.size else 0.0
            if normr < normrmin:
                ymin, ymin_lo = y, y_lo
                normrmin = normr
            if normr < restol:
                STOP = 1
            elif finew - fiprev < cgpars["stagtol"] * fiprev:
                STOP = 2
            elif k >= cgpars["maxiter"]:
                STOP = 2
        else:
            STOP = 1

    if STOP == 2:
        y, y_lo = ymin, ymin_lo

    if y is None:
        return None, k, None

    if k == 1:
        DAy = alpha * np.concatenate(
            [np.sqrt(d["l"]) * Ap[:Kl], asmDxq(d, Ap, K, ddotx=DApq), DAps]
        )
    else:
        DAy = _DAy_from_y(At, dense, d, K, y)
        if y_lo is not None:
            DAy = DAy + _DAy_from_y(At, dense, d, K, y_lo)

    return y, k, DAy


def wrapPcg(L, Lden, At, dense, d, DAt, K, rb, rv, cgpars, y0):
    """[y,dx,k,r] = wrapPcg(...): solve AP(d)A'*y = rb + AD(rv) (rb may
    be None/empty, matching the .m file's `if ~isempty(rb)` guard), with
    one unconditional CG-preconditioned step followed by refinement
    passes via loopPcg() as needed."""
    Kl = int(K["l"])
    restol = y0 * cgpars["restol"]

    dx = np.concatenate([np.sqrt(d["l"]) * rv[:Kl], asmDxq(d, rv, K), psdscale(d, rv, K, transp=True)])
    r = amul(At, dense, dx, transp=False)
    if rb is not None and rb.size:
        r = r + rb

    p = _native.fwdpr1(Lden, sparfwslv(L, r))
    y = p / L["d"]
    ssqrNew = float(p @ y)
    p = sparbwslv(L, _native.bwdpr1(Lden, y))
    x = _native.vecsym(amul(At, dense, p, transp=True), K)
    dx = np.concatenate([np.sqrt(d["l"]) * x[:Kl], asmDxq(d, x, K), psdscale(d, x, K)])
    ssqrdx = float(np.sum(dx**2))
    if ssqrdx <= 0.0:
        return np.zeros_like(r), rv, 0, r

    k = 1
    alpha = ssqrNew / ssqrdx
    y = alpha * p
    dx = rv - alpha * dx
    x = np.concatenate([np.sqrt(d["l"]) * dx[:Kl], asmDxq(d, dx, K), psdscale(d, dx, K, transp=True)])
    r = amul(At, dense, x, transp=False)
    if rb is not None and rb.size:
        r = r + rb
    normr = float(np.max(np.abs(r))) if r.size else 0.0
    if normr < restol:
        return y, dx, k, r

    trial = 0
    while True:
        dy, dk, x = loopPcg(L, Lden, At, dense, d, DAt, K, r, p, ssqrNew, cgpars, restol)
        if dy is None:
            return y, dx, k, r
        k = k + dk
        y = y + dy
        dx = dx - x
        x = np.concatenate(
            [np.sqrt(d["l"]) * dx[:Kl], asmDxq(d, dx, K), psdscale(d, dx, K, transp=True)]
        )
        r = amul(At, dense, x, transp=False)
        if rb is not None and rb.size:
            r = r + rb
        normr = float(np.max(np.abs(r))) if r.size else 0.0
        if normr < restol:
            return y, dx, k, r
        if trial >= cgpars["refine"]:
            return y, dx, k, r
        p = None
        trial += 1
