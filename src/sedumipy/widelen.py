"""Port of widelen.m: computes an approximate wide-region neighborhood
step length via a bounded bisection line search, doing the extensive
search only when it pays off (final rate at most twice the best
possible, step length at least half the best possible).

KNOWN UPSTREAM BUG NOT REPLICATED: widelen.m computes the Lorentz-block
eigenvalue term as

    tmp = halfxz.^2 - detxz;
    if all(tmp > 0)          % <- one global test for ALL blocks at once
        lab2q = halfxz + sqrt(tmp);
    else
        lab2q = halfxz;      % <- ...so ONE bad block degrades EVERY block
    end

`tmp` is not an arbitrary quantity that may legitimately go negative: the
two Lorentz eigenvalues this produces are `lab2q` and `detxz/lab2q`, whose
product is `detxz` and whose sum is `2*halfxz`, so `lab2q = halfxz +
sqrt(tmp)` means exactly `tmp = ((lab1-lab2)/2)^2` -- a perfect square,
hence >= 0 in exact arithmetic, dipping below zero only by rounding when
a block's two eigenvalues nearly coincide (and hitting *exactly* 0, which
`> 0` also rejects, whenever they coincide outright -- routine for
structured problems whose Lorentz blocks are copies of one another).
Upstream's `else` branch is a cheap safety net against handing sqrt() a
negative argument, but it is applied all-or-nothing: a single block
sitting on top of zero silently switches every other block onto the cruder
`lab2q = halfxz` formula too. This port instead clamps per block:

    lab2q = halfxz + sqrt(max(tmp, 0))

which is what upstream's own formula reduces to when `tmp <= 0` for that
one block (`sqrt(0) == 0`), keeps the accurate formula for every block
whose own discriminant is fine, and can never hand sqrt() a negative
argument -- i.e. it is strictly more accurate than either upstream branch
while preserving the safety property the `if` was there for.

Measured effect (DIMACS): nb_L2 goes from numerr=2 at iteration 10 to
numerr=0 at iteration 16, matching both the real Octave/MEX build's own
iteration count and the published objective -1.628972; qssp30old, which
the real build itself fails on with numerr=2, improves to numerr=1 with
cx and by agreeing to 7 digits. See CONTRIBUTING.md section 6.
"""

from __future__ import annotations

import numpy as np

from . import _native
from .cone import psdeig, psdfactor, psdscale, tdet
from .neighborhood import iswnbr


def _build_w(xM, zM, K: dict):
    Kl = int(K["l"])
    w: dict = {"tdetx": tdet(xM, K), "tdetz": tdet(zM, K)}
    detxz = w["tdetx"] * w["tdetz"] / 4

    lorN = len(K.get("q", []))
    if lorN == 0:
        lab2q = np.zeros(0)
    else:
        ix = np.asarray(K["mainblks"], dtype=np.int64).ravel()
        i1, i2, i3 = int(ix[0]), int(ix[1]), int(ix[2])
        halfxz = (
            xM[i1 - 1 : i2 - 1] * zM[i1 - 1 : i2 - 1]
            + _native.ddot(xM[i2 - 1 : i3 - 1], zM, K["qblkstart"])
        ) / 2
        tmp = halfxz**2 - detxz
        lab2q = halfxz + np.sqrt(np.maximum(tmp, 0.0))

    w["ux"], _ = psdfactor(xM, K)
    w["s"] = psdscale(w["ux"], zM, K)
    w["lab"] = np.concatenate([xM[:Kl] * zM[:Kl], detxz / lab2q, lab2q, psdeig(w["s"], K)])
    return w


def widelen(xc, zc, y0: float, dx, dz, dy0: float, d2y0: float, maxt: float, pars: dict, K: dict):
    """[t,wr,w] = widelen(xc,zc,y0,dx,dz,dy0,d2y0,maxt,pars,K)"""
    xc = np.asarray(xc, dtype=np.float64).ravel()
    zc = np.asarray(zc, dtype=np.float64).ravel()
    dx = np.asarray(dx, dtype=np.float64).ravel()
    dz = np.asarray(dz, dtype=np.float64).ravel()
    thetaSQR = pars["theta"] ** 2

    if dy0 < -1e-5 * y0:
        if d2y0 < 0:
            fullt = 2 * y0 / (-dy0 + np.sqrt(dy0**2 - 4 * y0 * d2y0))
        else:
            fullt = y0 / (-dy0)
        if fullt <= 0:
            raise AssertionError("widelen: fullt <= 0")
    else:
        fullt = 2 * maxt

    tR = min(maxt, fullt)
    if tR < 0:
        raise AssertionError("widelen: tR >= 0")

    t = 0.0
    ntry = 0
    w = wM = None
    deltaM = hM = alphaM = None

    while (t < 0.5 * tR) or ((fullt - tR) + (1e-7 * fullt) < (tR - t)) or ntry == 0:
        ntry = 1
        tM = 0.1 * t + 0.9 * tR if tR == maxt else 0.5 * (t + tR)
        xM = xc + tM * dx
        zM = zc + tM * dz
        wM = _build_w(xM, zM, K)
        deltaM, hM, alphaM = iswnbr(wM["lab"], thetaSQR)

        if (deltaM <= pars["beta"]) or ((tM < fullt / 10) and (deltaM < 1)):
            w = wM
            t = tM
            wr = {"h": hM, "alpha": alphaM, "delta": deltaM}
        else:
            tR = tM

    if t == 0:
        w = wM
        t = tM
        wr = {"h": hM, "alpha": alphaM, "delta": deltaM}

    wr["desc"] = 1  # always descent direction
    return t, wr, w
