"""Port of maxstep.m: computes the maximal step length to the boundary
of the cone K along direction dx from point x.

KNOWN UPSTREAM BUG NOT REPLICATED: maxstep.m guards its Lorentz square
root with the same global all-or-nothing test widelen.m uses (see
widelen.py's docstring), but with a fallback that is worse than
widelen.m's cruder-formula one:

    norm2 = reltr.^2 - tdet(dx,K).*auxx.tdet;
    if all(norm2 > 0)
        norm2 = sqrt(norm2);     % <- skipped entirely if ANY block fails
    end
    mindxq = min( (reltr - norm2)./auxx.tdet);

When the test fails, `norm2` is left *unrooted* and then subtracted from
`reltr` anyway -- mixing a squared quantity into a linear one, for every
block, not just the offending one. maxstep.m's own comment gives the
identity that makes this a bug rather than a judgement call:

    (lab2-lab1)^2 = (tr y)^2 - 4 det y = [(x'Jdx)^2 - tdetx*tdetdx]/detx^2

i.e. `norm2` is a perfect square, non-negative in exact arithmetic, and
non-positive only by rounding (or exactly 0 when a block's two
eigenvalues coincide, which `> 0` rejects too). This port therefore
clamps per block, `sqrt(max(norm2, 0))`, which is the same value upstream
computes whenever its own test passes.

This one is a safety fix, not just an accuracy one. For a discriminant
below 1 -- the common case once the problem is scaled -- `v < sqrt(v)`,
so leaving `norm2` unrooted makes `reltr - norm2` too large, and the step
length this function returns is an *over*estimate of the distance to the
cone boundary. Instrumented on DIMACS nb_L2, the fallback fires on 5 of
64 calls and overestimates the step every one of those 5 times.
"""

from __future__ import annotations

import numpy as np

from . import _native
from .cone import minpsdeig, psdinvscale, tdet


def maxstep(dx, x, auxx: dict, K: dict):
    """tp = maxstep(dx,x,auxx,K): auxx.tdet (Lorentz "t-determinant" of
    x, doubled -- see wregion.m's own `uxc.tdet = 2*vfrm.lab(...)`
    construction) and auxx.u (psdfactor(x,K)) are the per-block
    "reference point" data the step-length formulas are built around."""
    dx = np.asarray(dx, dtype=np.float64).ravel()
    x = np.asarray(x, dtype=np.float64).ravel()
    Kl = int(K["l"])

    mindx = float(np.min(dx[:Kl] / x[:Kl])) if Kl else np.inf

    if len(K.get("q", [])):
        ix = np.asarray(K["mainblks"], dtype=np.int64).ravel()
        i1, i2, i3 = int(ix[0]), int(ix[1]), int(ix[2])
        reltr = x[i1 - 1 : i2 - 1] * dx[i1 - 1 : i2 - 1] - _native.ddot(
            x[i2 - 1 : i3 - 1], dx, K["qblkstart"]
        )
        auxx_tdet = np.asarray(auxx["tdet"]).ravel()
        norm2 = reltr**2 - tdet(dx, K) * auxx_tdet
        norm2 = np.sqrt(np.maximum(norm2, 0.0))
        mindxq = float(np.min((reltr - norm2) / auxx_tdet))
        mindx = min(mindx, mindxq)

    if len(K.get("s", [])):
        reldx = psdinvscale(auxx["u"], dx, K)
        mindxs = minpsdeig(reldx, K)
        mindx = min(mindx, mindxs)

    return 1.0 / max(-mindx, 1e-16)
