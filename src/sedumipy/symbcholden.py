"""Port of symbcholden.m: symbolic forward Cholesky of the dense columns,
producing the `Lden` structure `deninfac.py`'s dpr1fact-based numeric
factorization will (eventually) consume.

`dense["A"]`'s columns are ordered `[LP-dense (nl), Lorentz-trace
placeholders (nq), Lorentz-norm-bound (nden)]` (see getdense.py's module
docstring for `dense.cols`'/`dense.q`'s differing index conventions,
which this same [nl, nq, nden] split mirrors). `DAt["denq"]` is the
`(m, nq)` sparse matrix upstream calls `DAt.denq` -- only its sparsity
*pattern* matters here (symbfwblk only reads pattern), so the very first
call (before any `adendotd()` numeric values exist) can pass
`getdense()`'s own `Adotdden` return value, exactly as sedumi.m does
(`[dense,DAt.denq] = getdense(...)` then later `symLden =
symbcholden(L,dense,DAt)`, with the real numeric `adendotd()` call
overwriting `DAt.denq`'s values only afterwards, inside getDAtm.m).

`L` is this port's usual symbchol()-produced dict (`L["perm"]`/
`L["xsuper"]` 0-indexed, matching symbchol.py/pcg.py's convention) --
_native.symbfwblk expects the raw 1-indexed Octave/MEX convention, so
this module adds the +1 offset at the call site rather than changing
symbfwblk's own already oracle-verified interface.

SEE ALSO symbfwblk (_native.py), incorder, finsymbden (_native.py).
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from . import _native
from .incorder import incorder


def symbcholden(L: dict, dense: dict, DAt: dict) -> dict:
    """Lden = symbcholden(L, dense, DAt): see module docstring."""
    nl = int(dense["l"])
    nq = int(np.asarray(dense["q"]).size)
    i1 = nl
    i2 = nl + nq

    A = dense["A"]
    A = A.tocsc() if sp.issparse(A) else sp.csc_matrix(A)

    L_1idx = {
        "L": L["L"],
        "perm": np.asarray(L["perm"], dtype=np.int64) + 1,
        "xsuper": np.asarray(L["xsuper"], dtype=np.int64) + 1,
    }

    LAD = sp.hstack(
        [
            _native.symbfwblk(L_1idx, A[:, :i1]),
            _native.symbfwblk(L_1idx, DAt["denq"]),
            _native.symbfwblk(L_1idx, A[:, i2:]),
            _native.symbfwblk(L_1idx, A[:, i1:i2]),
        ],
        format="csc",
    )

    nperm = int(np.asarray(dense["cols"]).size)
    perm, dz = incorder(LAD[:, :nperm])

    firstq = nl + 1
    return _native.finsymbden(LAD, perm, dz, firstq)
