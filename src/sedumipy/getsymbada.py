"""Port of getsymbada.m: builds the one-time (pre-main-loop) 0/1 sparsity
pattern of ADA = A*P(d)*A' used for the `K.s != 0` path's symbolic
Cholesky factorization (`symbchol.py`). Pure MATLAB (no MEX component),
built here from the already-bound `_native.extractA`/`_native.findblks`.

SEE ALSO sedumi, partitA, getada1, getada2.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from . import _native


def _spars(mat) -> float:
    size = mat.shape[0] * mat.shape[1]
    return (mat.nnz / size) if size else 0.0


def getsymbada(At, Ablkjc, q_pattern, psdblkstart):
    """SYMBADA = getsymbada(At, Ablkjc, q_pattern, psdblkstart): m x m
    0/1 sparse pattern of ADA. `Ablkjc` is `partitA(At, K["mainblks"])`'s
    output. `q_pattern` is the pre-loop Lorentz ddotA sparsity pattern
    (this port's `build_aord()` builds it as `Aord["q_pattern"]`) --
    NOT the same as `getdatm.getDAtm()`'s per-iteration numeric
    `DAt["q"]`, despite the real `sedumi.m` reusing (and later
    clobbering) that field name for both. `psdblkstart` is
    `K["sblkstart"]`.

    Falls back to a fully-dense m x m pattern (matching the real .m's
    `spars(...)==1` / `spars(...)>0.9` short-circuits) whenever any
    intermediate pattern is already that dense -- avoids computing a
    products that would be dense anyway.
    """
    m = At.shape[1]
    dense_pattern = sp.csc_matrix(np.ones((m, m)))

    Alpq = _native.extractA(At, Ablkjc, 0, 3, (1, int(psdblkstart[0])))
    Alpq.data[:] = 1.0
    Ablks = _native.findblks(At, Ablkjc, 3, None, psdblkstart)

    have_q = q_pattern is not None and q_pattern.shape[0] > 0
    if _spars(Ablks) == 1.0 or _spars(Alpq) == 1.0 or (have_q and _spars(q_pattern) == 1.0):
        return dense_pattern

    if q_pattern is None:
        symbada = sp.csc_matrix((m, m))
    else:
        symbada = (q_pattern.T @ q_pattern).tocsc()
    if _spars(symbada) > 0.9:
        return dense_pattern

    symbada = (symbada + Alpq.T @ Alpq).tocsc()
    if _spars(symbada) > 0.9:
        return dense_pattern

    symbada = (symbada + Ablks.T @ Ablks).tocsc()
    return symbada
