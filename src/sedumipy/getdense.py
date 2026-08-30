"""Port of getdense.m: the dense-column *detection* heuristic that feeds
the dense-columns optimization (see CONTRIBUTING.md's dense.cols/dense.q
notes). Pure `.m` logic (no MEX) built entirely on top of the already-bound
`_native.extractA`/`_native.findblks`.

Terminology note (easy to get backwards): `dense.cols`/`dense.l` are
1-indexed GLOBAL row-subscripts into `At` (i.e. subscripts of the LP and
Lorentz "primitive" variables, restricted to `1..K.lq` -- PSD variables can
never be flagged dense since PSD blocks aren't allowed to be split up).
`dense.q`, by contrast, holds 1-indexed indices LOCAL to the Lorentz-block
numbering `1..len(K.q)` (i.e. "the k-th Lorentz block", not "row k of
At") -- this matches `adendotd.c`'s own convention where `q` indexes
per-block arrays via `K.qblkstart`. Do not treat `dense.q` as a row
subscript the way `dense.cols` is.

SEE ALSO sedumi, getdatm (consumes dense.q via adendotd), deninfac
(consumes dense.cols via adenscale/dpr1fact), amul (consumes dense.cols/
dense.A directly).
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from . import _native

_NORMDEN = 5


def getdense(At, Ablkjc, K: dict, pars: dict):
    """dense, Adotdden = getdense(At, Ablkjc, K, pars): flags a small
    proportion of LP/Lorentz "columns" (rows of `At`) that appear in a
    large proportion of the primal constraints as dense, mirroring
    getdense.m exactly. `dense` is a dict with "l" (int), "cols", "q"
    (1-indexed int64 arrays, see module docstring for their differing
    index conventions). `Adotdden` is the `(m, len(dense.q))` sparse
    0/1 pattern `Ablkq(dense.q, :)'` upstream also returns.
    """
    N, m = At.shape
    lq = int(K["lq"])
    mainblks = np.asarray(K["mainblks"], dtype=np.int64).ravel()
    i1, i2 = int(mainblks[0]), int(mainblks[1])

    Apart = _native.extractA(At, Ablkjc, 0, 3, (1, lq + 1))
    colnz = Apart.getnnz(axis=1).astype(np.float64)

    sblkstart = np.asarray(K["sblkstart"], dtype=np.int64).ravel()
    blknz = np.asarray(_native.findblks(At, Ablkjc, 3, None, sblkstart).sum(axis=1)).ravel()
    h = max(_NORMDEN, int(blknz.max()) if blknz.size else 0)

    # Replace colnz-entries for the Lorentz-trace ("x1") variables by
    # nnz-constraints for that Lorentz block: a block can only be removed
    # as a whole (trace and norm-bound part together).
    Ablkq = _native.extractA(At, Ablkjc, 1, 2, (i1, i2))
    if i1 < i2:
        qblkstart = np.asarray(K["qblkstart"], dtype=np.int64).ravel()
        Ablkq2 = _native.findblks(At, Ablkjc, 2, 3, qblkstart)
        pattern = (Ablkq != 0).astype(np.float64) + Ablkq2
        pattern.data[:] = 1.0
        Ablkq = pattern.tocsc()
        colnz[i1 - 1 : i2 - 1] = np.asarray(Ablkq.sum(axis=1)).ravel()

    # Find the denq-quantile for dense columns: e.g. with denq=0.75,
    # denf=10, find the 75% quantile spquant, tag anything denser than
    # 10*spquant as dense. spquant is chosen so all columns with
    # nnz <= h are left of it (PSD blocks are never split up).
    bigcolnz = colnz[colnz > h]
    denqN = int(np.ceil(pars["denq"] * colnz.size)) - (N - bigcolnz.size)
    if denqN < 1:
        spquant = float(h)
    else:
        spquant = float(np.sort(bigcolnz)[denqN - 1])

    threshold = pars["denf"] * spquant
    dense_cols = np.nonzero(colnz > threshold)[0].astype(np.int64) + 1
    local_q = colnz[i1 - 1 : i2 - 1]
    dense_q = np.nonzero(local_q > threshold)[0].astype(np.int64) + 1
    dense_l = int(np.sum(dense_cols < i1))

    # Dense columns should be few relative to the number of constraints --
    # otherwise a single full Cholesky factorization is preferable.
    if dense_cols.size > m / 2:
        dense_l = 0
        dense_cols = np.zeros(0, dtype=np.int64)
        dense_q = np.zeros(0, dtype=np.int64)

    if dense_q.size == 0:
        Adotdden = sp.csc_matrix((m, 0))
    else:
        Adotdden = Ablkq[dense_q - 1, :].T.tocsc()

    dense = {"l": dense_l, "cols": dense_cols, "q": dense_q}
    return dense, Adotdden
