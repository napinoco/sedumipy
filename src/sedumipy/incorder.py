"""Port of incorder.c's real kernel (spPartTransp()+incorder() -- its
mexFunction is only argument marshalling). NOT wrapped via ctypes,
deliberately: incorder.c's greedy pivot loop has no qsort/UB dependency
(unlike sortnnz.c/iswnbr.c), but it is a small, easily-reimplemented
algorithm, and this project already prefers a direct Python port for such
deterministic kernels (see neighborhood.py, sortnnz() in _native.py) over
round-tripping through ctypes marshalling for a function this size.

SEE ALSO getada3 (this is `incorder`'s one real caller in sedumi.m's
main loop: `[Aord.sperm, Aord.dz] = incorder(A, Ablkjc(:,3),
K.mainblks(3))`), dpr1fact.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp


def incorder(At, Ajc1=None, ifirst: int = 1):
    """perm, dz = incorder(At, Ajc1=None, ifirst=1): greedily orders the
    columns of `At` (a sparse matrix, `lenfull` rows x `m` columns, in
    SeDuMi's usual per-constraint-column convention) by ascending number
    of not-yet-covered "PSD" nonzeros -- i.e. nonzeros in the row range
    `ifirst:lenfull` (1-indexed, inclusive of `ifirst`). Each pivot step
    picks the first not-yet-ordered column with the fewest remaining
    uncovered subscripts in that range, then marks its subscripts
    covered.

    `Ajc1` (length m, 0-indexed positions into `At`'s CSC `indices`/
    `data` arrays -- matching `_native.partitA()`'s own output
    convention, e.g. `Ablkjc[:, 2]` for the 3rd row-block boundary)
    restricts each column to `Ajc1[j]:end` before counting nonzeros;
    default is the whole column (equivalent to `ifirst=1` covering
    every row). `ifirst` is 1-indexed, matching `_native.extractA`'s/
    `_native.findblks`'s own `blkstart` convention.

    Returns:
      perm -- length-m int64 array, a 1-indexed permutation (SeDuMi's
        usual constraint-index convention, matching `_native.sortnnz`'s
        output): perm[k] is the (k+1)-th column processed.
      dz -- `lenfull x m` sparse 0/1 CSC matrix; column k (0-indexed,
        i.e. corresponding to constraint `perm[k]`) lists the row
        subscripts newly covered at step k, incremental w.r.t. all
        previous columns' subscripts.
    """
    A = At.tocsc()
    lenfull, m = A.shape
    first = int(ifirst) - 1
    lenud = lenfull - first

    indptr = A.indptr
    indices = A.indices
    atjc1 = (
        np.ascontiguousarray(indptr[:m], dtype=np.int64)
        if Ajc1 is None
        else np.ascontiguousarray(Ajc1, dtype=np.int64)
    )
    atjc2 = np.ascontiguousarray(indptr[1:], dtype=np.int64)

    # Adjacency of "PSD" row subscript -> columns touching it (spPartTransp).
    rows_list = []
    cols_list = []
    for j in range(m):
        seg = indices[atjc1[j] : atjc2[j]]
        rows_list.append(seg.astype(np.int64) - first)
        cols_list.append(np.full(seg.size, j, dtype=np.int64))
    rows = np.concatenate(rows_list) if rows_list else np.zeros(0, dtype=np.int64)
    cols = np.concatenate(cols_list) if cols_list else np.zeros(0, dtype=np.int64)
    order = np.argsort(rows, kind="stable")
    rows_sorted = rows[order]
    cols_sorted = cols[order]
    row_ptr = np.searchsorted(rows_sorted, np.arange(lenud + 1))

    # Greedy pivoting (incorder()).
    iwork = (atjc2 - atjc1).copy()
    discard = np.zeros(lenud, dtype=bool)
    perm = np.arange(m, dtype=np.int64)
    dzjc = np.zeros(m + 1, dtype=np.int64)
    dzir = np.zeros(lenud, dtype=np.int64)

    for k in range(m):
        kmin = k
        lenmin = iwork[perm[k]]
        for j in range(k + 1, m):
            if iwork[perm[j]] < lenmin:
                lenmin = iwork[perm[j]]
                kmin = j
        perm[k], perm[kmin] = perm[kmin], perm[k]
        permk = int(perm[k])

        jnz = dzjc[k]
        for inz in range(atjc1[permk], atjc2[permk]):
            i = int(indices[inz]) - first
            if not discard[i]:
                discard[i] = True
                dzir[jnz] = i
                jnz += 1
        dzjc[k + 1] = jnz

        for jj in range(dzjc[k], dzjc[k + 1]):
            i = int(dzir[jj])
            for inz in range(row_ptr[i], row_ptr[i + 1]):
                iwork[cols_sorted[inz]] -= 1

    nnz = int(dzjc[m])
    dz = sp.csc_matrix(
        (np.ones(nnz, dtype=np.float64), dzir[:nnz] + first, dzjc), shape=(lenfull, m)
    )
    return perm + 1, dz
