"""Port of getDAtm.m: computes DAt.q[k,j] = d[k]'*Aj[k] for each Lorentz
block k and constraint column j -- the per-iteration "scaled A" quantity
loopPcg.m/PopK.m/sddir.m need for the Lorentz part -- plus the
dense-Lorentz-block correction DAt.denq (via `_native.adendotd`), which
folds the dense.q blocks' contribution into a separate Woodbury-update
term and zeroes their rows out of DAt.q proper.

Implemented directly with NumPy/SciPy sparse arithmetic rather than by
wrapping extractA()/ddot()'s C kernels for the DAt.q part: ddot()'s
existing binding (_native.py, cluster 1) only wraps ddotxj() (the
dense-X path); getDAtm.m's own `ddot(d.q2, A, K.qblkstart, Ablkjc)` call
uses the sparse-X path (spddotxj(), not bound) precisely because A is
the full sparse constraint matrix here -- reimplementing the same
reduction via sparse matrix algebra avoids needing that binding at all.
The DAt.denq part instead wraps the real `adendotd()` C kernel directly
(_native.adendotd), since that computation is intricate enough (the
Sherman-Morrison-Woodbury dense-column bookkeeping) that reimplementing
it in NumPy would risk silently drifting from the reference algorithm.

`Ablkjc` is `_native.partitA(A, K["mainblks"])`'s output -- unused by
this port's DAt.q computation (which builds it via sparse slicing
instead), kept only for signature fidelity with getDAtm.m. `DAtdenq` is
the previous DAt.denq (or getdense()'s own Adotdden the very first
time) -- only its sparsity pattern (shape m x len(dense.q)) is used, as
a template `adendotd()` fills in-place; see _native.adendotd's
docstring.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from . import _native


def getDAtm(A, Ablkjc, dense: dict, DAtdenq, d: dict, K: dict) -> dict:
    """DAt = getDAtm(A,Ablkjc,dense,DAtdenq,d,K): A is the internal-format
    At (N x m, SeDuMi's own storage convention: rows are the N cone
    variables, columns are the m constraints)."""
    lorN = len(K.get("q", []))
    m = A.shape[1]
    A = A.tocsr() if sp.issparse(A) else sp.csr_matrix(A)

    ix = np.asarray(K["mainblks"], dtype=np.int64).ravel()
    i1, i2 = int(ix[0]), int(ix[1])
    trace_rows = np.arange(i1 - 1, i2 - 1)
    A_trace = A[trace_rows, :]
    q1 = np.asarray(d.get("q1", np.zeros(lorN))).ravel()
    DAt_q = sp.diags(q1, shape=(lorN, lorN)) @ A_trace

    qblkstart = np.asarray(K["qblkstart"], dtype=np.int64).ravel()
    vec_start = int(qblkstart[0]) - 1
    vec_end = int(qblkstart[-1]) - 1
    if lorN > 0 and vec_end > vec_start:
        vec_rows = np.arange(vec_start, vec_end)
        block_lens = np.diff(qblkstart)
        block_of = np.repeat(np.arange(lorN), block_lens)
        A_vec = A[vec_rows, :]
        W = sp.csr_matrix(
            (np.asarray(d["q2"]).ravel(), (block_of, np.arange(vec_rows.size))),
            shape=(lorN, vec_rows.size),
        )
        DAt_q = DAt_q + W @ A_vec

    DAt_q = DAt_q.tocsr()

    dense_q = np.asarray(dense.get("q", np.zeros(0, dtype=np.int64))).ravel().astype(np.int64)
    if dense_q.size:
        adotd_in = sp.csc_matrix(DAt_q[dense_q - 1, :].T)
    else:
        adotd_in = sp.csc_matrix((m, 0))
    denq = _native.adendotd(dense, d, adotd_in, DAtdenq, qblkstart)

    if dense_q.size:
        keep = np.ones(lorN)
        keep[dense_q - 1] = 0.0
        DAt_q = sp.diags(keep) @ DAt_q

    return {"q": DAt_q, "denq": denq}
