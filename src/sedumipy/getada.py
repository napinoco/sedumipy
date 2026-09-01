"""Port of getada.m: builds ADA = A*P(d)*A' for LP(+Lorentz) problems
directly via sparse linear algebra -- no C kernels, pure NumPy/SciPy.
Used both by optstep.m's LP-optimality-projection path (dead in
practice there, since sedumi.m only calls optstep.m when `lponly`,
which forces K.q and K.s empty too) AND, for real, by this port's own
sedumi.py main loop whenever `has_psd` is false -- i.e. every LP+SOCP
problem with no PSD blocks (K.s empty), where the Lorentz-block loop
below is very much live code.

DAt["q"] (lorN x m) and the ADA it feeds into can be either sparse or
dense, picked once per solve by sedumi.py's `is_dense` hint (from
getsymbada()'s structural ADA pattern density) and threaded through via
DAt["q"]'s own representation -- this function just follows whatever
form getDAtm() handed it. Sparse is required on large-lorN problems like
DIMACS's nql180/qssp180 (m ~ 1.3e5, lorN ~ 3.2e4 -- a dense DAt_q alone
is >30 GiB, and the dense ADA = DAt_q.T @ DAt_q this function used to
build unconditionally would be m x m, far bigger still); dense is
faster whenever ADA comes out dense/near-dense anyway (e.g. nb.mat,
m=123), where sparse @ sparse pays real bookkeeping overhead a plain
BLAS matmul doesn't. Both branches build bug-for-bug identical values;
ADA is always handed back as sp.csc_matrix regardless (numeric_cholesky
needs that either way).
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp


def getada(A, K: dict, d: dict, DAt: dict):
    """absd = getada(A,K,d,DAt): also returns the freshly built ADA
    matrix (unlike the .m file, which stores it into the global
    ADA_sedumi_ instead of returning it) -- callers here get it back as
    a normal return value: (ADA, absd)."""
    A = A.tocsc() if sp.issparse(A) else sp.csc_matrix(A)
    m = A.shape[1]

    ix3 = int(np.asarray(K["mainblks"]).ravel()[2])
    Alq = A[: ix3 - 1, :]

    lorN = len(K.get("q", []))
    ix2 = int(np.asarray(K["mainblks"]).ravel()[1]) if lorN else ix3
    Kl = int(K["l"])
    scalingvector = np.concatenate([d["l"], -d["det"], np.zeros(ix3 - ix2)])
    if lorN:
        qblkstart = np.asarray(K["qblkstart"], dtype=np.int64).ravel()
        for i in range(lorN):
            lo, hi = int(qblkstart[i]), int(qblkstart[i + 1])
            scalingvector[lo - 1 : hi - 1] = d["det"][i]

    DAt_q = DAt["q"]
    if sp.issparse(DAt_q):
        DAt_q = DAt_q.tocsr()
        ADA = DAt_q.T @ DAt_q
        ADA = ADA + Alq.T @ sp.diags(scalingvector) @ Alq
    else:
        DAt_q = np.asarray(DAt_q)
        ADA = DAt_q.T @ DAt_q
        ADA = ADA + (Alq.T @ sp.diags(scalingvector) @ Alq).toarray()
    ADA = sp.csc_matrix(ADA)
    absd = np.asarray(ADA.diagonal()).ravel()
    return ADA, absd
