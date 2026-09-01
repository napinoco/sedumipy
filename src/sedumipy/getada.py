"""Port of getada.m: builds ADA = A*P(d)*A' for LP(+Lorentz) problems
directly via sparse linear algebra -- no C kernels, pure NumPy/SciPy.
Used both by optstep.m's LP-optimality-projection path (dead in
practice there, since sedumi.m only calls optstep.m when `lponly`,
which forces K.q and K.s empty too) AND, for real, by this port's own
sedumi.py main loop whenever `has_psd` is false -- i.e. every LP+SOCP
problem with no PSD blocks (K.s empty), where the Lorentz-block loop
below is very much live code.

DAt["q"] (lorN x m) and the ADA it feeds into are kept sparse
throughout: densifying either is what made getdatm.py's old
`DAt_q.todense()` OOM on large-lorN problems like DIMACS's
nql180/qssp180 (m ~ 1.3e5, lorN ~ 3.2e4 -- a dense DAt_q alone is
>30 GiB, and the dense ADA = DAt_q.T @ DAt_q this function used to
build would be m x m, far bigger still). A*P(d)*A' stays about as
sparse as A itself for these problems, so sparse @ sparse is both the
correct fix and the fast path.
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
    DAt_q = DAt_q.tocsr() if sp.issparse(DAt_q) else sp.csr_matrix(DAt_q)
    ADA = DAt_q.T @ DAt_q
    ADA = ADA + Alq.T @ sp.diags(scalingvector) @ Alq
    ADA = sp.csc_matrix(ADA)
    absd = np.asarray(ADA.diagonal()).ravel()
    return ADA, absd
