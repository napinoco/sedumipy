"""PSD (K.s != 0) counterpart to getada.py: the one-time pre-main-loop
`Aord`/sparsity-pattern setup (`build_aord`, mirroring sedumi.m lines
~330-382) and the per-iteration ADA rebuild (`getada_psd`, mirroring
sedumi.m lines ~450-452's `getada1`->`getada2`->`getada3` sequence).

SEE ALSO sedumi, incorder, getsymbada, _native.getada1/getada2/getada3.
"""

from __future__ import annotations

import numpy as np

from . import _native
from .getsymbada import getsymbada
from .incorder import incorder


def build_aord(A, K: dict):
    """Ablkjc, Aord, symbada = build_aord(A, K): the one-time setup real
    sedumi.m does before its main loop starts (before `symbchol`).
    `dense.q`'s zeroing-out of `DAt.q` rows is not ported -- always a
    no-op in this port, since `sedumi.py` always uses an empty `dense`
    (see its module docstring).

    `Aord` has "lqperm", "qperm", "sperm" (1-indexed permutations),
    "dz" (incorder()'s covered-subscript output), and "q_pattern" (the
    pre-loop Lorentz ddotA sparsity pattern -- `None` when K has no
    Lorentz blocks; NOT the same as getdatm.getDAtm()'s per-iteration
    numeric `DAt["q"]`, see getsymbada.py's docstring).
    """
    A = A.tocsc()
    m = A.shape[1]
    mainblks = np.asarray(K["mainblks"], dtype=np.int64).ravel()
    sblkstart = np.asarray(K["sblkstart"], dtype=np.int64).ravel()

    Ablkjc = _native.partitA(A, mainblks)

    Aord: dict = {"lqperm": _native.sortnnz(A, None, Ablkjc[:, 2])}

    lorN = len(K.get("q", []))
    q_pattern = None
    if lorN:
        qblkstart = np.asarray(K["qblkstart"], dtype=np.int64).ravel()
        q_pattern = _native.findblks(A, Ablkjc, 2, 3, qblkstart)
        if q_pattern.nnz:
            Alp = _native.extractA(A, Ablkjc, 1, 2, (int(mainblks[0]), int(mainblks[1])))
            Alp.data[:] = 1.0
            q_pattern = (q_pattern + Alp).tocsc()
            Aord["qperm"] = _native.sortnnz(q_pattern, None, None)
        else:
            Aord["qperm"] = np.arange(1, m + 1, dtype=np.int64)
    else:
        Aord["qperm"] = np.arange(1, m + 1, dtype=np.int64)
    Aord["q_pattern"] = q_pattern

    Aord["sperm"], Aord["dz"] = incorder(A, Ablkjc[:, 2], int(mainblks[2]))

    symbada = getsymbada(A, Ablkjc, q_pattern, sblkstart)

    return Ablkjc, Aord, symbada


def getada_psd(ADA, A, Ablkjc, Aord: dict, DAt: dict, d: dict, K: dict):
    """(ADA, absd) = getada_psd(ADA, A, Ablkjc, Aord, DAt, d, K): rebuilds
    ADA = A*P(d)*A' for one main-loop iteration when K.s is nonempty,
    via getada1 (fresh LP+Lorentz-diagonal part) -> getada2 (+= Lorentz
    ddotA part) -> getada3 (+= PSD part, then symmetrize), exactly as
    sedumi.m's main loop calls them. `d` needs "l", "det", "u" (see
    sdinit.py/updtransfo.py); `d["perm"]` is passed through to
    invcholfac() -- but as `None` when empty (sdinit.py's own
    `d["perm"]` initialization on the very first iteration is a length-0
    array, not `None`, and invcholfac() treats "perm given" as "index
    each block's `nk`-length slice out of it", which would read past the
    end of a length-0 array).
    """
    Ajc = Ablkjc[:, 2]
    ADA = _native.getada1(ADA, A, Ajc, Aord["lqperm"], d, K["qblkstart"])
    ADA = _native.getada2(ADA, DAt, Aord, K)
    perm = d.get("perm")
    if perm is not None and np.asarray(perm).size == 0:
        perm = None
    udsqr = _native.invcholfac(d["u"], K, perm)
    ADA, absd = _native.getada3(ADA, A, Ajc, Aord, udsqr, K)
    return ADA, absd
