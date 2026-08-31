"""Port of deninfac.m: factors the dense columns into the PCG
preconditioner, producing the product-form rank-1 update `Lden` (via
dpr1fact) that lets wrapPcg/loopPcg apply `(L*diag(Ld)*L') + Ad*diag(smult)*Ad'`
as a preconditioner without ever densifying the sparse Cholesky factor
`L` itself.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from . import _native
from .pcg import sparfwslv


def deninfac(symLden, L: dict, dense: dict, DAt: dict, d: dict, absd, qblkstart, pars: dict):
    """[Lden,Ld] = deninfac(symLden,L,dense,DAt,d,absd,qblkstart,pars)"""
    cols = np.asarray(dense.get("cols", np.zeros(0))).ravel().astype(np.int64)

    if cols.size:
        nl = int(dense["l"])
        q_field = np.asarray(dense["q"], dtype=np.int64).ravel()
        nq = q_field.size
        i1 = nl
        i2 = nl + nq

        A = dense["A"]
        A = A.tocsc() if sp.issparse(A) else sp.csc_matrix(A)

        Ad = sp.hstack([A[:, :i1], DAt["denq"], A[:, i2:], A[:, i1:i2]], format="csc")

        detd = np.asarray(d["det"], dtype=np.float64).ravel()
        smult = np.concatenate([
            np.asarray(d["l"], dtype=np.float64).ravel()[cols[:nl] - 1],
            np.ones(nq, dtype=np.float64),
            np.asarray(_native.adenscale(dense, d, qblkstart), dtype=np.float64).ravel(),
            -detd[q_field - 1],
        ])

        LAD = sparfwslv(L, Ad, symLden["LAD"])
        Lden, Ld = _native.dpr1fact(LAD, L["d"], symLden, smult, pars.get("maxuden", 10))
        Lden["dz"] = symLden["dz"]
        Lden["first"] = symLden["first"]
        Lden["perm"] = symLden["perm"]
    else:
        Lden = {"betajc": np.array([0]), "rowperm": np.arange(1, L["d"].size + 1)}
        Ld = L["d"].copy()

    skip = np.asarray(L["skip"], dtype=np.int64).ravel()
    if skip.size:
        perm = np.asarray(L["perm"], dtype=np.int64).ravel()
        dtol = pars["canceltol"] * np.asarray(absd)[perm[skip]]
        dtol = np.maximum(dtol, pars["abstol"])
        mask = Ld[skip] <= dtol
        Ld[skip[mask]] = 1.0

    return Lden, Ld
