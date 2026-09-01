"""Verify symbchol.py's branch selection between the ordmmd()+
symbolic_cholesky() path (sparse ADA) and the symbolic_cholesky_dense()
identity-permutation shortcut (fully dense ADA).

This is the exact `spars(ADA)==1` branch symbchol.py's own module
docstring calls out as this port's first real oracle mismatch against
Octave -- worth pinning directly, independent of the much larger
end-to-end pipelines (test_sedumi.py's dense case, test_symbcholden.py's
sparse case) that only exercise it incidentally. No Octave oracle is
needed here: both branches are checked against this port's own,
independently-tested _native primitives (ordmmd/symbolic_cholesky are
verified against Octave in test_ordmmd.py/test_symfct.py), so this test
isolates symbchol()'s own routing/wiring logic instead.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import scipy.io
import scipy.sparse as sp

sedumipy = pytest.importorskip("sedumipy")
from sedumipy import _native  # noqa: E402
from sedumipy.symbchol import symbchol  # noqa: E402

FIXTURE_DIR = Path(__file__).parent / "fixtures" / "ordmmd"


def _fully_dense_symmetric(m):
    rng = np.random.default_rng(0)
    X = rng.standard_normal((m, m))
    A = X @ X.T + m * np.eye(m)
    assert np.all(A != 0.0)
    return sp.csc_matrix(A)


def test_symbchol_dense_ADA_uses_identity_permutation_and_one_supernode():
    m = 6
    ADA = _fully_dense_symmetric(m)

    L = symbchol(ADA)

    np.testing.assert_array_equal(L["perm"], np.arange(m, dtype=np.int64))
    np.testing.assert_array_equal(L["xsuper"], np.array([0, m], dtype=np.int64))
    np.testing.assert_array_equal(
        np.asarray(L["L"].todense()), np.tril(np.ones((m, m)))
    )


def test_symbchol_dense_ADA_matches_symbolic_cholesky_dense_directly():
    m = 8
    ADA = _fully_dense_symmetric(m)

    L = symbchol(ADA)
    expected = _native.symbolic_cholesky_dense(m)

    np.testing.assert_array_equal(L["perm"], expected["perm"])
    np.testing.assert_array_equal(L["xsuper"], expected["xsuper"])
    assert L["tmpsiz"] == expected["tmpsiz"]


@pytest.mark.skipif(not FIXTURE_DIR.exists(), reason="ordmmd fixtures not generated")
@pytest.mark.parametrize("name", ["grid5x5", "path5", "rand10", "rand25", "star5"])
def test_symbchol_sparse_ADA_routes_through_ordmmd_and_symbolic_cholesky(name):
    data = scipy.io.loadmat(FIXTURE_DIR / f"{name}.mat")
    A = data["A"].tocsc()
    m = A.shape[0]
    density = A.nnz / (m * m)
    assert density < 1.0, "fixture must be sparse to exercise the non-dense branch"

    L = symbchol(A)

    perm0 = _native.ordmmd(A)
    expected = _native.symbolic_cholesky(A, perm0)

    np.testing.assert_array_equal(L["perm"], expected["perm"])
    np.testing.assert_array_equal(L["xsuper"], expected["xsuper"])
    got_L, exp_L = L["L"], expected["L"]
    assert got_L.shape == exp_L.shape
    np.testing.assert_array_equal(got_L.indptr, exp_L.indptr)
    np.testing.assert_array_equal(got_L.indices, exp_L.indices)

    # A non-trivial sparse ADA should not merely coincide with the dense
    # shortcut's trivial identity/one-supernode structure -- otherwise this
    # test would not actually be distinguishing the two branches.
    assert L["xsuper"].size - 1 > 1 or not np.array_equal(
        L["perm"], np.arange(m, dtype=np.int64)
    )


def test_symbchol_accepts_dense_ndarray_input_not_just_sparse():
    m = 5
    ADA_sparse = _fully_dense_symmetric(m)
    ADA_dense_array = np.asarray(ADA_sparse.todense())

    L_from_sparse = symbchol(ADA_sparse)
    L_from_dense = symbchol(ADA_dense_array)

    np.testing.assert_array_equal(L_from_sparse["perm"], L_from_dense["perm"])
    np.testing.assert_array_equal(L_from_sparse["xsuper"], L_from_dense["xsuper"])
