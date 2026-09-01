"""Phase 4: .mat I/O (matio.py) -- not an oracle-comparison test like the
rest of this suite (there's no .m file to port here, see matio.py's own
docstring), just round-trip/parsing checks against fixtures already used
elsewhere in this suite."""

from pathlib import Path

import numpy as np
import pytest
import scipy.io

sedumipy = pytest.importorskip("sedumipy")

FIXTURE_DIR = Path(__file__).parent / "fixtures" / "sedumi"
EXAMPLES_DIR = (
    Path(__file__).parent.parent / "vendor" / "sedumi-upstream" / "examples"
)


@pytest.mark.skipif(not FIXTURE_DIR.exists(), reason="sedumi fixtures not generated")
def test_read_mat_matches_hand_unwrapped_fixture():
    """lp_socp_sdp_dense_feasible.mat has a "pars" field (pars.denf=3,
    see test_sedumi.py's own docstring) -- exercises that path too."""
    path = FIXTURE_DIR / "lp_socp_sdp_dense_feasible.mat"
    A, b, c, K, pars = sedumipy.read_mat(path)

    raw = scipy.io.loadmat(path)
    np.testing.assert_allclose(np.asarray(A.todense() if hasattr(A, "todense") else A), raw["At"])
    np.testing.assert_allclose(b, np.asarray(raw["b"]).ravel())
    np.testing.assert_allclose(c, np.asarray(raw["c"]).ravel())
    assert K["l"] == int(raw["K"][0, 0]["l"].item())
    assert pars["denf"] == 3


@pytest.mark.skipif(not FIXTURE_DIR.exists(), reason="sedumi fixtures not generated")
@pytest.mark.parametrize(
    "name",
    ["lp_feasible", "socp_feasible", "sdp_feasible", "sdp_mixed_cones_feasible"],
)
def test_read_mat_solves_correctly(name):
    A, b, c, K, pars = sedumipy.read_mat(FIXTURE_DIR / f"{name}.mat")
    x, y, info = sedumipy.sedumi(A, b, c, K, **pars)
    assert info["numerr"] in (0, 1)


@pytest.mark.skipif(not EXAMPLES_DIR.exists(), reason="vendor submodule not checked out")
def test_read_mat_handles_at_orientation_and_sparse_b_c():
    """vendor's own example problems store the matrix as "At" and b/c as
    sparse column vectors -- this is exactly what surfaced pretransfo.py's
    sparse-b/c bug (CONTRIBUTING.md section 6), so it's worth its own
    explicit check here rather than relying on it only indirectly through
    test_golden_end_to_end.py."""
    A, b, c, K, pars = sedumipy.read_mat(EXAMPLES_DIR / "nb.mat")
    assert b.ndim == 1 and c.ndim == 1
    assert pars == {}
    x, y, info = sedumipy.sedumi(A, b, c, K)
    assert info["numerr"] == 0


def test_read_mat_treats_empty_K_l_and_K_f_as_zero(tmp_path):
    """Some real-world .mat files (e.g. DIMACS's examples/dimacs/data/
    TRUSS/truss5.mat.gz) store K.l/K.f as an explicit-but-empty MATLAB
    array rather than omitting the field or storing 0 -- previously this
    unwrapped to an empty NumPy array, which crashed pretransfo.py's
    `K.get("l", 0) or 0` on the empty array's ambiguous truth value."""
    K_struct = np.zeros((1, 1), dtype=[("l", "O"), ("f", "O"), ("s", "O")])
    K_struct["l"][0, 0] = np.zeros((0, 0))
    K_struct["f"][0, 0] = np.zeros((0, 0))
    K_struct["s"][0, 0] = np.array([[2]])

    path = tmp_path / "empty_Kl.mat"
    scipy.io.savemat(
        path,
        {
            "At": np.eye(4),
            "b": np.array([1.0, 0.0, 0.0, 1.0]),
            "c": np.array([1.0, 0.0, 0.0, 1.0]),
            "K": K_struct,
        },
    )

    A, b, c, K, pars = sedumipy.read_mat(path)
    assert K["l"] == 0
    assert K["f"] == 0
    np.testing.assert_array_equal(K["s"], [2])

    x, y, info = sedumipy.sedumi(A, b, c, K)
    assert info["numerr"] in (0, 1)


def test_read_mat_transparently_gunzips_gz_paths(tmp_path):
    """DIMACS ships its .mat problem files gzip-compressed
    (examples/dimacs/data/*/*.mat.gz); scipy.io.loadmat's own gzip
    auto-detection does not reliably trigger on every real .mat.gz file
    from that collection, so read_mat() must gunzip a ".gz" path itself
    rather than delegating to scipy."""
    import gzip

    K_struct = np.zeros((1, 1), dtype=[("l", "O")])
    K_struct["l"][0, 0] = np.array([[2]])

    raw_path = tmp_path / "plain.mat"
    scipy.io.savemat(
        raw_path,
        {
            "At": np.eye(2),
            "b": np.array([1.0, 1.0]),
            "c": np.array([1.0, 2.0]),
            "K": K_struct,
        },
    )

    gz_path = tmp_path / "compressed.mat.gz"
    with open(raw_path, "rb") as src, gzip.open(gz_path, "wb") as dst:
        dst.write(src.read())

    A, b, c, K, pars = sedumipy.read_mat(gz_path)
    x, y, info = sedumipy.sedumi(A, b, c, K)
    assert info["numerr"] == 0
    # At = I is square and invertible, so x = b = [1, 1] is the only
    # feasible point -- this just confirms the gunzipped data round-tripped
    # correctly, not any actual optimization.
    np.testing.assert_allclose(np.asarray(x).ravel(), [1.0, 1.0])


def test_write_solution_mat_roundtrip(tmp_path):
    x = np.array([1.0, 2.0, 3.0])
    y = np.array([4.0, 5.0])
    info = {"iter": 7, "numerr": 0, "pinf": 0, "dinf": 0}

    out = tmp_path / "solution.mat"
    sedumipy.write_solution_mat(out, x, y, info)

    data = scipy.io.loadmat(out)
    np.testing.assert_allclose(np.asarray(data["x"]).ravel(), x)
    np.testing.assert_allclose(np.asarray(data["y"]).ravel(), y)
    assert int(data["info"][0, 0]["iter"].item()) == 7
