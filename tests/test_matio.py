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
