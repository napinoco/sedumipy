"""Phase 4: SDPA sparse (.dat-s) format I/O (sdpa.py).

read_sdpa is a port of conversion/fromsdpa.m -- test_read_sdpa_matches_
oracle checks it against a real Octave fromsdpa.m run (oracle generated
by tools/generate_sdpa_oracle.m, committed so this doesn't need Octave to
run). write_sdpa has no upstream equivalent (see sdpa.py's own docstring)
so it's checked by round-tripping through read_sdpa itself, plus (done
manually while writing this, not repeated here since it needs a live
Octave + built SeDuMi) a direct check that real Octave's own fromsdpa.m
correctly parses this port's write_sdpa output for a real problem
(vendor's arch0.mat) and reconstructs the identical (At,b,c) --
see CONTRIBUTING.md's Phase 4 note.
"""

from pathlib import Path

import numpy as np
import pytest
import scipy.io
import scipy.sparse as sp

sedumipy = pytest.importorskip("sedumipy")
from sedumipy.sdpa import read_sdpa, write_sdpa  # noqa: E402

FIXTURE_DIR = Path(__file__).parent / "fixtures" / "sdpa"
pytestmark = pytest.mark.skipif(
    not FIXTURE_DIR.exists(), reason="sdpa fixtures not generated"
)


def test_read_sdpa_matches_oracle():
    At, b, c, K = read_sdpa(FIXTURE_DIR / "test_problem.dat-s")
    oracle = scipy.io.loadmat(FIXTURE_DIR / "test_problem_oracle.mat")

    np.testing.assert_allclose(At.toarray(), sp.csc_matrix(oracle["At"]).toarray())
    np.testing.assert_allclose(b, np.asarray(oracle["b"]).ravel())
    np.testing.assert_allclose(c.toarray().ravel(), sp.csc_matrix(oracle["c"]).toarray().ravel())

    Ko = oracle["K"][0, 0]
    assert K["l"] == int(Ko["l"].item())
    np.testing.assert_allclose(np.sort(K["s"]), np.sort(Ko["s"].ravel()))


def test_read_sdpa_solves_correctly():
    At, b, c, K = read_sdpa(FIXTURE_DIR / "test_problem.dat-s")
    x, y, info = sedumipy.sedumi(At, b, c, K)
    assert info["numerr"] == 0
    assert info["pinf"] == 0
    assert info["dinf"] == 0


def test_write_sdpa_roundtrip(tmp_path):
    At, b, c, K = read_sdpa(FIXTURE_DIR / "test_problem.dat-s")

    out = tmp_path / "roundtrip.dat-s"
    write_sdpa(out, At, b, c, K)
    At2, b2, c2, K2 = read_sdpa(out)

    np.testing.assert_allclose(At.toarray(), At2.toarray())
    np.testing.assert_allclose(b, b2)
    np.testing.assert_allclose(c.toarray().ravel(), c2.toarray().ravel())
    assert K["l"] == K2["l"]
    np.testing.assert_allclose(np.sort(K["s"]), np.sort(K2["s"]))


def test_write_sdpa_rejects_socp_blocks():
    At = sp.csc_matrix(np.eye(3))
    b = np.ones(3)
    c = np.ones(3)
    K = {"l": 0, "q": np.array([3])}
    with pytest.raises(ValueError):
        write_sdpa("/dev/null", At, b, c, K)
