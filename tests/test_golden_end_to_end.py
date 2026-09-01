"""Phase 5: end-to-end regression of sedumipy.sedumi() against the Phase 0
golden reference, on the real (not synthetic) SDPLIB-derived problems used
to validate real SeDuMi itself.

This is the companion test test_golden_reference.py's own docstring
promised ("a companion test module will load the same golden files and
compare them against sedumipy's own output"): tests/golden/*_golden.mat
only has the *outputs* Octave's real sedumi.m produced (Phase 0); the
matching *inputs* (At, b, c, K) live in the vendor/sedumi-upstream
submodule's examples/ directory (not committed to this repo, since they're
upstream SeDuMi's own example data), so this test needs the submodule
checked out (`git submodule update --init --recursive`) to run -- it
skips cleanly if that hasn't been done, same as test_golden_reference.py
does for tests/golden/ itself.

quantum.mat is excluded: its K.scomplex/K.ycomplex fields mean it's a
genuinely complex-Hermitian-PSD problem, which is out of scope for this
port (see updtransfo.py's own "real-symmetric PSD blocks only" scope
note) -- attempting it raises inside minpsdeig's eigvalsh, not a
regression in anything this test could otherwise catch.

These are real, moderately large problems (up to ~400k rows for trto3),
so this module is slower than the rest of the suite (collectively under
a minute on a typical dev machine) -- that cost buys genuine validation
that small synthetic test fixtures elsewhere in this suite cannot: it was
exactly this test (run manually against a from-scratch Octave oracle,
before this file existed) that caught a real bug in the K.s==0 branch's
one-time ADA symbolic-pattern setup, only reproducible at this kind of
scale (see sedumi.py's own comment at the `d_symbolic` override for the
full explanation).
"""

from pathlib import Path

import numpy as np
import pytest
import scipy.io
import scipy.sparse as sp

sedumipy = pytest.importorskip("sedumipy")
from sedumipy.sedumi import sedumi  # noqa: E402

GOLDEN_DIR = Path(__file__).parent / "golden"
EXAMPLES_DIR = (
    Path(__file__).parent.parent / "vendor" / "sedumi-upstream" / "examples"
)

# name -> expected optimal value (same numbers as test_golden_reference.py
# and tools/generate_golden.m). quantum.mat is deliberately omitted (see
# module docstring).
EXPECTED_PROBLEMS = {
    "arch0": -5.665170e-01,
    "control07": -2.062510e01,
    "nb": -5.070309e-02,
    "OH_2Pi_STO-6GN9r12g1T2": 7.946708e01,
    "trto3": -1.279999e04,
}

TOL = 1e-5


def _available():
    return GOLDEN_DIR.exists() and EXAMPLES_DIR.exists()


pytestmark = pytest.mark.skipif(
    not _available(),
    reason=(
        "vendor/sedumi-upstream examples or tests/golden not present -- run "
        "`git submodule update --init --recursive` (and, if tests/golden is "
        "missing, the Phase 0 golden-generation script) from the repository root"
    ),
)


def _load_K(Kmat):
    K = {}
    for fld in Kmat.dtype.names:
        val = Kmat[fld]
        if val.size == 1 and fld in ("f", "l"):
            K[fld] = int(val.item())
        else:
            K[fld] = val.ravel()
    return K


def _scalar(value):
    if sp.issparse(value):
        value = value.toarray()
    return float(np.asarray(value).squeeze())


@pytest.mark.parametrize("name", sorted(EXPECTED_PROBLEMS))
def test_sedumi_matches_golden_on_real_problems(name):
    data = scipy.io.loadmat(EXAMPLES_DIR / f"{name}.mat")
    K = _load_K(data["K"][0, 0])

    x, y, info = sedumi(data["At"], data["b"], data["c"], K)

    cx = _scalar(data["c"].T @ x)
    by = _scalar(data["b"].T @ y)
    expected = EXPECTED_PROBLEMS[name]

    assert abs(cx - expected) / abs(expected) < TOL
    assert abs(by - expected) / abs(expected) < TOL
    assert info["numerr"] in (0, 1), f"solver reported serious numerical error ({info['numerr']})"
    assert info["pinf"] == 0
    assert info["dinf"] == 0
