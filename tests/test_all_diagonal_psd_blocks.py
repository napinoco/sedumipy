"""Regression test: an SDP whose PSD blocks are *all* diagonal used to
crash in pretransfo() with an IndexError before reaching the solver.

pretransfo() detects PSD blocks whose data (A and c together) touches
only the block's diagonal entries and rewrites them as nonnegative
scalars in K.l instead -- a real cone the interior-point loop never has
to factor. That mask is `sdiag`, and a size-1 block is always in it.

The branch that builds the *remaining* (genuinely matrix-valued) PSD
blocks reads its data through `sreal` (= ~sdiag, minus complex blocks)
but used to be guarded by `if K_rsdpN:` instead. Those two disagree:
K_rsdpN is len(K.s) on the no-complex path and count_nonzero(sreal) on
the complex one, so when every block was diagonal the guard was true
while `Ks[sreal]` was empty, and indexing the empty block-marker array
raised `IndexError: index 0 is out of bounds for axis 0 with size 0`.
Guarding on the data itself (`np.any(sreal)`, matching the `np.any(
sdiag)` branch just above it) is what makes the two agree.

Nothing here needs Octave or a committed oracle: every case has a
closed-form optimum, so the assertions are exact rather than a
comparison against recorded output. That also means this file, unlike
tests/test_sedumi.py and tests/test_pretransfo.py, is not skipped when
the .mat fixtures are absent -- a crash regression should not be able to
hide behind a missing fixture.

The `K.s=[1]` case is the one most likely to be hit by accident: a
one-by-one PSD block is just a nonnegative scalar, and it crashed
unconditionally, whatever A and c contained.
"""

import numpy as np
import pytest

sedumipy = pytest.importorskip("sedumipy")
from sedumipy.sedumi import sedumi  # noqa: E402


def _solve(A, b, c, K):
    A = np.asarray(A, dtype=float)
    b = np.asarray(b, dtype=float)
    c = np.asarray(c, dtype=float)
    x, y, info = sedumi(A, b, c, K)
    return x, y, info, float(c @ x)


# min <C,X> s.t. trace(X) = 1, X psd  is the smallest eigenvalue of C,
# and for a diagonal C that is its smallest diagonal entry -- so every
# `s` block below has a hand-checkable optimum.
@pytest.mark.parametrize(
    "name,A,b,c,K,expected",
    [
        # A 2x2 block whose constraint (trace) and objective (diag(2,5))
        # both touch only the diagonal: the whole block is `sdiag`.
        ("diagonal_only_2x2", [[1, 0, 0, 1]], [1], [2, 0, 0, 5], {"s": [2]}, 2.0),
        # A 1x1 PSD block: `sdiag` by size alone, no matter the data.
        ("psd_block_of_size_one", [[1.0]], [1], [3.0], {"s": [1]}, 3.0),
        # K.l alongside an all-diagonal K.s -- the LP part does not save
        # the SDP part, since the guard looked at K.s only.
        (
            "lp_plus_diagonal_only_sdp",
            [[1, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 1]],
            [1, 1],
            [1, 1, 2, 0, 0, 5],
            {"l": 2, "s": [2]},
            3.0,
        ),
    ],
)
def test_all_psd_blocks_diagonal(name, A, b, c, K, expected):
    x, _y, info, obj = _solve(A, b, c, K)
    assert info["numerr"] == 0
    assert obj == pytest.approx(expected, abs=1e-7)
    assert np.all(x > -1e-8)


def test_one_diagonal_and_one_matrix_block():
    """The mixed case, which never crashed: one all-diagonal block plus
    one with off-diagonal data leaves `sreal` non-empty, so the old
    guard happened to agree with the data. Kept as the control that the
    fix did not change behaviour where it was already correct."""
    A = np.zeros((2, 8))
    A[0, 0] = A[0, 3] = 1.0  # trace of block 1
    A[1, 4] = A[1, 7] = 1.0  # trace of block 2
    c = np.array([2, 0, 0, 5, 1, 0.5, 0.5, 1], dtype=float)
    _x, _y, info, obj = _solve(A, [1, 1], c, {"s": [2, 2]})
    assert info["numerr"] == 0
    # min diag entry of block 1 is 2; block 2's smallest eigenvalue is
    # 1 - 0.5 = 0.5; the two blocks are independent, so 2 + 0.5.
    assert obj == pytest.approx(2.5, abs=1e-7)
