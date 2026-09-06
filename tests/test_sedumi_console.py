"""Tests for sedumi.py's console progress printout (my_fprintf.m/pars.fid
port -- see sedumi.py's `_fprintf()`). Not an oracle-comparison test (this
is purely diagnostic/display, no effect on the returned (x,y,info) --
see CONTRIBUTING.md Sec. 5): just confirms pars.fid=0 stays silent,
pars.fid=1 (default) writes to stdout, a file-like object receives the
same text, and none of this changes the solve's actual result.
"""

import io

import numpy as np
import pytest

sedumipy = pytest.importorskip("sedumipy")
from sedumipy.sedumi import sedumi  # noqa: E402

# minimize 7x1+4x2+10x3 s.t. 3x1+x2+2x3=9, x1+2x2+4x3=8, x>=0 -- README's example.
A = np.array([[3.0, 1.0, 2.0], [1.0, 2.0, 4.0]])
B = np.array([9.0, 8.0])
C = np.array([7.0, 4.0, 10.0])
K = {"l": 3}


def test_fid_zero_is_silent(capsys):
    x, y, info = sedumi(A, B, C, K, fid=0)
    captured = capsys.readouterr()
    assert captured.out == ""
    np.testing.assert_allclose(x, [2.0, 3.0, 0.0], atol=1e-6)


def test_fid_default_writes_to_stdout(capsys):
    x, y, info = sedumi(A, B, C, K)
    captured = capsys.readouterr()
    assert "SeDuMi" in captured.out or "sedumipy" in captured.out
    assert "iter seconds" in captured.out
    np.testing.assert_allclose(x, [2.0, 3.0, 0.0], atol=1e-6)


def test_fid_file_like_object():
    buf = io.StringIO()
    x, y, info = sedumi(A, B, C, K, fid=buf)
    text = buf.getvalue()
    assert "iter seconds" in text
    assert "Detailed timing" in text
    np.testing.assert_allclose(x, [2.0, 3.0, 0.0], atol=1e-6)


def test_console_output_does_not_change_result(capsys):
    x_quiet, y_quiet, info_quiet = sedumi(A, B, C, K, fid=0)
    x_loud, y_loud, info_loud = sedumi(A, B, C, K, fid=1)
    capsys.readouterr()
    np.testing.assert_allclose(x_quiet, x_loud)
    np.testing.assert_allclose(y_quiet, y_loud)
    assert info_quiet == info_loud
