"""Phase 5: timing benchmark for sedumipy.sedumi() over the same
SDPLIB-derived real problems tests/test_golden_end_to_end.py validates
against (vendor/sedumi-upstream/examples/*.mat), so results are for
problems already confirmed to solve correctly, not just "fast but wrong".

quantum.mat is skipped (complex-Hermitian PSD, out of scope -- see
tests/test_golden_end_to_end.py's docstring).

This only times this port's own Python solve. To compare against the real
Octave/MEX sedumi.m, build it first (`git submodule update --init
--recursive`, then from vendor/sedumi-upstream:
`octave-cli --no-gui --eval "install_sedumi"`, which needs Octave +
mkoctfile + a Fortran BLAS dev package installed) and time
`[x,y,info] = sedumi(data.At, data.b, data.c, data.K, struct('fid', 0))`
the same way -- see CONTRIBUTING.md section 2's Phase 5 note for a
recorded comparison.

Run from the repository root:
    python3 tools/benchmark_examples.py
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import scipy.io

from sedumipy.sedumi import sedumi

EXAMPLES_DIR = Path(__file__).parent.parent / "vendor" / "sedumi-upstream" / "examples"

PROBLEMS = ["nb", "arch0", "control07", "trto3", "OH_2Pi_STO-6GN9r12g1T2"]


def _load_K(Kmat):
    K = {}
    for fld in Kmat.dtype.names:
        val = Kmat[fld]
        K[fld] = int(val.item()) if (val.size == 1 and fld in ("f", "l")) else val.ravel()
    return K


def main():
    if not EXAMPLES_DIR.exists():
        raise SystemExit(
            "vendor/sedumi-upstream/examples not found -- run "
            "`git submodule update --init --recursive` first"
        )

    print(f"{'problem':<25} {'m':>8} {'N':>8} {'seconds':>10} {'iter':>6} {'numerr':>7}")
    for name in PROBLEMS:
        data = scipy.io.loadmat(EXAMPLES_DIR / f"{name}.mat")
        K = _load_K(data["K"][0, 0])

        t0 = time.perf_counter()
        x, y, info = sedumi(data["At"], data["b"], data["c"], K)
        elapsed = time.perf_counter() - t0

        m = data["At"].shape[1]
        N = data["c"].shape[0]
        print(f"{name:<25} {m:>8} {N:>8} {elapsed:>10.2f} {info['iter']:>6} {info['numerr']:>7}")


if __name__ == "__main__":
    main()
