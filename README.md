# sedumipy

[![ci](https://github.com/napinoco/sedumipy/actions/workflows/ci.yml/badge.svg)](https://github.com/napinoco/sedumipy/actions/workflows/ci.yml)
[![wheels](https://github.com/napinoco/sedumipy/actions/workflows/wheels.yml/badge.svg)](https://github.com/napinoco/sedumipy/actions/workflows/wheels.yml)
[![docs](https://github.com/napinoco/sedumipy/actions/workflows/docs.yml/badge.svg)](https://napinoco.github.io/sedumipy/)
[![License: GPL v2](https://img.shields.io/badge/License-GPLv2-blue.svg)](LICENSE)

**A MATLAB/Octave-free Python port of [SeDuMi](https://github.com/sqlp/sedumi)** —
an interior-point solver for linear (LP), second-order cone (SOCP), and
semidefinite (SDP) programs — as a standalone C kernel library plus a
Python (NumPy/SciPy) package. No MATLAB or GNU Octave runtime is needed
to install or run it.

```python
import numpy as np
import sedumipy

# minimize x1 + x2  s.t.  x1 = x2 = 1,  x >= 0
x, y, info = sedumipy.sedumi(np.eye(2), np.array([1.0, 1.0]), np.array([1.0, 1.0]), {"l": 2})
```

This mirrors real SeDuMi's own `[x, y, info] = sedumi(A, b, c, K)` call.

**New contributor?** Read [`CONTRIBUTING.md`](CONTRIBUTING.md) first — it
has the current phase-by-phase status, the porting workflow this project
follows, known scope limitations, and the prioritized list of remaining
work. (It's in Japanese; ask if you'd like an English translation.)

## Status

LP, second-order-cone (SOCP), and semidefinite (SDP, `K.s`) problems are
all fully ported and verified against real Octave/SeDuMi output, to tight
numerical tolerances, including on published
[SDPLIB](https://github.com/vsdp/SDPLIB)
and [DIMACS](https://github.com/vsdp/DIMACS) benchmark problems (see
[Benchmarks](#benchmarks) below). Dense-column preconditioning is also
implemented. `pip`-installable wheels build for Linux (manylinux), macOS,
and Windows via [`wheels.yml`](.github/workflows/wheels.yml)'s
cibuildwheel job (not yet published to PyPI; the Windows build uses an
MSYS2/MinGW toolchain and has only been exercised on GitHub Actions'
hosted runner, not hand-verified on a real Windows machine) — see
[`CONTRIBUTING.md`](CONTRIBUTING.md) for the full phase-by-phase status
and known limitations.

## Documentation

Full documentation (installation, the problem/solver API, and the
internals reference) is published at
**[napinoco.github.io/sedumipy](https://napinoco.github.io/sedumipy/)**
(built with Sphinx from [`docs/`](docs/), auto-published to GitHub Pages
by [`.github/workflows/docs.yml`](.github/workflows/docs.yml) on every
push to `main`). To build it locally instead:

```sh
pip install -r docs/requirements.txt
sphinx-build -b html docs docs/_build/html
```

## Repository layout

The original MATLAB/Octave/MEX implementation this project ports from is
kept as a reference-only git submodule, pinned to the commit the port
started from:

```
sedumipy/
  vendor/sedumi-upstream/   # submodule: sqlp/sedumi (reference only, not built by default)
  examples/sdplib/          # submodule: vsdp/SDPLIB (published benchmark problems + optimal values)
  examples/dimacs/          # submodule: vsdp/DIMACS (published benchmark problems + optimal values)
  csrc/                     # forked, MEX-free standalone C kernels (source for libsedumi.so)
  src/sedumipy/             # the Python package
  tests/                    # test suite + committed Octave-generated oracle fixtures
  tools/                    # libsedumi build script + oracle/golden-reference generators
  docs/                     # Sphinx documentation source
```

## Getting started

```sh
git clone --recurse-submodules <this-repo-url>
cd sedumipy
python -m venv .venv
.venv/bin/pip install -e .[test]
.venv/bin/python -m pytest tests/ -q
```

If you already cloned without `--recurse-submodules`, run
`git submodule update --init --recursive` first. Building `libsedumi.so`
(the compiled C kernel library) requires a C compiler and a BLAS
development package (e.g. `apt install build-essential libblas-dev` on
Debian/Ubuntu, or nothing extra at all on macOS — it links against the
system Accelerate framework there); it's then built automatically the
first time `sedumipy` is imported, via `tools/build_libsedumi.sh`. On
Windows, install [MSYS2](https://www.msys2.org/) and its
`mingw-w64-x86_64-gcc`/`mingw-w64-x86_64-openblas` packages first — see
the Windows note in [`docs/installation.rst`](docs/installation.rst) or
[`CONTRIBUTING.md`](CONTRIBUTING.md). The Octave submodule is only
needed to regenerate oracle/golden-reference data, not to run the
existing test suite.

## Benchmarks

`tests/test_benchmarks.py` solves published [SDPLIB](https://github.com/vsdp/SDPLIB)
and [DIMACS](https://github.com/vsdp/DIMACS) problems (added as git
submodules under `examples/`) and checks the result against each
collection's own official optimal-value table -- a correctness/regression
check against real reference numbers, not a synthetic self-check. It also
prints a timing/iteration-count summary and writes it to
`benchmark_results.csv`.

Coverage is essentially the full published set: all 92 SDPLIB problems
and every DIMACS problem sedumi() can read directly, minus a documented
handful this port can't yet solve in bounded time/memory or that hit a
real solver limitation -- see the module's own docstring for the full
exclusion list and why each one is excluded.

```sh
git submodule update --init --recursive   # if not already done
.venv/bin/python -m pytest tests/test_benchmarks.py -v          # everything, ~101 problems (~10 min)
.venv/bin/python -m pytest tests/test_benchmarks.py -v -m mini  # fastest subset only (~35s)
```

Problems are grouped by `pytest.mark.mini` (<2s each) / `timing` (2-20s)
/ `extended` (20s-130s) by measured solve time -- see the module's own
docstring for the sign conventions each collection/family needs.

## License

SeDuMi is licensed under the GNU General Public License v2 (see
[`LICENSE`](LICENSE)); this port, being a derivative work, is licensed
the same way.
