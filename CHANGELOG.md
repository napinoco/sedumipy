# Changelog

All notable changes to this project are documented here. Format loosely
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). See
[`CONTRIBUTING.md`](CONTRIBUTING.md) §2 for the full phase-by-phase
status and history behind these entries, and
[`RELEASING.md`](RELEASING.md) for how a release is cut.

## [Unreleased]

## [0.0.2] - 2026-09-05

### Added

- `sedumipy.cvxpy_interface`: sedumipy as a
  [cvxpy](https://www.cvxpy.org/) solver, via cvxpy's custom-solver hook
  (`problem.solve(solver=SEDUMIPY())`) -- so it works against a stock
  cvxpy install, with nothing to merge into cvxpy itself. Covers the
  same LP/SOCP/SDP scope as the solver (cvxpy refuses exponential/power
  cone and mixed-integer problems for it rather than mis-solving them),
  returns primal and dual values and cvxpy's own solver statuses, and
  passes `solve()`'s extra keyword arguments through as `pars`. Install
  with `pip install sedumipy[cvxpy]`; see `docs/usage.rst`.

### Fixed

- `sedumi()` raised `IndexError: index 0 is out of bounds for axis 0
  with size 0` from `pretransfo()` on any problem whose PSD blocks were
  *all* diagonal -- either size 1 (`K.s=[1]`, which is just a
  nonnegative scalar and crashed whatever the data) or with `A` and `c`
  touching only the block's diagonal entries. `pretransfo()` rewrites
  such blocks into `K.l`, and the branch building the remaining
  matrix-valued blocks read its data through `sreal` while being
  guarded by `K_rsdpN` -- which counts the diagonal blocks too on the
  no-complex path, so the guard passed with no data behind it. It now
  guards on the data (`np.any(sreal)`), matching the `np.any(sdiag)`
  branch beside it. A mix of diagonal and matrix-valued blocks was
  never affected.

### Changed

- CI (`.github/workflows/ci.yml`) now runs the full test suite on
  Windows against both BLAS choices `tools/build_libsedumi.sh` can make
  there, not just the MSYS2 one: scipy-openblas64, which is what the
  published wheels actually link and is ILP64 (a different `blasint` --
  see `csrc/sedumi_platform.h`), and MSYS2's OpenBLAS, the documented
  fallback. Until now the only Windows exercise of the ILP64 build was
  `wheels.yml`'s one-line `test-command`, which solves a 2x2 LP and
  never reaches the sparse Cholesky, dense-column or PSD cone paths
  where an integer-width mistake would show. Linux already covered
  ILP64 this way.

## [0.0.1] - 2026-09-05

First release published to [PyPI](https://pypi.org/project/sedumipy/):
`pip install sedumipy`. Everything below is the work that got it there.

### Added

- Ported `sedumi()`, SeDuMi's top-level LP/SOCP/SDP interior-point
  solver driver, to pure Python (NumPy/SciPy) plus a standalone C kernel
  library (`libsedumi.so`/`.dylib`, no MATLAB/Octave/MEX dependency),
  loaded via `ctypes`. LP, second-order-cone (SOCP), and semidefinite
  (SDP, `K.s`) problems are in scope and verified against
  original Octave/SeDuMi output (to tight numerical tolerances -- see
  `docs/status.rst`), including on published
  [SDPLIB](https://github.com/vsdp/SDPLIB) and
  [DIMACS](https://github.com/vsdp/DIMACS) benchmark problems. Dense-
  column preconditioning is implemented.
- `.mat` and sparse SDPA (`.dat-s`) problem/solution file I/O
  (`sedumipy.read_mat`/`write_solution_mat`, `read_sdpa`/`write_sdpa`).
- Sphinx documentation (`docs/`), published to GitHub Pages.
- CI (`.github/workflows/ci.yml`): the test suite (Octave-generated
  oracle fixtures, no Octave needed to run it) plus the fastest SDPLIB/
  DIMACS benchmark subset, both on every push/PR, on Linux, macOS and
  Windows. The BLAS implementation is the axis that actually matters for
  these numerical comparisons, so the Linux job runs it against three:
  reference Netlib and OpenBLAS (both a system package), and
  scipy-openblas64 (pip-installed, ILP64 -- see below); Windows runs
  OpenBLAS via MSYS2 and macOS runs Accelerate.
- `tools/build_libsedumi.sh` links [scipy-openblas64](
  https://pypi.org/project/scipy-openblas64/) -- a pip-installable,
  prebuilt ILP64 OpenBLAS, pinned to an exact version, the same package
  numpy/scipy themselves build against -- as `libsedumi`'s BLAS on
  Linux and Windows, whenever it's importable by the Python building it,
  falling back to that OS's own system-BLAS story (`-lblas`/
  `-lopenblas` on Linux, MSYS2's `-lopenblas` on Windows) only when it
  isn't. It's a build-time-only dependency -- never installed at runtime
  for end users -- since a wheel build vendors the resulting shared
  library into the wheel itself, same as it already did for
  `-lopenblas`. macOS deliberately keeps linking the system Accelerate
  framework unconditionally instead: it never had Linux/Windows's
  build/install pain (Accelerate needs no install step at all, unlike
  `libblas-dev`/MSYS2), so there's nothing here for scipy-openblas64 to
  save, only wheel bytes and a possible speed regression on Apple
  Silicon to risk. See `csrc/sedumi_platform.h`'s `SEDUMI_BLAS_ILP64`/
  `BLAS_SYMBOL_PREFIX`/`BLAS_SYMBOL_SUFFIX` for how the ILP64 integer
  width and the symbol-mangling scipy-openblas64 needs (to coexist with
  other OpenBLAS copies in one process) are handled.
- Wheel builds via cibuildwheel (`.github/workflows/wheels.yml`) for
  Linux (manylinux), macOS, and Windows (MSYS2/MinGW toolchain, gcc
  only now), uploaded to PyPI on release via Trusted Publishing (see
  [`RELEASING.md`](RELEASING.md)). Linux and Windows now link
  scipy-openblas64 (see above) as their BLAS, vendored into the wheel by
  `auditwheel` on Linux and `delvewheel repair` (via
  `tools/repair_windows_wheel.py`) on Windows; macOS is unchanged --
  still the system Accelerate framework, nothing to vendor. The wheel
  jobs print each repaired wheel's bundled libraries so this stays
  verifiable rather than assumed.
- Windows support: `tools/build_libsedumi.sh` builds `libsedumi.dll` via
  an MSYS2 MinGW64 toolchain (`mingw-w64-x86_64-gcc`) instead of
  requiring MSVC -- `libsedumi.dll` is a plain ctypes-loaded DLL, not a
  CPython extension, so it doesn't need to match whatever compiler built
  Python itself. Exercised in CI on GitHub Actions' hosted Windows
  runner; not yet hand-verified on a real Windows machine.

### Known limitations

- Complex Hermitian PSD problems (`K.scomplex`/`K.ycomplex`) are out of
  scope.
- No wheels for Linux `aarch64`, Intel macOS, Alpine/musl, or 32-bit/ARM
  Windows: those platforms install from the source distribution, which
  needs a C compiler and a BLAS (see `docs/installation.rst`).
