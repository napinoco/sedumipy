# Changelog

All notable changes to this project are documented here. Format loosely
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/); this
project has not yet made a numbered release, so entries so far are
grouped by phase (see [`CONTRIBUTING.md`](CONTRIBUTING.md) §2 for the
full phase-by-phase status and history).

## [Unreleased]

### Added

- Ported `sedumi()`, SeDuMi's top-level LP/SOCP/SDP interior-point
  solver driver, to pure Python (NumPy/SciPy) plus a standalone C kernel
  library (`libsedumi.so`/`.dylib`, no MATLAB/Octave/MEX dependency),
  loaded via `ctypes`. LP, second-order-cone (SOCP), and semidefinite
  (SDP, `K.s`) problems are in scope and verified against
  real Octave/SeDuMi output (to tight numerical tolerances -- see
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
  only now); not yet published to PyPI. Linux and Windows now link
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
- Not yet published to PyPI.
