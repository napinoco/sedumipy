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
  (SDP, `K.s`) problems are in scope and verified bit-for-bit against
  real Octave/SeDuMi output, including on published
  [SDPLIB](https://github.com/vsdp/SDPLIB) and
  [DIMACS](https://github.com/vsdp/DIMACS) benchmark problems. Dense-
  column preconditioning is implemented.
- `.mat` and sparse SDPA (`.dat-s`) problem/solution file I/O
  (`sedumipy.read_mat`/`write_solution_mat`, `read_sdpa`/`write_sdpa`).
- Sphinx documentation (`docs/`), published to GitHub Pages.
- CI (`.github/workflows/ci.yml`): the test suite (Octave-generated
  oracle fixtures, no Octave needed to run it) plus the fastest SDPLIB/
  DIMACS benchmark subset, both on every push/PR.
- Wheel builds via cibuildwheel (`.github/workflows/wheels.yml`) for
  Linux (manylinux) and macOS; not yet published to PyPI.

### Known limitations

- Windows is not supported: `tools/build_libsedumi.sh` (which compiles
  `libsedumi.so`) is a bash/gcc script with no `cl.exe`/MSVC equivalent
  yet.
- Complex Hermitian PSD problems (`K.scomplex`/`K.ycomplex`) are out of
  scope.
- Not yet published to PyPI, and manylinux wheels are not yet
  `auditwheel repair`'d to bundle their own BLAS (currently dynamically
  linked against the build host's `libblas`).
