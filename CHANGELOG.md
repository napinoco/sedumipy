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
  these numerical comparisons, so the suite runs against four of them:
  reference Netlib and OpenBLAS on Linux, OpenBLAS on Windows, and
  Accelerate on macOS.
- Wheel builds via cibuildwheel (`.github/workflows/wheels.yml`) for
  Linux (manylinux), macOS, and Windows (MSYS2/MinGW toolchain); not yet
  published to PyPI. Each wheel carries whatever BLAS it needs, by a
  different route per platform: `auditwheel` vendors `libblas`/
  `libopenblas` into `sedumipy.libs/` on Linux, `delvewheel repair
  --analyze-existing` vendors `libopenblas.dll` and the mingw runtime
  DLLs on Windows, and macOS needs no vendoring at all because
  Accelerate is a system framework. The wheel jobs print each repaired
  wheel's bundled libraries so this stays verifiable rather than
  assumed.
- Windows support: `tools/build_libsedumi.sh` now builds
  `libsedumi.dll` via an MSYS2 MinGW64 toolchain
  (`mingw-w64-x86_64-gcc`/`-openblas`) instead of requiring MSVC --
  `libsedumi.dll` is a plain ctypes-loaded DLL, not a CPython extension,
  so it doesn't need to match whatever compiler built Python itself.
  Exercised in CI on GitHub Actions' hosted Windows runner; not yet
  hand-verified on a real Windows machine.

### Known limitations

- Complex Hermitian PSD problems (`K.scomplex`/`K.ycomplex`) are out of
  scope.
- Not yet published to PyPI.
