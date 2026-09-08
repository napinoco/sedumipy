# Changelog

All notable changes to this project are documented here. Format loosely
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). See
[`CONTRIBUTING.md`](CONTRIBUTING.md) §2 for the full phase-by-phase
status and history behind these entries, and
[`RELEASING.md`](RELEASING.md) for how a release is cut.

## [Unreleased]

## [0.0.5] - 2026-09-08

### Added

- Wheels now also cover Linux `aarch64` and Intel macOS, via GitHub's
  native `ubuntu-24.04-arm` and `macos-15-intel` hosted runners
  (`.github/workflows/wheels.yml`) -- no cross-compilation involved,
  cibuildwheel just autodetects the runner's own architecture. Linux
  aarch64 links scipy-openblas64 the same way x86_64 already did (it
  ships wheels for that platform too); Intel macOS needs no changes at
  all, since it already linked the system Accelerate framework the
  same way Apple silicon does. `macos-13`, the previous Intel label,
  was not an option: GitHub retired it in December 2025.
- Wheels now also cover Windows ARM64 (CPython 3.11-3.13; 3.10 stays
  skipped -- numpy/scipy publish no win_arm64 wheel for it, confirmed
  by trying in CI), via a new `build-windows-arm64` job on the
  `windows-11-arm` hosted runner (GA for public repos since August
  2025). MSYS2 ships no ARM64 build of its own gcc, only the
  CLANGARM64 environment's clang-based cross toolchain, so this job
  installs `mingw-w64-clang-aarch64-gcc-compat` instead of
  `mingw-w64-x86_64-gcc` -- a compatibility layer providing a `gcc.exe`
  wrapper around clang that accepts the same GNU-style flags
  `tools/build_libsedumi.sh` already passed for x64. That script now
  self-discovers which MSYS2 subdirectory (`mingw64` vs `clangarm64`)
  actually exists rather than switching on `uname -m`: a first attempt
  at the latter failed in real CI, since the bash running the script
  is Git for Windows' own bundled bash (found via `shutil.which`, not
  the MSYS2 install), which runs under emulation on this runner and
  reports its own x86_64 regardless of the ARM64 host underneath it.
  `tools/repair_windows_wheel.py` picks whichever toolchain directory
  actually exists the same way. scipy-openblas64 ships a win_arm64
  wheel same as win_amd64, so no BLAS-side changes were needed.
- Wheels now also cover Linux musllinux (Alpine), on both x86_64 and
  aarch64 -- removed from pyproject.toml's skip list once verified that
  scipy-openblas64, numpy and scipy all ship musllinux_1_2 wheels for
  every CPython version this project builds (no cp310-style gap like
  win_arm64's). No wheels.yml changes needed: cibuildwheel already
  builds both the manylinux and musllinux variant of each wheel
  identifier from the same Linux runner by default, and cibuildwheel
  2.21 already defaults musllinux to musllinux_1_2 -- current enough
  that, unlike manylinux2014, no musllinux-\\*-image override was
  needed either.
- Still no wheels for 32-bit Windows: win32 has no scipy-openblas64
  build to link.

## [0.0.4] - 2026-09-06

Re-release of 0.0.3, whose PyPI upload failed because
`pyproject.toml`'s `version` was not bumped for it. No code changes.

## [0.0.3] - 2026-09-06

### Added

- `sedumi()` now prints original SeDuMi's console progress output:
  `sedumi.py` gained an `_fprintf()` helper wired in at every point
  `sedumi.m` calls `my_fprintf(pars.fid, ...)` -- the welcome banner,
  alg/theta/beta, the preprocessing summary, pre-loop eqs/nnz stats,
  the per-iteration progress table, "Run into numerical problems"/
  "Maximum number of iterations reached", the final result summary
  (both the feasible and infeasible paths), detailed timing, and
  max-norms/Cholesky stats. `pars["fid"]` keeps its upstream default of
  `1` (stdout); `fid=0` stays silent, and any file-like object with
  `.write()` is also accepted. Purely diagnostic -- confirmed not to
  affect the returned `(x, y, info)` (`tests/test_sedumi_console.py`).
  The welcome banner no longer echoes the original authors' credit line
  from `sedumi.m` (see README's "A note on citation and attribution");
  it now points at `fid=0` for quiet mode instead.

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
  loaded via `ctypes`. LP (`K.l`), second-order-cone (SOCP, `K.q`/`K.r`),
  and semidefinite (SDP, `K.s`) problems are in scope and verified against
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
