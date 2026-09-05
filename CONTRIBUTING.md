# sedumipy Development Guide (Working Notes)

This document is a handoff resource for anyone new joining this project
to port SeDuMi off MATLAB/Octave — human or AI agent alike — so they can
keep working without getting lost. It summarizes what's done, what's
left, and the procedures/policies followed so far. For the detailed,
session-by-session bug-hunting narrative and the full history of how
each piece of work got done, see [`DEVLOG.md`](DEVLOG.md) — this file
stays focused on the current state and how to work on the project.

**On the repository layout:** this port originally lived under
`python_port/` inside `napinoco/sedumi` (a fork of `sqlp/sedumi`),
alongside the original MATLAB/C implementation. This repository
(`sedumipy`) is a clean spin-off containing only the code the port
needs. The original `.m`/MEX `.c` sources are kept for reference under
`vendor/sedumi-upstream/` (a submodule of `sqlp/sedumi`, pinned to the
commit the fork started from); oracle regeneration (producing comparison
data on real Octave) uses this submodule.

**Development environment prerequisites:** `pip install -e .[test]`
builds `libsedumi.so` on first import (`_native.py`'s `_ensure_built()`),
so a C compiler (`gcc`) must already be installed (`apt install
build-essential` on Debian/Ubuntu). For BLAS: **macOS always uses the
system Accelerate framework** (nothing to install — `tools/
build_libsedumi.sh`'s Darwin branch does this unconditionally, since
macOS never had Linux/Windows's "install/build your own BLAS" pain to
begin with, so scipy-openblas64 buys nothing there, and there's no
reason to give up Accelerate's Apple Silicon tuning just to match the
other platforms). On Linux/Windows, running `pip install
scipy-openblas64==0.3.34.106.0` before `pip install -e .[test]
--no-build-isolation` uses the same pip-installable, prebuilt ILP64
OpenBLAS (`scipy-openblas64`) that the published wheels themselves link
(`tools/build_libsedumi.sh` auto-detects whether it's importable; see
that file and `csrc/sedumi_platform.h`'s `SEDUMI_BLAS_ILP64` for
details). The version is pinned deliberately (so a breaking change
upstream doesn't suddenly break the build) — if you ever bump it, update
the same string in `.github/workflows/ci.yml`/`README.md`/
`docs/installation.rst` per the comment in `pyproject.toml`, and re-run
the tests and benchmarks. Without that package, it falls back to each
OS's system-build BLAS (on Linux, the `libblas`/`libopenblas` dev
headers, e.g. `apt install libopenblas-dev`; LAPACK is actually unused —
linking is BLAS-only). A bare container/CI environment may have neither
installed, causing the build to fail — if you hit a wall in a new
environment, suspect this first. (`.github/workflows/ci.yml` exercises
both paths on every run.)

**Windows development environment:** `tools/build_libsedumi.sh` is a
bash script and can't be run directly on Windows (`setup.py`/
`_native.py` explicitly invoke it via `bash tools/build_libsedumi.sh
...` only on Windows (`sys.platform == "win32"`)). Install
[MSYS2](https://www.msys2.org/), then from the MINGW64 shell run

```
pacman -S mingw-w64-x86_64-gcc
```

to put `bash`/`gcc` on Windows's `PATH` (`C:\msys64\mingw64\bin` and
`C:\msys64\usr\bin`); after that, building and testing follow the same
steps as on other OSes (BLAS uses the `scipy-openblas64` path above; if
you skip that, you'll also need to install
`mingw-w64-x86_64-openblas`). `libsedumi.dll` is dynamically linked
against the BLAS DLL and mingw runtime DLLs (it's a plain ctypes-loaded
DLL, not a CPython extension, so it was never tied to a particular
compiler), so producing a distributable wheel requires bundling those
with `delvewheel repair` (see `.github/workflows/wheels.yml`/
`tools/repair_windows_wheel.py`/`pyproject.toml`'s
`[tool.cibuildwheel.windows]`). This session had no real Windows machine
to verify behavior on (development on this repository happens on
Linux), so CI's results are the practical verification.

`tools/build_libsedumi.sh`'s Windows branch calls `gcc` neither via
`PATH` nor via MSYS2's `/mingw64` pseudo-path, but by a **Windows-style
absolute path**, `C:\msys64\mingw64\bin\gcc.exe` (overridable via the
`MSYS2_ROOT` environment variable). Two reasons, both confirmed on real
CI: (1) GitHub Actions' Windows runner ships with an unrelated,
pre-installed MinGW toolchain at `C:\mingw64`, and a bare `gcc` call
picks that one up instead (link fails for lack of `-lopenblas`). (2)
MSYS2's own POSIX-style `/mingw64` mount also doesn't resolve for a
non-interactive/non-login `bash script.sh` invocation (confirmed on CI:
even with `mingw-w64-x86_64-gcc` freshly installed, this produced
`/mingw64/bin/gcc.exe: No such file or directory`).

## 1. Project goal

Port SeDuMi (an interior-point solver for SDP/SOCP problems that runs on
MATLAB/Octave — roughly 96 `.m` files plus about 54 MEX `.c` files) to a
form with **zero dependency on either MATLAB or Octave**. Architecture
adopted (decided in an earlier session, "Option B"):

- The algorithm itself (interior-point iteration logic) → **Python
  (NumPy/SciPy)**
- Performance-critical low-level kernels (Cholesky factorization, cone
  arithmetic, etc.) → build the existing `.c` files with their MEX
  dependency stripped out as a **standalone C library**
  (`libsedumi.so`), called from Python via **ctypes**

The end goal is a `pip install`-able Python package with no
MATLAB/Octave dependency.

## 2. Overall phases and current progress

Per-phase tasks are tracked in the task-management tool (TaskList). As
of 2026-08-31:

| Phase | What | Status |
|---|---|---|
| Phase 0 | Verification harness (golden reference on real Octave) | **Done** |
| Phase 1 | Strip MEX from C kernels → standalone C library | **Done** |
| Phase 2 | Python bindings (ctypes) (clusters 1-5) | **Done** |
| Phase 3-a | Turn thin MEX-wrapper `.m` files into the public Python API | **Done** |
| Phase 3-b | Port cone-math utilities (eigK, psdeig, psdscale, etc.) | **Done** |
| Phase 3-c | Port interior-point iteration control logic (sdinit..optstep) | **Done** |
| Phase 3-d | Port `sedumi.m` itself + full verification against golden reference | **Done (LP+SOCP+PSD scope)** |
| Phase 4 | High-level API and I/O compatibility layer (.mat/SDPA) | **Done** |
| Phase 5 | Verification and benchmarking | **Done** |
| Phase 6 | Packaging and release | **Done** (Linux/macOS/Windows cibuildwheel builds in `.github/workflows/wheels.yml`, BLAS bundled into the wheels by `auditwheel`/`delvewheel`, source distribution completed by `MANIFEST.in` and install-tested in CI, published to PyPI on release via Trusted Publishing -- see [`RELEASING.md`](RELEASING.md)) |

Phase 3 (porting the interior-point algorithm itself) is **complete for
LP + SOCP (second-order cone) + PSD (positive semidefinite cone)
problems**: calling `sedumipy.sedumi.sedumi(A,b,c,K)` has been confirmed,
via real-oracle comparison, to return solutions that exactly match
Octave SeDuMi (`tests/test_sedumi.py`. The PSD-cone main-loop wiring is
in `getada_psd.py` (`build_aord`/`getada_psd`, using `incorder.py`/
`getsymbada.py`/`_native.getada1`/`getada2`/`getada3`), verified by the
same file's `sdp_feasible`/`sdp_mixed_cones_feasible` cases).

**Dense-column optimization is also ported** (`getdense.py`/
`symbcholden.py`/`deninfac.py`/`pcg.py`'s product-form preconditioning
correction, orchestrating `_native.symbfwblk`/`adendotd`/`adenscale`/
`dpr1fact`). Even on problems that actually contain dense columns
(`tests/test_sedumi.py::test_sedumi_dense_matches_octave`, deliberately
lowering `getdense.m`'s detection threshold via `pars.denf=3`), `iter`/
`numerr`/`pinf`/`dinf`/`x` all match the Octave version (see
[`DEVLOG.md`](DEVLOG.md) on "dual solution non-uniqueness" for `y`).

**Phase 5 (verification on real problems) is also complete**
(`tests/test_golden_end_to_end.py`): running `sedumipy.sedumi()` on the
real problems the Phase 0 golden reference targeted (from SDPLIB, under
`vendor/sedumi-upstream/examples/`: `nb`/`arch0`/`control07`/`trto3`/
`OH_2Pi_STO-6GN9r12g1T2` — `quantum` is excluded as out of scope, since
it's a complex-Hermitian PSD problem using `K.scomplex`/`K.ycomplex`)
confirmed the objective value matches the Octave golden reference. This
verification effort found and fixed two real bugs along the way — see
[`DEVLOG.md`](DEVLOG.md) for the details:
- The `K.s==0` (LP+SOCP only) path's ADA symbolic-Cholesky ordering was
  missing part of the sparsity pattern that depends on the Lorentz
  cone's arrow term (`d.q2`) — as `d.q2` grows, the Cholesky
  factorization becomes inaccurate and PCG diverges (surfaced on
  `nb.mat`, which has 396 SOCP blocks).
- `cpspdiag` (the diagonal-extraction helper the `K.s==0` branch of
  `getada3` calls for diagnostics) used `bsearch()` via the `ibsearch`
  macro, whose comparator hit the same qsort/bsearch
  undefined-behavior class of bug as `sortnnz.c`/`iswnbr.c`. However,
  since the actual `sedumi.py` call path never invokes `getada3` when
  `K.s==0`, the real-world impact was limited to a test
  (`test_getada_no_psd_blocks`).

## 3. Directory layout

```
sedumipy/                    # repository root
  vendor/
    sedumi-upstream/          # submodule: sqlp/sedumi (pinned to the fork's starting commit, reference only)
  csrc/                       # standalone C kernel sources with mex.h stripped out (source for libsedumi.so)
    *.c / *.h                  # forked to build MEX-free via sedumi_platform.h
    sedumi_platform.c / .h
    kernel_smoke/smoke_test.c  # Phase 1 smoke test
  src/
    sedumipy/                 # the ported Python package itself
      _native.py               # all ctypes bindings collected here (every C-kernel call goes through this)
      libsedumi.so              # built shared library (produced by tools/build_libsedumi.sh, gitignored)
      cone.py                   # cone-math utilities (eigK, psdeig, psdscale, frameit, ...)
      pretransfo.py / posttransfo.py   # external <-> internal format conversion
      sdinit.py                 # initial-point generation
      sdfactor.py / sddir.py    # self-dual embedding factorization and direction computation
      pcg.py                    # preconditioned conjugate gradient (loopPcg/wrapPcg)
      wregion.py                # one iteration's predictor(+corrector) step
      updtransfo.py             # scaling-point update
      maxstep.py / widelen.py / stepdif.py / trydif.py  # step-length computation
      getada.py / getdatm.py / deninfac.py  # ADA matrix construction/factorization (K.s==0)
      incorder.py / getsymbada.py / getada_psd.py  # ADA matrix construction (K.s!=0)
      symbchol.py                # symbolic Cholesky of ADA (run once)
      optstep.py                 # LP early-optimality check (optstep.m)
      amul.py / checkpars.py     # auxiliary utilities
      sedumi.py                  # top-level driver (wires everything together)
      matio.py                   # Phase 4: read/write .mat problem/solution files
      sdpa.py                    # Phase 4: read/write SDPA sparse (.dat-s) format
  tests/
    test_*.py                  # per-module verification tests (oracle comparison)
    fixtures/                  # oracle data generated on real Octave (committed)
    golden/                    # Phase 0 golden reference
  tools/
    generate_*_oracle.m        # per-test scripts that generate oracles from vendor/sedumi-upstream's Octave/MEX build
    build_libsedumi.sh         # builds libsedumi.so from csrc/
  pyproject.toml
  setup.py                      # Phase 6: build_ext hook that compiles libsedumi.so at wheel-build time
  README.md                    # (somewhat stale — this CONTRIBUTING.md has the more current phase overview)
```

## 4. Development workflow (steps for porting a single `.m` file)

This is the procedure followed consistently so far. Follow it when
porting a new function.

1. **Read the target `.m` file's entire real source.** Understand the
   implementation logic line by line, not just the comments. If it's a
   thin MEX wrapper (a stub that just calls
   `sedumi_binary_error()`), read the corresponding `.c` file instead.
2. **Write a faithful Python port.** Match variable names and
   processing order to the `.m`/`.c` source as closely as possible, so
   the two are easy to compare later. In the docstring, spell out what
   was ported, which lines of the `.m` file it corresponds to, and any
   parts deliberately omitted and why.
3. **Generate an oracle on real Octave.** Write
   `tools/generate_<name>_oracle.m`, call the actual `.m`/MEX build
   under `vendor/sedumi-upstream/` (with `install_sedumi` already run)
   to save the input data and output to a `.mat` file.
   ```
   octave-cli --no-gui --eval "cd tools; generate_<name>_oracle"
   ```
   Reuse an existing fixture (e.g. `pretransfo`'s `K2`/`prep`) whenever
   possible, rather than paying the cost of recomputing from scratch.
4. **Write `test_<name>.py` and compare.** Load the `.mat` file and
   check the Python output against it with
   `np.testing.assert_allclose`.
5. **When something doesn't match, debug by comparing intermediate
   values one at a time.** Don't "fix it by feel" — print the same
   variable from both the Octave side and the Python side and pin down
   exactly where they first diverge. Every real bug found so far (see
   [`DEVLOG.md`](DEVLOG.md)) was found this way.
6. **Run all of `python_port/tests/ -q` to check for regressions before
   committing.** Commit messages should state what was ported, any bugs
   found and fixed, and any deliberate scope limitations and why.

### A note on writing oracle scripts (random streams)

Octave's `rand('seed', N)` only resets the state of the uniform random
generator. `randn` is a separate stream, so any case using `randn` must
also call `randn('seed', N)` (forgetting this means a different dataset
is generated on every run, making the fixture unreproducible — this
actually happened and cost time in this session).

Also, inserting new code **in the middle** of an existing oracle-
generation script shifts every `rand()` call after it, silently
rewriting unrelated existing fixtures too. Always append new cases at
the **end** of the script, and check `git diff --stat` to confirm no
existing fixture file's size changed.

## 5. Scope limitations (deliberately unimplemented parts)

The following are not "gaps we didn't notice" but constraints
**deliberately rejected with `NotImplementedError`**, each documented
with its reasoning in the relevant module's docstring.

- **Rotated second-order cones (`K.r`)** are converted to standard
  second-order cones (`K.q`) at the `pretransfo.py` stage, so no code
  downstream needs to be aware of them (`sedumi.py`'s verification
  tests include one rotated-cone case).
- **Console output, the v-plot, `pars.stopat`'s debug break, the
  pre-solve rank diagnostic, and the DIMACS error metrics
  (`info.err`)** are not ported. All of these are purely
  diagnostic/display and have no effect on the returned `(x, y, info)`
  values, so they were deprioritized.

## 6. Running the tests (and regenerating oracle fixtures)

```
cd sedumipy
git submodule update --init --recursive   # first time only: fetches vendor/sedumi-upstream
.venv/bin/pip install -e .[test]           # libsedumi.so is built automatically on first import
.venv/bin/python -m pytest tests/ -q
```

The test suite passes even without Octave installed (oracle data is
pre-generated and committed as `.mat` files, and the tests never call
Octave at runtime). Octave is only needed to **regenerate** an oracle
(with `install_sedumi` already run under `vendor/sedumi-upstream`):

```
octave-cli --no-gui --eval "cd tools; generate_<name>_oracle"
```

## 7. Coding conventions and naming

- **`K` (the cone-structure struct) is represented as a Python `dict`.**
  Field names match the `.m` version's `K.f`/`K.l`/`K.q`/... as closely
  as possible (`K["l"]`, `K["mainblks"]`, etc.). Internal-only fields
  computed by `pretransfo.py` (`mainblks`, `qblkstart`, `sblkstart`,
  `lq`, `N`, `rsdpN`, etc.) are threaded through unchanged to every
  downstream function.
- **Indexing is 0-indexed (Python-style) by default.** However, the
  contents of `K` (e.g. `mainblks`) often still carry the original
  `.m` file's 1-indexed values as-is, and callers convert with something
  like `int(x) - 1` at each use site (there's no unified conversion
  layer yet). When writing new code, mimic the indexing-conversion
  pattern of the nearest existing similar function.
- **`scipy.sparse.csc_matrix` is the default sparse-matrix format.**
- **Every C-kernel call is collected in `_native.py`**; other modules
  only ever touch C code through `_native.py`.
- **Every module's docstring states what it implements and what it
  deliberately does not** — this project's consistent culture. A reader
  should be able to tell the scope just from reading it.

## 8. Bugs found and detailed work history

See [`DEVLOG.md`](DEVLOG.md) for the full log of real bugs found during
porting (with root-cause analysis and measured before/after numbers) and
the session-by-session history of how packaging, benchmarking, and
performance work got done.
