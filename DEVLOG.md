# sedumipy Development Log

This is the detailed, session-by-session record of bugs found while
porting SeDuMi, and the narrative history of how each piece of work got
done (including false starts and revised judgments along the way). It's
kept separate from [`CONTRIBUTING.md`](CONTRIBUTING.md) because that file
is meant to stay a compact, current-state guide for someone starting
work on this project — goal, architecture, directory layout, workflow,
scope, and conventions — while this file is the archival trail behind
it, useful mainly when you want to know *why* something is implemented
the way it is, or want the full context before touching one of the
trickier corners of the algorithm (Lorentz-cone step-length fallbacks,
ADA symbolic Cholesky reuse, dense-column PCG, ...). It grows with every
session and isn't meant to be read start to finish — search it for the
module or symptom you care about.

## Bugs and gotchas found so far (lessons learned)

- **Real SeDuMi's "all-or-nothing across every block" branch is not
  ported as-is — clamp block by block instead.** `widelen.m`/
  `trydif.m`/`maxstep.m` all contain, for the Lorentz-cone blocks'
  discriminant,

  ```matlab
  tmp = halfxz.^2 - detxz;
  if all(tmp > 0)          % <- decided once, across ALL blocks together
      lab2q = halfxz + sqrt(tmp);
  else
      lab2q = halfxz;      % <- one bad block degrades every block
  end
  ```

  a **global all-or-nothing branch** (`maxstep.m` has an isomorphic one
  for `norm2`). At first glance this looks like "a safeguard against
  passing a negative value to `sqrt`," but the discriminant is
  **always non-negative in exact arithmetic**: the two eigenvalues this
  expression produces are `lab2q` and `detxz/lab2q`, whose product is
  `detxz` and whose sum is `2*halfxz`, so `tmp` is identically
  `((lab1-lab2)/2)^2` — a perfect square. It only goes negative from
  rounding error when a single block's two eigenvalues nearly coincide,
  and when they coincide *exactly* it becomes exactly 0 — which real
  SeDuMi's strict `> 0` also rejects (a routine occurrence in problems
  structured from duplicates of a single Lorentz block). In practice,
  the observed triggering value wasn't a tiny negative number but
  **exactly `-0.0`**.

  So the port **clamps per block** instead:

  ```python
  lab2q = halfxz + np.sqrt(np.maximum(tmp, 0.0))
  ```

  For the one block whose discriminant is non-positive, this is exactly
  equivalent to real SeDuMi's fallback (`sqrt(0)==0`); for every other
  block, it keeps the exact expression real SeDuMi was discarding, while
  still satisfying real SeDuMi's own safety goal of never passing a
  negative value to `sqrt` — in other words, this isn't a tradeoff, it
  is **strictly more accurate than either of real SeDuMi's branches**.

  Note that `maxstep.m`'s version is a **safety bug, not just an
  accuracy one**: in the fallback, `norm2` is used unsquared in
  `reltr - norm2` without ever taking a square root (a dimensional
  mismatch, subtracting a linear quantity from a squared one). When the
  discriminant is below 1 — the normal case after scaling — `v <
  sqrt(v)`, so real SeDuMi **overestimates** the step length to the cone
  boundary. Measured on nb_L2, this fired 5 times out of 64, and all 5
  were overestimates.

  Effect (measured on DIMACS, toggling all three sites between real
  SeDuMi's branch and the clamp):

  | Problem | Real SeDuMi's branch | Clamped | Real Octave/MEX |
  |---|---|---|---|
  | nb_L2 | numerr=2, iter=10 | **numerr=0, iter=16** | numerr=0, iter=16 |
  | nql180old | numerr=2, iter=12 (cx=18.08 vs by=7.08) | **numerr=1, iter=42** (cx≈by to 8 digits) | numerr=1, iter=54 |
  | qssp30old | numerr=2 (cx=6.6017 vs by=6.3582) | **numerr=1** (cx=6.496695, published value 6.4966749) | numerr=2 (real SeDuMi also fails) |

  The lesson: **when a branch in real SeDuMi's code exists as a
  safeguard against a "mathematically impossible" case, porting it
  verbatim can make that safeguard fire far too broadly.** "It's written
  that way in real SeDuMi" is the right default porting policy, but
  checking the invariant the condition is trying to protect (here, "the
  discriminant is a perfect square, hence non-negative") can uniquely
  determine an implementation that is strictly better than the
  original.
- **MATLAB's pass-by-value semantics vs. ctypes' in-place mutation.**
  `fwsolve`/`bwsolve` are direct ctypes bindings of C kernels
  (`fwblkslv.c`/`bwblkslv.c`) that mutate their buffers in place, so
  porting MATLAB code that reuses a buffer across calls **breaks
  silently**. Solved in `pcg.py`'s `sparfwslv`/`sparbwslv` by
  consistently copying the buffer before every call. The same trap can
  apply to other in-place native functions (`fwdpr1`/`bwdpr1`, etc.), so
  always suspect it when binding something new.
- **mexFunction's automatic slicing.** Some MEX kernels, such as
  `psdframeit.c`/`psdinvjmul.c`, accept either a "short array covering
  only the PSD part" or a "full-length L+Q+PSD array," and when given
  the full-length form, automatically skip the offset
  (`x += cK.lpN + 2*cK.lorN`). Phase 2's ctypes bindings originally
  forgot to reproduce this auto-slicing, which surfaced as a real bug
  in `psdinvjmul` (found while porting `wregion.py`). Watch for this
  kind of branch in the original `mexFunction` when adding a new
  binding.
- **Undefined behavior in a qsort comparator.** `sortnnz.c`/`iswnbr.c`
  cast a comparator returning `signed char` to
  `int(*)(const void*,const void*)` before passing it to `qsort()`,
  which was confirmed on real hardware to cause non-deterministic
  behavior. Functions like this are **not ctypes-bound**; instead, the
  algorithm described in the comments is rewritten directly in Python
  (`neighborhood.py`'s `iswnbr` is one example).
- **`symbchol.m`'s fully-dense-matrix branch.** When ADA is completely
  dense (every entry nonzero), the real `symbchol.m` skips minimum-
  degree ordering (MMD) and uses identity ordering plus a single giant
  supernode directly. Not reproducing this branch and always calling
  `ordmmd` still converges, but **the iteration count drifts from the
  Octave version** (confirmed on a small dense test problem).
  `symbchol.py` reproduces this branch exactly
  (`_native.symbolic_cholesky_dense`).
- **MATLAB's `'` is the conjugate transpose.** Even for real-valued
  arrays, in expressions involving complex numbers (e.g.
  `posttransfo.m`'s `(x'*prep.QR)'`), the difference between a plain
  transpose and a conjugate transpose affects the result. Always watch
  for this when porting.
- **`optstep.m`'s `sum(K.s)!=0` branch is effectively dead code.**
  `sedumi.m` only calls `optstep` when `lponly = (K.l==length(c))`,
  which forces `K.q`/`K.s` to both be empty. So even though
  `optstep.m` itself has a branch for the PSD case, the real call path
  can never reach it. This kind of "dead code whose unreachability can
  be proven from the caller's own condition" is fine to leave unported,
  behind a `NotImplementedError` — but only after actually confirming
  the caller's condition (never on a guess of "probably dead code").
- **`sparfwslv`/`sparbwslv` (`pcg.py`) must gather/scatter by `L.perm`
  internally.** The real `fwblkslv.c`/`bwblkslv.c` do the gather/scatter
  by `L.perm` *inside* the call itself — forward substitution is
  `y = L\b(L.perm)`, back substitution is `y(L.perm) = L'\b`. The
  version of this port before dense-column optimization was added
  omitted this and instead had every caller of `loopPcg`/`wrapPcg`
  (tests included) apply/undo `L.perm` manually and consistently from
  the outside — since PCG converges to the correct solution as long as
  "the internal indexing convention is self-consistent across the whole
  call," this was numerically correct as long as there were no dense
  columns (i.e. `Lden` was the identity). But once `deninfac.py` started
  building a **genuine** (non-identity) `Lden` that actually permutes
  `Ad` by `L.perm`, combining the `Lden` term with the `sparfwslv` term
  within the same `loopPcg` iteration (as in
  `fwdpr1(Lden, sparfwslv(L,r))`) hit an indexing-convention mismatch
  that broke the iteration (this surfaced as `iter` still matching but
  the intermediate residual diverging). Fixed by making `sparfwslv`/
  `sparbwslv` themselves do the gather/scatter exactly like the real C
  kernels, while restoring `wrapPcg`/`sdfactor.py` to never permute
  anything, exactly as in the real `.m` files — a textbook case where a
  "make the caller compensate" workaround can look fine until a
  non-trivial later feature (here, dense columns) breaks it.
- **Dual solutions `y` are non-unique on test fixtures whose `At` is
  rank-deficient.** `tests/fixtures/sedumi/lp_socp_sdp_dense_feasible.mat`
  (used by the dense-column-optimization end-to-end test) has a 23-row,
  20-column `At` with `rank(At)=16` (pinning the dense-column detection
  floor `h` to `NORMDEN=5` required lowering the PSD block's background
  row density to `0.02`, which made some rows structurally all-zero; an
  exhaustive search over seeds 1-400 never produced a full-rank
  (`rank(At)=20`) case). With an equality-constraint system whose `At`
  isn't full rank, the optimal dual solution `y` is non-unique, since
  shifting it by any `delta` with `At@delta=0` gives the same objective
  value and the same dual slack (`s=c-At@y`). Even though the Python and
  Octave versions match exactly on `iter`/`numerr`/`pinf`/`dinf` and `x`
  agrees to floating-point tolerance (i.e. the algorithm behaves
  identically), `y` differs by about 5 in 2-norm along the degenerate
  direction (confirmed by both `At@(y_py-y_oct)` and
  `b@(y_py-y_oct)` being essentially zero). Because of this,
  `test_sedumi_dense_matches_octave` compares dual feasibility/
  optimality (`c-At@y` and `b@y`) rather than requiring exact agreement
  on `y`. **For the same reason, when creating a new randomly generated
  test fixture whose `At` isn't guaranteed full rank, check
  `rank(At)==m` before comparing `y` exactly.**
- **The `K.s==0` path's one-time ADA symbolic Cholesky ordering was
  missing part of the sparsity pattern that depends on the Lorentz
  cone's arrow term.** `sedinit.py` always initializes the scaling
  point's `d["q2"]` (each Lorentz cone's arrow-part scalar) to exactly
  0, matching `sdinit.m`'s `d.q2 = zeros(...)`. The `K.s==0` branch of
  `sedumi.py` used to compute `ADA = getada(A,K,d,DAt)` once from this
  initial `d` (with `d.q2=0`), and pass its **numeric** sparsity pattern
  straight to `symbchol(ADA)`, reusing that pattern for every subsequent
  iteration. But in `getada.py`'s Lorentz term (`DAt_q.T @ DAt_q`, where
  `DAt.q[k,j] = d.q1[k]*Aj[k] + d.q2[k]*(...)`), when `d.q2=0` the
  contribution from any pair of constraints `(i,j)` that share the same
  Lorentz block but don't overlap in a row structurally vanishes
  entirely — meaning the ADA at iteration 1 has an
  under-sized sparsity pattern that's missing positions that actually
  become nonzero once `d.q2` grows in later iterations. Since
  `numeric_cholesky` cannot write outside this (fixed) symbolic pattern,
  the Cholesky factorization becomes progressively inaccurate starting
  from the iteration where `d.q2` grows (in practice, iteration 3-5
  onward), PCG's preconditioning degrades, the iteration count pins at
  its cap, and the solver ultimately converges to the wrong answer.
  Confirmed on `vendor/sedumi-upstream/examples/nb.mat` (LP + 396 SOCP
  blocks): real Octave SeDuMi converges in 20 iterations with
  `numerr=0`, while the unfixed Python version reported `numerr=2`
  (a serious numerical error) after 9 iterations, returning a
  completely different value. Real `sedumi.m` avoids this by always
  building the structural (value-independent) pattern via
  `getsymbada.m` before `sdinit` runs, regardless of `sum(K.s)`'s value
  (this file previously noted this, and `sedumi.py`'s own SCOPE
  docstring flagged it, as a "known difference from a simplified
  implementation"). **Fix**: in `sedumi.py`'s `K.s==0` branch, build the
  one-time ADA passed to `symbchol()` from a local, throwaway copy of
  `d` whose `d["q2"]` is forced nonzero, without touching the `d`/`DAt`/
  `ADA` actually used in each subsequent iteration.
  `tests/test_golden_end_to_end.py` (Phase 5, see item 3 in the work log
  below) confirms all five problems — `nb`/`arch0`/`control07`/`trto3`/
  `OH_2Pi_STO-6GN9r12g1T2` — now match the real Octave results.
- **`cpspdiag`, called by `getada3`'s `K.s==0` branch, hit the same
  qsort/bsearch undefined-behavior bug as `sortnnz.c`/`iswnbr.c`.**
  `cpspdiag` looks up ADA's diagonal entries via the `ibsearch` macro in
  `blksdp.h` (i.e. the standard library's `bsearch()`). `ibsearch`
  passes it a comparator, `icmp()`, that returns `char` cast to
  `COMPFUN` (`int(*)(const void*,const void*)`) — also undefined
  behavior. In practice, this port's build of `bsearch()` never
  successfully found the diagonal entries, so `absd` was always all
  0.0 even though the diagonal entries were present and properly sorted
  (found via `tests/test_getada.py::test_getada_no_psd_blocks`).
  However, since `getada3` itself is only ever called when
  `has_psd=True` (i.e. `K.s` is non-empty, so `sdpN>0` internally), this
  `sdpN==0` branch is unreachable dead code in real usage and had no
  practical impact. Resolved the same way as `sortnnz`/`iswnbr`:
  `cpspdiag` is no longer ctypes-bound, and instead reimplemented in
  Python reading the diagonal directly via `scipy.sparse`'s
  `.diagonal()` (matching what `cpspdiag.c`'s own doc comment says it
  was meant to do).

## Work log: investigation history and fixes

Numbered roughly in the order these were tackled. Most items below are
struck through and marked **Done**/**Resolved**/**Fixed** — this is a
log of how the project got to its current state (see
[`CONTRIBUTING.md`](CONTRIBUTING.md) §2 and
[`docs/status.rst`](docs/status.rst) for that current state itself), not
a live task list. As of this writing the only item that's still
genuinely open is part of item 4 (Phase 6): publishing to PyPI.

1. ~~**Phase 3-a: finish the public API for thin MEX-wrapper `.m`
   files.**~~ **Done.** Cross-checking `install_sedumi.m`'s list of MEX
   build targets against `_native.py`'s bindings found: every real MEX
   kernel that's actually used is already collected in `_native.py`,
   each called as `_native.xxx()` from the appropriate higher-level
   module (`getdense.py`/`getdatm.py`/`pcg.py`/`cone.py`/
   `updtransfo.py`/`wregion.py`/`sdinit.py`/`getada_psd.py`/
   `symbchol.py`/`symbcholden.py`, etc.) — no "bound but not wired up"
   gaps were found (`incorder`/`iswnbr` are the only two deliberately
   *not* ctypes-bound, to avoid qsort undefined behavior, and instead
   exist as pure-Python implementations in `incorder.py`/
   `neighborhood.py` — a known, intentional design choice). Conversely,
   7 bindings in `_native.py` (`realdot`/`realssqr`/`scalarmul`/
   `addscalarmul`/`blkmul`/`mJdetd`/`cholsplit`) turned out to be called
   from nowhere else, but each was confirmed **unused in real SeDuMi
   itself** too (`blkmul.c`/`mJdetd.c` aren't even in
   `install_sedumi.m`'s MEX build target list — dead code in the
   original already; `cholsplit()`'s output, `L.split`, doesn't appear
   in `blkchol.c`'s mex argument list and is never read even in real
   SeDuMi; `realdot` etc. are BLAS-style helpers bound only for Phase
   1's smoke test, with no independent MEX target of their own).
   Concluded no further work is needed, and documented this inventory
   itself, including why each is unused, in `_native.py`'s module
   docstring.
2. ~~**Phase 4: high-level API and I/O compatibility layer.**~~
   **Done.**
   - Top-level API: `import sedumipy; sedumipy.sedumi(A,b,c,K)` now
     works (previously only reachable via the `sedumipy.sedumi`
     submodule). `__init__.py` does `from .sedumi import sedumi`, and
     it's confirmed that importing the `sedumipy.sedumi` submodule
     before or after doesn't change which object the name refers to
     (Python's `sys.modules` cache means the parent package's attribute
     is only overwritten the first time the submodule is imported).
     `sedumi()` also now accepts individual options via `**kwargs` in
     addition to the `pars` dict (e.g. `sedumi(A,b,c,K,eps=1e-9)`).
   - `.mat` I/O: `matio.py` (`read_mat`/`write_solution_mat`). SeDuMi
     problem files are just plain MATLAB structs (there's no `.m` file
     to port this from), so unlike other modules this is new code
     specific to this port, not a port. Handles both `A`/`At`
     orientations and cases where `b`/`c` are stored sparse (confirmed
     against `vendor/sedumi-upstream/examples/*.mat`).
   - SDPA sparse format (`.dat-s`) read/write: `sdpa.py`
     (`read_sdpa`/`write_sdpa`). `read_sdpa` is a faithful port of
     `conversion/fromsdpa.m` (confirmed against a real Octave oracle,
     `tools/generate_sdpa_oracle.m`/`tests/fixtures/sdpa/`).
     `write_sdpa` has no counterpart in real SeDuMi (real SeDuMi's
     `conversion/writesdp.m` writes a different, unrelated format,
     SDPpack, not SDPA) and is new code, but manually confirmed correct
     by writing out `vendor/sedumi-upstream/examples/arch0.mat` with
     `write_sdpa` and reading it back with real Octave's `fromsdpa.m`,
     matching the original `(At,b,c)` exactly (`write_sdpa` explicitly
     rejects `K.q`/`K.r` with a `ValueError`, since SDPA format can't
     represent them).
3. ~~**Phase 5: verification and benchmarking.**~~ **Done.**
   `tests/test_golden_end_to_end.py` runs `sedumipy.sedumi()` against
   the Phase 0 golden-reference problems
   (`vendor/sedumi-upstream/examples/`) and verifies the result matches
   real Octave (see CONTRIBUTING.md §2 "Phase 5" and the two bugs fixed
   above). Performance benchmarks live in `tools/benchmark_examples.py`
   (see that script's own docstring for how to run it). Measured on
   this environment (Octave built locally for comparison; absolute
   numbers are environment-dependent on CPU/core count, treat as
   indicative):

   | problem | m | N (=length(c)) | Python (s) | Octave/MEX (s) | iter |
   |---|---:|---:|---:|---:|---:|
   | nb | 123 | 2383 | 3.0 | 0.9 | 20 |
   | arch0 | 174 | 56197 | 2.5 | 2.4 | 31-32 |
   | control07 | 666 | 6125 | 9.3 | 9.2 | 40 |
   | trto3 | 544 | 398977 | 18.2 | 19.8 | 60 |
   | OH_2Pi_STO-6GN9r12g1T2 | 948 | 240720 | 34.4 | 34.8 | 20 |

   On the smallest problem (`nb`), Python-side overhead (function calls,
   NumPy array allocation, the cost of crossing the ctypes boundary)
   dominates, running about 3x slower than the Octave/MEX version; as
   problems grow, the native C-kernel compute time dominates instead,
   and on medium-to-large problems (`arch0` and up) the two are roughly
   on par or the port is slightly faster. `arch0`'s `iter` differs from
   Octave's by exactly one (31 vs. 32) due to accumulated floating-point
   rounding differences at this scale (as noted above, `test_sedumi_
   matches_octave`'s exact `iter`-match requirement is specific to small
   synthetic fixtures; a one-iteration difference like this is expected
   at real-problem scale), and both `cx`/`by` match the expected values,
   so there's no practical impact.
4. **Phase 6: packaging.** ~~Investigate how to bundle `libsedumi.so`
   (currently the prebuilt binary is committed directly to the
   repository)~~ — that description was itself stale: in reality,
   `libsedumi.so`/`.dylib` are gitignored and **not committed to the
   repository**; `_native.py`'s `_ensure_built()` automatically calls
   `tools/build_libsedumi.sh` on first import to build it on the spot
   (development's `pip install -e .[test]` depends on this). This
   approach works for an editable install, but **breaks for a real
   installed wheel** (`csrc`/`tools` aren't bundled as part of the
   `sedumipy` package, so there's nothing to recompile from on import,
   and there's no guarantee the end user's environment even has gcc/BLAS
   dev headers installed).

   **What was done:** added a custom `build_ext`-overriding step
   (`BuildLibsedumi`) to `setup.py` that runs the same compile command
   as `tools/build_libsedumi.sh` exactly once during `pip install`/
   `python -m build --wheel`, building `libsedumi.so` directly into
   `build_lib/sedumipy/` so it gets bundled into the wheel
   (`_native.py`'s `_ensure_built()` is unchanged, so a wheel install
   finds the file already present and does nothing, while an editable
   install still builds it on first import as before — both paths
   coexist). Registering one sourceless `Extension` is purely a trick to
   make setuptools correctly mark the wheel as platform-specific (not
   `py3-none-any`).

   **Confirmed in this environment:** `python -m build --wheel` builds
   a wheel tagged `cp311-cp311-linux_x86_64` with `libsedumi.so` bundled
   inside it, and installing that wheel with `pip install` into an
   isolated virtualenv with no access to this repository's `csrc`/
   `tools` at all, `import sedumipy; sedumipy.sedumi(...)` works
   correctly (Linux only, in this environment).

   **Not yet verified (no Docker daemon available in this
   environment):** actually running `cibuildwheel` itself (the
   `[tool.cibuildwheel]` settings were added to `pyproject.toml`, but an
   actual build inside a manylinux container is unverified), building on
   macOS/Windows (`tools/build_libsedumi.sh` assumes gcc and doesn't
   support Windows's `cl.exe`), and distributable portability of the
   dynamic link against `libblas` (`ldd` shows it dynamically linked
   against `libblas.so.3`/`libopenblas.so.0`; a real PyPI-distributable
   manylinux wheel would need `auditwheel repair` to bundle these, or a
   switch to static linking — not done in this pass).

   **Additional work done in a later session (CI setup, macOS/Windows
   support):** added `.github/workflows/{ci,wheels,docs}.yml`, moving
   the Docker/macOS/Windows runs this environment couldn't do onto real
   GitHub Actions runners.
   - **macOS:** `before-all` was yum/apt-only, and neither exists on
     macOS, so the original config was guaranteed to break there.
     Fixed `tools/build_libsedumi.sh` to link the OS-standard Accelerate
     framework (`-framework Accelerate`), which needs no Homebrew at
     all (also switched the compiler to `cc` instead of `gcc`, since
     macOS's `gcc` is normally just a clang alias).
   - **Windows:** rather than porting to `cl.exe`, chose to use MSYS2's
     MinGW64 toolchain (`mingw-w64-x86_64-gcc`/
     `mingw-w64-x86_64-openblas`) — since `libsedumi.dll` is a plain
     ctypes-loaded DLL rather than a `PyInit_*`-exporting CPython
     extension, it never needed to match whatever compiler built
     Python, and `sedumi_platform.h`'s BLAS symbol-naming convention
     (`FORT(x) = x##_`) is the same across Linux/macOS/Windows and
     matches OpenBLAS, so this was judged far lower-risk than an MSVC
     port. `setup.py`/`_native.py` were changed to explicitly invoke
     `bash tools/build_libsedumi.sh ...` via `bash` when
     `sys.platform == "win32"` (Windows doesn't interpret shebangs).
     The resulting `libsedumi.dll` dynamically links `libopenblas.dll`
     and mingw runtime DLLs, so a distributable wheel needs
     `delvewheel repair` to bundle them (the Windows analog of
     `auditwheel`/`delocate`; not included in cibuildwheel's Windows
     defaults, so configured explicitly in `pyproject.toml`'s
     `[tool.cibuildwheel.windows]`).
   - **On real-hardware verification:** development on this repository
     happens on Linux, with no macOS/Windows machine available. The
     changes above are practically verified by actual macOS/Windows
     runners on GitHub Actions running `ci.yml`/`wheels.yml` (the
     intended workflow is to open a PR and check the CI results).

   **Still open:** publishing to PyPI has not been done yet.
5. **`getdatm.py`'s OOM fix (always building `DAt.q` sparse) was, in
   turn, slowing down small-to-medium problems where dense is actually
   faster. Now fixed.** Making the `has_psd=False` (LP+SOCP only,
   `K.s==0`) path always build `DAt.q`/`ADA` sparse (correct as an OOM
   fix) ended up slowing down problems where `m` is small and ADA is
   effectively dense (e.g. `nb.mat`, m=123) by about 36% (confirmed with
   `cProfile`: the sparse-sparse product `csr_matmat` accounted for 56%
   of total time). Fixed by determining `is_dense` once, from the
   density of the structural ADA pattern `getsymbada()` computes once
   (reusing the existing 0.9 threshold), and wiring `getDAtm()`/
   `getada()` from `sedumi.py` to switch between dense (numpy, BLAS
   matmul) and sparse (scipy, avoids OOM) accordingly. Both branches are
   bug-for-bug identical in the values they produce (confirmed no
   regressions across the full existing test suite and benchmarks).
   Measured on `nb.mat`: 2.12s with sparse forced (right after the OOM
   fix) → 1.75s after hybridizing (close to the original 1.56s with
   dense forced). Large problems (`nql180`/`qssp180`, m~1.3e5) still go
   through the sparse path, so the OOM fix still holds for them.
6. ~~**DIMACS `nb_L2`'s numerr=2: root cause identified (fix
   deferred).**~~ **Resolved. Fixed** (see the "all-or-nothing branch"
   entry earlier in this file). An earlier session had fully traced the
   cause to `widelen.py`'s `all(tmp>0)` global branch, but deliberately
   deferred fixing it, reasoning it was "the algorithm's own chaotic
   sensitivity, and changing the branch's numerics has unpredictable
   effects" — **that judgment was overly conservative**. `tmp` isn't an
   arbitrary quantity; it's `((lab1-lab2)/2)^2`, a **perfect square**, so
   it's always non-negative in exact arithmetic, and the only way it
   goes negative is rounding error — meaning "which branch to take" is
   actually uniquely determined. Remeasuring this time also showed the
   value wasn't a tiny negative number like the previously-assumed
   `±1.78e-15`, but **exactly `-0.0`**: real SeDuMi's strict `> 0` was
   simply rejecting zero (which happens routinely on problems built from
   duplicates of a single Lorentz block, where two eigenvalues coincide
   exactly). Clamping with `np.sqrt(np.maximum(tmp, 0.0))` block by
   block leaves the affected block exactly equivalent to real SeDuMi's
   fallback, keeps the exact expression real SeDuMi discarded for the
   other 838 blocks, and still never passes a negative value to `sqrt`.
   Result: **nb_L2 improved from numerr=2/iter=10 to numerr=0/iter=16**,
   matching both the real Octave/MEX build's iteration count (16) and
   the published value (`-1.62897198`, this port now gets
   `-1.628971959`). This also removed the need for the earlier symptomatic
   workaround of forcing `stepdif=1` (the `pars` defaults are unchanged).
   The same branch also existed in `trydif.m` (a verbatim copy) and
   `maxstep.m` (isomorphic, with an even worse fallback), so all three
   sites were fixed. What follows is the previous session's record of
   how the root cause was tracked down (kept as-is):

   Dumping `ADA`/`d`/`DAt.q` at iterations 1-3 from both the real
   Octave/MEX build (`vendor/sedumi-upstream`, with `octave`/
   `liboctave-dev`/`libopenblas-dev` installed in this environment and
   built via `install_sedumi -rebuild`) and this port, and comparing
   them directly: `d.l`/`d.det` (the LP/trace part) match to floating-
   point tolerance (~1e-13) through iteration 3, and `d.q1`/`d.q2` (the
   Lorentz-cone scaling point) also match the same way up through the
   start of iteration 2, but the `d` produced by iteration 2's step
   (i.e. the `d` used at iteration 3) has `d.q1`'s largest component
   diverging by about 15% relative error. This lines up exactly with the
   iteration where `err["kcg"]`/`Lsd["kcg"]` jumps from the real
   hardware's 1/1 to 6/5. That's as far as the earlier session's record
   went.

   **This time, the investigation was carried further to a complete root
   cause.** The earlier session had concluded "`updtransfo.py` was
   audited line-by-line against `updtransfo.m` with no differences
   found," but this time, beyond the line-by-line audit, it was also
   **verified by actually running it**: a temporary debug `save()` was
   inserted into the real Octave `wregion.m` (an uncommitted, throwaway
   patch), dumping iteration 2's `wregion` outputs `xscl`/`zscl`/`w`
   (and the preceding `d`/`K`) to a `.mat` file, then feeding that
   **directly** into Python's `updtransfo()` — which reproduced real
   iteration 3's `d.q1`/`d.q2` bit-for-bit (matching to the 10th
   decimal digit of `max|q1|`). In other words, `updtransfo.py` really
   is innocent — the divergence happens before `updtransfo`, in the
   computation of `xscl`/`zscl`/`w` themselves.

   Next, dumping this port's own computed `xscl`/`zscl`/`w` at iteration
   2 the same way and comparing directly against real Octave's values:
   `xscl`/`zscl` differ by an absolute error of ~1.5e-13 (floating-point
   noise level — an excellent match for two independent implementations
   agreeing on a 4196-dimensional vector), and `w["tdetx"]`/
   `w["tdetz"]` similarly agree to ~1e-13 to ~1e-12, yet **`w["lab"]`
   alone diverged by as much as 7.6 in absolute error**. Since
   `w["lab"]` should be computed almost directly from `tdetx`/`tdetz`,
   this disproportion was the decisive clue.

   Looking at `widelen.py`'s `_build_w()` (a direct port of the
   corresponding part of `widelen.m`), the Lorentz-cone eigenvalue-like
   quantity `lab2q` is computed as

   ```python
   tmp = halfxz**2 - detxz
   if np.all(tmp > 0):        # widelen.m: if all(tmp > 0)
       lab2q = halfxz + np.sqrt(tmp)
   else:
       lab2q = halfxz          # all 839 blocks fall back here
   ```

   a **single global all-or-nothing branch across all 839 Lorentz-cone
   blocks** (if even one block's discriminant `tmp` is non-positive,
   every other block, all 838 of them, also gets pushed to the less
   accurate fallback expression). Actually computing `tmp` from
   iteration 2's `xscl`/`zscl` showed block 838 of 839 (0-indexed block
   396)'s `tmp` was `+1.78e-15` on real hardware but `-1.78e-15` in this
   port — **a value on the exact knife-edge where only the sign
   flips** (every other block's `tmp` was comfortably positive on both
   sides). Even though `xscl`/`zscl` themselves agree to ~1e-13, this
   one block happened to sit exactly on the zero crossing, so a tiny
   rounding difference between two independent floating-point pipelines
   (NumPy/SciPy plus this port's own C kernels, vs. Octave plus real
   hardware BLAS, differing in internal summation order and BLAS
   implementation) alone flips `all(tmp>0)`'s truth value, causing
   `lab2q` (and hence all of `w["lab"]`) to be computed by a completely
   different expression, which then makes the scaling-point update in
   the following `updtransfo` diverge substantially — this is the
   complete causal chain identified this time. (Confirmed directly:
   feeding real hardware's `xscl`/`zscl` into this port's own
   `_build_w()` reproduces real hardware's `w["lab"]` bit-for-bit, and
   conversely, feeding this port's own `xscl`/`zscl` back in reproduces
   `tmp[396]` flipping negative and `all(tmp>0)` becoming False.)

   **This `all(tmp>0)` global branch itself exists as-is in real
   `widelen.m`** (readable as a deliberate safety measure: if even one
   block's discriminant could go negative, fall back to a conservative
   expression across all blocks at once, to avoid ever passing a
   NaN/complex value to `sqrt`) — this is not an error this port
   introduced on its own. The situation of two independent
   implementations landing on opposite sides of rounding error for a
   quantity that sits exactly on a zero crossing is, in itself
   (reproducible with the right random seed), inherent chaotic
   sensitivity in the algorithm, and was judged as something no single
   line of `updtransfo.py`/`widelen.py`/`tdet`/`ddot` etc. could resolve.
   So **the decision at the time was to deliberately defer a fix** (the
   same reasoning as the item noted earlier in this log, "forcing
   `stepdif=1` would solve nb_L2, but changing a `pars` default for
   every problem is too large a side effect to accept" — this session's
   progress was in fully identifying the cause and concluding "the
   current implementation is fine as is," rather than a band-aid fix for
   the symptom).
7. ~~**`nql180`/`qssp180`'s numerr=2 — re-checked and it turned out to
   already be fixed (`nql180old` remained an open, real robustness gap;
   `qssp180old` is now verified, confirmed not a porting bug).**~~
   **The `nql180old` gap is also resolved** (see the "addendum" at the
   end of this item, resulting from the "all-or-nothing branch" fix
   above). What follows is the history:
   The earlier session's record ("still numerr=2 after a few iterations
   even after the OOM fix") was stale. Running these with the item-5
   dense/sparse hybridization in place, loaded via `matio.read_mat()`:
   - `nql180` (m=226,802, no PSD blocks): converges cleanly with
     **numerr=0, iter=16, ~39s** (confirmed via internal consistency —
     `cx`≈`by`, `feasratio`→1, `r0`=1e-8 — since the DIMACS README's
     reference value is "N/A").
   - `qssp180`: converges cleanly with **numerr=0, iter=42, ~249s**
     (likewise confirmed via internal consistency, no reference value
     available).
   On the other hand, `nql180old`, part of the same "old (legacy,
   deprecated) formulation" family as `nql30old`/`qssp30old`, was
   compared against a real Octave/MEX build in this environment
   (`install_sedumi -rebuild`, as above): **the two do not fail the same
   way** — real hardware grinds on to `iter=54` and ends with `numerr=1`
   (accuracy capped at `pars.bigeps` but not a hard failure, though
   numerically demanding enough that the console prints `skip=5361`
   worth of skipped Cholesky pivots), while this port returns
   `numerr=2` (a complete failure) at `iter=27` (`feasratio=0.90`,
   `r0=0.53`) — unlike `nql30old`/`qssp30old`, where "real SeDuMi fails
   the same way too, so it's not a porting bug," this was a **genuine
   robustness gap**, failing earlier and worse than real SeDuMi.
   `qssp180old` (the largest problem in this family, ~36MB) didn't
   finish within the earlier session's time budget (550s for both the
   Python and real-hardware versions) and was left unverified, but
   **this session gave it a larger time budget and let both run to
   completion, settling the question**: the real Octave/MEX build fails
   with **iter=30, numerr=2** (1705s measured with `tic`/`toc` on an
   environment with `install_sedumi -rebuild` already run). Running the
   same file through this port via `matio.read_mat()` also failed with
   **iter=30, numerr=2** (3557s measured; monkey-patching `wregion()` to
   log per-iteration elapsed time confirmed the iterations were
   proceeding steadily rather than hanging partway through). **The
   failing iteration number (30) matches exactly**, with none of the
   "fails earlier and worse than real SeDuMi" robustness gap seen on
   `nql180old` — confirming this is the same genre of problem as
   `qssp30old`/`nql30old`, where real SeDuMi fails the same way too (not
   a porting bug).
   **Addendum: the `nql180old` robustness gap is also resolved** (via
   the `all(tmp>0)` fix in item 6 above). This had been the one open
   item this section still listed as "fails earlier and worse than real
   SeDuMi = a genuine robustness gap," but remeasuring with all three
   sites clamped block-by-block:

   | nql180old | numerr | iter | cx vs by |
   |---|---|---|---|
   | Real SeDuMi's branch (all 3 sites) | **2** (complete failure) | 12 | 18.08 vs 7.08 |
   | Block-wise clamp | **1** | 42 | match to 8 digits |
   | Real Octave/MEX build | 1 | 54 | - |

   Where real SeDuMi's behavior breaks down at iteration 12 with `cx`/
   `by` still 2.5x apart, the clamped version converges to `cx=
   0.9311428505`/`by=0.9311428684` (matching to 8 digits) and ends with
   `numerr=1` — **reaching the same `numerr=1` in fewer iterations (42)
   than the real hardware build (54)** — closing the "fails earlier and
   worse than real SeDuMi" gap (confirmed via internal consistency,
   since the DIMACS README's reference value is "N/A" for this
   problem). This problem triggers the fallback at an unusually high
   rate (`widelen`: 6 of 12 calls = 50%; `maxstep`: 5 of 50), showing
   real SeDuMi's branch was firing constantly here.
   The same fix also turns `qssp30old` from `numerr=2` into `numerr=1`,
   and moreover now returns **a solution matching the DIMACS README's
   published value of `6.4966749`** (`cx=6.496695`) — since **the real
   Octave/MEX build itself fails with `numerr=2`** on this problem, this
   case isn't "on par with real SeDuMi" but **better than real SeDuMi**.
   Together with `nql30old` (published value `0.9460`), both were
   promoted from `tests/test_benchmarks.py`'s exclusion list to
   parametrized tests checked against their published values.
8. **Fixed two Python-level performance bugs found via `qssp180old`'s
   cProfile output (the main reason this port was about 2.1x slower
   than real SeDuMi; now fixed).**
   While confirming `qssp180old`'s `numerr=2` match in item 7, an
   additional question came up: "why is this port about 2.1x slower
   than real SeDuMi (1705s on real hardware vs. 3557s in this port)?"
   Profiling just the first 5 iterations with `cProfile`
   (`pars["maxiter"]=5` to cut it short, 462s) found two issues, both
   pure Python-level overhead with zero effect on the computed result:

   - **`_native.fwsolve()`/`bwsolve()` were re-converting
     `L_csc.indptr`/`.indices`/`.data` via
     `np.ascontiguousarray(..., dtype=np.uintp)` on every call.** The
     `L_csc` (`scipy.sparse.csc_matrix`) returned by
     `numeric_cholesky()` is the same object reused across an entire
     outer iteration's PCG loop, with only the right-hand-side vector
     changing — but since `scipy.sparse.csc_matrix` actually normalizes
     any `uintp`-typed indices passed to its constructor back to
     int32/int64 (confirmed directly), `_as_index_array()`'s recast back
     to `uintp` was faithfully happening on every single fwsolve/bwsolve
     call. Since `qssp180old` has a very large `nnz(L)~8.3e7`,
     `numpy.ascontiguousarray` alone accounted for 130 of the profiled
     462 seconds (28%). **Fix**: cache the converted arrays as
     attributes on the `L_csc` object itself
     (`_cached_csc_solve_arrays()`). Since a new `L_csc` object is
     created for every outer iteration, the cache can never go stale,
     and it was confirmed that `deninfac.py`/`sdfactor.py` never mutate
     `L_csc.data`/`.indices`.
   - **`_native.qblkmul()` (per-Lorentz-block scalar multiplication,
     ported directly to NumPy since `qblkmul.c`'s `mexFunction` has no
     standalone C function counterpart) processed blocks one at a time
     in a Python `for k in range(nblk)` loop.** `qssp180old` has as many
     as 65,341 Lorentz-cone blocks, and this loop body alone accounted
     for 49.7 of the profiled 462 seconds (across 334 calls). **Fix**:
     vectorized as `np.repeat(mu, block_sizes) * d` (replicating each
     block's `mu[k]` across that block's width, then multiplying in one
     shot — the same computation, same order, same values, just
     rewritten). While vectorizing, this briefly introduced a regression
     by missing a pattern where a caller passes a longer `d` than a
     single block's length (the original loop silently ignored the
     extra tail) — but the full test suite (247 tests) caught it
     immediately, and it was fixed by explicitly trimming with
     `d = d[:span]`.

   **Effect**: `qssp180old`'s first-5-iteration cProfile went from
   **462s to 153s (3.0x faster)**. `numpy.ascontiguousarray`'s own time
   dropped from 130s to 4.2s, and `qblkmul` disappeared from the top of
   the profile entirely. The largest remaining cost is
   `numeric_cholesky` (the real C-kernel Cholesky factorization, 75.5s)
   — which is actually the healthy outcome, since that's the dominant
   cost in real SeDuMi too. The full test suite (247 tests),
   `pytest -m mini` (46 problems checked against published SDPLIB/DIMACS
   reference values), and `pytest -m extended` (16 more, including
   problems with a high LP+SOCP ratio) all pass with no regressions.
   This has no effect on `qssp180old`'s `numerr=2` match itself (item
   7) — the computed result is bug-for-bug identical.
9. **Continuing item 8: investigated the remaining performance
   bottleneck — rebuilding ADA/DAt.q as sparse matrices turns out to be
   real SeDuMi's own design, and the `numeric_cholesky` uintp/int64
   round-trip can only be partially reduced.**

   Further analyzing the post-item-8 `qssp180old` first-5-iteration
   cProfile (153s) to characterize the remaining major costs:

   - **The cost of rebuilding ADA/DAt.q as scipy sparse matrices from
     scratch every iteration (8.8s inside scipy internals via
     `numpy.array`) turns out to be real SeDuMi's own design.**
     Reading the real `.m`-only (not MEX) `vendor/sedumi-upstream/
     getada.m` used by the `sum(K.s)==0` path (`qssp180old`'s path)
     shows it **creates a brand-new empty sparse matrix every time**
     via `ADA_sedumi_ = sparse([],[],[],m,m,nnz(ADA_sedumi_))` on a
     `global ADA_sedumi_`, and rebuilds it via the sparse-times-sparse
     operations `ADA_sedumi_ + DAt.q'*DAt.q`/`ADA_sedumi_ +
     Alq'*diag(sparse(scalingvector))*Alq` — exactly the same pattern
     as this port's `getada.py` (the equivalent sparse-times-sparse
     operations in `scipy.sparse`). `getDAtm.m` similarly rebuilds from
     scratch every time, after extracting via `extractA` (MEX):
     `spdiags(d.q1,0,nq,nq) * DAt.q` (creating a fresh sparse diagonal
     matrix and multiplying, every call). So this cost isn't a porting
     inefficiency but **the result of a faithful port**, and is judged
     out of scope for further fixing. (Aside: by contrast, the
     PSD-block-bearing `sum(K.s)!=0` path's `getada1.c`/`getada2.c`/
     `getada3.c` reuse the sparse pattern that `getsymbada` fixes once,
     writing **only the values** in place into the same global array
     `ADA_sedumi_` — real SeDuMi has an asymmetric design here, reusing
     when PSD blocks are present and rebuilding from scratch when they
     aren't; the latter path, which `qssp180old` uses, appears to simply
     never have been optimized in real SeDuMi itself.)
   - **`numeric_cholesky`'s `Lir.astype(np.int64)`/
     `Ljc.astype(np.int64)` calls (nearly all of item 8's 18s of
     `astype` time) are an unavoidable cost of using this container
     type**, stemming from `scipy.sparse.csc_matrix` never actually
     keeping index arrays as `uintp` — it always normalizes them to
     int32/int64 (confirmed directly) — so item 8's estimate (that this
     entire `astype` round-trip could be eliminated) was wrong. What
     could actually be eliminated was only the portion where
     `fwsolve()`/`bwsolve()`'s cache (`_cached_csc_solve_arrays()`) was
     re-converting `L_csc.indptr`/`.indices` (int64) back to `uintp` on
     its first access each iteration (once per iteration, 5 times
     total). **Fix**: have `numeric_cholesky` inject the `uintp`-typed
     Ljc/Lir it already holds right before constructing `L_csc` directly
     into the cache itself (`numeric_cholesky` populates
     `L_csc._sedumipy_solve_cache` itself), skipping this last
     re-conversion. **Effect: 153s to 138s (about 10% reduction)** —
     short of the originally estimated 16%. Cutting further would
     require replacing `Lnum["L"]`'s scipy sparse matrix with a custom
     lightweight structure that keeps `uintp` natively (a larger type-
     level refactor whose blast radius could reach the test code too —
     deferred for now).

   Confirmed no regressions on the full test suite (247 tests) and
   `pytest -m mini` (46 problems).
