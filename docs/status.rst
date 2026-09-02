Project status
==============

sedumipy ports SeDuMi's ~96 ``.m`` files and ~54 MEX ``.c`` files to a
MATLAB/Octave-free C library + Python package. As of this writing:

.. list-table::
   :header-rows: 1

   * - Phase
     - What
     - Status
   * - Phase 0
     - Octave-based golden reference harness
     - Done
   * - Phase 1
     - C kernels: drop MEX, build standalone libsedumi.so
     - Done
   * - Phase 2
     - ctypes Python bindings over libsedumi.so
     - Done
   * - Phase 3
     - Interior-point solver logic (LP + SOCP + PSD)
     - Done
   * - Phase 4
     - High-level API + .mat/SDPA I/O
     - Done
   * - Phase 5
     - Verification against published benchmarks
     - Done
   * - Phase 6
     - Packaging: Linux/macOS/Windows wheels build in CI; not yet
       published to PyPI
     - Partial

**Scope.** LP, second-order cone (SOCP, ``K.q``/``K.r``), and
semidefinite (SDP, ``K.s``) problems are fully ported and verified
bit-for-bit against real Octave/SeDuMi output on both synthetic fixtures
and published SDPLIB/DIMACS benchmark problems (:doc:`usage`'s
Benchmarks section). Dense-column preconditioning is implemented.

**Not ported** (deliberately out of scope, no effect on the returned
``(x, y, info)``): the console progress printout, ``pars.vplot``'s
v-plot, ``pars.stopat``'s interactive debug break, the optional pre-solve
rank/infeasibility diagnostic, and the DIMACS error-measures block
(``info.err``). Complex Hermitian PSD problems (``K.scomplex``/
``K.ycomplex``) are also out of scope.

**Packaging.** A wheel with ``libsedumi.so`` bundled builds and installs
correctly (verified in an isolated virtualenv with no access to the
source tree). ``cibuildwheel`` builds run in CI
(``.github/workflows/wheels.yml``) on all three platforms:

* **Linux**: manylinux containers, linked against the container's own
  ``libblas``.
* **macOS**: linked against the system Accelerate framework, no
  Homebrew/external BLAS dependency.
* **Windows**: built with an MSYS2 MinGW64 toolchain
  (``mingw-w64-x86_64-gcc``/``-openblas``) rather than MSVC --
  ``libsedumi.dll`` is a plain ctypes-loaded DLL, not a CPython
  extension, so it doesn't need to be built with the same compiler as
  Python itself. ``delvewheel repair`` bundles the resulting
  ``libopenblas.dll``/mingw-runtime dependencies into the wheel. This
  path has only been exercised on GitHub Actions' hosted Windows
  runner, not hand-verified on a real Windows machine.

Not yet done: publishing to PyPI, and manylinux-compliant BLAS bundling
(currently dynamically linked to the manylinux container's own
``libblas``, which would need ``auditwheel repair`` or static linking to
be redistributable outside that container).

For the full phase-by-phase history, the porting workflow, known bugs
found and fixed along the way, and the prioritized list of remaining
work, see :doc:`contributing`.
