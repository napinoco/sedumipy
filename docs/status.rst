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
against real Octave/SeDuMi output on both synthetic fixtures
and published SDPLIB/DIMACS benchmark problems (:doc:`usage`'s
Benchmarks section). Dense-column preconditioning is implemented.

The test suite compares against Octave-generated oracle fixtures with
``numpy.testing.assert_allclose`` -- ``rtol`` between 1e-12 and 1e-10 for
the individual C kernels, and 1e-5 for whole end-to-end solves, where the
interior-point iteration accumulates rounding. Integer-valued results
(orderings, permutations, sparsity patterns) are compared exactly. Several
individual kernels were confirmed bit-identical to the MEX build during
development -- see ``CONTRIBUTING.md`` -- but that is a development
observation, not what CI asserts: the tolerances above are what has to
hold across the four BLAS implementations CI covers (reference Netlib,
OpenBLAS, and scipy-openblas64 -- a pip-installable, prebuilt ILP64
OpenBLAS, see below -- on Linux; OpenBLAS on Windows; Accelerate on
macOS).

**Not ported** (deliberately out of scope, no effect on the returned
``(x, y, info)``): the console progress printout, ``pars.vplot``'s
v-plot, ``pars.stopat``'s interactive debug break, the optional pre-solve
rank/infeasibility diagnostic, and the DIMACS error-measures block
(``info.err``). Complex Hermitian PSD problems (``K.scomplex``/
``K.ycomplex``) are also out of scope.

**Packaging.** A wheel with ``libsedumi.so`` bundled builds and installs
correctly (verified in an isolated virtualenv with no access to the
source tree). ``cibuildwheel`` builds run in CI
(``.github/workflows/wheels.yml``) on all three platforms. Linux and
Windows link `scipy-openblas64
<https://pypi.org/project/scipy-openblas64/>`_ -- a pip-installable,
prebuilt ILP64 OpenBLAS, pinned to an exact version, the same package
numpy/scipy themselves build against -- as their BLAS; macOS keeps
linking the system Accelerate framework unconditionally, since it never
had the other two's build/install problem for scipy-openblas64 to solve
(see ``tools/build_libsedumi.sh``'s Darwin case for the full reasoning):

* **Linux**: manylinux containers; ``auditwheel`` vendors the resulting
  ``libscipy_openblas64_.so`` into the wheel.
* **macOS**: no vendoring needed -- Accelerate is a system framework
  present on every Mac, not a bundled shared library, and no Homebrew
  dependency either.
* **Windows**: built with an MSYS2 MinGW64 toolchain
  (``mingw-w64-x86_64-gcc``) rather than MSVC -- ``libsedumi.dll`` is a
  plain ctypes-loaded DLL, not a CPython extension, so it doesn't need
  to be built with the same compiler as Python itself.
  ``tools/repair_windows_wheel.py`` (``delvewheel repair
  --analyze-existing`` under the hood) bundles the resulting
  ``libscipy_openblas64_.dll`` and mingw-runtime dependencies into the
  wheel. This path has only been exercised on GitHub Actions' hosted
  Windows runner, not hand-verified on a real Windows machine.

Not yet done: publishing to PyPI.

For the full phase-by-phase history, the porting workflow, known bugs
found and fixed along the way, and the prioritized list of remaining
work, see :doc:`contributing`.
