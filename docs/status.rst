Project status
==============

sedumipy ports SeDuMi's ~96 ``.m`` files and ~54 MEX ``.c`` files to a
MATLAB/Octave-free C library + Python package. As of this writing:

============ ================================================= =========
Phase        What                                              Status
============ ================================================= =========
Phase 0      Octave-based golden reference harness              Done
Phase 1      C kernels: drop MEX, build standalone libsedumi.so  Done
Phase 2      ctypes Python bindings over libsedumi.so            Done
Phase 3      Interior-point solver logic (LP + SOCP + PSD)       Done
Phase 4      High-level API + .mat/SDPA I/O                      Done
Phase 5      Verification against published benchmarks           Done
Phase 6      Packaging (wheel build validated; CI/manylinux not) Partial
============ ================================================= =========

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
correctly on this project's Linux development environment (verified in
an isolated virtualenv with no access to the source tree). Not yet
verified: ``cibuildwheel``/manylinux builds (no Docker available in the
environment this was developed in), macOS/Windows builds (the build
script assumes ``gcc``; there is no Windows equivalent yet), and
manylinux-compliant BLAS/LAPACK linking (currently dynamically linked to
the build host's ``libblas``/``libopenblas``, which would need
``auditwheel repair`` or static linking for real PyPI distribution).

For the full phase-by-phase history, the porting workflow, known bugs
found and fixed along the way, and the prioritized list of remaining
work, see :doc:`contributing`.
