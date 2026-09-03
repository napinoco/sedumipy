Installation
============

Prebuilt wheels
---------------

Where a wheel matches your platform, nothing below is needed: the
compiled kernel library and the BLAS it needs are already inside it, so
no compiler and no BLAS install are required.

Wheels are published for CPython 3.10-3.13 on:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Platform
     - BLAS carried by the wheel
   * - Linux x86_64 (manylinux)
     - OpenBLAS, vendored by ``auditwheel``
   * - Windows x64
     - OpenBLAS, vendored by ``delvewheel``
   * - macOS (Apple silicon)
     - none needed -- the system Accelerate framework

Anything else -- 32-bit or ARM Windows, Alpine/musl, Linux ``aarch64``,
Intel macOS -- has no wheel, so ``pip`` falls back to the source
distribution and compiles on your machine, which needs the toolchain in
`Requirements`_ below. On Windows that means MSYS2 specifically, which
is a real obstacle rather than a formality; ``win_arm64`` is skipped
because MSYS2's OpenBLAS package is x86_64-only.

Requirements
------------

These apply when building from source -- either from a checkout, or
because ``pip`` fell back to the source distribution as described above.

* Python >= 3.10
* A C compiler and BLAS, to build ``libsedumi.so``/``.dylib``/``.dll``
  (the compiled kernel library) the first time ``sedumipy`` is imported:

  * **Linux**: ``gcc`` plus a BLAS development package (e.g. ``apt
    install build-essential libopenblas-dev``). Either OpenBLAS or the
    reference Netlib BLAS works; the build links OpenBLAS when it can
    find it, since it runs the kernels this library calls roughly 3x
    faster, and falls back to ``-lblas`` otherwise.
  * **macOS**: nothing extra -- links against the system Accelerate
    framework.
  * **Windows**: `MSYS2 <https://www.msys2.org/>`_ with the MINGW64
    ``mingw-w64-x86_64-gcc`` and ``mingw-w64-x86_64-openblas`` packages
    installed (``pacman -S ...``), with ``C:\msys64\mingw64\bin`` and
    ``C:\msys64\usr\bin`` on ``PATH``. Not MSVC -- ``libsedumi.dll`` is
    a plain ctypes-loaded DLL, not a CPython extension, so it doesn't
    need to match whatever compiler built Python itself. This path is
    exercised in CI (see :doc:`status`) but has not been hand-verified
    on a real Windows machine.

From source
-----------

.. code-block:: sh

   git clone --recurse-submodules https://github.com/napinoco/sedumipy.git
   cd sedumipy
   python -m venv .venv
   .venv/bin/pip install -e .[test]

If you already cloned without ``--recurse-submodules``, run
``git submodule update --init --recursive`` first -- this pulls in the
`SDPLIB <https://github.com/vsdp/SDPLIB>`_ and
`DIMACS <https://github.com/vsdp/DIMACS>`_ benchmark submodules used by
the test suite (``examples/``) and the original MATLAB/MEX
implementation this project ports from, kept for reference
(``vendor/sedumi-upstream/``). Neither submodule is required just to
``import sedumipy`` and solve a problem.

The compiled kernel library is built automatically the first time
``sedumipy`` is imported (via ``tools/build_libsedumi.sh``, invoked
through ``bash`` on Windows); no separate build step is needed. It is
not committed to the repository (see ``.gitignore``), so a fresh
checkout always builds its own.

Verifying the install
----------------------

.. code-block:: sh

   .venv/bin/python -m pytest tests/ -q

Octave is *not* required to run the existing test suite: the Octave/MEX
oracle data the tests check against is committed as ``.mat`` fixtures. It
is only needed to *regenerate* that oracle data (see :doc:`contributing`).
