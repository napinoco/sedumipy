Installation
============

Requirements
------------

* Python >= 3.10
* A C compiler and BLAS, to build ``libsedumi.so``/``.dylib``/``.dll``
  (the compiled kernel library) the first time ``sedumipy`` is imported:

  * **Linux**: ``gcc`` plus a BLAS development package (e.g. ``apt
    install build-essential libblas-dev``).
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
