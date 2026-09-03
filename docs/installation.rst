Installation
============

Prebuilt wheels
---------------

Where a wheel matches your platform, nothing below is needed: the
compiled kernel library and the BLAS it needs are already inside it, so
no compiler and no BLAS install are required.

Wheels are published for CPython 3.10-3.13 on Linux x86_64 (manylinux),
Windows x64, and macOS (Apple silicon). All three carry
`scipy-openblas64 <https://pypi.org/project/scipy-openblas64/>`_ (a
prebuilt, ILP64 OpenBLAS build) as their BLAS, vendored into the wheel
at build time (``auditwheel``/``delocate``/``delvewheel`` -- see
`Requirements`_ below for what building it yourself needs).

Anything else -- 32-bit or ARM Windows, Alpine/musl, Linux ``aarch64``,
Intel macOS -- has no wheel, so ``pip`` falls back to the source
distribution and compiles on your machine, which needs the toolchain in
`Requirements`_ below. On Windows that means MSYS2 specifically, which
is a real obstacle rather than a formality.

Requirements
------------

These apply when building from source -- either from a checkout, or
because ``pip`` fell back to the source distribution as described above.

* Python >= 3.10
* A C compiler, to build ``libsedumi.so``/``.dylib``/``.dll`` (the
  compiled kernel library) the first time ``sedumipy`` is imported:

  * **Linux**: ``gcc`` (e.g. ``apt install build-essential``).
  * **macOS**: nothing extra -- Xcode's command line tools provide ``cc``.
  * **Windows**: `MSYS2 <https://www.msys2.org/>`_ with the MINGW64
    ``mingw-w64-x86_64-gcc`` package installed (``pacman -S
    mingw-w64-x86_64-gcc``), with ``C:\msys64\mingw64\bin`` and
    ``C:\msys64\usr\bin`` on ``PATH``. Not MSVC -- ``libsedumi.dll`` is
    a plain ctypes-loaded DLL, not a CPython extension, so it doesn't
    need to match whatever compiler built Python itself. This path is
    exercised in CI (see :doc:`status`) but has not been hand-verified
    on a real Windows machine.

* A BLAS, either of:

  * **Recommended, same on every OS**: ``pip install
    scipy-openblas64==0.3.34.106.0`` into the environment you're
    installing sedumipy into, *before* running ``pip install -e .`` --
    and pass ``--no-build-isolation`` to that install, since
    scipy-openblas64 isn't declared as one of this project's own build
    requirements (see ``pyproject.toml``'s ``[tool.cibuildwheel]``
    comment for why) and a normal isolated build wouldn't see it
    otherwise. This is what the published wheels themselves link, so it
    exercises the exact same code path as a wheel install and needs
    nothing beyond ``pip``, on Linux, macOS or Windows alike. Pinned to
    an exact version deliberately -- see ``pyproject.toml``'s comment on
    why and where else to update it if you ever bump it.
  * **Fallback, no network access needed** (used automatically when
    scipy-openblas64 isn't importable by the Python running the build):
    a system BLAS per OS, exactly as before -- **Linux**: a BLAS
    development package (e.g. ``apt install libopenblas-dev``; either
    OpenBLAS or the reference Netlib BLAS works, OpenBLAS preferred when
    both are present since it runs the kernels this library calls
    roughly 3x faster). **macOS**: nothing extra -- the system
    Accelerate framework. **Windows**: also install
    ``mingw-w64-x86_64-openblas`` alongside the compiler above.

  See ``tools/build_libsedumi.sh`` for exactly how the choice between
  the two is made.

From source
-----------

.. code-block:: sh

   git clone --recurse-submodules https://github.com/napinoco/sedumipy.git
   cd sedumipy
   python -m venv .venv
   .venv/bin/pip install -e .[test]

This uses whatever system BLAS you have, per `Requirements`_ above. To
build against scipy-openblas64 instead (recommended -- see
`Requirements`_ for why), install it first and pass
``--no-build-isolation``:

.. code-block:: sh

   .venv/bin/pip install scipy-openblas64==0.3.34.106.0
   .venv/bin/pip install -e .[test] --no-build-isolation

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
