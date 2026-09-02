Installation
============

Requirements
------------

* Python >= 3.10
* A C compiler (``gcc`` on Linux/macOS) and BLAS/LAPACK development
  headers, to build ``libsedumi.so``/``.dylib`` (the compiled kernel
  library) the first time ``sedumipy`` is imported. Windows is not yet
  supported -- see :doc:`status`.

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

``libsedumi.so`` is built automatically the first time ``sedumipy`` is
imported (via ``tools/build_libsedumi.sh``); no separate build step is
needed. This compiled library is not committed to the repository (see
``.gitignore``), so a fresh checkout always builds its own.

Verifying the install
----------------------

.. code-block:: sh

   .venv/bin/python -m pytest tests/ -q

Octave is *not* required to run the existing test suite: the Octave/MEX
oracle data the tests check against is committed as ``.mat`` fixtures. It
is only needed to *regenerate* that oracle data (see :doc:`contributing`).
