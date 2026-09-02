Contributing
============

New contributor? Read `CONTRIBUTING.md
<https://github.com/napinoco/sedumipy/blob/main/CONTRIBUTING.md>`_ in the
repository root first. It has the current phase-by-phase status, the
porting workflow this project follows (how a single ``.m`` file gets
ported and verified against a real Octave oracle), known scope
limitations, coding/naming conventions, and the prioritized list of
remaining work.

It is written in Japanese; ask in an issue if an English translation
would help.

Development setup
------------------

.. code-block:: sh

   git clone --recurse-submodules https://github.com/napinoco/sedumipy.git
   cd sedumipy
   python -m venv .venv
   .venv/bin/pip install -e .[test]
   .venv/bin/python -m pytest tests/ -q

Regenerating oracle fixtures (only needed when changing what a module is
checked against, not for routine development) requires a real Octave
build of upstream SeDuMi under ``vendor/sedumi-upstream/`` -- see
``CONTRIBUTING.md`` section 8.

Building the documentation
----------------------------

Building the docs alone does not require a C compiler or BLAS/LAPACK,
``libsedumi.so``, or a full ``sedumipy`` install (``docs/conf.py`` puts
``src/`` on ``sys.path`` directly and mocks the ``sedumipy._native``
import); ``docs/requirements.txt`` lists exactly what building the docs
needs:

.. code-block:: sh

   .venv/bin/pip install -r docs/requirements.txt
   .venv/bin/sphinx-build -b html docs docs/_build/html

``pip install -e .[docs]`` also works if you already have a full dev
install (:doc:`installation`), and additionally makes the ``docs`` extra
available for other tooling that expects it.
