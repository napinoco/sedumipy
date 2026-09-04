Contributing
============

New contributor? Read `CONTRIBUTING.md
<https://github.com/napinoco/sedumipy/blob/main/CONTRIBUTING.md>`_ in the
repository root first. It has the current phase-by-phase status, the
porting workflow this project follows (how a single ``.m`` file gets
ported and verified against a real Octave oracle), known scope
limitations, and coding/naming conventions. The detailed, session-by-
session history of bugs found and fixed during porting -- and how
packaging/benchmarking/performance work actually got done -- lives
separately in `DEVLOG.md
<https://github.com/napinoco/sedumipy/blob/main/DEVLOG.md>`_.

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
``CONTRIBUTING.md`` section 6.

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

Docs are published to GitHub Pages automatically on every push to
``main`` that touches ``docs/``, ``src/sedumipy/``, or
``pyproject.toml`` (see ``.github/workflows/docs.yml``); there's no
manual publish step.

A handful of ported docstrings use ``*``/``'`` in ways that read as
unterminated RST emphasis to docutils (e.g. a bare apostrophe in
"A*x=b"-style math or a MATLAB-style ``K.q`` cross-reference). These
show up as harmless ``WARNING: Inline emphasis/strong start-string
without end-string`` lines during the build; the docs build (both
locally and in CI) does not use ``-W``, so they don't fail it. Fixing a
docstring's RST escaping is welcome but optional -- don't reflow a
carefully-written docstring's wording just to silence one.
