"""Sphinx configuration for sedumipy's documentation.

Building these docs does not require a working `libsedumi.so` (autodoc
mocks `sedumipy._native` -- see `autodoc_mock_imports` below), so `pip
install -r docs/requirements.txt` alone is enough; no compiler/
BLAS/LAPACK toolchain or `git submodule update` is needed just to build
the docs. Published to GitHub Pages by .github/workflows/docs.yml.
"""

from __future__ import annotations

import sys
from importlib.metadata import PackageNotFoundError, version as _dist_version
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

project = "sedumipy"
copyright = "The sedumipy contributors"
author = "The sedumipy contributors"
# Read from installed package metadata rather than `import sedumipy`
# itself: sedumipy/__init__.py imports _native, which builds/loads
# libsedumi.so as a side effect of being imported (see _native.py) --
# reading metadata instead means this file never needs that toolchain
# just to know the version.
try:
    release = _dist_version("sedumipy")
except PackageNotFoundError:
    release = "0.0.0+unknown"
version = release

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "myst_parser",
]

# sedumipy/__init__.py does `from . import _native`, and _native.py builds
# (or loads) libsedumi.so as a side effect of being imported -- mocking it
# here means the docs build never needs a compiler or BLAS/LAPACK, and
# works identically on a fresh checkout, CI, or Read the Docs.
autodoc_mock_imports = ["sedumipy._native"]
autodoc_member_order = "bysource"
autosummary_generate = True
napoleon_numpy_docstring = False
napoleon_google_docstring = False

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}

html_theme = "furo"
html_static_path = ["_static"]
