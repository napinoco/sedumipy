"""Sphinx configuration for sedumipy's documentation.

Building these docs does not require a working `libsedumi.so` (autodoc
mocks `sedumipy._native` -- see `autodoc_mock_imports` below), so `pip
install -r docs/requirements.txt` alone is enough; no compiler/
BLAS/LAPACK toolchain or `git submodule update` is needed just to build
the docs. Published to GitHub Pages by .github/workflows/docs.yml.
"""

from __future__ import annotations

import re
import sys
from importlib.metadata import PackageNotFoundError, version as _dist_version
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

project = "sedumipy"
copyright = "The sedumipy contributors"
author = "The sedumipy contributors"


def _version_from_pyproject() -> str | None:
    """The version straight out of pyproject.toml, the file that
    actually defines it.

    Needed because the published docs are built *without* installing
    sedumipy: .github/workflows/docs.yml runs `pip install -r
    docs/requirements.txt` only, deliberately, since installing the
    package would trigger setup.py's libsedumi build and drag a C
    compiler plus BLAS/LAPACK into a job that just renders HTML. That
    left the metadata lookup below with nothing to find, so every
    published page was titled "sedumipy 0.0.0+unknown documentation".

    tomllib is stdlib from 3.11; on 3.10 (which this project still
    supports) fall back to a regex over the [project] table's own
    version line -- enough for one well-known field in a file we
    control, and it keeps `sphinx-build` working on every Python the
    package itself runs on.
    """
    pyproject = REPO_ROOT / "pyproject.toml"
    if not pyproject.is_file():
        return None
    try:
        import tomllib
    except ImportError:
        match = re.search(
            r'^version\s*=\s*"([^"]+)"', pyproject.read_text(encoding="utf-8"), re.MULTILINE
        )
        return match.group(1) if match else None
    with pyproject.open("rb") as fh:
        return tomllib.load(fh).get("project", {}).get("version")


# pyproject.toml first, since it is the source of truth and is always
# there in a checkout; installed metadata second, for a build from an
# installed sedumipy with no repository around it. Deliberately never
# `import sedumipy` to ask: sedumipy/__init__.py imports _native, which
# builds/loads libsedumi.so as a side effect (see _native.py), so that
# would reintroduce the very toolchain requirement this file avoids.
release = _version_from_pyproject()
if release is None:
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
