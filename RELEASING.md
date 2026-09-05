# Releasing to PyPI

`sedumipy` is published to [PyPI](https://pypi.org/project/sedumipy/) by
[`.github/workflows/wheels.yml`](.github/workflows/wheels.yml), which
already builds the artifacts a release consists of: CPython 3.10-3.13
wheels for Linux x86_64 (manylinux), Windows x64 and macOS (Apple
silicon), plus a source distribution for every platform without one.

Publishing uses PyPI's **Trusted Publishing** (OIDC). No API token is
stored in this repository, in a secret or anywhere else: the upload
credential is minted per workflow run, is valid for minutes, and is
bound to this repository, this workflow file and the GitHub environment
the job declares. That binding is the security property, so the
one-time setup below has to name all four.

## One-time setup (per index)

Do this once on TestPyPI and once on PyPI. Both accept a *pending*
publisher, i.e. one created before the project exists on that index --
which is the case here, since nothing has been uploaded yet.

1. Sign in, go to **Your account -> Publishing -> Add a new pending
   publisher**, choose GitHub, and enter:

   | Field | Value |
   | --- | --- |
   | PyPI project name | `sedumipy` |
   | Owner | `napinoco` |
   | Repository name | `sedumipy` |
   | Workflow name | `wheels.yml` |
   | Environment name | `testpypi` (on TestPyPI) / `pypi` (on PyPI) |

   The environment name is not optional here: both publish jobs declare
   an `environment:`, and an upload from a run that did not go through
   that environment will be rejected.

2. Nothing to do on the GitHub side. Both environments are created
   automatically the first time a job referencing them runs. Optionally
   add a required reviewer to the `pypi` environment
   (**Settings -> Environments -> pypi**) if you want a human approval
   step in front of every real upload.

## Rehearsal: TestPyPI

Always do this first for a version number you have not published before.
PyPI accepts a given version exactly once, forever -- a mistake caught
here is free, the same mistake on PyPI costs a version number.

1. **Actions -> wheels -> Run workflow**, set **publish** to `testpypi`,
   run it from `main`.
2. Wait for `publish-testpypi` to go green, then install what it
   uploaded, from a clean virtualenv, with the real dependencies coming
   from real PyPI:

   ```sh
   python -m venv /tmp/tp && /tmp/tp/bin/pip install \
     --index-url https://test.pypi.org/simple/ \
     --extra-index-url https://pypi.org/simple/ \
     sedumipy
   cd /tmp && /tmp/tp/bin/python -c "
   import sedumipy, numpy as np
   x, y, info = sedumipy.sedumi(np.eye(2), np.array([1.0, 1.0]),
                                np.array([1.0, 1.0]), {'l': 2})
   print(sedumipy.__version__, info['numerr'], x)
   "
   ```

   Check that pip took a *wheel* rather than falling back to the sdist
   (the install log names the file), and read the rendered project page
   at <https://test.pypi.org/project/sedumipy/> -- this is the last
   chance to fix how the README, metadata and links look.

## The real release

1. Bump `version` in `pyproject.toml`, move the `CHANGELOG.md`
   `[Unreleased]` entries under the new version heading, and merge that
   to `main`.
2. Confirm `ci` and `wheels` are green on the merge commit.
3. **Releases -> Draft a new release**: create the tag `vX.Y.Z` on
   `main`, title it `vX.Y.Z`, paste the changelog section as the notes,
   and **Publish release**.
4. That fires `wheels.yml` a final time on the `release: published`
   event, which rebuilds everything from the tag and runs both
   `publish-to-release` (attaches the files to the GitHub release) and
   `publish-pypi` (uploads them to PyPI).
5. Verify from a clean virtualenv: `pip install sedumipy`, then the same
   smoke test as above.

## If an upload goes wrong

A published version cannot be replaced -- deleting a release from PyPI
frees neither the version number nor the filenames. Publish `X.Y.Z+1`
with the fix instead. `yank` (on the PyPI project page) is the right
tool for a release that is actively broken: it stays installable for
anyone who pinned it exactly, but pip stops resolving to it.
