# sedumipy — MATLAB/Octave-free port of SeDuMi

`sedumipy` is an in-progress port of [SeDuMi](https://github.com/sqlp/sedumi)
(a linear/quadratic/semidefinite programming solver originally written for
MATLAB/Octave) to a standalone C library + Python (NumPy/SciPy) package,
with no MATLAB or GNU Octave runtime dependency.

**New contributor?** Read [`CONTRIBUTING.md`](CONTRIBUTING.md) first — it
has the current phase-by-phase status, the porting workflow this project
follows, known scope limitations, and the prioritized list of remaining
work. (It's in Japanese; ask if you'd like an English translation.)

## Repository layout

The original MATLAB/Octave/MEX implementation this project ports from is
kept as a reference-only git submodule, pinned to the commit the port
started from:

```
sedumipy/
  vendor/sedumi-upstream/   # submodule: sqlp/sedumi (reference only, not built by default)
  examples/sdplib/          # submodule: vsdp/SDPLIB (published benchmark problems + optimal values)
  examples/dimacs/          # submodule: vsdp/DIMACS (published benchmark problems + optimal values)
  csrc/                     # forked, MEX-free standalone C kernels (source for libsedumi.so)
  src/sedumipy/             # the Python package
  tests/                    # test suite + committed Octave-generated oracle fixtures
  tools/                    # libsedumi build script + oracle/golden-reference generators
```

## Getting started

```sh
git clone --recurse-submodules <this-repo-url>
cd sedumipy
python -m venv .venv
.venv/bin/pip install -e .[test]
.venv/bin/python -m pytest tests/ -q
```

If you already cloned without `--recurse-submodules`, run
`git submodule update --init --recursive` first. `libsedumi.so` (the
compiled C kernel library) is built automatically the first time
`sedumipy` is imported, via `tools/build_libsedumi.sh`; the Octave
submodule is only needed to regenerate oracle/golden-reference data, not
to run the existing test suite.

## Status

LP and SOCP (second-order cone) problems are fully ported and verified
bit-for-bit against real Octave/SeDuMi output. PSD cones (`K.s`) and the
dense-columns optimization are not yet implemented — see
[`CONTRIBUTING.md`](CONTRIBUTING.md) for details and the current
priority list.

## Benchmarks

`tests/test_benchmarks.py` solves published [SDPLIB](https://github.com/vsdp/SDPLIB)
and [DIMACS](https://github.com/vsdp/DIMACS) problems (added as git
submodules under `examples/`) and checks the result against each
collection's own official optimal-value table -- a correctness/regression
check against real reference numbers, not a synthetic self-check. It also
prints a timing/iteration-count summary and writes it to
`benchmark_results.csv`.

Coverage is essentially the full published set: all 92 SDPLIB problems
and every DIMACS problem sedumi() can read directly, minus a documented
handful this port can't yet solve in bounded time/memory or that hit a
real solver limitation -- see the module's own docstring for the full
exclusion list and why each one is excluded.

```sh
git submodule update --init --recursive   # if not already done
.venv/bin/python -m pytest tests/test_benchmarks.py -v          # everything, ~101 problems (~10 min)
.venv/bin/python -m pytest tests/test_benchmarks.py -v -m mini  # fastest subset only (~35s)
```

Problems are grouped by `pytest.mark.mini` (<2s each) / `timing` (2-20s)
/ `extended` (20s-130s) by measured solve time -- see the module's own
docstring for the sign conventions each collection/family needs.

## License

SeDuMi is licensed under the GNU General Public License v2 (see
[`LICENSE`](LICENSE)); this port, being a derivative work, is licensed
the same way.
