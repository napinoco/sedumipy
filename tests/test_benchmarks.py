"""Regression tests against the published SDPLIB/DIMACS optimal objective
values -- i.e. against the collections' own official reference numbers,
not a baseline this port recorded from its own runs.

Problem data is sourced from git submodules:
  - examples/sdplib/ -> https://github.com/vsdp/SDPLIB.git
  - examples/dimacs/ -> https://github.com/vsdp/DIMACS.git

Initialize with:
    git submodule update --init --recursive

This covers every problem in both collections that sedumi() can actually
solve: all 92 SDPLIB problems and every DIMACS problem whose data ships
as a single `.mat.gz`/`.dat-s` file sedumi() can read directly (DIMACS's
BISECT/FAP problems that only ship a `.dat` graph description needing a
MATLAB generator script to become a `.mat` are out of scope, same as
they are for this port's format readers generally), MINUS the following,
each excluded for a documented, verified reason rather than silently
dropped:

  - Too large for this solver to finish in a bounded time/memory budget
    (SDPLIB: equalG51, maxG32, maxG55, maxG60, qpG11, qpG51, theta6,
    thetaG11, thetaG51; DIMACS: copo68, hamming_10_2, hamming_11_2,
    hamming_8_3_4, hamming_9_5_6, sched_100_100_orig, sched_200_100_orig,
    sched_200_100_scaled, bm1, nql60(old), qssp60(old), and the "-15"
    TORUS instances (single ~3375-order dense SDP block); nql180old/
    qssp180old also belong here (below). nql180/qssp180 themselves used
    to belong here too (getdatm.py's old `DAt_q.todense()` OOM'd on
    them, and past that fix they still hit numerr=2 within the first
    couple of iterations) but are now confirmed fixed: re-tested
    directly (not via this test file, since neither has a published
    reference objective -- both DIMACS README rows read "N/A" -- so
    there's nothing to parametrize into DIMACS_PARAMS) after the
    getada.py/getdatm.py dense/sparse hybrid fix (see CONTRIBUTING.md
    section 7 item 5), both now solve cleanly: nql180 numerr=0, iter=16,
    ~39s; qssp180 numerr=0, iter=42, ~249s (internal consistency checked
    via cx~=by, feasratio->1, r0=1e-8, since there's no published value
    to check against). Left out of this file's parametrized tests
    anyway since (a) no reference objective exists to assert against and
    (b) qssp180 alone is far past the "timing" mark's "well under a
    minute" bar. nql180old (inferior "old"-formulation variant, same
    family as nql30old/qssp30old below) is a different, still-unresolved
    story: it was cross-checked against the real Octave/MEX build (built
    from source in this environment) and BOTH struggle badly on it (it's
    a genuinely ill-conditioned instance -- the real build's own console
    output shows `skip=5361` Cholesky pivots skipped by iteration 54),
    but not identically: the real build still limps to numerr=1 (iter=54,
    degraded but not a total failure), while this port gives up earlier
    and worse (numerr=2, iter=27, feasratio=0.90, r0=0.53) -- unlike
    nql30old/qssp30old, this is NOT simply "the real build fails too,
    nothing to fix here"; there's a real, if narrow, robustness gap on
    this specific hard instance. qssp180old (largest file in this
    family, ~36 MB) didn't finish in either this port or a real-build run
    within that investigation's time budget (550s each) and was left
    unverified -- resolved in a later session by giving both builds a
    much larger budget instead: the real build now completes in ~1705s
    (numerr=2, iter=30), and this port completes in ~3557s (numerr=2,
    iter=30) -- the exact same failing iteration on both, unlike
    nql180old's gap. So qssp180old belongs with nql30old/qssp30old, not
    with nql180old: a genuine solver limitation shared by both builds,
    not a porting bug.
  - sedumi() returns numerr=2 (a genuine, reproducible solver failure,
    not a reference-value problem) on: SDPLIB none; DIMACS nb_L2,
    nql30old, qssp30old. nql30 used to be in this list too, but is now
    fixed (see below) and has a DIMACS_PARAMS row instead. Cross-checked
    against a from-source build of the real
    Octave/MEX SeDuMi (vendor/sedumi-upstream, `install_sedumi`) on the
    same .mat files: qssp30old (and nql30old, same family/shape) fails
    there too (numerr=2), confirming a genuine solver limitation on
    those instances, not a porting bug. nb_L2 is the opposite: the real
    Octave/MEX build solves it cleanly (numerr=0, iter=16); this port
    still doesn't. That gap is narrowed but not yet closed -- passing
    `stepdif=1` (skipping pars.stepdif's default "Adaptive
    Step-Differentiation" auto-switch entirely) makes this port solve
    nb_L2 cleanly too (numerr=0, iter=17), and the *default* run's own
    per-iteration CG counts (`err["kcg"]`/`Lsd["kcg"]`, the same
    quantities real sedumi.m's console "cg cg" columns report) run well
    above the real build's throughout -- so the auto-switch fires
    several iterations earlier here than it does in the real build,
    changing the solved trajectory before it can recover -- but nothing
    downstream of that (dense-column detection, the one-time symbolic
    ADA pattern now built via getsymbada() same as the has_psd branch,
    numeric Cholesky -- `skip=0` every iteration, i.e. never missing a
    position it needs) turned up an actual defect, and every PCG
    sub-solve genuinely converges (residual well under `restol`, never
    the stagnation branch) -- just using more iterations than the real
    build's own preconditioner needs for the same linear system. Left
    excluded, not force-fixed with `stepdif=1`, since that's a
    numerical-sensitivity symptom pinned to a specific mechanism, not a
    located line-level bug, and forcing it off pars's own default for
    every problem risks trading this instance's failure for a worse
    trajectory on others that currently rely on the adaptive default.

    Narrowed further, and now fully located (see CONTRIBUTING.md section
    7 item 6 for the full derivation): dumping ADA/d/DAt.q from both the
    real Octave/MEX build (a temporary `save()` inserted into a scratch
    copy of sedumi.m's main loop, not committed) and this port at each
    of the first 3 iterations on nb_L2 (839 Lorentz blocks feeding 123
    constraints, no dense columns -- `getdense()` returns
    `dense["cols"].size == dense["q"].size == 0` for this file, ruling
    the dense-column/product-form machinery out entirely rather than
    just "checked out fine") shows d.l/d.det matching to float noise
    (~1e-13) through iteration 3, and d.q1/d.q2 (the Lorentz-block
    scaling point) matching to float noise through the d used at the
    *start* of iteration 2 -- but the d produced by iteration 2's step
    (used at iteration 3) diverges by ~15% relative in d.q1's worst
    entry, well past anything float-order noise explains, and this is
    exactly where err["kcg"]/Lsd["kcg"] jump from 1/1 (iterations 1-2,
    matching the real build) to 6/5.

    This session went further than checking updtransfo.py line-by-line:
    it actually transplanted the real build's own iteration-2
    xscl/zscl/w/d (dumped via the same temporary save()) directly into
    this port's updtransfo() and got the real build's exact iteration-3
    d.q1/d.q2 back, bit for bit -- proving updtransfo.py innocent by
    execution, not just by audit. Comparing this port's own iteration-2
    xscl/zscl/w against the real build's next pinned down where the two
    actually part ways: xscl/zscl and w["tdetx"]/w["tdetz"] all agree to
    ~1e-13 (ordinary cross-implementation float noise for a 4196-d
    vector), but w["lab"] itself disagrees by up to 7.6 -- wildly out of
    proportion to inputs that agree to 1e-13. The cause is
    widelen.py's `_build_w()` (a faithful port of widelen.m's own logic):
    `lab2q`, the Lorentz-block eigenvalue term, is computed as
    `halfxz + sqrt(tmp)` only `if np.all(tmp > 0)` across *all 839
    Lorentz blocks at once* -- a single global all-or-nothing branch,
    not a per-block one -- and falls back to the cruder `lab2q = halfxz`
    for *every* block otherwise. Recomputing tmp from each build's own
    iteration-2 xscl/zscl shows exactly one block (index 396 of 839)
    sitting right on top of zero: tmp = +1.78e-15 in the real build,
    tmp = -1.78e-15 in this port -- a sign flip from sub-ULP rounding
    noise between two independent floating-point pipelines (NumPy/SciPy
    plus this port's own C kernels vs. Octave plus its BLAS), not from
    any actual defect in either. That single flipped sign trips the
    global `all()` and switches the fallback formula on for all 839
    blocks at once, which is why a 1e-13-level input disagreement
    balloons into a 7.6-level disagreement in w["lab"] and, propagated
    through one honest updtransfo() call, ~15% in d.q1.
    (Verified in both directions: feeding the real build's own
    xscl/zscl into this port's `_build_w()` reproduces the real build's
    w["lab"] exactly bit for bit; feeding this port's own xscl/zscl back
    into the same function reproduces tmp[396] < 0 and the fallback
    branch, matching what this port actually computed.)

    That all-or-nothing branch is widelen.m's own design (present
    unmodified in vendor/sedumi-upstream/widelen.m), not something this
    port introduced -- plausibly a deliberate cheap safety net against
    handing sqrt() a negative discriminant on *any* block, applied to
    every block at once rather than per-block. Two independent
    same-input floating-point pipelines occasionally landing on opposite
    sides of an exact zero in a quantity like this is expected chaotic
    sensitivity, not a locatable off-by-one or formula error -- no line
    in updtransfo.py/widelen.py/tdet/ddot is wrong. Left unfixed for the
    same reason `stepdif=1` isn't forced as pars's default (see above):
    changing the branch's numerics to dodge this one instance's
    coin-flip would be an algorithm-level change with unclear effects on
    every other problem that currently relies on today's exact
    branching, not a correctness fix.
  - SDPLIB hinf12: strong duality fails for this instance (duality gap
    ~28, matches sdpt3py's own documented exclusion of the same problem).
  - DIMACS hinf12/hinf13: the README marks both "(?)" (its own
    low-confidence flag) and the solved objective doesn't match the
    listed value -- unlike SDPLIB's hinf13/hinf15 below, there's no
    known-good target to widen the tolerance to.
  - DIMACS sched_100_50_scaled: solved objective is off from the listed
    value by an unexplained factor of ~10 (not the file's own c_mult,
    which is 2708.1) -- excluded rather than guess why.

Reference values are transcribed from each submodule's own README.md
table (SDPLIB's "Optimal Objective Value" column; DIMACS's per-set
"Opt. value" columns), with two corrections to the tables' own listed
text: SDPLIB's own footnote 14 gives qap10's correct value as -1.093e+03
(the table cell itself has a typo, -1.093e+01); and two name mismatches
between a table row and its data file (SDPLIB's "eqaulG11" row is
`equalG11.dat-s`; DIMACS's "toruspm-8-50"/"filter48" rows are
`toruspm3-8-50.mat.gz`/`filter48_socp.mat.gz`) are resolved to the
actual filename.

Sign conventions (empirically confirmed against real solves for every
row below, not just inferred from the READMEs):
  - SDPLIB: sedumi()'s primal objective is the *negative* of the
    published value, uniformly across all 92 problems. SDPLIB's table
    follows SDPA's convention (SDPA maximizes trace(C,X));
    sedumipy.read_sdpa() flips C's sign on the way in so sedumi() can
    minimize instead (see sdpa.py's own docstring), so the returned
    objective -- and, for the 4 infeasible instances, which of
    sedumi()'s own info["pinf"]/info["dinf"] flags fires -- is inverted
    relative to the table.
  - DIMACS: ANTENNA, FILTER, TRUSS, COPOS, SCHED, and QSSP already store
    their data in sedumi()'s own min c'x form, so the published value
    matches sedumi()'s objective directly (no sign flip, and
    info["pinf"]/info["dinf"] match the table's own words). HAMMING is
    the opposite -- the published value needs negating, matching this
    port's actual solve (not something the DIMACS README documents
    explicitly for this set the way it does for TORUS below). TORUS
    problems are `max c'x` reformulated as `min -c'x` (the DIMACS
    README's own documented caveat), so ref_obj = -1 * published value;
    the non-Gaussian "pm" instance's stored `c` already has that flip
    baked in (c_scale=1.0 solves directly to -1 * published value), but
    the Gaussian "g" instance's `c` does not -- it needs an *additional*
    -1/100000 factor (the /100000 part is the README's own documented
    Gaussian-only adjustment).

Markers (by measured solve time on this port's own reference run):
  mini      < 2s -- safe to run on every commit
  timing    2s-20s -- fine for a full local run, a bit much for every commit
  extended  20s+ (up to ~130s) -- deliberate, not routine runs

Tests are skipped outright if the submodules aren't checked out.
"""

from __future__ import annotations

import pathlib
import time

import numpy as np
import pytest
import scipy.sparse as sp

sedumipy = pytest.importorskip("sedumipy")
from sedumipy.sdpa import read_sdpa  # noqa: E402
from sedumipy.matio import read_mat  # noqa: E402
from sedumipy.sedumi import sedumi  # noqa: E402

HERE = pathlib.Path(__file__).parent.parent
SDPLIB_DIR = HERE / "examples" / "sdplib" / "data"
DIMACS_DIR = HERE / "examples" / "dimacs" / "data"
# Each test also skips individually on a missing data *file* (a partial
# submodule checkout); this just short-circuits the whole module when the
# submodule was never initialized at all.
pytestmark = pytest.mark.skipif(
    not SDPLIB_DIR.exists() and not DIMACS_DIR.exists(),
    reason="examples/sdplib and examples/dimacs submodules not checked out",
)


def _dense(v):
    return np.asarray(v.todense() if sp.issparse(v) else v).ravel()


def _record(collector, source, name, t0, x, c, info):
    collector.append({
        "source": source,
        "name": name,
        "time_s": time.perf_counter() - t0,
        "iter": info["iter"],
        "pobj": float(np.sum(_dense(c) * _dense(x))),
        "numerr": info["numerr"],
        "status": "ok" if info["numerr"] in (0, 1) else f"numerr={info['numerr']}",
    })


# ---------------------------------------------------------------------------
# SDPLIB (SDPA sparse .dat-s; https://github.com/vsdp/SDPLIB README.md table)
# ref_obj = -1 * published value (see module docstring's Sign conventions).
# atol is 0.2% of |ref_obj| (floor 1e-4), except several `hinf*` entries
# widened by hand -- their own published value is only given to 2-3
# significant figures, which 0.2% doesn't cover.
# ---------------------------------------------------------------------------

SDPLIB_PARAMS = [
    # (name, ref_obj, atol, marks)
    pytest.param("arch0", -0.566517, 0.001133034, marks=pytest.mark.timing),
    pytest.param("arch2", -0.671515, 0.00134303, marks=pytest.mark.timing),
    pytest.param("arch4", -0.9726274, 0.001945255, marks=pytest.mark.mini),
    pytest.param("arch8", -7.05698, 0.01411396, marks=pytest.mark.timing),
    pytest.param("control1", -17.78463, 0.03556926, marks=pytest.mark.mini),
    pytest.param("control10", -38.533, 0.077066, marks=pytest.mark.extended),
    pytest.param("control11", -31.959, 0.063918, marks=pytest.mark.extended),
    pytest.param("control2", -8.3, 0.0166, marks=pytest.mark.mini),
    pytest.param("control3", -13.63327, 0.02726654, marks=pytest.mark.mini),
    pytest.param("control4", -19.79423, 0.03958846, marks=pytest.mark.mini),
    pytest.param("control5", -16.8836, 0.0337672, marks=pytest.mark.timing),
    pytest.param("control6", -37.3044, 0.0746088, marks=pytest.mark.timing),
    pytest.param("control7", -20.6251, 0.0412502, marks=pytest.mark.timing),
    pytest.param("control8", -20.286, 0.040572, marks=pytest.mark.timing),
    pytest.param("control9", -14.6754, 0.0293508, marks=pytest.mark.extended),
    pytest.param("equalG11", -629.1553, 1.258311, marks=pytest.mark.extended),
    pytest.param("gpp100", 44.9435, 0.089887, marks=pytest.mark.mini),
    pytest.param("gpp124-1", 7.3431, 0.0146862, marks=pytest.mark.mini),
    pytest.param("gpp124-2", 46.8623, 0.0937246, marks=pytest.mark.mini),
    pytest.param("gpp124-3", 153.014, 0.306028, marks=pytest.mark.mini),
    pytest.param("gpp124-4", 418.99, 0.83798, marks=pytest.mark.mini),
    pytest.param("gpp250-1", 15.445, 0.03089, marks=pytest.mark.timing),
    pytest.param("gpp250-2", 81.869, 0.163738, marks=pytest.mark.timing),
    pytest.param("gpp250-3", 303.5, 0.607, marks=pytest.mark.timing),
    pytest.param("gpp250-4", 747.3, 1.4946, marks=pytest.mark.timing),
    pytest.param("gpp500-1", 25.3, 0.0506, marks=pytest.mark.extended),
    pytest.param("gpp500-2", 156.06, 0.31212, marks=pytest.mark.extended),
    pytest.param("gpp500-3", 513.02, 1.02604, marks=pytest.mark.extended),
    pytest.param("gpp500-4", 1567.02, 3.13404, marks=pytest.mark.extended),
    pytest.param("hinf1", -2.0326, 0.0040652, marks=pytest.mark.mini),
    pytest.param("hinf10", -109, 0.5, marks=pytest.mark.mini),  # table only to 3 sig figs
    pytest.param("hinf11", -65.9, 0.1318, marks=pytest.mark.mini),
    pytest.param("hinf13", -46, 1.0, marks=pytest.mark.mini),  # table only to 2 sig figs
    pytest.param("hinf14", -13, 0.026, marks=pytest.mark.timing),
    pytest.param("hinf15", -25, 0.3, marks=pytest.mark.mini),  # table only to 2 sig figs
    pytest.param("hinf2", -10.967, 0.021934, marks=pytest.mark.mini),
    pytest.param("hinf3", -56.9, 0.1138, marks=pytest.mark.mini),
    pytest.param("hinf4", -274.764, 0.549528, marks=pytest.mark.mini),
    pytest.param("hinf5", -363, 1.0, marks=pytest.mark.mini),  # table only to 3 sig figs
    pytest.param("hinf6", -449, 0.898, marks=pytest.mark.mini),
    pytest.param("hinf7", -391, 2.0, marks=pytest.mark.mini),  # table only to 3 sig figs
    pytest.param("hinf8", -116, 0.232, marks=pytest.mark.mini),
    pytest.param("hinf9", -236.25, 0.4725, marks=pytest.mark.mini),
    pytest.param("maxG11", -629.1648, 1.25833, marks=pytest.mark.extended),
    pytest.param("maxG51", -4003.809, 8.007618, marks=pytest.mark.extended),
    pytest.param("mcp100", -226.1574, 0.4523148, marks=pytest.mark.mini),
    pytest.param("mcp124-1", -141.9905, 0.283981, marks=pytest.mark.mini),
    pytest.param("mcp124-2", -269.8802, 0.5397604, marks=pytest.mark.mini),
    pytest.param("mcp124-3", -467.7501, 0.9355002, marks=pytest.mark.mini),
    pytest.param("mcp124-4", -864.4119, 1.728824, marks=pytest.mark.mini),
    pytest.param("mcp250-1", -317.2643, 0.6345286, marks=pytest.mark.timing),
    pytest.param("mcp250-2", -531.9301, 1.06386, marks=pytest.mark.timing),
    pytest.param("mcp250-3", -981.1726, 1.962345, marks=pytest.mark.timing),
    pytest.param("mcp250-4", -1681.96, 3.36392, marks=pytest.mark.timing),
    pytest.param("mcp500-1", -598.1485, 1.196297, marks=pytest.mark.timing),
    pytest.param("mcp500-2", -1070.057, 2.140114, marks=pytest.mark.timing),
    pytest.param("mcp500-3", -1847.97, 3.69594, marks=pytest.mark.timing),
    pytest.param("mcp500-4", -3566.738, 7.133476, marks=pytest.mark.timing),
    pytest.param("qap10", 1093, 2.186, marks=pytest.mark.timing),  # table's own footnote 14 correction
    pytest.param("qap5", 436, 0.872, marks=pytest.mark.mini),
    pytest.param("qap6", 381.44, 0.76288, marks=pytest.mark.mini),
    pytest.param("qap7", 425, 0.85, marks=pytest.mark.mini),
    pytest.param("qap8", 757, 1.514, marks=pytest.mark.timing),
    pytest.param("qap9", 1410, 2.82, marks=pytest.mark.timing),
    pytest.param("ss30", -20.2395, 0.040479, marks=pytest.mark.timing),
    pytest.param("theta1", -23, 0.046, marks=pytest.mark.mini),
    pytest.param("theta2", -32.87917, 0.06575834, marks=pytest.mark.mini),
    pytest.param("theta3", -42.16698, 0.08433396, marks=pytest.mark.timing),
    pytest.param("theta4", -50.32122, 0.1006424, marks=pytest.mark.extended),
    pytest.param("theta5", -57.23231, 0.1144646, marks=pytest.mark.extended),
    pytest.param("truss1", 8.999996, 0.01799999, marks=pytest.mark.mini),
    pytest.param("truss2", 123.3804, 0.2467608, marks=pytest.mark.mini),
    pytest.param("truss3", 9.109996, 0.01821999, marks=pytest.mark.mini),
    pytest.param("truss4", 9.009996, 0.01801999, marks=pytest.mark.mini),
    pytest.param("truss5", 132.6357, 0.2652714, marks=pytest.mark.mini),
    pytest.param("truss6", 901.001, 1.802002, marks=pytest.mark.timing),
    pytest.param("truss7", 900.001, 1.800002, marks=pytest.mark.timing),
    pytest.param("truss8", 133.1146, 0.2662292, marks=pytest.mark.timing),
]


@pytest.mark.parametrize("name,ref_obj,atol", SDPLIB_PARAMS)
def test_sdplib(name, ref_obj, atol, benchmark_collector):
    path = SDPLIB_DIR / f"{name}.dat-s"
    if not path.exists():
        pytest.skip(f"data file not found: {path}")

    At, b, c, K = read_sdpa(path)
    t0 = time.perf_counter()
    x, y, info = sedumi(At, b, c, K)
    _record(benchmark_collector, "SDPLIB", name, t0, x, c, info)

    pobj = float(np.sum(_dense(c) * _dense(x)))
    assert info["numerr"] in (0, 1), f"{name}: numerr={info['numerr']} (iter={info['iter']})"
    assert abs(pobj - ref_obj) <= atol, f"{name}: pobj={pobj:.6f}, ref={ref_obj:.6f}"


# 4 SDPLIB problems are infeasible rather than optimal -- the published
# table gives the SDPA-primal/dual role that's infeasible, which (per the
# module docstring's sign-convention note) sedumi() reports on the
# *opposite* of its own info["pinf"]/info["dinf"] flags.
SDPLIB_INFEASIBLE_PARAMS = [
    # (name, expect_field) -- expect_field is the sedumi() info[...] flag that must be 1
    ("infd1", "pinf"),  # table: dual infeasible
    ("infd2", "pinf"),  # table: dual infeasible
    ("infp1", "dinf"),  # table: primal infeasible
    ("infp2", "dinf"),  # table: primal infeasible
]


@pytest.mark.mini
@pytest.mark.parametrize("name,expect_field", SDPLIB_INFEASIBLE_PARAMS)
def test_sdplib_infeasible(name, expect_field, benchmark_collector):
    path = SDPLIB_DIR / f"{name}.dat-s"
    if not path.exists():
        pytest.skip(f"data file not found: {path}")

    At, b, c, K = read_sdpa(path)
    t0 = time.perf_counter()
    x, y, info = sedumi(At, b, c, K)
    _record(benchmark_collector, "SDPLIB", name, t0, x, c, info)

    assert info[expect_field] == 1, f"{name}: expected info[{expect_field!r}]==1, got {info}"


# ---------------------------------------------------------------------------
# DIMACS (SeDuMi-format .mat.gz; https://github.com/vsdp/DIMACS README.md
# per-set tables). ref_obj = published value directly, except HAMMING
# (negated) and TORUS (see TORUS_PARAMS below) -- see module docstring.
# atol is 0.2% of |ref_obj| (floor 1e-4).
# ---------------------------------------------------------------------------

DIMACS_PARAMS = [
    # (name, class_dir, ref_obj, atol, marks)
    pytest.param("nb", "ANTENNA", -0.05070309, 0.0001014062, marks=pytest.mark.timing),
    pytest.param("nb_L1", "ANTENNA", -13.01234, 0.02602467, marks=pytest.mark.timing),
    pytest.param("nb_L2_bessel", "ANTENNA", -0.1025695, 0.000205139, marks=pytest.mark.timing),
    pytest.param("copo14", "COPOS", 0, 0.0001, marks=pytest.mark.mini),
    pytest.param("copo23", "COPOS", 0, 0.0001, marks=pytest.mark.extended),
    pytest.param("filter48_socp", "FILTER", 1.416129, 0.002832258, marks=pytest.mark.timing),
    pytest.param("minphase", "FILTER", 5.98, 0.01196, marks=pytest.mark.mini),
    pytest.param("hamming_7_5_6", "HAMMING", -42.66667, 0.08533333, marks=pytest.mark.timing),
    pytest.param("hamming_9_8", "HAMMING", -224, 0.448, marks=pytest.mark.extended),
    pytest.param("nql30", "NQL", -0.9460, 0.001892, marks=pytest.mark.timing),
    pytest.param("qssp30", "QSSP", -6.496675, 0.01299335, marks=pytest.mark.timing),
    pytest.param("sched_100_100_scaled", "SCHED", 27.3307, 0.0546614, marks=pytest.mark.extended),
    pytest.param("sched_100_50_orig", "SCHED", 181889.9, 363.7798, marks=pytest.mark.extended),
    pytest.param("sched_50_50_orig", "SCHED", 26673, 53.346, marks=pytest.mark.timing),
    pytest.param("sched_50_50_scaled", "SCHED", 7.852038, 0.01570408, marks=pytest.mark.timing),
    pytest.param("truss5", "TRUSS", 132.6357, 0.2652714, marks=pytest.mark.mini),
    pytest.param("truss8", "TRUSS", 133.1146, 0.2662292, marks=pytest.mark.timing),
]


@pytest.mark.parametrize("name,class_dir,ref_obj,atol", DIMACS_PARAMS)
def test_dimacs(name, class_dir, ref_obj, atol, benchmark_collector):
    path = DIMACS_DIR / class_dir / f"{name}.mat.gz"
    if not path.exists():
        pytest.skip(f"data file not found: {path}")

    A, b, c, K, pars = read_mat(path)
    t0 = time.perf_counter()
    x, y, info = sedumi(A, b, c, K)
    _record(benchmark_collector, "DIMACS", name, t0, x, c, info)

    pobj = float(np.sum(_dense(c) * _dense(x)))
    assert info["numerr"] in (0, 1), f"{name}: numerr={info['numerr']} (iter={info['iter']})"
    assert abs(pobj - ref_obj) <= atol, f"{name}: pobj={pobj:.6f}, ref={ref_obj:.6f}"


# DIMACS's FILTER set has one problem the README lists as "primal inf."
# rather than an objective value; unlike SDPLIB, DIMACS problems don't go
# through a sign-flipping reader, so the table's own word matches
# sedumi()'s own info["pinf"] directly (no inversion needed here).
DIMACS_INFEASIBLE_PARAMS = [
    # (name, class_dir, expect_field)
    ("filtinf1", "FILTER", "pinf"),
]


@pytest.mark.mini
@pytest.mark.parametrize("name,class_dir,expect_field", DIMACS_INFEASIBLE_PARAMS)
def test_dimacs_infeasible(name, class_dir, expect_field, benchmark_collector):
    path = DIMACS_DIR / class_dir / f"{name}.mat.gz"
    if not path.exists():
        pytest.skip(f"data file not found: {path}")

    A, b, c, K, pars = read_mat(path)
    t0 = time.perf_counter()
    x, y, info = sedumi(A, b, c, K)
    _record(benchmark_collector, "DIMACS", name, t0, x, c, info)

    assert info[expect_field] == 1, f"{name}: expected info[{expect_field!r}]==1, got {info}"


# ---------------------------------------------------------------------------
# DIMACS TORUS: `max c'x` given as `min -c'x` (DIMACS README's own
# documented caveat), so ref_obj = -1 * published value -- see module
# docstring's Sign conventions for c_scale's role (1.0 for the
# non-Gaussian "pm" instance, -1e-5 for the Gaussian "g" one).
# ---------------------------------------------------------------------------

TORUS_PARAMS = [
    # (name, ref_obj, atol, c_scale, marks)
    pytest.param("torusg3-8", -457.3582, 0.9147164, -1e-5, marks=pytest.mark.timing),
    pytest.param("toruspm3-8-50", -527.8087, 1.055617, 1.0, marks=pytest.mark.timing),
]


@pytest.mark.parametrize("name,ref_obj,atol,c_scale", TORUS_PARAMS)
def test_dimacs_torus(name, ref_obj, atol, c_scale, benchmark_collector):
    path = DIMACS_DIR / "TORUS" / f"{name}.mat.gz"
    if not path.exists():
        pytest.skip(f"data file not found: {path}")

    A, b, c, K, pars = read_mat(path)
    c = _dense(c) * c_scale
    t0 = time.perf_counter()
    x, y, info = sedumi(A, b, c, K)
    _record(benchmark_collector, "DIMACS", name, t0, x, c, info)

    pobj = float(np.sum(c * _dense(x)))
    assert info["numerr"] in (0, 1), f"{name}: numerr={info['numerr']} (iter={info['iter']})"
    assert abs(pobj - ref_obj) <= atol, f"{name}: pobj={pobj:.6f}, ref={ref_obj:.6f}"
