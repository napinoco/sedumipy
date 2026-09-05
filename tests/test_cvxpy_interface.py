"""cvxpy interface (cvxpy_interface.py): sedumipy driven as a cvxpy
custom solver.

Every solve here is checked against cvxpy's own default conic solver
(CLARABEL, a hard dependency of cvxpy itself) on the same problem --
objective value, variable values, and dual values -- since what this
module has to get right is not the interior-point math (that is what the
rest of the suite covers, against Octave) but the two translations
around it: cvxpy's (A, b, c, cones) canonical form into SeDuMi's own
(At, b, c, K) dual form, and SeDuMi's (x, y, info) back into cvxpy's
primal/dual variables and status. A sign or ordering slip in either
direction shows up immediately as a mismatch against CLARABEL.

Both solvers are asked for a tighter-than-default tolerance and the
comparison is then made at 1e-5. That is not a statement about how
accurately sedumipy solves anything (docs/status.rst and the Octave-
oracle tests cover that): at each solver's default tolerance the
optimal *point* of a flat problem like the SOCP below is only pinned
down to ~1e-4 by either solver, which says nothing about the interface,
and this way a genuine translation bug -- a sign flip, a reordered cone
block -- cannot hide inside the tolerance either.
"""

from __future__ import annotations

import numpy as np
import pytest

sedumipy = pytest.importorskip("sedumipy")
cp = pytest.importorskip("cvxpy")

from sedumipy.cvxpy_interface import SEDUMIPY, dims_to_solver_dict  # noqa: E402

TOL = 1e-5

# Tighter than either solver's default stopping tolerance -- see this
# module's docstring.
SEDUMIPY_OPTS = {"eps": 1e-12}

# The reference solver every problem below is cross-checked against.
REFERENCE = cp.CLARABEL
REFERENCE_OPTS = {"tol_gap_abs": 1e-11, "tol_gap_rel": 1e-11, "tol_feas": 1e-11}


def _flat(value) -> np.ndarray:
    """Variable/dual values as one flat float array -- an SOC constraint's
    dual_value is a list of (scalar, vector) rather than one array."""
    if isinstance(value, (list, tuple)):
        return np.concatenate([np.atleast_1d(np.asarray(v, dtype=float)).ravel() for v in value])
    return np.atleast_1d(np.asarray(value, dtype=float)).ravel()


def _solve(problem, solver, opts, variables, constraints):
    """(value, variable values, dual values) from one solve."""
    value = problem.solve(solver=solver, **opts)
    return (
        value,
        [_flat(v.value) for v in variables],
        [_flat(c.dual_value) for c in constraints],
    )


def _assert_matches_reference(problem, variables, constraints) -> None:
    ours = _solve(problem, SEDUMIPY(), SEDUMIPY_OPTS, variables, constraints)
    assert problem.status == cp.settings.OPTIMAL
    ref = _solve(problem, REFERENCE, REFERENCE_OPTS, variables, constraints)

    assert ours[0] == pytest.approx(ref[0], abs=TOL, rel=TOL)
    for got, want in zip(ours[1], ref[1]):
        np.testing.assert_allclose(got, want, atol=TOL, rtol=TOL)
    for got, want in zip(ours[2], ref[2]):
        np.testing.assert_allclose(got, want, atol=TOL, rtol=TOL)


def test_lp() -> None:
    """The README's own LP, as a cvxpy problem."""
    A = np.array([[3.0, 1.0, 2.0], [1.0, 2.0, 4.0]])
    b = np.array([9.0, 8.0])
    c = np.array([7.0, 4.0, 10.0])
    x = cp.Variable(3)
    constraints = [A @ x == b, x >= 0]
    problem = cp.Problem(cp.Minimize(c @ x), constraints)

    _assert_matches_reference(problem, [x], constraints)
    np.testing.assert_allclose(x.value, [2.0, 3.0, 0.0], atol=TOL)
    assert problem.value == pytest.approx(26.0, abs=TOL)


def test_lp_with_free_variables() -> None:
    """A K.f (free) block: variables with no sign constraint at all."""
    rng = np.random.default_rng(0)
    A = rng.standard_normal((3, 5))
    c = rng.standard_normal(5)
    x = cp.Variable(5)
    constraints = [A @ x == np.ones(3), x <= 4, x >= -4]
    problem = cp.Problem(cp.Minimize(c @ x), constraints)

    _assert_matches_reference(problem, [x], constraints)


def test_socp() -> None:
    """K.q blocks, from an explicit cp.SOC and from a norm atom."""
    rng = np.random.default_rng(1)
    A = rng.standard_normal((6, 3))
    b = rng.standard_normal(6)
    c = rng.standard_normal(3)
    x = cp.Variable(3)
    t = cp.Variable()
    constraints = [cp.SOC(t, A @ x - b), cp.norm(x, 2) <= 2]
    problem = cp.Problem(cp.Minimize(t + 0.5 * (c @ x)), constraints)

    _assert_matches_reference(problem, [x, t], constraints)


def test_sdp() -> None:
    """A K.s block: minimize <C, X> over the spectraplex."""
    rng = np.random.default_rng(2)
    C = rng.standard_normal((4, 4))
    C = (C + C.T) / 2
    X = cp.Variable((4, 4), symmetric=True)
    constraints = [X >> 0, cp.trace(X) == 1]
    problem = cp.Problem(cp.Minimize(cp.trace(C @ X)), constraints)

    _assert_matches_reference(problem, [X], constraints)


def test_mixed_cones_with_objective_offset() -> None:
    """All four cone kinds at once (f/l/q/s), plus a constant term in the
    objective -- the offset cvxpy keeps out of the solver data and adds
    back in invert()."""
    rng = np.random.default_rng(3)
    C = rng.standard_normal((3, 3))
    C = (C + C.T) / 2
    X = cp.Variable((3, 3), symmetric=True)
    y = cp.Variable(2)
    z = cp.Variable(4)
    constraints = [
        X >> 0,
        cp.trace(X) == 1,
        cp.norm(y, 2) <= 1,
        y[0] + z[0] == 0.5,
        z >= 0,
        z <= 3,
    ]
    objective = cp.Minimize(cp.trace(C @ X) + y[1] + cp.sum(z) + 3.5)
    problem = cp.Problem(objective, constraints)

    _assert_matches_reference(problem, [X, y, z], constraints)


def test_parameters_resolve() -> None:
    """A parametrized problem, solved twice: cvxpy re-applies the
    parameter values to the same canonicalized data between solves."""
    A = np.array([[3.0, 1.0, 2.0], [1.0, 2.0, 4.0]])
    c = np.array([7.0, 4.0, 10.0])
    rhs = cp.Parameter(2)
    x = cp.Variable(3)
    problem = cp.Problem(cp.Minimize(c @ x), [A @ x == rhs, x >= 0])

    rhs.value = np.array([9.0, 8.0])
    assert problem.solve(solver=SEDUMIPY()) == pytest.approx(26.0, abs=TOL)

    rhs.value = np.array([9.0, 7.0])
    ours = problem.solve(solver=SEDUMIPY())
    assert ours == pytest.approx(problem.solve(solver=REFERENCE, **REFERENCE_OPTS), abs=TOL)


def test_infeasible() -> None:
    """cvxpy's problem is SeDuMi's *dual*, so an infeasible cvxpy problem
    is the info["dinf"] certificate, not info["pinf"]."""
    x = cp.Variable(2)
    problem = cp.Problem(cp.Minimize(cp.sum(x)), [x >= 1, x <= 0])

    problem.solve(solver=SEDUMIPY())
    assert problem.status == cp.settings.INFEASIBLE
    assert problem.value == np.inf


def test_unbounded() -> None:
    """...and an unbounded one is info["pinf"] (an improving direction)."""
    x = cp.Variable(2)
    problem = cp.Problem(cp.Minimize(cp.sum(x)), [x <= 1])

    problem.solve(solver=SEDUMIPY())
    assert problem.status == cp.settings.UNBOUNDED
    assert problem.value == -np.inf


def test_solver_stats() -> None:
    """iteration count, solve time and the raw `info` dict come back
    through cvxpy's own solver_stats."""
    x = cp.Variable(3)
    A = np.array([[3.0, 1.0, 2.0], [1.0, 2.0, 4.0]])
    problem = cp.Problem(cp.Minimize(cp.sum(x)), [A @ x == np.array([9.0, 8.0]), x >= 0])
    problem.solve(solver=SEDUMIPY())

    stats = problem.solver_stats
    assert stats.solver_name == "SEDUMIPY"
    assert stats.num_iters > 0
    assert stats.solve_time > 0
    assert stats.extra_stats["numerr"] == 0
    assert set(stats.extra_stats) >= {"iter", "pinf", "dinf", "numerr", "r0"}


def test_solver_options_are_passed_through() -> None:
    """`solve()`'s extra keyword arguments reach sedumi() as `pars`: two
    iterations is not enough for this problem, so it comes back as a
    solver error rather than an answer."""
    x = cp.Variable(3)
    A = np.array([[3.0, 1.0, 2.0], [1.0, 2.0, 4.0]])
    problem = cp.Problem(cp.Minimize(cp.sum(x)), [A @ x == np.array([9.0, 8.0]), x >= 0])

    with pytest.raises(cp.error.SolverError):
        problem.solve(solver=SEDUMIPY(), maxiter=2)

    assert problem.solve(solver=SEDUMIPY(), maxiter=100) == pytest.approx(
        problem.solve(solver=REFERENCE, **REFERENCE_OPTS), abs=TOL
    )


def test_unsupported_cone_is_refused() -> None:
    """No exponential cone in sedumipy, so cvxpy must not route an
    exp/log problem here (rather than getting a wrong answer)."""
    x = cp.Variable()
    problem = cp.Problem(cp.Minimize(cp.exp(x)), [x >= 1])

    with pytest.raises(cp.error.SolverError):
        problem.solve(solver=SEDUMIPY())


def test_mixed_integer_is_refused() -> None:
    x = cp.Variable(2, integer=True)
    problem = cp.Problem(cp.Minimize(cp.sum(x)), [x >= 0.5])

    with pytest.raises(cp.error.SolverError):
        problem.solve(solver=SEDUMIPY())


def test_name_does_not_collide_with_a_built_in_solver() -> None:
    """cvxpy rejects a custom solver named like one of its own, so this
    name must stay out of cvxpy's table -- including in future cvxpy
    versions that might add a built-in SeDuMi interface."""
    assert SEDUMIPY().name() == "SEDUMIPY"
    assert SEDUMIPY().name() not in cp.settings.SOLVERS
    assert SEDUMIPY().is_installed()


def test_dims_to_solver_dict() -> None:
    """cvxpy's ConeDims -> SeDuMi's K, in cvxpy's f/l/q/s order."""
    x = cp.Variable(3)
    X = cp.Variable((2, 2), symmetric=True)
    t = cp.Variable()
    problem = cp.Problem(
        cp.Minimize(cp.sum(x) + t + cp.trace(X)),
        [cp.sum(x) == 1, x >= 0, cp.SOC(t, x), X >> 0],
    )
    data = problem.get_problem_data(solver=SEDUMIPY())[0]

    K = dims_to_solver_dict(data["dims"])
    assert K["f"] == 1
    assert K["l"] == 3
    assert K["q"] == [4]
    assert K["s"] == [2]
    # SeDuMi's K covers exactly the rows cvxpy stuffed into the data.
    assert K["f"] + K["l"] + sum(K["q"]) + sum(n * n for n in K["s"]) == data["A"].shape[0]
