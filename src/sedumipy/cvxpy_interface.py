"""cvxpy interface: lets sedumipy be used as a solver from cvxpy,
without anything having to be merged into cvxpy itself.

cvxpy resolves ``problem.solve(solver=...)`` either by name (a string,
looked up in its own table of built-in solvers) or by being handed a
``cvxpy.reductions.solvers.solver.Solver`` instance directly -- the
latter is its "custom solver" path: ``construct_solving_chain()`` in
cvxpy/reductions/solvers/solving_chain.py returns any ``Solver``
instance it is given, after checking that the instance's ``name()``
does not collide with a built-in solver's. This module takes that
path::

    import cvxpy as cp
    from sedumipy.cvxpy_interface import SEDUMIPY

    x = cp.Variable(3)
    prob = cp.Problem(cp.Minimize(c @ x), [A @ x == b, x >= 0])
    prob.solve(solver=SEDUMIPY())

so it works against a stock cvxpy install (>= 1.3), with no patched or
forked cvxpy anywhere.

PROBLEM FORM. After canonicalization cvxpy hands a conic solver a
triple (A, b, c) in SCS's convention,

    minimize  c'z  such that  b - A*z in K,

where K is a product of cones in cvxpy's own fixed order: zero (the
equality constraints), nonnegative orthant, second-order, PSD. SeDuMi
solves the primal-dual pair

    (P) min c_s'x s.t. A_s x = b_s, x in K
    (D) max b_s'y s.t. c_s - A_s'y in K,

so cvxpy's problem *is* SeDuMi's (D) -- take y := z, A_s' := A,
c_s := b, and b_s := -c, and (D) becomes "maximize -c'z such that
b - A*z in K", i.e. cvxpy's problem with the objective negated (SeDuMi
maximizes where cvxpy minimizes; the optimal value itself is recomputed
as c'z below rather than read off SeDuMi's objective).

sedumipy.sedumi() takes its first argument in SeDuMi's own "At"
orientation -- ``length(c) x length(b)``, i.e. A_s' -- so cvxpy's A goes
in untransposed (pretransfo.py's shape check picks that orientation
from the row count, unambiguously, since A has as many rows as K has
dimensions):

    x, y, info = sedumi(A, -c, b, K)

and then

  * ``y`` (one entry per cvxpy variable) is the primal solution z;
  * ``x`` (one entry per cone dimension, x in K) holds the cvxpy dual
    variables. No sign flip is needed: (P)'s own constraint A_s x = b_s
    reads A'x = -c, which is exactly the dual feasibility condition
    cvxpy's duals satisfy under the same convention (cf. SCS's own
    interface, cvxpy/reductions/solvers/conic_solvers/scs_conif.py).
    It is split at ``dims.zero`` into the equality- and inequality-
    constraint duals.

Because cvxpy's problem is SeDuMi's dual, the two infeasibility flags
swap on the way back: SeDuMi reporting *its* dual infeasible
(``info["dinf"]``) means cvxpy's problem is infeasible, and SeDuMi
reporting its primal infeasible (``info["pinf"]``, i.e. an improving
direction for cvxpy's problem) means cvxpy's problem is unbounded.

SCOPE. Zero, nonnegative, second-order and PSD cones -- the cones
sedumipy itself supports. Exponential/power cones are not in scope
(``log``, ``exp``, ``entr``, ``power`` with a fractional exponent, ...),
and neither are mixed-integer problems; cvxpy raises its own
``SolverError`` when a problem needs a cone this solver does not list.
PSD constraints are taken in cvxpy's ``PSD`` (not ``SvecPSD``) form,
whose rows are the full ``n**2`` entries of the matrix in column-major
order -- exactly SeDuMi's own ``K["s"]`` block layout, so no triangle
packing and no sqrt(2) scaling of off-diagonals is involved, and
``PSD_TRIANGLE_KIND``/``PSD_SQRT2_SCALING`` stay unset.

ON THE NAME. The solver calls itself ``SEDUMIPY``, not ``SEDUMI``:
cvxpy rejects a custom solver whose name collides with one of its
built-in solvers, so were cvxpy ever to gain a built-in ``SEDUMI``
interface, a name of "SEDUMI" here would stop working overnight. It
also names the package a user actually installs.
"""

from __future__ import annotations

import time

import numpy as np

# cvxpy is an optional dependency: this module is importable only where
# it is installed (`pip install sedumipy[cvxpy]`), and nothing in the
# sedumipy package itself imports it.
import cvxpy.settings as s
from cvxpy.constraints import PSD, SOC
from cvxpy.reductions.solution import Solution, failure_solution
from cvxpy.reductions.solvers import utilities
from cvxpy.reductions.solvers.conic_solvers.conic_solver import ConicSolver

#: The name cvxpy knows this solver by, as returned by ``SEDUMIPY.name()``.
SOLVER_NAME = "SEDUMIPY"


def dims_to_solver_dict(cone_dims) -> dict:
    """K = dims_to_solver_dict(cone_dims): cvxpy's ConeDims -> SeDuMi's
    cone-structure dict. The two orderings already agree (zero/f,
    nonneg/l, soc/q, psd/s), so this is a rename, not a reordering."""
    return {
        "f": int(cone_dims.zero),
        "l": int(cone_dims.nonneg),
        "q": [int(v) for v in cone_dims.soc],
        "s": [int(v) for v in cone_dims.psd],
    }


def _status(info: dict) -> str:
    """cvxpy status string from sedumi()'s `info` dict. `info["pinf"]`/
    `info["dinf"]` are SeDuMi's certificates for *its* primal/dual, which
    are cvxpy's dual/primal respectively (see this module's docstring),
    and `info["numerr"]` grades the accuracy actually reached: 0 = to
    `pars["eps"]`, 1 = to `pars["bigeps"]` only, 2 = neither (sedumi.m's
    own "numerical problems" outcome, which leaves no usable answer)."""
    pinf, dinf = int(info["pinf"]), int(info["dinf"])
    numerr = int(info["numerr"])
    if numerr == 2 or (pinf and dinf):
        return s.SOLVER_ERROR
    inaccurate = numerr == 1
    if dinf:  # SeDuMi's dual == cvxpy's problem
        return s.INFEASIBLE_INACCURATE if inaccurate else s.INFEASIBLE
    if pinf:
        return s.UNBOUNDED_INACCURATE if inaccurate else s.UNBOUNDED
    return s.OPTIMAL_INACCURATE if inaccurate else s.OPTIMAL


class SEDUMIPY(ConicSolver):
    """sedumipy as a cvxpy conic solver: ``prob.solve(solver=SEDUMIPY())``.

    An instance is what cvxpy wants (not the class), and one instance is
    reusable across solves -- it holds no per-problem state.
    """

    # Solver capabilities.
    MIP_CAPABLE = False
    SUPPORTED_CONSTRAINTS = ConicSolver.SUPPORTED_CONSTRAINTS + [SOC, PSD]
    # sedumi() needs a cone to work over; an unconstrained problem has none.
    REQUIRES_CONSTR = True

    def name(self) -> str:
        """The name of the solver."""
        return SOLVER_NAME

    def import_solver(self) -> None:
        """Imports the solver."""
        import sedumipy  # noqa: F401

    def apply(self, problem):
        """Returns a new problem and data for inverting the new solution.

        ConicSolver's own implementation already produces (A, b, c) in
        the ``b - A*z in K`` convention this interface wants, with the
        cones in the order SeDuMi's K expects.
        """
        return super().apply(problem)

    def solve_via_data(
        self, data, warm_start: bool, verbose: bool, solver_opts, solver_cache=None
    ):
        """Solves the canonicalized problem via sedumipy.sedumi().

        `solver_opts` is passed straight through as sedumi()'s `pars`
        (`eps`, `maxiter`, `bigeps`, ... -- see checkpars.py for the
        full list). `warm_start` and `verbose` are both ignored: SeDuMi
        is an interior-point method with no warm-start interface, and
        this port does not implement sedumi.m's console progress
        printout (see sedumi.py's own docstring).
        """
        from sedumipy import sedumi

        A = data[s.A]
        b = np.asarray(data[s.B], dtype=np.float64).ravel()
        c = np.asarray(data[s.C], dtype=np.float64).ravel()
        K = dims_to_solver_dict(data[self.DIMS])

        start = time.perf_counter()
        x, y, info = sedumi(A, -c, b, K, dict(solver_opts))
        solve_time = time.perf_counter() - start

        x = np.asarray(x, dtype=np.float64).ravel()
        y = np.asarray(y, dtype=np.float64).ravel()

        solution = {
            s.STATUS: _status(info),
            s.PRIMAL: y,
            s.EQ_DUAL: x[: K["f"]],
            s.INEQ_DUAL: x[K["f"]:],
            s.VALUE: float(c @ y),
            s.SOLVE_TIME: solve_time,
            s.NUM_ITERS: int(info["iter"]),
            s.EXTRA_STATS: info,
        }
        return solution

    def invert(self, solution, inverse_data):
        """Returns the solution to the original problem given the
        inverse_data."""
        status = solution[s.STATUS]
        attr = {
            s.SOLVE_TIME: solution[s.SOLVE_TIME],
            s.NUM_ITERS: solution[s.NUM_ITERS],
            s.EXTRA_STATS: solution[s.EXTRA_STATS],
        }

        if status not in s.SOLUTION_PRESENT:
            return failure_solution(status, attr)

        primal_vars = {inverse_data[self.VAR_ID]: solution[s.PRIMAL]}
        dual_vars = utilities.get_dual_values(
            solution[s.EQ_DUAL],
            utilities.extract_dual_value,
            inverse_data[self.EQ_CONSTR],
        )
        dual_vars.update(
            utilities.get_dual_values(
                solution[s.INEQ_DUAL],
                utilities.extract_dual_value,
                inverse_data[self.NEQ_CONSTR],
            )
        )
        opt_val = solution[s.VALUE] + inverse_data[s.OFFSET]
        return Solution(status, opt_val, primal_vars, dual_vars, attr)

    def cite(self, data) -> str:
        """Returns bibtex citation for the solver."""
        return """
@article{sedumi,
    author = {Sturm, Jos F.},
    title = {Using {SeDuMi} 1.02, a {MATLAB} toolbox for optimization
             over symmetric cones},
    journal = {Optimization Methods and Software},
    volume = {11},
    number = {1-4},
    pages = {625--653},
    year = {1999},
    doi = {10.1080/10556789908805766},
}
"""
