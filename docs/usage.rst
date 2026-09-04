Usage
=====

The public entry point is :func:`sedumipy.sedumi`, matching real SeDuMi's
own signature:

.. code-block:: python

   x, y, info = sedumipy.sedumi(A, b, c, K, pars=None, **pars_kwargs)

which solves the primal problem

.. math::

   \text{(P)} \qquad \min_x\ c^T x \quad \text{such that } Ax = b,\ x \in K

together with its dual

.. math::

   \text{(D)} \qquad \max_{y,s}\ b^T y \quad \text{such that } A^T y + s = c,\ s \in K^*

One call solves both: ``sedumi()`` returns the primal optimum ``x`` and
the dual optimum ``y``. The dual slack ``s = c - A^T y`` is not returned,
but is one line to recompute. Every cone ``sedumipy`` supports is
self-dual (:math:`K^* = K`) except the free block ``K.f``, whose dual is
:math:`\{0\}` -- so its dual-slack entries are pinned to zero. This
primal-dual pair, and the meaning of each ``K`` field below, follow real
SeDuMi's own convention; see the `Addendum to the SeDuMi User Guide
<https://sedumi.ie.lehigh.edu/sedumi/files/sedumi-downloads/SeDuMi_Guide_11.pdf>`_
(Pólik, 2005), Section 2.

``A`` may be given either as :math:`m \times n` or its transpose
:math:`n \times m` (SeDuMi disambiguates by shape, matching real SeDuMi);
dense NumPy arrays and SciPy sparse matrices are both accepted.

The cone structure ``K``
-------------------------

``K`` is a plain ``dict`` describing how the columns of ``x`` are carved
up into cones, in the same field-name convention as real SeDuMi:

============ ===================================================
``K["f"]``   number of free (unrestricted) variables
``K["l"]``   number of variables constrained to the nonnegative
             orthant (``x >= 0``)
``K["q"]``   sizes of Lorentz (second-order) cone blocks
``K["r"]``   sizes of rotated Lorentz cone blocks
``K["s"]``   sizes of positive semidefinite (PSD) cone blocks
============ ===================================================

Only the fields you need have to be present; omitted fields default to
none of that cone type. Blocks appear in ``x`` in ``f, l, q, r, s``
order, each block occupying that many consecutive entries (``s`` blocks
occupy ``size**2`` entries, stored column-major/vec).

.. code-block:: python

   import numpy as np
   import sedumipy

   # LP: minimize x1 + x2 s.t. x1 = x2 = 1, x >= 0
   x, y, info = sedumipy.sedumi(np.eye(2), np.array([1.0, 1.0]), np.array([1.0, 1.0]), {"l": 2})

Worked example: turning a written model into ``A``, ``b``, ``c``, ``K``
------------------------------------------------------------------------

Models rarely arrive already in the shape of (P). Far more often you have
a handful of decision variables and a list of constraints, each asking
some *affine function of those variables* to be nonnegative, or to have
bounded norm, or to be positive semidefinite. That is precisely the shape
of the dual (D), so that is the side to map onto. Take this problem,
which mixes an LP block, two second-order-cone blocks, and a 3x3 PSD
block:

.. math::

   \begin{aligned}
   \max_{y \in \mathbb{R}^3} \quad & 6y_1 + 4y_2 + 5y_3 \\[4pt]
   \text{s.t.} \quad
   & 16y_1 - 14y_2 + 5y_3 \le -3, \qquad 7y_1 + 2y_2 \le 5, \\[6pt]
   & \left\|
       \begin{pmatrix} 8y_1 + 13y_2 - 12y_3 - 2 \\
                       -8y_1 + 18y_2 + 6y_3 - 14 \\
                       y_1 - 3y_2 - 17y_3 - 13 \end{pmatrix}
     \right\| \le -24y_1 - 7y_2 + 15y_3 + 12,
     \qquad
     \left\| \begin{pmatrix} y_1 \\ y_2 \\ y_3 \end{pmatrix} \right\| \le 10, \\[6pt]
   & \begin{pmatrix}
       7y_1+3y_2+9y_3 & -5y_1+13y_2+6y_3 & y_1-6y_2-6y_3 \\
       -5y_1+13y_2+6y_3 & y_1+12y_2-7y_3 & -7y_1-10y_2-7y_3 \\
       y_1-6y_2-6y_3 & -7y_1-10y_2-7y_3 & -4y_1-28y_2-11y_3
     \end{pmatrix}
     \preceq
     \begin{pmatrix} 68 & -30 & -19 \\ -30 & 99 & 23 \\ -19 & 23 & 10 \end{pmatrix}
   \end{aligned}

The cone and the objective vector can be read straight off: there are
:math:`m = 3` variables, so ``b`` is the objective's coefficient vector,
and the five constraints group into four cone blocks.

.. code-block:: python

   b = np.array([6.0, 4.0, 5.0])
   K = {"l": 2, "q": [4, 4], "s": [3]}

That leaves ``A`` and ``c``.

Step 1: rewrite every constraint as "constant minus linear"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

(D)'s constraint is :math:`s = c - A^T y \in K^*`, so each block's job is
to express its own slack :math:`s` as a *constant minus something linear
in* :math:`y`. Those constants become that block's slice of ``c``, and
the coefficients of :math:`y_k` become its slice of column :math:`k` of
:math:`A^T`. Working block by block:

**The LP block** (``K.l = 2``). Move each ``<=`` to a nonnegativity:

.. math::

   16y_1 - 14y_2 + 5y_3 \le -3
   \;\Longleftrightarrow\;
   \underbrace{(-3)}_{c_1} - \underbrace{(16y_1 - 14y_2 + 5y_3)}_{\text{row of } A^T} \ge 0

and likewise :math:`5 - (7y_1 + 2y_2) \ge 0`, giving
:math:`c = (-3, 5)` and rows :math:`(16, -14, 5)`, :math:`(7, 2, 0)`.

**The Lorentz blocks** (``K.q = [4, 4]``). A ``K.q`` block of size 4 is
one 4-vector :math:`s = (s_0; \bar{s})` constrained by
:math:`s_0 \ge \|\bar{s}\|` -- **the bound comes first**, then the three
components whose norm it bounds. So stack the right-hand side on top of
the norm's argument, then split each of those four entries the same way
as above:

.. math::

   s = \begin{pmatrix}
     \;12 - (24y_1 + 7y_2 - 15y_3) \\
     -2 - (-8y_1 - 13y_2 + 12y_3) \\
     -14 - (8y_1 - 18y_2 - 6y_3) \\
     -13 - (-y_1 + 3y_2 + 17y_3)
   \end{pmatrix}
   \quad\Longrightarrow\quad
   c = \begin{pmatrix} 12 \\ -2 \\ -14 \\ -13 \end{pmatrix},\qquad
   \text{rows of } A^T =
   \begin{pmatrix} 24 & 7 & -15 \\ -8 & -13 & 12 \\
                   8 & -18 & -6 \\ -1 & 3 & 17 \end{pmatrix}

The second block, :math:`\|(y_1, y_2, y_3)\| \le 10`, is the same idea
with nothing to rearrange: :math:`c = (10, 0, 0, 0)` and the rows
:math:`(0,0,0)`, :math:`(-1,0,0)`, :math:`(0,-1,0)`, :math:`(0,0,-1)`.

**The PSD block** (``K.s = [3]``). Rearranged the same way,
:math:`M(y) \preceq C` is :math:`C - M(y) \succeq 0`, where the matrix on
the left splits into one coefficient matrix per variable,
:math:`M(y) = y_1 A_1 + y_2 A_2 + y_3 A_3`:

.. math::

   A_1 = \begin{pmatrix} 7 & -5 & 1 \\ -5 & 1 & -7 \\ 1 & -7 & -4 \end{pmatrix},\;
   A_2 = \begin{pmatrix} 3 & 13 & -6 \\ 13 & 12 & -10 \\ -6 & -10 & -28 \end{pmatrix},\;
   A_3 = \begin{pmatrix} 9 & 6 & -6 \\ 6 & -7 & -7 \\ -6 & -7 & -11 \end{pmatrix}

So the "constant" is the matrix :math:`C` and the "coefficient of
:math:`y_k`" is the matrix :math:`A_k`. Both are matrices, and neither
``c`` nor :math:`A^T` has anywhere to put a matrix -- which brings us to
the one genuinely fiddly part.

Step 2: flatten the PSD block's matrices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``sedumipy`` has no matrix-valued slots. An ``n``-by-``n`` PSD block is
stored as ``n**2`` consecutive *scalar* entries of ``x``, ``c``, and
``s``, in column-major (``vec``) order:

.. math::

   \operatorname{vec}
   \begin{pmatrix} S_{11} & S_{12} & S_{13} \\
                   S_{21} & S_{22} & S_{23} \\
                   S_{31} & S_{32} & S_{33} \end{pmatrix}
   = (S_{11}, S_{21}, S_{31},\; S_{12}, S_{22}, S_{32},\; S_{13}, S_{23}, S_{33})^T

which is exactly NumPy's ``M.flatten(order="F")``, inverted by
``v.reshape(3, 3, order="F")``. Note that this is the *full* ``n**2``
vectorization: all nine entries, with the off-diagonal ones genuinely
written twice, and no :math:`\sqrt{2}` scaling anywhere. (Some solvers
instead take a half-vectorized, :math:`\sqrt{2}`-scaled ``svec`` of
length ``n(n+1)/2``; SeDuMi does not.)

Apply ``vec`` to the matrices from step 1 and the PSD block stops being
special: ``vec(C)`` is that block's slice of ``c``, and
``vec(A_k)`` is that block's slice of column ``k`` of :math:`A^T` -- nine
rows, just as the LP block contributed two and each Lorentz block four.

Step 3: stack the blocks
~~~~~~~~~~~~~~~~~~~~~~~~~

Blocks occupy consecutive entries in ``f, l, q, r, s`` order, so ``c`` is
the blocks' constants concatenated, and :math:`A^T` is the blocks'
coefficients stacked vertically -- one row per cone coordinate, one
column per variable :math:`y_k`:

.. list-table::
   :header-rows: 1
   :widths: 20 18 62

   * - block
     - rows of :math:`A^T`
     - what those rows hold
   * - ``K.l = 2``
     - ``0:2``
     - coefficients of :math:`y_k` in the two LP rows
   * - ``K.q[0] = 4``
     - ``2:6``
     - norm bound first, then its three components
   * - ``K.q[1] = 4``
     - ``6:10``
     - norm bound first, then its three components
   * - ``K.s = [3]``
     - ``10:19``
     - ``vec(A_k)``, column-major, ``9 = 3**2`` entries

That totals :math:`n = 2 + 4 + 4 + 9 = 19` rows against :math:`m = 3`
columns. Since ``sedumi()`` accepts ``A`` in either orientation and
:math:`19 \neq 3`, the stacked :math:`A^T` can be passed as-is:

.. code-block:: python

   import numpy as np
   import sedumipy

   # K.l block: 2 rows
   At_l = np.array([[16.0, -14.0, 5.0],
                    [7.0, 2.0, 0.0]])
   c_l = np.array([-3.0, 5.0])

   # K.q blocks: 4 rows each, norm bound first
   At_q1 = np.array([[24.0, 7.0, -15.0],
                     [-8.0, -13.0, 12.0],
                     [8.0, -18.0, -6.0],
                     [-1.0, 3.0, 17.0]])
   c_q1 = np.array([12.0, -2.0, -14.0, -13.0])

   At_q2 = np.array([[0.0, 0.0, 0.0],
                     [-1.0, 0.0, 0.0],
                     [0.0, -1.0, 0.0],
                     [0.0, 0.0, -1.0]])
   c_q2 = np.array([10.0, 0.0, 0.0, 0.0])

   # K.s block: 9 rows, one vec() per variable
   A1 = np.array([[7.0, -5.0, 1.0], [-5.0, 1.0, -7.0], [1.0, -7.0, -4.0]])
   A2 = np.array([[3.0, 13.0, -6.0], [13.0, 12.0, -10.0], [-6.0, -10.0, -28.0]])
   A3 = np.array([[9.0, 6.0, -6.0], [6.0, -7.0, -7.0], [-6.0, -7.0, -11.0]])
   C = np.array([[68.0, -30.0, -19.0], [-30.0, 99.0, 23.0], [-19.0, 23.0, 10.0]])

   At_s = np.column_stack([M.flatten(order="F") for M in (A1, A2, A3)])
   c_s = C.flatten(order="F")

   At = np.vstack([At_l, At_q1, At_q2, At_s])        # (19, 3)
   c = np.concatenate([c_l, c_q1, c_q2, c_s])        # (19,)
   b = np.array([6.0, 4.0, 5.0])
   K = {"l": 2, "q": [4, 4], "s": [3]}

   x, y, info = sedumipy.sedumi(At, b, c, K, fid=0)

This converges (``info["numerr"] == 0``) to
``y = [-1.2209, 0.0966, 3.5775]`` with optimal value ``b @ y = 10.9485``.

Reading the answer back out
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``y`` is the model's :math:`(y_1, y_2, y_3)` directly. The slacks come
from the same ``s = c - A^T y``, sliced by the table above -- and
un-flattened, for the PSD block, by the inverse of step 2:

.. code-block:: python

   s = c - At @ y

   s_l = s[0:2]                            # K.l -- nonnegative
   q1, q2 = s[2:6], s[6:10]                # K.q -- each is (bound, *components)
   S = s[10:19].reshape(3, 3, order="F")   # K.s -- back to a 3x3 matrix

   q1[0] >= np.linalg.norm(q1[1:])         # the Lorentz condition, per block
   np.linalg.eigvalsh(S)                   # all >= 0

Here that gives ``s_l = [0, 13.35]``, both Lorentz blocks within their
norm bounds (the first one tight, ``94.288`` on both sides), and
``eigvals(S) = [0, 50.16, 165.18]`` -- so the LP block's first row, the
first cone constraint, and the PSD constraint are all active at the
optimum.

Solver options (``pars``)
--------------------------

Options are passed either as a ``pars`` dict or as individual keyword
arguments, which override the same key in ``pars`` when both are given:

.. code-block:: python

   x, y, info = sedumipy.sedumi(A, b, c, K, eps=1e-9)

Commonly used options (see ``checkpars.py`` for the full list and
defaults, ported line-for-line from ``checkpars.m``):

``eps``
   Desired accuracy (default ``1e-8``).
``bigeps``
   Accuracy still considered a usable ("numerically OK, but not to
   ``eps``") solution.
``maxiter``
   Maximum interior-point iterations (default ``150``).
``fid``
   Set to ``0`` to suppress solver progress output.

The ``info`` dict
------------------

``info`` reports the outcome: ``info["iter"]`` (iteration count),
``info["pinf"]``/``info["dinf"]`` (primal/dual infeasibility flags), and
``info["numerr"]`` (``0`` = solved to ``eps``, ``1`` = solved only to
``bigeps``, ``2`` = failed).

Reading and writing problem files
-----------------------------------

Two file formats are supported for problem/solution I/O, both under
``sedumipy`` directly:

.. code-block:: python

   # SeDuMi-style .mat problem/solution files
   At, b, c, K, pars = sedumipy.read_mat("problem.mat")
   sedumipy.write_solution_mat("solution.mat", x, y, info)

   # sparse SDPA format (.dat-s), as used by SDPLIB/DIMACS
   At, b, c, K = sedumipy.read_sdpa("problem.dat-s")
   sedumipy.write_sdpa("problem.dat-s", At, b, c, K)  # K.q/K.r not supported

See :doc:`api` for the full signatures.

Benchmarks
----------

``tests/test_benchmarks.py`` solves published `SDPLIB
<https://github.com/vsdp/SDPLIB>`_ and `DIMACS
<https://github.com/vsdp/DIMACS>`_ problems and checks the result against
each collection's own official optimal-value table:

.. code-block:: sh

   git submodule update --init --recursive   # if not already done
   .venv/bin/python -m pytest tests/test_benchmarks.py -v          # everything, ~101 problems (~10 min)
   .venv/bin/python -m pytest tests/test_benchmarks.py -v -m mini  # fastest subset only (~35s)

Measured comparison against real Octave/MEX SeDuMi
----------------------------------------------------

The table below is a **timing and optimal-value comparison of this port
against a from-source build of the real Octave/MEX SeDuMi**
(``vendor/sedumi-upstream``, built via ``install_sedumi -rebuild``),
solving all 105 problems ``tests/test_benchmarks.py`` covers (both
solvable and infeasible SDPLIB/DIMACS/TORUS instances), measured
2026-09-04 on Octave 8.4.0 / a 4-core Linux container. **Absolute
numbers are environment-dependent** (CPU/core count/thermal state) --
an earlier pass in this same environment, with the two suites run under
different background load (Octave/MEX measured right after a heavy
``mkoctfile`` rebuild, Python measured on an otherwise-idle machine),
showed Octave/MEX ~22% faster in aggregate purely from that noise. The
numbers below are from both suites run **back-to-back on an otherwise-
idle machine** (Python immediately followed by Octave/MEX, no rebuild or
install in between) to remove that skew; treat the ratio column as
indicative, not a tight bound. See :doc:`contributing` (DEVLOG's Phase 5
entry) for an earlier, 5-problem version of this same comparison.

**Timing, by problem size** (bucketed by the real Octave/MEX solve time):

.. list-table::
   :header-rows: 1

   * - Bucket
     - Problems
     - Python total
     - Octave/MEX total
     - Ratio (Python / Octave)
   * - < 2 s
     - 61
     - 51.1 s
     - 44.5 s
     - 1.15
   * - 2-20 s
     - 34
     - 261.1 s
     - 278.1 s
     - 0.94
   * - 20 s+
     - 10
     - 501.4 s
     - 478.3 s
     - 1.05
   * - **All 105**
     - 105
     - **813.6 s**
     - **800.9 s**
     - **1.02**

Under matched conditions, this port and real Octave/MEX SeDuMi run
**within a few percent of each other in aggregate** -- close enough that
the remaining gap is within this environment's own run-to-run noise, not
a clear, reproducible slowdown in either direction. Individual problems
still vary more (see the full table below): DIMACS's small ANTENNA-family
SOCP problems (``nb``, ``nb_L2``, ``nql30old``, ...) consistently show
this port running 2-3x slower, plausibly Python interpreter/NumPy-
allocation/ctypes-crossing overhead dominating on problems this small,
while several others (``qssp30``, ``hinf13``, ``truss6``, ...) run faster
than real Octave/MEX -- see the ``ratio_py_over_oct`` column.

**Optimal value (pobj) agreement.** Both solvers' final primal objective
values, formatted in scientific notation to the same number of
significant digits so they line up for a direct compare:

.. csv-table:: pobj (7 significant digits), all 105 problems
   :file: _static/benchmark_pobj_comparison.csv
   :header-rows: 1
   :widths: 10, 24, 22, 22, 12

``reldiff`` is ``|python_pobj - octave_pobj| / max(|python_pobj|,
|octave_pobj|)``. Across all 105 problems: 59 agree to a relative
difference under ``1e-8``, 81 under ``1e-6``, and 99 under ``1e-4``. The
handful of larger outliers are already-documented cases, not porting
bugs:

- ``hinf7`` (``reldiff`` ~1.5e-2): SDPLIB's own published reference value
  for the ``hinf*`` family is only given to 2-3 significant figures, and
  ``tests/test_benchmarks.py``'s tolerance for this family is widened by
  hand accordingly (see that file's ``SDPLIB_PARAMS``).
- ``qssp30old`` (``reldiff`` ~1.4e-2): the real Octave/MEX build itself
  returns ``numerr=2`` on this instance (a genuine solver limitation, not
  something this port introduced) -- see :doc:`contributing`'s DEVLOG
  reference on ``nb_L2``/``nql30old``/``qssp30old`` for the full story.

Raw per-run CSVs (:download:`Python <_static/benchmark_results_python.csv>`,
:download:`Octave/MEX <_static/benchmark_results_octave.csv>`) and the
merged :download:`comparison table <_static/benchmark_comparison.csv>`
are downloadable for anyone who wants to slice the data differently.
