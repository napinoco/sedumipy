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

where :math:`K^*` is the dual cone of :math:`K`. ``sedumipy.sedumi``
returns the primal optimum ``x`` and the dual optimum ``y``; the dual
slack ``s = c - A^T y`` is not returned directly, but is one line to
recover (see the worked examples below). Every cone ``sedumipy`` supports
is self-dual (:math:`K^* = K`) except the free block ``K.f``, whose dual
is :math:`\{0\}` -- i.e. the corresponding dual-slack entries are forced
to zero. This primal-dual pair, and the field-by-field meaning of ``K``
below, follow real SeDuMi's own convention -- see the `Addendum to the
SeDuMi User Guide
<https://sedumi.ie.lehigh.edu/sedumi/files/sedumi-downloads/SeDuMi_Guide_11.pdf>`_
(Pólik, 2005), Section 2, for the original reference.

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

Combining cone blocks in one problem
--------------------------------------

A single problem can mix any of the block types above -- ``x`` is just one
vector partitioned ``[x_f | x_l | x_q | x_r | x_s]`` in that fixed order,
and the columns of ``A`` (and entries of ``c``) follow the same layout.
The example below uses all of ``K.f``, ``K.l``, ``K.q``, and ``K.s`` in a
single solve, on a problem small enough to check by hand:

.. code-block:: python

   import numpy as np
   import sedumipy

   # x = [ x1 (free) | x2, x3 (l >= 0) | x4, x5, x6 (SOC) | s1, s2, s3, s4 (2x2 PSD) ]
   n = 1 + 2 + 3 + 4
   m = 4

   A = np.zeros((m, n))
   A[0, 0] = 1.0                    # x1 = 1
   A[1, 1] = 1.0; A[1, 2] = 1.0     # x2 + x3 = 2
   A[2, 3] = 1.0                    # x4 = 2   (the SOC block's norm bound)
   A[3, 6] = 1.0; A[3, 9] = 1.0     # s1 + s4 = 2  (trace of the 2x2 PSD block)
   b = np.array([1.0, 2.0, 2.0, 2.0])

   c = np.zeros(n)
   c[0] = 1.0                        # free-variable cost
   c[1], c[2] = 1.0, 2.0             # l costs -- picks x2 = 2, x3 = 0
   c[4] = -1.0                       # q cost -- pushes x5 to the SOC boundary
   c[7], c[8] = -0.5, -0.5           # s cost -- maximizes the PSD block's off-diagonal

   K = {"f": 1, "l": 2, "q": [3], "s": [2]}
   x, y, info = sedumipy.sedumi(A, b, c, K, fid=0)

which returns ``x = [1, 2, 0, 2, 2, 0, 1, 1, 1, 1]``. The ``K.s`` block
(``x[6:10]``, reshaped column-major as ``mat(x[6:10], 2)``) is
``[[1, 1], [1, 1]]`` -- PSD, with eigenvalues ``[0, 2]``, sitting exactly
on the boundary the trace-2 constraint allows.

Worked example: LP, SOCP, and SDP together (from the literature)
--------------------------------------------------------------------

Example 5 of Ito, *A Study on the Algorithm and Implementation of SDPT3*
(`arXiv:2512.24623 <https://arxiv.org/abs/2512.24623>`_), is a realistic
mixed-cone problem -- one ``K.l`` block, two ``K.q`` blocks, and a 3x3
``K.s`` block -- given in SDPT3's own ``[blk, At, C, b]`` input format:

.. math::

   \begin{aligned}
   \max_{y \in \mathbb{R}^3} \quad & 6y_1 + 4y_2 + 5y_3 \\[4pt]
   \text{s.t.} \quad
   & 16y_1 - 14y_2 + 5y_3 \le -3, \\
   & 7y_1 + 2y_2 \le 5, \\[4pt]
   & \left\|
       \begin{pmatrix} 8y_1 + 13y_2 - 12y_3 - 2 \\
                        -8y_1 + 18y_2 + 6y_3 - 14 \\
                        y_1 - 3y_2 - 17y_3 - 13 \end{pmatrix}
     \right\| \le -24y_1 - 7y_2 + 15y_3 + 12, \\[4pt]
   & \left\| \begin{pmatrix} y_1 \\ y_2 \\ y_3 \end{pmatrix} \right\| \le 10, \\[4pt]
   & \begin{pmatrix}
       7y_1+3y_2+9y_3 & -5y_1+13y_2+6y_3 & y_1-6y_2-6y_3 \\
       -5y_1+13y_2+6y_3 & y_1+12y_2-7y_3 & -7y_1-10y_2-7y_3 \\
       y_1-6y_2-6y_3 & -7y_1-10y_2-7y_3 & -4y_1-28y_2-11y_3
     \end{pmatrix}
     \preceq
     \begin{pmatrix} 68 & -30 & -19 \\ -30 & 99 & 23 \\ -19 & 23 & 10 \end{pmatrix}
   \end{aligned}

This is a maximization over ``y`` subject to cone-membership constraints
on affine expressions of ``y`` -- exactly SeDuMi's dual problem (D) above,
and every one of its cones (``K.l``, ``K.q``, ``K.s``) is self-dual. SDPT3's
``At{p}`` cell array is literally :math:`(A^p)^T` -- one column per
constraint -- so stacking SDPT3's ``At`` blocks vertically gives
``sedumipy``'s ``A`` directly, passed in its :math:`n \times m` transposed
form (which ``sedumi`` auto-detects from its shape). The only real
translation needed is ``svec`` -> full ``vec`` for the PSD block:
``sedumipy`` (like real SeDuMi) stores ``K.s`` blocks as the full
``n**2``-entry column-major matrix, not SDPT3's ``n(n+1)/2``-entry
``svec``:

.. code-block:: python

   import numpy as np
   import sedumipy

   # K.l block (n1 = 2)
   At1 = np.array([[16.0, -14.0, 5.0], [7.0, 2.0, 0.0]])
   c1 = np.array([-3.0, 5.0])

   # K.q block 1 (n2 = 4)
   At2 = np.array([[24.0, 7.0, -15.0], [-8.0, -13.0, 12.0],
                    [8.0, -18.0, -6.0], [-1.0, 3.0, 17.0]])
   c2 = np.array([12.0, -2.0, -14.0, -13.0])

   # K.q block 2 (n3 = 4)
   At3 = np.array([[0.0, 0.0, 0.0], [-1.0, 0.0, 0.0],
                    [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]])
   c3 = np.array([10.0, 0.0, 0.0, 0.0])

   # K.s block (n4 = 3, i.e. 9 vec entries); SeDuMi wants vec, not svec
   A1_sdp = np.array([[7.0, -5.0, 1.0], [-5.0, 1.0, -7.0], [1.0, -7.0, -4.0]])
   A2_sdp = np.array([[3.0, 13.0, -6.0], [13.0, 12.0, -10.0], [-6.0, -10.0, -28.0]])
   A3_sdp = np.array([[9.0, 6.0, -6.0], [6.0, -7.0, -7.0], [-6.0, -7.0, -11.0]])
   C4 = np.array([[68.0, -30.0, -19.0], [-30.0, 99.0, 23.0], [-19.0, 23.0, 10.0]])
   At4 = np.column_stack([m.flatten(order="F") for m in (A1_sdp, A2_sdp, A3_sdp)])
   c4 = C4.flatten(order="F")

   A = np.vstack([At1, At2, At3, At4])  # n x m -- sedumi auto-detects the transpose
   c = np.concatenate([c1, c2, c3, c4])
   b = np.array([6.0, 4.0, 5.0])
   K = {"l": 2, "q": [4, 4], "s": [3]}

   x, y, info = sedumipy.sedumi(A, b, c, K, fid=0)
   print(y, b @ y, info["numerr"])   # y is (y1, y2, y3) above; b @ y is the optimal value

which converges (``info["numerr"] == 0``) to dual objective
``b @ y ≈ 10.9485`` at ``y ≈ [-1.2209, 0.0966, 3.5775]``. Reconstructing
the dual slack ``s = c - A @ y`` block by block confirms each piece lands
in its cone: ``s[:2] >= 0`` (the ``K.l`` block), both ``K.q`` blocks
satisfy ``s[0] >= norm(s[1:])``, and the 3x3 ``K.s`` block's eigenvalues
are all non-negative.

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
