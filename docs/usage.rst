Usage
=====

The public entry point is :func:`sedumipy.sedumi`, matching real SeDuMi's
own signature:

.. code-block:: python

   x, y, info = sedumipy.sedumi(A, b, c, K, pars=None, **pars_kwargs)

which solves

.. math::

   \text{minimize } c^T x \quad \text{such that } Ax = b,\ x \in K

and its dual. ``A`` may be given either as :math:`m \times n` or its
transpose :math:`n \times m` (SeDuMi disambiguates by shape, matching
real SeDuMi); dense NumPy arrays and SciPy sparse matrices are both
accepted.

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
numbers are environment-dependent** (CPU/core count/thermal state --
back-to-back runs in this same environment varied by up to ~25% on the
largest problems), so treat the ratio column as indicative, not a tight
bound; see :doc:`contributing` (DEVLOG's Phase 5 entry) for an earlier,
5-problem version of this same comparison.

**Timing, by problem size** (bucketed by the real Octave/MEX solve time):

.. list-table::
   :header-rows: 1

   * - Bucket
     - Problems
     - Python total
     - Octave/MEX total
     - Ratio (Python / Octave)
   * - < 2 s
     - 59
     - 62.7 s
     - 41.8 s
     - 1.50
   * - 2-20 s
     - 38
     - 373.0 s
     - 321.1 s
     - 1.16
   * - 20 s+
     - 8
     - 547.5 s
     - 442.4 s
     - 1.24
   * - **All 105**
     - 105
     - **983.1 s**
     - **805.3 s**
     - **1.22**

Small problems (Python interpreter/NumPy-allocation/ctypes-crossing
overhead per iteration) still show this port running visibly slower, as
in the 5-problem Phase 5 measurement; on larger problems, where the
native C-kernel Cholesky factorization dominates total time, the two are
much closer (within ~20-25%, the same order as this environment's own
run-to-run noise).

**Optimal value (pobj) agreement.** Both solvers' final primal objective
values, formatted in scientific notation to the same number of
significant digits so they line up for a direct compare:

.. csv-table:: pobj (7 significant digits) and timing, all 105 problems
   :file: _static/benchmark_comparison.csv
   :header-rows: 1
   :widths: 8, 16, 12, 12, 8, 6, 6, 6, 6, 8, 8, 8

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
