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
