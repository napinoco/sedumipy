sedumipy
========

**sedumipy** is a MATLAB/Octave-free port of `SeDuMi
<https://github.com/sqlp/sedumi>`_, an interior-point solver for linear
(LP), second-order cone (SOCP), and semidefinite (SDP) programs, to a
standalone C library plus a Python (NumPy/SciPy) package. No MATLAB or
GNU Octave runtime is required.

It solves conic programs over a symmetric cone :math:`K`, returning both
sides of the primal-dual pair from one call:

.. math::

   \begin{aligned}
   \text{(P)} \quad \min_{x} \;\; & c^T x
     & \qquad \text{(D)} \quad \max_{y,\,s} \;\; & b^T y \\
   \text{s.t.} \;\; & Ax = b, \;\; x \in K
     & \qquad \text{s.t.} \;\; & A^T y + s = c, \;\; s \in K^{*}
   \end{aligned}

:math:`K` is a product of blocks of five kinds -- free, nonnegative
orthant, second-order (Lorentz), rotated Lorentz, and positive
semidefinite -- so LP, SOCP, and SDP are the cases where every block is
of one kind, and a single problem may mix them freely. See :doc:`usage`
for each cone's definition and for a worked mixed-cone example.

To solve, say,

.. math::

   \min_x\ 7x_1 + 4x_2 + 10x_3 \quad \text{such that}\quad
   \begin{aligned}
   3x_1 + x_2 + 2x_3 &= 9, \\
   x_1 + 2x_2 + 4x_3 &= 8,
   \end{aligned}
   \quad x_1, x_2, x_3 \ge 0

feed the objective coefficients as ``c``, the constraint rows as ``A``,
their right-hand sides as ``b``, and say that all three variables are
nonnegative with ``K``:

.. code-block:: python

   import numpy as np
   import sedumipy

   A = np.array([[3.0, 1.0, 2.0],      # one row per equality constraint
                 [1.0, 2.0, 4.0]])
   b = np.array([9.0, 8.0])            # their right-hand sides
   c = np.array([7.0, 4.0, 10.0])      # objective coefficients
   K = {"l": 3}                        # all 3 variables are >= 0

   x, y, info = sedumipy.sedumi(A, b, c, K)
   # x = [2., 3., 0.]   the optimum, costing c @ x = 26

This mirrors original SeDuMi's own ``[x, y, info] = sedumi(A, b, c, K)``
call signature and semantics -- see :doc:`usage` for the problem format
(the ``K`` cone-structure dict, and how to mix LP, second-order-cone, and
PSD blocks in one problem) and :doc:`api` for the full public API.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   installation
   usage
   api
   internals
   status
   contributing

Project status
--------------

LP and SOCP problems are fully ported and verified against
original Octave/SeDuMi output; PSD (``K.s``) cones are also implemented and
verified against real reference solves. See :doc:`status` for the current
scope and :doc:`contributing` for the full phase-by-phase project history.

A note on citation and attribution
-----------------------------------

sedumipy is an independent, unofficial re-implementation of SeDuMi,
created without involvement from the original SeDuMi authors or
maintainers. Although it aims to reproduce SeDuMi's numerical behavior
faithfully, it is a from-scratch port and may differ from the original in
ways not yet identified -- in short, **it may not always behave
identically to the original SeDuMi.**

If you use this software in research, please cite the original SeDuMi
paper to give credit where it is due:

    Sturm, J.F. (1999). Using SeDuMi 1.02, a MATLAB toolbox for
    optimization over symmetric cones. *Optimization Methods and
    Software*, 11(1-4), 625-653.

At the same time, please make clear in your own work that results were
produced using **sedumipy, an unofficial Python port** -- not the
original SeDuMi. Any discrepancy, bug, or unexpected numerical behavior
you observe here is a property of this port, not of SeDuMi itself, and
should not be attributed to the original project or its authors.

License
-------

SeDuMi is licensed under the GNU General Public License v2 (see
`LICENSE <https://github.com/napinoco/sedumipy/blob/main/LICENSE>`_);
this port, being a derivative work, is licensed the same way.

Indices
-------

* :ref:`genindex`
* :ref:`modindex`
