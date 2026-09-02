sedumipy
========

**sedumipy** is a MATLAB/Octave-free port of `SeDuMi
<https://github.com/sqlp/sedumi>`_, an interior-point solver for linear
(LP), second-order cone (SOCP), and semidefinite (SDP) programs, to a
standalone C library plus a Python (NumPy/SciPy) package. No MATLAB or
GNU Octave runtime is required.

.. code-block:: python

   import numpy as np
   import sedumipy

   # minimize c'x s.t. Ax = b, x >= 0
   A = np.eye(2)
   b = np.array([1.0, 1.0])
   c = np.array([1.0, 1.0])
   x, y, info = sedumipy.sedumi(A, b, c, {"l": 2})

This mirrors real SeDuMi's own ``[x, y, info] = sedumi(A, b, c, K)`` call
signature and semantics -- see :doc:`usage` for the problem format (the
``K`` cone-structure dict) and :doc:`api` for the full public API.

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

LP and SOCP problems are fully ported and verified bit-for-bit against
real Octave/SeDuMi output; PSD (``K.s``) cones are also implemented and
verified against real reference solves. See :doc:`status` for the current
scope and :doc:`contributing` for the full phase-by-phase project history.

License
-------

SeDuMi is licensed under the GNU General Public License v2 (see
`LICENSE <https://github.com/napinoco/sedumipy/blob/main/LICENSE>`_);
this port, being a derivative work, is licensed the same way.

Indices
-------

* :ref:`genindex`
* :ref:`modindex`
