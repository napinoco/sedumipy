API reference
=============

Top-level package
------------------

.. automodule:: sedumipy
   :members:
   :undoc-members:

Solver driver
-------------

.. autofunction:: sedumipy.sedumi.sedumi

File I/O
--------

.. automodule:: sedumipy.matio
   :members:

.. automodule:: sedumipy.sdpa
   :members: read_sdpa, write_sdpa

Solver options
--------------

.. automodule:: sedumipy.checkpars
   :members:

Internals
---------

Everything documented on the :doc:`internals` page is a direct,
module-per-``.m``-file port of SeDuMi's own internals (interior-point
iteration, cone math, symbolic/numeric factorization, ...) and is not a
stable public API -- it is there for contributors reading or extending
the port, not for end users. See :doc:`contributing` for how these
modules fit together.
