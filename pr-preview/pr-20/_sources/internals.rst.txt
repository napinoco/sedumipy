Internals
=========

Everything on this page is a direct, module-per-``.m``-file port of
original SeDuMi's own internals (interior-point iteration, cone math,
symbolic/numeric factorization, ...) and is not a stable public API --
this exists for contributors reading or extending the port, not for end
users. Each module's own docstring names the ``.m``/``.c`` file it ports
and any deliberate deviations. See :doc:`contributing` for how these
modules fit together into :func:`sedumipy.sedumi.sedumi`'s main loop.

.. autosummary::
   :toctree: generated

   sedumipy.amul
   sedumipy.cone
   sedumipy.deninfac
   sedumipy.getada
   sedumipy.getada_psd
   sedumipy.getdatm
   sedumipy.getdense
   sedumipy.getsymbada
   sedumipy.incorder
   sedumipy.maxstep
   sedumipy.neighborhood
   sedumipy.optstep
   sedumipy.pcg
   sedumipy.posttransfo
   sedumipy.pretransfo
   sedumipy.sddir
   sedumipy.sdfactor
   sedumipy.sdinit
   sedumipy.stepdif
   sedumipy.symbchol
   sedumipy.symbcholden
   sedumipy.trydif
   sedumipy.updtransfo
   sedumipy.widelen
   sedumipy.wregion
