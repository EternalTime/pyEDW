Guide: metrics
==============

:mod:`pyEDW.metrics` turns an endpoint ensemble into viability and the
information measures of Section 4 of the paper\ :footcite:`sowinski2025exo`.
The pipeline is: histogram the joint distribution
:math:`p_{AE}(f_B, f_W, T, L)` at each luminosity, tabulate the Shannon
entropies\ :footcite:`shannon1948` of every marginal, then combine them.

Entropy table
-------------

.. code-block:: python

   from pyEDW import metrics

   H  = metrics.entropy_table(bio)       # (n_lum, 16)
   He = metrics.entropy_table_env(env)   # (n_lum, 4), agent-free

Each row of ``H`` holds the fifteen marginal entropies of :math:`p_{AE}` plus
the bin count. The column layout is documented in the module and referenced by
name (``metrics.H_A1A2E1E2`` and friends) so the derived measures read like
their equations. Binning uses the square-root rule
:math:`N_{bins} = 1 + \lceil\sqrt{n}\rceil` over the living instances
(:math:`f_B f_W > 10^{-7}`); luminosities with fewer than two living instances
are ``NaN``.

Derived measures
----------------

.. code-block:: python

   V   = metrics.viability(bio, p.f)        # Eq. 9
   E   = metrics.efficacy(bio)              # <T> - 1
   IAE = metrics.mutual_information(H)       # I(A:E)
   dI  = metrics.delta_I(H, He)              # Eq. 12
   C   = metrics.cooperation(H)              # Eq. 14

All return one value per luminosity. :func:`~pyEDW.metrics.mutual_information`,
:func:`~pyEDW.metrics.cooperation`, and
:func:`~pyEDW.metrics.intra_environment_information` also accept a 2-D table of
shape ``(n_bandwidth, n_lum, 16)`` and broadcast over the leading axis, so a
full bandwidth :math:`\times` luminosity sweep is a single call.

Interpreting the signs
----------------------

:math:`I(A{:}E) \ge 0` always. :math:`\Delta I > 0` means the biome tightens the
temperature--luminosity correlation; :math:`\Delta I < 0` means it loosens it.
Negative :math:`C(a_1{:}a_2\|E)` indicates the environment synergistically
enhances the correlation between the two daisy species.

References
^^^^^^^^^^

.. footbibliography::
