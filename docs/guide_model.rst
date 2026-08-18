Guide: the model
================

:class:`~pyEDW.model.Parameters` holds the physical constants of Table 1 and
builds the dimensionless ``theta`` vector consumed by the integrator.
:class:`~pyEDW.model.ExoDaisyWorld` wraps the compiled SRK1
step\ :footcite:`roberts2012` behind a small stateful interface that keeps the
MATLAB usage.

Parameters
----------

.. code-block:: python

   from pyEDW import Parameters

   p = Parameters()             # paper defaults (Table 1)
   p = Parameters(dT=15.0)      # narrower growth-rate bandwidth
   theta = p.to_theta()         # 11-vector in the documented layout

``with_bandwidth`` and ``with_lambda`` return modified copies, handy when
sweeping the growth-rate bandwidth :math:`\Delta T` or the luminosity offset
:math:`\lambda`.

Stepping and running
--------------------

.. code-block:: python

   import numpy as np
   from pyEDW import ExoDaisyWorld

   env = ExoDaisyWorld(p, x0=[0.1, 0.1, 1.0, 1.0], rng=np.random.default_rng(0))
   env.evolve()                 # one SRK1 step, updates env.x
   end  = env.run(2000)         # advance in place, return endpoint
   path = env.run(2000, record=True)   # full (2001, 4) trajectory

The state is ``[f_B, f_W, T, L]``. Daisy fractions are clamped non-negative
after every step, as in the MATLAB.

Ensembles
---------

:meth:`~pyEDW.model.ExoDaisyWorld.ensemble` is a parallel (Numba ``prange``)
sweep over a luminosity axis. Each column sets :math:`\lambda = L - 1` so the
Ornstein--Uhlenbeck mean tracks the axis, samples ``N`` initial conditions the
same way ``eDW_BHsim.m`` does, integrates ``nsteps`` steps, and records the
endpoints:

.. code-block:: python

   Ls  = 1.0 + np.linspace(-0.7, 1.4, 400)
   bio = ExoDaisyWorld.ensemble(p, Ls=Ls, nsteps=2000, N=500, seed=1)

Passing an integer ``seed`` makes the result reproducible regardless of how the
parallel loop is scheduled: each luminosity column seeds its own RNG stream
with ``seed + column_index``.

Reproducing the original runs
-----------------------------

The original MATLAB drivers and figure scripts ship under ``pyEDW/matlab/`` for
reference, and ``validate_against_matlab.py`` checks a Python ensemble against a
saved ``data_*.mat``.

References
^^^^^^^^^^

.. footbibliography::
