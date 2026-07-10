Getting started
===============

Installation
------------

.. code-block:: bash

   git clone https://github.com/EternalTime/pyEDW.git
   cd pyEDW
   python3 -m venv .venv
   source .venv/bin/activate
   pip install -e .

Requires Python 3.8+. ``numpy``, ``scipy``, ``numba``, and ``matplotlib`` are
installed automatically. The integrator hot loop is JIT-compiled by Numba on
first call, so the first ``run`` or ``ensemble`` in a session pays a one-time
compilation cost.

A single planet
---------------

Build the parameters, seed a planet, and integrate it forward:

.. code-block:: python

   import numpy as np
   from pyEDW import Parameters, ExoDaisyWorld

   p = Parameters(dT=30.0)                    # growth-rate bandwidth ΔT
   env = ExoDaisyWorld(p, rng=np.random.default_rng(0))
   fB, fW, T, L = env.run(2000)               # endpoint [f_B, f_W, T, L]

Pass ``record=True`` to keep the whole trajectory:

.. code-block:: python

   path = env.run(2000, record=True)          # shape (2001, 4)

An ensemble over luminosity
---------------------------

The paper's experiments sweep the mean stellar luminosity and run many
independent instances per luminosity. :meth:`~pyEDW.model.ExoDaisyWorld.ensemble`
does this in parallel:

.. code-block:: python

   Ls  = 1.0 + np.linspace(-0.7, 1.4, 400)
   bio = ExoDaisyWorld.ensemble(p, Ls=Ls, nsteps=2000, N=500, seed=1)
   # bio.shape == (400, 500, 4)

Set ``agent_free=True`` to seed both daisy fractions at zero and integrate the
bare planet, the reference used for the agent-induced correlation change
:math:`\Delta I`.

Reading the information architecture
------------------------------------

.. code-block:: python

   from pyEDW import metrics

   env = ExoDaisyWorld.ensemble(p, Ls=Ls, nsteps=200, N=500,
                                agent_free=True, seed=2)
   H  = metrics.entropy_table(bio)
   He = metrics.entropy_table_env(env)

   V   = metrics.viability(bio, p.f)
   IAE = metrics.mutual_information(H)
   dI  = metrics.delta_I(H, He)
   C   = metrics.cooperation(H)

See :doc:`guide_metrics` for what each measure means and :doc:`theory` for the
equations behind them.
