Getting started
===============

Installation
------------

.. code-block:: bash

   git clone https://github.com/EternalTime/pyEDW.git
   cd pyEDW
   python3 -m venv .venv
   source .venv/bin/activate
   python -m pip install --upgrade pip
   pip install -e .

Requires Python 3.8+; ``numpy``, ``scipy``, ``numba``, and ``matplotlib`` come
along with it. Numba JIT-compiles the integrator hot loop, so the first ``run``
or ``ensemble`` in a session pays a one-time cost.

A single planet
---------------

Seed a planet and integrate it forward:

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

The paper sweeps the mean stellar luminosity with many independent instances at
each. :meth:`~pyEDW.model.ExoDaisyWorld.ensemble` does this in parallel:

.. code-block:: python

   Ls  = 1.0 + np.linspace(-0.7, 1.4, 400)
   bio = ExoDaisyWorld.ensemble(p, Ls=Ls, nsteps=2000, N=500, seed=1)
   # bio.shape == (400, 500, 4)

``agent_free=True`` seeds both daisy fractions at zero and integrates the bare
planet, the reference for the agent-induced correlation change
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

See :doc:`guide_metrics` for what each measure means, and :doc:`theory` for the
equations.

Running the tests
-----------------

The suite needs ``pytest``, which ships in the ``test`` extra:

.. code-block:: bash

   source .venv/bin/activate
   pip install -e ".[test]"
   pytest

Building the documentation
--------------------------

Sphinx builds this site from ``docs/``, using the ``docs`` extra:

.. code-block:: bash

   source .venv/bin/activate
   pip install -e ".[docs]"
   cd docs
   make html

The HTML lands in ``docs/_build/html``.
