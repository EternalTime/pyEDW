# pyEDW

pyEDW runs Exo-Daisy World, a stochastic generalization of Watson & Lovelock's
Daisy World tuned to M-dwarf exoplanets, with the measures that read off its
*information architecture*. Black and white daisies share a planet of habitable
fraction `f`, their albedos reining the surface temperature, while the star's
luminosity drifts as an Ornstein–Uhlenbeck process. The thermal timescale is
comparable to the stellar one, so the classic model's instantaneous-equilibrium
constraint breaks and the system becomes a genuine stochastic differential
equation.

`ExoDaisyWorld` integrates that agent (daisies) plus environment (temperature,
luminosity) SDE with a strong-order-1 stochastic Runge–Kutta step, and
`pyEDW.metrics` turns an endpoint ensemble into viability and the measures
I(A:E), ΔI, and cooperation C(a₁:a₂‖E). The model and its informational reading
come from Sowinski et al, *Planet. Sci. J.* **6**, 176 (2025), where the
derivations live. The library ports the original MATLAB classes and keeps their
interface; the `.m` files ship under `pyEDW/matlab/`.

## Installation

```
git clone https://github.com/EternalTime/pyEDW.git
cd pyEDW
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -e .
```

Requires Python 3.8+; numpy, scipy, numba, and matplotlib come along with it.

## Quick start

Seed a planet at one luminosity and integrate it:

```python
import numpy as np
from pyEDW import Parameters, ExoDaisyWorld

p = Parameters(dT=30.0)                       # growth-rate bandwidth ΔT
env = ExoDaisyWorld(p, rng=np.random.default_rng(0))
fB, fW, T, L = env.run(2000)                  # endpoint [f_B, f_W, T, L]
```

Sweep the luminosity for an ensemble, then read its information architecture:

```python
from pyEDW import metrics

Ls  = 1.0 + np.linspace(-0.7, 1.4, 400)       # mean luminosities
bio = ExoDaisyWorld.ensemble(p, Ls=Ls, nsteps=2000, N=500, seed=1)
env = ExoDaisyWorld.ensemble(p, Ls=Ls, nsteps=200,  N=500,
                             agent_free=True, seed=2)

H  = metrics.entropy_table(bio)               # per-luminosity entropy table
He = metrics.entropy_table_env(env)           # agent-free environment

V   = metrics.viability(bio, p.f)             # V = E[(f_B+f_W)/f]   (Eq. 9)
IAE = metrics.mutual_information(H)            # I(A:E)
dI  = metrics.delta_I(H, He)                  # agent-induced ΔI     (Eq. 12)
C   = metrics.cooperation(H)                  # C(a₁:a₂‖E)           (Eq. 14)
```

Regenerate the paper's figures:

```python
from pyEDW import figures
fig, ax = figures.ensemble_scatter(bio, p, Ls)   # Fig. 1
fig.savefig("fig1.png", dpi=150)
```

## Validation

`validate_against_matlab.py` reruns a Python ensemble against a `data_*.mat`
from the original `eDW_BHsim.m`, comparing the mean temperature, luminosity, and
viability curves. The pytest suite also checks the SRK1 step bit-for-bit against
an independent re-derivation of the MATLAB update; it needs the `test` extra:

```
source .venv/bin/activate
pip install -e ".[test]"
pytest
```

## Documentation

Documentation lives at <https://damiansowinski.com/pyEDW/>; its
[Getting started](https://damiansowinski.com/pyEDW/getting_started.html) page is
the authority on installation. To build it locally:

```
source .venv/bin/activate
pip install -e ".[docs]"
cd docs
make html
```

## License

MIT. See `LICENSE`.
