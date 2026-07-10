# pyEDW

pyEDW runs Exo-Daisy World: a stochastic generalization of Watson &
Lovelock's Daisy World, tuned to M-dwarf exoplanets, together with the
information-theoretic measures used to read off its *information architecture*.
A planet of habitable fraction `f` is shared by black and white daisies whose
albedos rein the surface temperature; the star's luminosity drifts as an
Ornstein–Uhlenbeck process, and because the thermal timescale is comparable to
the stellar one, the classic model's instantaneous-equilibrium constraint is
broken and the system becomes a genuine stochastic differential equation.

The library ships one environment class and a metrics module. `ExoDaisyWorld`
is the coupled agent (daisies) + environment (temperature, luminosity) SDE,
integrated with a strong-order-1 stochastic Runge–Kutta step; `pyEDW.metrics`
turns an endpoint ensemble into viability and the information measures I(A:E),
ΔI, and cooperation C(a₁:a₂‖E). The model and its informational reading are in
Sowinski, Ghoshal & Frank, *Planet. Sci. J.* **6**, 176 (2025).

This is a Python port of the original MATLAB classes, whose interface it keeps.
The original `.m` files ship under `pyEDW/matlab/` for reference.

## Installation

```
git clone https://github.com/EternalTime/pyEDW.git
cd pyEDW
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

Requires Python 3.8+; numpy, scipy, numba, and matplotlib are installed
automatically. The integrator hot loop is JIT-compiled by Numba on first call.

## Quick start

Seed a planet with a few daisies at one luminosity and integrate it:

```python
import numpy as np
from pyEDW import Parameters, ExoDaisyWorld

p = Parameters(dT=30.0)                       # growth-rate bandwidth ΔT
env = ExoDaisyWorld(p, rng=np.random.default_rng(0))
fB, fW, T, L = env.run(2000)                  # endpoint [f_B, f_W, T, L]
```

Sweep the stellar luminosity to build an ensemble, then read its information
architecture:

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

`validate_against_matlab.py` loads a `data_*.mat` produced by the original
`eDW_BHsim.m`, reruns the Python ensemble at the same bandwidth and luminosity
axis, and compares the mean temperature, luminosity, and viability curves. The
pytest suite (`pytest`) additionally checks the SRK1 step bit-for-bit against
an independent re-derivation of the MATLAB update.

## License

MIT. See `LICENSE`.
