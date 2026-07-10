"""pyEDW -- Exo-Daisy World.

A Python port of the MATLAB Exo-Daisy World (eDW) model: a stochastic
generalization of Watson & Lovelock's Daisy World tailored to M-dwarf
exoplanets, together with the information-theoretic measures used to read off
its "information architecture."

The model is a 4-dof stochastic differential equation for the black- and
white-daisy fractions (the biotic *agent*) coupled to the planetary
temperature and stellar luminosity (the *environment*). Fluctuations enter
through the luminosity as an Ornstein-Uhlenbeck process; the temperature
timescale is comparable to the stellar one, breaking the instantaneous
thermal-equilibrium constraint of the classic model.

Vertical slice
--------------
- :class:`~pyEDW.model.Parameters` -- physical parameters and the dimensionless
  ``theta`` vector.
- :class:`~pyEDW.model.ExoDaisyWorld` -- the coupled agent+environment SDE,
  integrated with a strong-order-1 stochastic Runge-Kutta step (ports
  ``updateExoDaisyWorld.m``).
- :mod:`~pyEDW.metrics` -- viability, efficacy, and the information measures
  I(A:E), delta I, and cooperation C(a1:a2||E).

Reference: D. R. Sowinski, G. Ghoshal & A. Frank, "Exo-Daisy World:
Revisiting Gaia Theory through an Informational Architecture Perspective,"
*Planet. Sci. J.* **6**, 176 (2025).
"""


def docs():
    """Open the online pyEDW documentation in a web browser."""
    import webbrowser
    webbrowser.open("https://damiansowinski.com/pyEDW/")


from pyEDW.model import Parameters, ExoDaisyWorld
from pyEDW import metrics

__all__ = ["Parameters", "ExoDaisyWorld", "metrics", "docs"]
