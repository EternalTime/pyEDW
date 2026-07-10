Theory
======

The dimensionless model
-----------------------

Exo-Daisy World is a four-degree-of-freedom stochastic differential equation
for the black- and white-daisy fractions :math:`f_B, f_W` (the biotic *agent*)
coupled to the planetary temperature :math:`T` and stellar luminosity :math:`L`
(the *environment*). In the dimensionless form of Appendix B of Sowinski,
Ghoshal & Frank (2025), with temperature in units of the optimal temperature
:math:`T_{opt}` and luminosity in units of :math:`L_{opt}`,

.. math::

   \frac{df_B}{dt} &= w(T_B-1)\,(f - f_B - f_W)\,f_B - \frac{\gamma_D}{\gamma_G} f_B, \\
   \frac{df_W}{dt} &= w(T_W-1)\,(f - f_B - f_W)\,f_W - \frac{\gamma_D}{\gamma_G} f_W, \\
   \frac{dT}{dt}   &= \frac{1}{\gamma_G \tau_E}\!\left[\left(1 - \frac{\delta A_B f_B + \delta A_W f_W}{1-A_G}\right)L - T^4\right], \\
   \frac{dL}{dt}   &= \frac{1}{\gamma_G \tau_S}(1 + \lambda - L) + \sqrt{\frac{2}{\gamma_G \tau_S}}\,\delta\,(1+\lambda)\,\eta,

with the smooth growth window

.. math::

   w(x) = e^{-\alpha x^4}, \qquad \alpha = 8\left(\frac{T_{opt}}{\Delta T}\right)^4,

and the scaled daisy temperatures

.. math::

   T_\alpha^4 = T^4 + Q\,\delta A_B f_B + Q\,\delta A_W f_W - Q\,\delta A_\alpha,
   \qquad \delta A_\alpha = A_\alpha - A_G.

The luminosity is the only noisy dof, an Ornstein--Uhlenbeck process with mean
:math:`1+\lambda`, so the driving Wiener process is scalar. This is what lets
:class:`~pyEDW.model.ExoDaisyWorld` use the strong-order-1 stochastic
Runge--Kutta scheme of A. Roberts (2012); the step is a line-for-line port of
``updateExoDaisyWorld.m``.

Rein control
------------

Because :math:`A_B < A_G < A_W`, black daisies warm their surroundings and
white daisies cool them. As the star brightens, the population mix shifts from
black to white, adjusting the planetary albedo :math:`A = A_G + \sum_\alpha
(A_\alpha - A_G) f_\alpha` so that the surface temperature is held near
:math:`T_{opt}` across a wide band of luminosities -- the Gaian *rein control*
mechanism.

The information architecture
----------------------------

Partitioning the system into agent :math:`a = (f_B, f_W)` and environment
:math:`e = (T, L)`, the coevolution is characterized by information measures
estimated from the joint distribution :math:`p_{AE}` (Section 4):

- **Viability** :math:`V = \mathbb{E}^A[(f_B+f_W)/f]` (Eq. 9), the expected
  occupied fraction of the habitable area.
- **Mutual information** :math:`I(A{:}E)`, the agent--environment correlation.
- **Correlation change** :math:`\Delta I_{\varnothing\to A} = I(e_1{:}e_2) -
  I_0(e_1{:}e_2)` (Eq. 12), how much the biome strengthens or weakens the
  intrinsic temperature--luminosity correlation.
- **Cooperation** :math:`C(a_1{:}a_2\|E)` (Eq. 14), the interaction information
  among the two species and the environment; negative values signal a
  synergy mediated by the environment.

All entropies are estimated by histogramming with the square-root binning rule
:math:`N_{bins} = 1 + \lceil\sqrt{n}\rceil`, matching Appendix B.3.
