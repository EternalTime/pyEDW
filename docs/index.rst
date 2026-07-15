pyEDW
=====

pyEDW runs Exo-Daisy World: a stochastic generalization of Watson &
Lovelock's Daisy World\ :footcite:`watson1983`, tuned to M-dwarf exoplanets,
together with the information-theoretic measures used to read off its
*information architecture*.
A planet of habitable fraction ``f`` is shared by black and white daisies whose
albedos rein the surface temperature; the star's luminosity drifts as an
Ornstein--Uhlenbeck process\ :footcite:`uhlenbeck1930`, and because the
thermal timescale is comparable to
the stellar one the classic model's instantaneous-equilibrium constraint is
broken and the system becomes a genuine stochastic differential equation.

The library ships one environment class and a metrics module.
:class:`~pyEDW.model.ExoDaisyWorld` is the coupled agent (daisies) +
environment (temperature, luminosity) SDE, integrated with a strong-order-1
stochastic Runge--Kutta step; :mod:`pyEDW.metrics` turns an endpoint ensemble
into viability and the information measures I(A:E), :math:`\Delta I`, and
cooperation :math:`C(a_1{:}a_2\|E)`.

The model and its informational reading are in Sowinski, Ghoshal & Frank,
`Exo-Daisy World: Revisiting Gaia Theory through an Informational Architecture
Perspective <https://damiansowinski.com/assets/docs/papers/ExoDaisy_Sowinski_2025.pdf>`_,
*Planet. Sci. J.* **6**, 176 (2025), `doi:10.3847/PSJ/ade310
<https://doi.org/10.3847/PSJ/ade310>`_\ :footcite:`sowinski2025exo`. If
you're new here, start with
:doc:`getting_started`, take the model for a spin in the :doc:`applet`, then
read :doc:`theory` for the equations. The library ports a set of MATLAB
classes, whose interface it keeps.

Guide
^^^^^

.. toctree::
   :maxdepth: 1

   getting_started
   applet
   guide_model
   guide_metrics
   theory

Reference
^^^^^^^^^

.. toctree::
   :maxdepth: 2

   api/pyEDW
   license

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

References
^^^^^^^^^^

.. footbibliography::
