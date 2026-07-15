Applet
======

.. raw:: html

   <link rel="stylesheet" href="_static/edw_applet.css">
   <div id="edw-app"></div>
   <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
   <script src="_static/edw_applet.js"></script>

Exo-Daisy World, running live in your browser. The physics is the same one
pyEDW integrates in Python: the dimensionless equations of motion of
:doc:`theory`, stepped with the same strong-order-1 stochastic Runge--Kutta
scheme. Nothing is precomputed — the daisies are the state of a stochastic
differential equation being solved as you watch. Drag the planet to spin it.

The applet is yours to take apart. It ships with the library as
`docs/_static/edw_applet.js
<https://github.com/EternalTime/pyEDW/blob/main/docs/_static/edw_applet.js>`_:
535 lines with no build step and no dependency beyond three.js, of which the
first hundred or so are the model itself — the growth window, the equations of
motion, and the stochastic Runge--Kutta step, ported line for line from
:mod:`pyEDW.model`. The rest is scenery.

Brighten the star and the black daisies, which warm the ground beneath them,
give way to white ones, which cool it. The shifting mix raises the planetary
albedo and reins the surface temperature back toward :math:`T_{opt}` — this is
the Gaian feedback\ :footcite:`lovelock1974` the model was built to exhibit. In the plot, the green track
pulls away from the bare-planet Stefan--Boltzmann curve exactly where the biome
is alive; drive the luminosity far enough in either direction and the biome
collapses, the track snapping back onto the curve of a dead world.

The growth bandwidth :math:`\Delta T / T_{opt}` sets how fussy the daisies are
about temperature. Narrow it and the habitable band of luminosities shrinks to
a sliver; widen it and the biome clings on across the whole range.

Move the luminosity slowly and you will find the model is bistable — an
established biome survives at luminosities where a freshly seeded one dies out,
because a large standing population can hold the temperature where it needs it
while a sparse one cannot. This is why the ensembles in the
paper\ :footcite:`sowinski2025exo`, which sample a broad range of initial
daisy fractions, show life across a much wider band of luminosities than a
single small seeding would suggest. One liberty is taken
here: since the growth terms are multiplicative in :math:`f_B` and :math:`f_W`,
zero is an absorbing state, and a biome that collapsed could never return. The
applet therefore keeps a dormant seed bank of :math:`10^{-4}` — too small to
draw or to perturb a living biome, but enough to recolonize a planet that
becomes habitable again.

References
^^^^^^^^^^

.. footbibliography::
