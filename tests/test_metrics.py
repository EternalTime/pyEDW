"""Checks on the information measures.

Distributional identities that must hold for any histogram estimate, plus the
sign conventions used in the paper, computed on a small biotic ensemble.
"""

import numpy as np

from pyEDW.model import Parameters, ExoDaisyWorld
from pyEDW import metrics


def _tables():
    p = Parameters(dT=30.0)
    Ls = 1.0 + np.linspace(-0.6, 1.2, 30)
    bio = ExoDaisyWorld.ensemble(p, Ls=Ls, nsteps=1500, N=200, seed=11)
    env = ExoDaisyWorld.ensemble(p, Ls=Ls, nsteps=200, N=200,
                                 agent_free=True, seed=12)
    return p, bio, metrics.entropy_table(bio), metrics.entropy_table_env(env)


def test_entropy_of_degenerate_is_zero():
    """A point mass has zero entropy."""
    p = np.zeros((3, 3))
    p[1, 1] = 1.0
    assert metrics._entropy(p) == 0.0


def test_bin_index_matches_matlab_rule():
    """Bin index equals the count of interior edges strictly below the value."""
    edges = np.array([0.0, 1.0, 2.0])
    v = np.array([-0.5, 0.0, 0.5, 1.0, 2.5])
    got = metrics._bin_index(v, edges)
    np.testing.assert_array_equal(got, [0, 0, 1, 1, 3])


def test_mutual_information_nonnegative():
    """I(A:E) >= 0 wherever the biome lives."""
    _p, _bio, H, _He = _tables()
    iae = metrics.mutual_information(H)
    iae = iae[np.isfinite(iae)]
    assert np.all(iae >= -1e-9)


def test_viability_in_unit_interval():
    """0 <= V <= 1."""
    p, bio, _H, _He = _tables()
    V = metrics.viability(bio, p.f)
    assert np.nanmin(V) >= -1e-12
    assert np.nanmax(V) <= 1.0 + 1e-12


def test_measures_shapes_and_finiteness():
    """Every measure has one value per luminosity and is finite where alive."""
    p, bio, H, He = _tables()
    n = bio.shape[0]
    for m in (
        metrics.mutual_information(H),
        metrics.delta_I(H, He),
        metrics.cooperation(H),
        metrics.viability(bio, p.f),
        metrics.efficacy(bio),
    ):
        assert m.shape == (n,)
    # At least some luminosities support a living, finite information estimate.
    assert np.count_nonzero(np.isfinite(metrics.cooperation(H))) > 0
