"""Validation of the pyEDW dynamics against the MATLAB formulas.

Two kinds of check, mirroring the pyGD test suite:

1. Exact single step: with the noise switched off (``delta = 0`` makes the
   diffusion vanish) the SRK1 update is deterministic. We recompute the
   ``updateExoDaisyWorld.m`` increment independently, including the complex
   fourth-root branch, and assert machine-precision agreement.
2. Statistical: the biotic ensemble reproduces the Daisy World behavior --
   near-unit viability and a surface temperature reined toward ``T_opt`` where
   the bare planet would run hot/cold.
"""

import numpy as np

from pyEDW.model import Parameters, ExoDaisyWorld, _srk1_step


def _eom(x, th):
    """Independent re-implementation of the MATLAB ``EoM`` subfunction."""
    f, dA_G, dA_B, dA_W, Q, gammaR, kE, kS, alpha, _delta, lam = th
    fB, fW, T, L = x
    dAf = dA_B * fB + dA_W * fW
    Tg4 = T ** 4 + Q * dAf
    df = f - fB - fW
    TB = (Tg4 - Q * dA_B + 0j) ** 0.25
    TW = (Tg4 - Q * dA_W + 0j) ** 0.25
    wB = np.exp(-alpha * (TB - 1.0) ** 4)
    wW = np.exp(-alpha * (TW - 1.0) ** 4)
    a0 = df * wB * fB - gammaR * fB
    a1 = df * wW * fW - gammaR * fW
    a2 = kE * ((1.0 - dAf / dA_G) * L - T ** 4)
    a3 = kS * (1.0 + lam - L)
    return np.array([a0, a1, a2, a3])


def _step_reference(x, th, dt):
    """Deterministic (delta=0) SRK1 step recomputed from scratch."""
    a = _eom(x, th)
    k1 = np.real(a * dt)
    b = _eom(x + k1, th)
    k2 = np.real(b * dt)
    xn = x + 0.5 * (k1 + k2)
    xn[0] = max(xn[0], 0.0)
    xn[1] = max(xn[1], 0.0)
    return xn


def test_to_theta_layout():
    """theta packs the physical constants in the documented order."""
    p = Parameters()
    th = p.to_theta()
    assert th.shape == (11,)
    np.testing.assert_allclose(
        th,
        [
            0.88, 1 - 0.3, 0.1 - 0.3, 0.6 - 0.3, 0.1, 0.2 / 1.0,
            1 / (1.0 * 5.0), 1 / (1.0 * 3.0), 8 * (300.0 / 30.0) ** 4,
            0.05, 0.0,
        ],
        rtol=0, atol=1e-12,
    )


def test_single_step_exact():
    """One deterministic SRK1 step matches the hand-computed MATLAB formula."""
    p = Parameters(delta=0.0)
    th = p.to_theta()
    dt = 0.1
    rng = np.random.default_rng(3)
    for _ in range(50):
        x = np.array([
            rng.uniform(0, 0.5), rng.uniform(0, 0.5),
            rng.uniform(0.6, 1.4), rng.uniform(0.3, 2.4),
        ])
        got = _srk1_step(x, th, dt, 0.0, 1.0)
        expected = _step_reference(x, th, dt)
        np.testing.assert_allclose(got, expected, rtol=0, atol=1e-12)


def test_daisy_clamp_nonnegative():
    """Daisy fractions never go negative after a step."""
    p = Parameters()
    th = p.to_theta()
    x = np.array([1e-6, 1e-6, 1.6, 2.4])  # hot, dying daisies
    for _ in range(200):
        x = _srk1_step(x, th, 0.1, 0.0, 1.0)
        assert x[0] >= 0.0 and x[1] >= 0.0


def test_ensemble_determinism():
    """A fixed seed reproduces the ensemble bitwise."""
    p = Parameters()
    Ls = 1.0 + np.linspace(-0.5, 1.0, 8)
    d1 = ExoDaisyWorld.ensemble(p, Ls=Ls, nsteps=300, N=40, seed=7)
    d2 = ExoDaisyWorld.ensemble(p, Ls=Ls, nsteps=300, N=40, seed=7)
    np.testing.assert_array_equal(d1, d2)


def test_agent_free_stays_dead():
    """With both daisies seeded at zero the biome never appears."""
    p = Parameters()
    Ls = 1.0 + np.linspace(-0.5, 1.0, 8)
    d = ExoDaisyWorld.ensemble(p, Ls=Ls, nsteps=200, N=40, agent_free=True,
                               seed=1)
    assert np.allclose(d[:, :, :2], 0.0)


def test_rein_control_regulates_temperature():
    """The biome pulls the surface toward T_opt relative to the bare planet.

    At a luminosity where the biome thrives, the mean temperature deviation
    from optimal is smaller than the agent-free planet's deviation.
    """
    p = Parameters(dT=30.0)
    Ls = np.array([1.2])
    Ts = Ls ** 0.25
    bio = ExoDaisyWorld.ensemble(p, Ls=Ls, nsteps=2000, N=200, seed=2)
    bare = ExoDaisyWorld.ensemble(p, Ls=Ls, nsteps=400, N=200,
                                  agent_free=True, seed=3)
    dev_bio = abs(bio[0, :, 2].mean() - 1.0)
    dev_bare = abs(bare[0, :, 2].mean() - 1.0)
    assert bio[0, :, :2].mean() > 0.0        # daisies present
    assert dev_bio < dev_bare                 # regulation tightens toward T_opt
