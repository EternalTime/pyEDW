"""The Exo-Daisy World model and its stochastic integrator.

This module ports ``updateExoDaisyWorld.m`` (the single time step) and
``eDW_BHsim.m`` / ``exoDaisyWorld.m`` (the ensemble drivers) to Python. The
dynamics are the dimensionless system of Appendix B of Sowinski, Ghoshal &
Frank, *Planet. Sci. J.* **6**, 176 (2025):

- Eq. (B8)  ``d f_B / dt = w(T_B - 1)(f - f_B - f_W) f_B - (gamma_D/gamma_G) f_B``
- Eq. (B9)  ``d f_W / dt = w(T_W - 1)(f - f_B - f_W) f_W - (gamma_D/gamma_G) f_W``
- Eq. (B10) ``d T   / dt = (1/gamma_G tau_E)[(1 - (dA_B f_B + dA_W f_W)/(1-A_G)) L - T^4]``
- Eq. (B11) ``d L   / dt = (1/gamma_G tau_S)(1 + lambda - L) + sqrt(2/gamma_G tau_S) delta (1+lambda) eta``

with growth window ``w(x) = exp(-alpha x^4)`` (Eq. B7), scaled daisy
temperatures ``T_alpha^4 = T^4 + Q dA_B f_B + Q dA_W f_W - Q dA_alpha`` (Eq.
B12), and ``dA_alpha = A_alpha - A_G``. Temperatures are in units of the
optimal temperature ``T_opt`` and luminosity in units of ``L_opt``.

The noise is a single scalar Wiener process entering only the luminosity, so
the SDE is scalar-driven and the strong-order-1 stochastic Runge-Kutta scheme
of A. Roberts (2012) applies. :class:`Parameters` packs the physical constants
into the dimensionless ``theta`` vector consumed by the integrator, whose
layout mirrors the MATLAB comment block exactly:

===  ============================
idx  theta component
===  ============================
0    f
1    1 - A_G
2    A_B - A_G
3    A_W - A_G
4    Q
5    gamma_D / gamma_G
6    1 / (gamma_G tau_E)
7    1 / (gamma_G tau_S)
8    8 (T_opt / Delta T)^4
9    delta
10   lambda
===  ============================
"""

from dataclasses import dataclass, replace

import numpy as np
from numba import njit, prange

# State vector layout: x = [f_B, f_W, T, L].
_FB, _FW, _T, _L = 0, 1, 2, 3


@dataclass
class Parameters:
    """Physical parameters of Exo-Daisy World.

    All temperatures are expressed in units of the optimal temperature
    ``T_opt`` and luminosities in units of ``L_opt`` once dimensionalized; the
    fields below are the dimensional constants of Table 1, from which
    :meth:`to_theta` builds the dimensionless ``theta`` vector.

    Parameters
    ----------
    f : float
        Habitable land fraction.
    A_G, A_B, A_W : float
        Ground, black-daisy, and white-daisy albedos (``A_B < A_G < A_W``).
    T_opt : float
        Optimal growth temperature (sets the temperature scale).
    Q : float
        Temperature-shift heat constant coupling daisy albedo to local
        temperature.
    gamma_G, gamma_D : float
        Maximal growth rate and daisy decay rate.
    tau_S, tau_E : float
        Stellar-fluctuation and environmental (thermal) timescales, in units
        of ``1 / gamma_G``.
    dT : float
        Growth-rate bandwidth ``Delta T`` (range of bearable temperatures).
    delta : float
        Relative amplitude of the stellar-luminosity fluctuations.
    lam : float
        Luminosity offset ``lambda``; the mean luminosity is
        ``<L> = (1 + lambda) L_opt``.
    """

    f: float = 0.88
    A_G: float = 0.3
    A_B: float = 0.1
    A_W: float = 0.6
    T_opt: float = 300.0
    Q: float = 0.1
    gamma_G: float = 1.0
    gamma_D: float = 0.2
    tau_S: float = 3.0
    tau_E: float = 5.0
    dT: float = 30.0
    delta: float = 0.05
    lam: float = 0.0

    def to_theta(self):
        """Return the 11-element dimensionless ``theta`` vector.

        Returns
        -------
        numpy.ndarray
            ``theta`` in the layout documented in the module docstring,
            matching ``eDW_BHsim.m``.
        """
        return np.array(
            [
                self.f,
                1.0 - self.A_G,
                self.A_B - self.A_G,
                self.A_W - self.A_G,
                self.Q,
                self.gamma_D / self.gamma_G,
                1.0 / (self.gamma_G * self.tau_E),
                1.0 / (self.gamma_G * self.tau_S),
                8.0 * (self.T_opt / self.dT) ** 4,
                self.delta,
                self.lam,
            ],
            dtype=np.float64,
        )

    def with_bandwidth(self, dT):
        """Return a copy with a new growth-rate bandwidth ``Delta T``."""
        return replace(self, dT=dT)

    def with_lambda(self, lam):
        """Return a copy with a new luminosity offset ``lambda``."""
        return replace(self, lam=lam)


# ---------------------------------------------------------------------------
# Compiled core (ports updateExoDaisyWorld.m).
# ---------------------------------------------------------------------------


@njit(cache=True)
def _srk1_step(x, theta, dt, z, s):
    """One strong-order-1 stochastic Runge-Kutta step, in place semantics.

    Faithful port of ``updateExoDaisyWorld.m``: ``z`` is a standard normal and
    ``s`` is +/-1. The two daisy equations can involve a fourth root of a
    negative argument at low temperature; MATLAB carries this as a complex
    number and takes the real part of the increment, which we reproduce with
    complex arithmetic. Daisy fractions are clamped non-negative after the
    step.
    """
    f = theta[0]
    dA_G = theta[1]           # 1 - A_G
    dA_B = theta[2]           # A_B - A_G
    dA_W = theta[3]           # A_W - A_G
    Q = theta[4]
    gammaR = theta[5]         # gamma_D / gamma_G
    kE = theta[6]             # 1 / (gamma_G tau_E)
    kS = theta[7]             # 1 / (gamma_G tau_S)
    alpha = theta[8]          # 8 (T_opt / dT)^4
    delta = theta[9]
    lam = theta[10]

    bL = np.sqrt(2.0 * kS) * delta * (1.0 + lam)

    # --- drift/diffusion at x ------------------------------------------------
    a0, a1, a2, a3 = _eom(x, f, dA_G, dA_B, dA_W, Q, gammaR, kE, kS, alpha, lam)
    k1_0 = (a0 * dt).real
    k1_1 = (a1 * dt).real
    k1_2 = a2 * dt
    k1_3 = a3 * dt + (z - s) * bL * np.sqrt(dt)

    # --- drift/diffusion at x + k1 ------------------------------------------
    x1_0 = x[0] + k1_0
    x1_1 = x[1] + k1_1
    x1_2 = x[2] + k1_2
    x1_3 = x[3] + k1_3
    x1 = np.array([x1_0, x1_1, x1_2, x1_3])
    b0, b1, b2, b3 = _eom(x1, f, dA_G, dA_B, dA_W, Q, gammaR, kE, kS, alpha, lam)
    k2_0 = (b0 * dt).real
    k2_1 = (b1 * dt).real
    k2_2 = b2 * dt
    k2_3 = b3 * dt + (z + s) * bL * np.sqrt(dt)

    out = np.empty(4)
    out[0] = x[0] + 0.5 * (k1_0 + k2_0)
    out[1] = x[1] + 0.5 * (k1_1 + k2_1)
    out[2] = x[2] + 0.5 * (k1_2 + k2_2)
    out[3] = x[3] + 0.5 * (k1_3 + k2_3)
    if out[0] < 0.0:
        out[0] = 0.0
    if out[1] < 0.0:
        out[1] = 0.0
    return out


@njit(cache=True)
def _eom(x, f, dA_G, dA_B, dA_W, Q, gammaR, kE, kS, alpha, lam):
    """Deterministic drift of the four dof (ports the ``EoM`` subfunction).

    Returns the drift components; the first two are complex (real part is the
    physical increment), the last two are real.
    """
    fB = x[0]
    fW = x[1]
    T = x[2]
    L = x[3]

    dAf = dA_B * fB + dA_W * fW          # dA_B f_B + dA_W f_W
    Tg4 = T ** 4 + Q * dAf               # T^4 + Q dA_f
    df = f - fB - fW                     # f - f_B - f_W

    # Scaled daisy temperatures T_alpha = (Tg4 - Q dA_alpha)^{1/4}.
    baseB = Tg4 - Q * dA_B + 0j
    baseW = Tg4 - Q * dA_W + 0j
    TB = np.exp(0.25 * np.log(baseB))
    TW = np.exp(0.25 * np.log(baseW))

    wB = np.exp(-alpha * (TB - 1.0) ** 4)
    wW = np.exp(-alpha * (TW - 1.0) ** 4)

    a0 = df * wB * fB - gammaR * fB
    a1 = df * wW * fW - gammaR * fW
    a2 = kE * ((1.0 - dAf / dA_G) * L - T ** 4)
    a3 = kS * (1.0 + lam - L)
    return a0, a1, a2, a3


@njit(cache=True)
def _run_endpoint(x0, theta, dt, nsteps):
    """Integrate one trajectory and return the ``(4,)`` endpoint."""
    x = x0.copy()
    for _ in range(nsteps):
        z = np.random.standard_normal()
        s = 1.0 if np.random.random() < 0.5 else -1.0
        x = _srk1_step(x, theta, dt, z, s)
    return x


@njit(cache=True)
def _run_path(x0, theta, dt, nsteps):
    """Integrate one trajectory and return the ``(nsteps + 1, 4)`` path."""
    x = x0.copy()
    path = np.empty((nsteps + 1, 4))
    path[0] = x
    for t in range(nsteps):
        z = np.random.standard_normal()
        s = 1.0 if np.random.random() < 0.5 else -1.0
        x = _srk1_step(x, theta, dt, z, s)
        path[t + 1] = x
    return path


# ---------------------------------------------------------------------------
# Object interface (mirrors the MATLAB usage).
# ---------------------------------------------------------------------------


class ExoDaisyWorld:
    """The coupled agent+environment stochastic Daisy World.

    Wraps the compiled integrator behind a small stateful interface that keeps
    the MATLAB usage: build with :class:`Parameters`, advance with
    :meth:`evolve` (one step) or :meth:`run` (many), and reach for the static
    :meth:`ensemble` to sweep a luminosity axis in parallel.

    Parameters
    ----------
    params : Parameters
        Physical parameters.
    x0 : array_like, optional
        Initial ``[f_B, f_W, T, L]``. Defaults to a small equal daisy seed at
        the equilibrium temperature ``T = L**0.25`` for ``L = 1 + lambda``.
    dt : float, optional
        Integration time step (default ``0.1``, as in the MATLAB).
    rng : numpy.random.Generator, optional
        Used only to seed NumPy's global state for the compiled loops (Numba's
        random stream). Pass an integer-seeded generator for repeatability.
    """

    def __init__(self, params, x0=None, dt=0.1, rng=None):
        self.params = params
        self.dt = float(dt)
        self.theta = params.to_theta()
        if x0 is None:
            L0 = 1.0 + params.lam
            x0 = [0.05 * params.f, 0.05 * params.f, L0 ** 0.25, L0]
        self.x = np.asarray(x0, dtype=np.float64).copy()
        if rng is not None:
            seed = int(rng.integers(0, 2 ** 31 - 1))
            np.random.seed(seed)

    def evolve(self):
        """Advance the state by one SRK1 step and return it."""
        z = np.random.standard_normal()
        s = 1.0 if np.random.random() < 0.5 else -1.0
        self.x = _srk1_step(self.x, self.theta, self.dt, z, s)
        return self.x

    def run(self, nsteps, record=False):
        """Integrate ``nsteps`` steps from the current state.

        Parameters
        ----------
        nsteps : int
            Number of SRK1 steps.
        record : bool, optional
            If True, return the full ``(nsteps + 1, 4)`` trajectory including
            the current state; otherwise advance in place and return the
            ``(4,)`` endpoint.

        Returns
        -------
        numpy.ndarray
            Endpoint or trajectory, per ``record``.
        """
        if record:
            path = _run_path(self.x, self.theta, self.dt, int(nsteps))
            self.x = path[-1].copy()
            return path
        end = _run_endpoint(self.x, self.theta, self.dt, int(nsteps))
        self.x = end.copy()
        return end

    @staticmethod
    def ensemble(params, Ls=None, nsteps=2000, N=500, dt=0.1,
                 agent_free=False, seed=None):
        """Parallel endpoint ensemble over a luminosity axis.

        Ports ``eDW_BHsim.m`` (biotic) and the agent-free driver of
        ``exoDaisyWorld.m``. For each luminosity, ``N`` independent instances
        are integrated for ``nsteps`` steps from freshly sampled initial
        conditions and their endpoints recorded.

        Parameters
        ----------
        params : Parameters
            Physical parameters; ``params.lam`` is overridden per luminosity so
            that ``<L> = 1 + lambda`` tracks the axis.
        Ls : array_like, optional
            Mean luminosities. Defaults to ``1 + linspace(-0.7, 1.4, 400)``.
        nsteps : int, optional
            Steps per instance (default 2000; the agent-free run used 200).
        N : int, optional
            Instances per luminosity (default 500).
        dt : float, optional
            Time step (default 0.1).
        agent_free : bool, optional
            If True, seed both daisy fractions at zero (environment only).
        seed : int, optional
            Seeds NumPy's global RNG before the parallel loop.

        Returns
        -------
        numpy.ndarray
            Endpoints of shape ``(len(Ls), N, 4)``.
        """
        if Ls is None:
            Ls = 1.0 + np.linspace(-0.7, 1.4, 400)
        Ls = np.asarray(Ls, dtype=np.float64)
        Ts = Ls ** 0.25
        theta = params.to_theta()
        # lambda per luminosity is (L - 1); it enters the OU mean and noise, so
        # each column integrates with its own theta[10]. A base seed makes each
        # column reproducible independently of thread scheduling; -1 disables.
        base_seed = -1 if seed is None else int(seed)
        return _ensemble_axis(theta, dt, int(nsteps), Ls, Ts, int(N),
                              agent_free, base_seed)


@njit(parallel=True, cache=True)
def _ensemble_axis(theta, dt, nsteps, Ls, Ts, N, agent_free, base_seed):
    """Endpoint ensemble over a luminosity axis, ``lambda = L - 1`` per column.

    The luminosity offset ``lambda`` enters the OU mean and noise amplitude, so
    each luminosity column integrates with its own ``theta[10] = Ls[ll] - 1``.
    When ``base_seed >= 0`` the column's RNG stream is seeded with
    ``base_seed + ll``, making the result reproducible regardless of how the
    parallel loop is scheduled across threads.
    """
    f = theta[0]
    nL = Ls.shape[0]
    data = np.empty((nL, N, 4))
    for ll in prange(nL):
        if base_seed >= 0:
            np.random.seed(base_seed + ll)
        th = theta.copy()
        th[10] = Ls[ll] - 1.0
        for nn in range(N):
            if agent_free:
                fB0 = 0.0
                fW0 = 0.0
            else:
                r = np.random.random()
                y = np.random.random()
                f1 = r * f
                fB0 = y * f1
                fW0 = y * (1.0 - f1)
            x = np.array([fB0, fW0, Ts[ll], Ls[ll]])
            for _ in range(nsteps):
                z = np.random.standard_normal()
                s = 1.0 if np.random.random() < 0.5 else -1.0
                x = _srk1_step(x, th, dt, z, s)
            data[ll, nn, 0] = x[0]
            data[ll, nn, 1] = x[1]
            data[ll, nn, 2] = x[2]
            data[ll, nn, 3] = x[3]
    return data
