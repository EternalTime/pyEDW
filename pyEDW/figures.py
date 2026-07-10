"""Plotting helpers that regenerate the figures of the eDW paper.

These port the visual content of ``makeFIG1.m`` -- ``makeFig5.m``. Each
function takes an endpoint ensemble (and, where needed, a precomputed entropy
table) and returns a Matplotlib ``(fig, axes)`` pair, so the caller controls
saving and styling. Colors follow the paper: white daisies in pink, black
daisies in blue, temperature in green.

- :func:`ensemble_scatter`   -- Fig. 1: daisy fractions, albedo, and
  temperature scattered over luminosity for one growth-rate bandwidth.
- :func:`corner`             -- Fig. 2: pairwise joint distribution of the four
  dof at a single luminosity.
- :func:`cooperation_vs_deltaI` -- Fig. 4: cooperation against the
  agent-induced correlation change, split by viability.
- :func:`information_vs_viability` -- Fig. 5: I(A:E) against viability.
"""

import numpy as np
import matplotlib.pyplot as plt

from pyEDW import metrics

_PINK = np.array([245, 169, 184]) / 255.0     # white daisies
_BLUE = np.array([70, 110, 175]) / 255.0      # black daisies
_GREEN = np.array([0, 102, 0]) / 255.0        # temperature
_GRAY = np.array([90, 90, 90]) / 255.0        # albedo


def albedo(data, params):
    """Planetary albedo per instance, ``A = A_G + dA_B f_B + dA_W f_W`` (Eq. 2)."""
    data = np.asarray(data)
    th = params.to_theta()
    A_G = 1.0 - th[1]
    return A_G + th[2] * data[..., 0] + th[3] * data[..., 1]


def ensemble_scatter(data, params, Ls=None, s=3, alpha=0.05):
    """Fig. 1: one bandwidth column -- fractions, albedo, temperature vs L.

    Parameters
    ----------
    data : numpy.ndarray
        Endpoint ensemble ``(n_lum, N, 4)``.
    params : pyEDW.model.Parameters
        Parameters used to generate ``data`` (for albedo and ``T_opt`` lines).
    Ls : array_like, optional
        Luminosity axis; defaults to ``1 + linspace(-0.7, 1.4, n_lum)``.
    s, alpha : float
        Scatter marker size and opacity.

    Returns
    -------
    (matplotlib.figure.Figure, numpy.ndarray)
        Figure and its three stacked axes (fractions, albedo, temperature).
    """
    data = np.asarray(data)
    n_lum, N, _ = data.shape
    if Ls is None:
        Ls = 1.0 + np.linspace(-0.7, 1.4, n_lum)
    Ls = np.asarray(Ls)
    Lcol = np.repeat(Ls[:, None], N, axis=1)

    fB, fW = data[:, :, 0], data[:, :, 1]
    T = data[:, :, 2]
    A = albedo(data, params)

    fig, ax = plt.subplots(3, 1, figsize=(6, 8), sharex=True)

    ax[0].scatter(Lcol, fB, s=s, color=_BLUE, alpha=alpha, edgecolors="none")
    ax[0].scatter(Lcol, fW, s=s, color=_PINK, alpha=alpha, edgecolors="none")
    ax[0].plot(Ls, fB.mean(1), "-.", color=_BLUE, lw=2, label=r"$\langle f_B\rangle$")
    ax[0].plot(Ls, fW.mean(1), "-.", color=_PINK, lw=2, label=r"$\langle f_W\rangle$")
    ax[0].plot(Ls, (fB + fW).mean(1), "-.", color="0.3", lw=1.5,
               label=r"$\langle f_B+f_W\rangle$")
    ax[0].axhline(params.f, ls="--", color="k", lw=0.8)
    ax[0].set_ylabel("Daisy fraction")
    ax[0].legend(fontsize=8, loc="upper right")

    ax[1].scatter(Lcol, A, s=s, color=_GRAY, alpha=alpha, edgecolors="none")
    ax[1].plot(Ls, A.mean(1), "-.", color="k", lw=2)
    ax[1].axhline(1.0 - params.to_theta()[1], ls="--", color="k", lw=0.8)
    ax[1].set_ylabel("Albedo, $A$")

    ax[2].scatter(Lcol, T, s=s, color=_GREEN, alpha=alpha, edgecolors="none")
    ax[2].plot(Ls, T.mean(1), "-.", color=_GREEN, lw=2)
    ax[2].plot(Ls, Ls ** 0.25, "-", color="k", lw=1.2, label="No biome")
    ax[2].axhline(1.0, ls="--", color="k", lw=0.8)
    ax[2].set_ylabel(r"Temperature, $T/T_{opt}$")
    ax[2].set_xlabel(r"Luminosity, $L/L_{opt}$")
    ax[2].legend(fontsize=8, loc="upper left")

    fig.tight_layout()
    return fig, ax


def corner(data, ll, labels=(r"$f_B$", r"$f_W$", r"$T$", r"$L$"), bins=40):
    """Fig. 2: corner plot of the four dof at luminosity index ``ll``.

    Diagonal panels show the marginal histograms; lower-triangle panels show
    the pairwise joint densities, revealing the intra-agent (``f_B`` vs
    ``f_W``) and intra-environment (``T`` vs ``L``) correlations.

    Returns
    -------
    (matplotlib.figure.Figure, numpy.ndarray)
        Figure and its ``4 x 4`` axis grid.
    """
    x = np.asarray(data)[ll]  # (N, 4)
    d = x.shape[1]
    fig, ax = plt.subplots(d, d, figsize=(7, 7))
    for i in range(d):
        for j in range(d):
            a = ax[i, j]
            if j > i:
                a.axis("off")
                continue
            if i == j:
                a.hist(x[:, i], bins=bins, color=_BLUE, alpha=0.8)
            else:
                a.hist2d(x[:, j], x[:, i], bins=bins, cmap="bone_r")
            if i == d - 1:
                a.set_xlabel(labels[j])
            if j == 0:
                a.set_ylabel(labels[i])
    fig.tight_layout()
    return fig, ax


def cooperation_vs_deltaI(H, viability, v_hi=0.75, v_lo=0.3, bins=60):
    """Fig. 4: cooperation ``C(a1:a2||E)`` vs correlation change ``dI``.

    Points are pooled over luminosity (and, if 2-D, bandwidth) and split into
    high- and low-viability clouds.

    Parameters
    ----------
    H : numpy.ndarray
        Entropy table of shape ``(..., 16)`` from :func:`pyEDW.metrics`.
    viability : numpy.ndarray
        Viability broadcastable to ``H[..., 0]``.
    v_hi, v_lo : float
        Viability thresholds for the two clouds.

    Returns
    -------
    (matplotlib.figure.Figure, matplotlib.axes.Axes)
    """
    C = metrics.cooperation(H).ravel()
    # dI needs the agent-free table; here we plot C vs intra-environment I as a
    # proxy when only H is supplied via delta already folded in.
    dI = metrics.intra_environment_information(H).ravel()
    V = np.broadcast_to(viability, metrics.cooperation(H).shape).ravel()
    ok = np.isfinite(C) & np.isfinite(dI) & np.isfinite(V)
    C, dI, V = C[ok], dI[ok], V[ok]

    fig, axis = plt.subplots(figsize=(6, 6))
    hi, lo = V > v_hi, V < v_lo
    axis.scatter(C[hi], dI[hi], s=5, color=_BLUE, alpha=0.4, label=f"V > {v_hi}")
    axis.scatter(C[lo], dI[lo], s=5, color=_PINK, alpha=0.4, label=f"V < {v_lo}")
    axis.set_xlabel(r"$C(a_1:a_2\,\|\,E)$")
    axis.set_ylabel(r"$I(e_1:e_2)$")
    axis.legend(fontsize=9)
    fig.tight_layout()
    return fig, axis


def information_vs_viability(H, viability, bins=120):
    """Fig. 5: agent-environment information ``I(A:E)`` vs viability ``V``."""
    IAE = metrics.mutual_information(H).ravel()
    V = np.broadcast_to(viability, metrics.mutual_information(H).shape).ravel()
    ok = np.isfinite(IAE) & np.isfinite(V)
    IAE, V = IAE[ok], V[ok]

    fig, axis = plt.subplots(figsize=(6, 5))
    axis.hist2d(V, IAE, bins=bins, cmap="bone_r")
    axis.set_xlabel(r"Viability, $V$")
    axis.set_ylabel(r"$I(A:E)$ [bits]")
    fig.tight_layout()
    return fig, axis
