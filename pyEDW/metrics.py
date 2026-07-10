"""Viability and the information architecture of Exo-Daisy World.

This module ports the entropy pipeline of ``Daisy_World_Paper_Figures.m``
(lines 307-392) and the measure assembly of ``makeFIG4.m`` / ``makeFig5.m``.
Given an endpoint ensemble ``data`` of shape ``(n_lum, N, 4)`` -- the output
of :meth:`pyEDW.model.ExoDaisyWorld.ensemble` -- it estimates, at each
luminosity, the joint distribution ``p_AE(f_B, f_W, T, L)`` by histogramming
and reads off the Shannon information measures of Section 4.

Binning follows the paper's square-root rule: with ``n`` living instances the
number of bins is ``N_bins = 1 + ceil(sqrt(n))`` and each variable is split by
``N_bins - 1`` equally spaced interior edges spanning its observed range, so a
value lands in bin ``sum(v > edges)`` (0-indexed) -- exactly the MATLAB
``1 + sum(v > x)`` up to the index base.

Entropy table
-------------
:func:`entropy_table` returns, per luminosity, a 16-column array ``H`` whose
columns are the entropies of every marginal of ``p_AE`` needed downstream
(agent dof ``a_1 = f_B``, ``a_2 = f_W``; environment dof ``e_1 = T``,
``e_2 = L``):

===  =====================  ===  =====================
col  entropy                col  entropy
===  =====================  ===  =====================
0    H[a1]                  8    H[a1 e1]
1    H[a2]                  9    H[a1 a2]
2    H[e1]                  10   H[a2 e1 e2]
3    H[e2]                  11   H[a1 e1 e2]
4    H[e1 e2]               12   H[a1 a2 e2]
5    H[a2 e2]               13   H[a1 a2 e1]
6    H[a2 e1]               14   H[a1 a2 e1 e2]
7    H[a1 e2]               15   N_bins
===  =====================  ===  =====================

The information measures are then

- ``I(A:E)  = H[e1 e2] + H[a1 a2] - H[a1 a2 e1 e2]``          (:func:`mutual_information`)
- ``dI      = I(e1:e2) - I0(e1:e2)``                          (:func:`delta_I`)
- ``C(a1:a2||E) = H[a1] + H[a2] + H[e1 e2] + H[a1 a2 e1 e2]``
  ``           - H[a1 a2] - H[a1 e1 e2] - H[a2 e1 e2]``       (:func:`cooperation`)

Reference: Sowinski, Ghoshal & Frank, *Planet. Sci. J.* **6**, 176 (2025),
Section 4 and Appendix B.3.
"""

import numpy as np

# Column indices into the entropy table H (see module docstring).
H_A1, H_A2, H_E1, H_E2 = 0, 1, 2, 3
H_E1E2, H_A2E2, H_A2E1, H_A1E2, H_A1E1, H_A1A2 = 4, 5, 6, 7, 8, 9
H_A2E1E2, H_A1E1E2, H_A1A2E2, H_A1A2E1, H_A1A2E1E2 = 10, 11, 12, 13, 14
NBINS = 15

_LIVE_TOL = 1e-7


def _entropy(p):
    """Shannon entropy in bits of a (sub)distribution array (ports ``getH``)."""
    p = p[p > 0.0]
    return -np.sum(p * np.log2(p))


def _bin_index(v, edges):
    """0-indexed bin of each value: ``sum(v > edges)`` (MATLAB convention)."""
    return np.searchsorted(edges, v, side="left")


def _joint_pmf(fB, fW, T, L):
    """4-D joint pmf ``p_AE`` over already-filtered living instances.

    Returns ``(p, n_bins)`` with ``p`` of shape ``(n_bins,) * 4`` on axes
    ``(a1, a2, e1, e2) = (f_B, f_W, T, L)``.
    """
    n = fB.shape[0]
    n_bins = 1 + int(np.ceil(np.sqrt(n)))
    cols = (fB, fW, T, L)
    idx = np.empty((4, n), dtype=np.intp)
    for a, v in enumerate(cols):
        edges = np.linspace(v.min(), v.max(), n_bins - 1)
        idx[a] = _bin_index(v, edges)
    p = np.zeros((n_bins, n_bins, n_bins, n_bins))
    np.add.at(p, (idx[0], idx[1], idx[2], idx[3]), 1.0)
    p /= p.sum()
    return p, n_bins


def _marginal_entropies(p, n_bins):
    """The 16-column entropy row for one luminosity's joint pmf ``p``.

    Axes of ``p`` are ``(a1, a2, e1, e2)``; marginals are taken by summing out
    the complementary axes, mirroring the MATLAB ``sum(pAE, [...])`` calls.
    """
    H = np.empty(16)
    H[H_A1] = _entropy(p.sum(axis=(1, 2, 3)))
    H[H_A2] = _entropy(p.sum(axis=(0, 2, 3)))
    H[H_E1] = _entropy(p.sum(axis=(0, 1, 3)))
    H[H_E2] = _entropy(p.sum(axis=(0, 1, 2)))
    H[H_E1E2] = _entropy(p.sum(axis=(0, 1)))
    H[H_A2E2] = _entropy(p.sum(axis=(0, 2)))
    H[H_A2E1] = _entropy(p.sum(axis=(0, 3)))
    H[H_A1E2] = _entropy(p.sum(axis=(1, 2)))
    H[H_A1E1] = _entropy(p.sum(axis=(1, 3)))
    H[H_A1A2] = _entropy(p.sum(axis=(2, 3)))
    H[H_A2E1E2] = _entropy(p.sum(axis=0))
    H[H_A1E1E2] = _entropy(p.sum(axis=1))
    H[H_A1A2E2] = _entropy(p.sum(axis=2))
    H[H_A1A2E1] = _entropy(p.sum(axis=3))
    H[H_A1A2E1E2] = _entropy(p)
    H[NBINS] = n_bins
    return H


def entropy_table(data):
    """Per-luminosity entropy table for a biotic endpoint ensemble.

    Parameters
    ----------
    data : numpy.ndarray
        Ensemble endpoints of shape ``(n_lum, N, 4)`` with the last axis
        ``[f_B, f_W, T, L]`` (from
        :meth:`pyEDW.model.ExoDaisyWorld.ensemble`).

    Returns
    -------
    numpy.ndarray
        Array ``H`` of shape ``(n_lum, 16)``; see the module docstring for the
        column layout. Luminosities with fewer than two living instances
        (``f_B f_W > 1e-7``) are left as ``NaN`` except the bin count.
    """
    data = np.asarray(data)
    n_lum = data.shape[0]
    H = np.full((n_lum, 16), np.nan)
    for ll in range(n_lum):
        fB = data[ll, :, 0]
        fW = data[ll, :, 1]
        live = fB * fW > _LIVE_TOL
        if np.count_nonzero(live) > 1:
            p, n_bins = _joint_pmf(
                fB[live], fW[live], data[ll, live, 2], data[ll, live, 3]
            )
            H[ll] = _marginal_entropies(p, n_bins)
    return H


def entropy_table_env(data_env):
    """Environment-only entropy table for the agent-free ensemble.

    Ports the ``data_entropyE`` block: with the biome removed, only ``(T, L)``
    vary. The bin count is fixed by the full instance count ``N`` (not a living
    subset), matching ``Nbins = 1 + ceil(sqrt(N))`` in the MATLAB.

    Parameters
    ----------
    data_env : numpy.ndarray
        Agent-free endpoints of shape ``(n_lum, N, 4)``.

    Returns
    -------
    numpy.ndarray
        Array of shape ``(n_lum, 4)`` with columns ``[H[e1], H[e2], H[e1 e2],
        N_bins]``.
    """
    data_env = np.asarray(data_env)
    n_lum, N = data_env.shape[0], data_env.shape[1]
    n_bins = 1 + int(np.ceil(np.sqrt(N)))
    H = np.zeros((n_lum, 4))
    for ll in range(n_lum):
        T = data_env[ll, :, 2]
        L = data_env[ll, :, 3]
        eT = np.linspace(T.min(), T.max(), n_bins - 1)
        eL = np.linspace(L.min(), L.max(), n_bins - 1)
        p = np.zeros((n_bins, n_bins))
        np.add.at(p, (_bin_index(T, eT), _bin_index(L, eL)), 1.0)
        p /= p.sum()
        H[ll, 0] = _entropy(p.sum(axis=1))
        H[ll, 1] = _entropy(p.sum(axis=0))
        H[ll, 2] = _entropy(p)
        H[ll, 3] = n_bins
    return H


# ---------------------------------------------------------------------------
# Derived measures.
# ---------------------------------------------------------------------------


def viability(data, f):
    """Viability ``V = E^A[(f_B + f_W)/f]`` per luminosity (Eq. 9).

    Parameters
    ----------
    data : numpy.ndarray
        Ensemble endpoints ``(n_lum, N, 4)``.
    f : float
        Habitable land fraction.

    Returns
    -------
    numpy.ndarray
        ``V`` of shape ``(n_lum,)``, the ensemble-mean occupied fraction of the
        habitable area (``0 <= V <= 1``).
    """
    data = np.asarray(data)
    return (data[:, :, 0] + data[:, :, 1]).mean(axis=1) / f


def efficacy(data):
    """Thermoregulatory efficacy ``E = <T> - 1`` per luminosity.

    The mean planetary temperature's deviation from the optimum (in units of
    ``T_opt``); zero means the biome holds the surface exactly at ``T_opt``.
    """
    data = np.asarray(data)
    return data[:, :, 2].mean(axis=1) - 1.0


def mutual_information(H):
    """Agent-environment mutual information ``I(A:E)`` (bits).

    ``I(A:E) = H[e1 e2] + H[a1 a2] - H[a1 a2 e1 e2]``.
    """
    H = np.asarray(H)
    return H[..., H_E1E2] + H[..., H_A1A2] - H[..., H_A1A2E1E2]


def intra_environment_information(H):
    """Intra-environment correlation ``I(e1:e2) = H[e1] + H[e2] - H[e1 e2]``."""
    H = np.asarray(H)
    return H[..., H_E1] + H[..., H_E2] - H[..., H_E1E2]


def delta_I(H, H_env):
    """Agent-induced change in intra-environment correlation, ``dI``.

    ``dI = I(e1:e2) - I0(e1:e2)`` (Eq. 12): the intra-environment mutual
    information with the biome present minus its agent-free value. Positive
    means the biome strengthens the temperature-luminosity correlation.

    Parameters
    ----------
    H : numpy.ndarray
        Biotic entropy table, shape ``(..., 16)`` from :func:`entropy_table`.
    H_env : numpy.ndarray
        Agent-free table, shape ``(n_lum, 4)`` from :func:`entropy_table_env`,
        broadcast against ``H`` along the luminosity axis.

    Returns
    -------
    numpy.ndarray
        ``dI`` in bits.
    """
    H = np.asarray(H)
    H_env = np.asarray(H_env)
    I0 = H_env[..., 0] + H_env[..., 1] - H_env[..., 2]
    return intra_environment_information(H) - I0


def cooperation(H):
    """Intra-agent cooperation given the environment, ``C(a1:a2||E)`` (Eq. 14).

    The interaction information

    ``C = H[a1] + H[a2] + H[e1 e2] + H[a1 a2 e1 e2]``
    ``  - H[a1 a2] - H[a1 e1 e2] - H[a2 e1 e2]``,

    equal to ``I(a1:a2) - I(a1:a2 | E)``. Negative values indicate a synergistic
    effect where the environment enhances the correlation between the two daisy
    species.
    """
    H = np.asarray(H)
    return (
        H[..., H_A1] + H[..., H_A2] + H[..., H_E1E2] + H[..., H_A1A2E1E2]
        - H[..., H_A1A2] - H[..., H_A1E1E2] - H[..., H_A2E1E2]
    )
