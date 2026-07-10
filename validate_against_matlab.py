"""Validate the pyEDW ensemble against the original MATLAB ``.mat`` output.

Run this locally, where the MATLAB ``data_*.mat`` files are readable (the
sandbox mount used during development could not read them). It loads one
``data_XXX.mat`` produced by ``eDW_BHsim.m``, reruns the Python ensemble at the
same bandwidth and luminosity axis, and compares the ensemble-mean temperature,
luminosity, and viability curves. Agreement is statistical, not bit-for-bit:
the two codes use different RNG streams.

Usage
-----
    python validate_against_matlab.py /path/to/data_064.mat
"""

import sys

import numpy as np
from scipy.io import loadmat

from pyEDW.model import Parameters
from pyEDW.model import ExoDaisyWorld
from pyEDW import metrics

T_OPT = 300.0  # temperature scale used throughout eDW_BHsim.m


def main(path):
    m = loadmat(path)
    theta = m["theta"].ravel()
    Ls = m["Ls"].ravel()
    Ts = m["Ts"].ravel()
    data_ml = m["data"]  # (n_lum, N, 4)
    n_lum, N, _ = data_ml.shape

    # Recover the physical parameters from theta (layout in model.py).
    dT = T_OPT * (8.0 / theta[8]) ** 0.25
    p = Parameters(
        f=theta[0],
        A_G=1.0 - theta[1],
        A_B=theta[2] + (1.0 - theta[1]),
        A_W=theta[3] + (1.0 - theta[1]),
        T_opt=T_OPT,
        Q=theta[4],
        gamma_D=theta[5],
        gamma_G=1.0,
        tau_E=1.0 / theta[6],
        tau_S=1.0 / theta[7],
        dT=dT,
        delta=0.05,
    )
    print(f"loaded {path}: {n_lum} luminosities x {N} instances, dT={dT:.3f}")

    data_py = ExoDaisyWorld.ensemble(p, Ls=Ls, nsteps=2000, N=N, seed=0)

    for name, ml, py in [
        ("<T>", data_ml[:, :, 2].mean(1), data_py[:, :, 2].mean(1)),
        ("<L>", data_ml[:, :, 3].mean(1), data_py[:, :, 3].mean(1)),
        ("V", (data_ml[:, :, 0] + data_ml[:, :, 1]).mean(1) / p.f,
              metrics.viability(data_py, p.f)),
    ]:
        rms = np.sqrt(np.nanmean((ml - py) ** 2))
        print(f"  {name:>3}: RMS(MATLAB - Python) over luminosity = {rms:.4f}")

    print("Done. Curves should track within ensemble sampling error.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
