# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "numpy",
#     "scipy",
# ]
# ///

import numpy as np
from scipy.fft import ifft
from scipy.io import savemat
from pathlib import Path

alpha = np.array([0.7 - 0.5j, 0.7 + 0.5j, 0.9])  # mirror image poles
c = np.array([1 + 2j, 1 - 2j, 1])  # residues
n_ir = 512  # length of the impulse response
n_input = 5000  # length of the input signal


def main() -> None:
    # compute ir from poles and residues
    h = np.sum(c[:, np.newaxis] * np.power.outer(alpha.conj(), np.arange(n_ir)), axis=0)
    # construct u with a geometric sequence
    u = 0.9 ** np.arange(n_input)
    y = np.convolve(u, h, mode="same")
    savemat(
        Path(__file__).parent / "data" / "toy_example.mat",
        {
            "u": np.squeeze(u),
            "y": np.squeeze(y),
            "h": np.squeeze(h),
            "fs": 1,
        },
    )


if __name__ == "__main__":
    main()
