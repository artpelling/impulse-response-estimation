# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "irdl",
#     "pyfar",
#     "scipy",
# ]
#
# [tool.uv.sources]
# irdl = { git = "https://github.com/artpelling/ssm-tools", subdirectory = "packages/irdl" }
# ///


def main():
    from irdl import get_fabian
    import numpy as np
    from pathlib import Path
    import pyfar as pf
    from scipy.io import savemat

    ir = get_fabian(kind="measured", hato=0)["impulse_response"][0, 0]
    sweep = pf.signals.exponential_sweep_time(ir.n_samples, (20, 20000))
    sweep = pf.dsp.pad_zeros(sweep, ir.n_samples)
    out = pf.dsp.convolve(ir, sweep, mode="cut")
    savemat(
        Path(__file__).parent / "data" / "fabian.mat",
        {
            "u": np.squeeze(sweep.time),
            "y": np.squeeze(out.time),
            "h": np.squeeze(ir.time),
            "fs": ir.sampling_rate,
        },
    )


if __name__ == "__main__":
    main()
