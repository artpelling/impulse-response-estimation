# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "irdl",
#     "scipy",
# ]
#
# [tool.uv.sources]
# irdl = { git = "https://github.com/artpelling/ssm-tools", subdirectory = "packages/irdl" }
# ///


def main():
    from irdl import get_fabian
    from pathlib import Path
    from scipy.io import savemat

    ir = get_fabian(kind="measured", hato=0)["impulse_response"]
    savemat(
        Path(__file__).parent / "data" / "fabian.mat",
        {
            "h": ir.time[:, 0, 0],
            "fs": ir.sampling_rate,
        },
    )


if __name__ == "__main__":
    main()
