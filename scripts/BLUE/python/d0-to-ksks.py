#!/usr/bin/env python3
"""Plot the values of all existing measurements of ACP(D0 -> KS KS).

To check what options are available, run:

   python d0-to-ksks.py -h

"""

import argparse
import logging
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt

from charm_fitter.blue import Measurement, blue_parser, get_units_label
from charm_fitter.utils import setup_matplotlib


def get_measures(comb: str) -> list[Measurement]:
    """Get the list of measurements to be shown for a given summary plot.

    :param comb: Identifier for a given summary plot.
    :type comb: str
    """
    measures = [
        Measurement("CLEO", "hep-ex/0012054", -23e-2, 19e-2),
        Measurement("LHCb 2015", "1508.06087", -2.9e-2, 5.2e-2, 2.2e-2),
        # Measurement("Belle", "1705.05966", -0.02e-2, 1.53e-2, 0.02e-2, 0.17e-2),
        Measurement("LHCb 2021", "2105.01565", -3.1e-2, 1.2e-2, 0.4e-2, 0.1e-2),
        Measurement("CMS", "2405.11606", 6.3e-2, 3.0e-2, 0.2e-2),  # Neglect ACP(KS pi pi)
    ]
    if comb == "2024":
        measures.extend(
            [
                Measurement("Belle + Belle II 2024", "2411.00306", -1.4e-2, 1.3e-2, 0.1e-2),
                Measurement("World average 2024", None, -1.63e-2, 0.84e-2, 0.20e-2),  # p-value 4.67%
            ]
        )
    else:
        measures.extend(
            [
                Measurement("Belle + Belle II 2025", "2504.15881", -0.6e-2, 1.1e-2, 0.1e-2),
                # Measurement("World average March 2025", None, -1.19e-2, 0.77e-2, 0.17e-2),  # p-value 3.76%
                Measurement("LHCb 2025", "2510.14732", 1.86e-2, 1.04e-2, 0.41e-2),
                Measurement("LHCb average 2025", None, -0.37e-2, 0.78e-2, 0.29e-2),
                Measurement("World average 2025", None, -0.17e-2, 0.62e-2, 0.18e-2),  # p-value 0.99%
            ]
        )
    return measures


def plot_average(date: str, out_dir: Path) -> None:
    measures = get_measures(date)
    n_meas = len(measures)
    if n_meas == 0:
        logging.warning(f"No measurements available for combination {date}. Will not produce the plot")
        return

    # Axes and labels
    fig_ysize = 6
    fig, ax = plt.subplots(figsize=(6, fig_ysize))
    (x_min, x_max), units = (-0.45, 0.6), 1e-2
    y_min, y_max = -0.5, n_meas - 0.5
    plt.xlim(x_min / units, x_max / units)
    plt.ylim(y_min, y_max)
    plt.tick_params(axis="y", which="both", right=False, left=False, labelleft=False)
    plt.xlabel(rf"$A_{{CP}}(D^0 \to K^0_S K^0_S)$ $[{get_units_label(units)}]$", fontsize=24, ha="center")

    # Plot the measures and their numerical values
    x_text = max([meas.val + meas.err() for meas in measures]) + 0.05 * (x_max - x_min) / units
    for i in range(n_meas):
        meas = measures[i]
        plt.errorbar(
            meas.val / units,
            n_meas - 1 - i,
            xerr=meas.err() / units,
            fmt=".",
            markersize=8,
            capsize=7,
            color=meas.color,
        )
        plt.errorbar(
            meas.val / units,
            n_meas - 1 - i,
            xerr=meas.stat / units,
            capsize=5,
            color=meas.color,
        )

        ylim = ax.get_ylim()
        height_px = ax.bbox.height
        fontsize_yscale = fig.dpi / 72 * (ylim[1] - ylim[0]) / height_px

        fontsize = matplotlib.rcParams["font.size"] * 8.0 / len(measures)
        y_text = n_meas - i - 1 - 0.25 * (n_meas / 6)

        if args.arxiv and meas.arxiv is not None:
            arxiv_fontsize = fontsize * 0.8
            y_text += 0.5 * fontsize_yscale * arxiv_fontsize
            plt.text(
                x_text,
                y_text - fontsize_yscale * arxiv_fontsize * 0.8,
                rf"\href{{https://arxiv.org/abs/{meas.arxiv}}}{{arXiv:{meas.arxiv}}}"
                if args.latex
                else f"arXiv:{meas.arxiv}",
                fontsize=arxiv_fontsize,
                color="deepskyblue",
            )
        plt.text(
            x_text,
            y_text,
            f"{meas.label}\n{meas.result_str(units)}",
            fontsize=fontsize,
            color=meas.color,
        )

    # Vertical and horizontal lines for world average
    if "average" in measures[-2].label:
        plt.plot(
            [x_min / units, x_max / units],
            [1.5, 1.5],
            linestyle="-",
            linewidth=1,
            color="k",
        )
    if "World" in measures[-1].label:
        plt.plot(
            [measures[-1].val / units, measures[-1].val / units],
            [y_min, y_max],
            linestyle="--",
            linewidth=1,
            color="k",
        )
        plt.plot(
            [x_min / units, x_max / units],
            [0.5, 0.5],
            linestyle="-",
            linewidth=1,
            color="k",
        )

    # Save figure
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_name = f"acp-d0-to-ksks-{date}"
    for ext in ["pdf"]:
        plt.savefig(out_dir / f"{fig_name}.{ext}")


def parse_args() -> argparse.Namespace:
    parser = blue_parser("d0-to-ksks")
    parser.add_argument(
        "-c",
        "--date",
        type=str,
        default="2025",
        help="Date of the combination",
        choices=["2024", "2025"],
    )
    return parser.parse_args()


if __name__ == "__main__":
    logging.basicConfig(encoding="utf-8", level=logging.INFO)

    args = parse_args()
    setup_matplotlib(usetex=args.latex)
    plot_average(args.date, args.outdir)
