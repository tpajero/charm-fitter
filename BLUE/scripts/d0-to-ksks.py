#!/usr/bin/env python3
"""Plot the values of all existing measurements of ACP(D0 -> KS KS).

To check what options are available, run:

   python d0-to-ksks.py -h

"""

import argparse
from pathlib import Path

from charm_fitter.blue import Measurement, blue_parser, get_units_label, plot_measurements
from charm_fitter.utils import setup_matplotlib


def get_measures(comb: str) -> list[Measurement]:
    """Get the list of measurements to be shown for a given summary plot.

    :param comb: Identifier for a given summary plot.
    :type comb: str
    """
    measures = [
        Measurement("CLEO", -23e-2, 19e-2, arxiv="hep-ex/0012054"),
        Measurement("LHCb 2015", -2.9e-2, 5.2e-2, 2.2e-2, arxiv="1508.06087"),
        # Measurement("Belle", -0.02e-2, 1.53e-2, 0.02e-2, 0.17e-2, arxiv="1705.05966"),
        Measurement("LHCb 2021", -3.1e-2, 1.2e-2, 0.4e-2, 0.1e-2, arxiv="2105.01565"),
        Measurement("CMS", 6.3e-2, 3.0e-2, 0.2e-2, arxiv="2405.11606"),  # Neglect ACP(KS pi pi)
    ]
    if comb == "2024":
        measures.extend(
            [
                Measurement("Belle + Belle II 2024", -1.4e-2, 1.3e-2, 0.1e-2, arxiv="2411.00306"),
                Measurement("World average 2024", -1.63e-2, 0.84e-2, 0.20e-2, is_average=True),  # p-value 4.67%
            ]
        )
    else:
        measures.extend(
            [
                Measurement("Belle + Belle II 2025", -0.6e-2, 1.1e-2, 0.1e-2, arxiv="2504.15881"),
                # Measurement("World average March 2025", -1.19e-2, 0.77e-2, 0.17e-2, is_average=True),  # p-value 3.76%
                Measurement("LHCb 2025", 1.86e-2, 1.04e-2, 0.41e-2, arxiv="2510.14732"),
                Measurement("LHCb average 2025", -0.37e-2, 0.78e-2, 0.29e-2, is_average=True),
                Measurement("World average 2025", -0.17e-2, 0.62e-2, 0.18e-2, is_average=True),  # p-value 0.99%
            ]
        )
    return measures


def plot_average(date: str, out_dir: Path, *, arxiv: bool, pub: bool) -> None:
    measures = get_measures(date)
    units = 1e-2
    plot_measurements(
        measures,
        rf"$A_{{CP}}(D^0 \to K^0_S K^0_S)$ $[{get_units_label(units)}]$",
        out_dir / f"acp-d0-to-ksks-{date}.pdf",
        units=units,
        arxiv=arxiv,
        pub=pub,
    )


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
    args = parse_args()
    setup_matplotlib(usetex=args.latex)
    plot_average(args.date, args.outdir, arxiv=args.arxiv, pub=args.pub)
