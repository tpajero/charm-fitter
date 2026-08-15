#!/usr/bin/env python3
"""Plot the values of all existing measurements of yCP (CP-even final states).

To check what options are available, run:

   python ycp.py -h

"""

import argparse
from pathlib import Path

from charm_fitter.blue import Measurement, blue_parser, get_units_label, plot_measurements
from charm_fitter.utils import setup_matplotlib

AVG_MEASURES = {
    "wa-2015": Measurement("World average 2015", -3.70e-3, 5.56e-3, 4.32e-3, is_average=True),
}


def get_measures(comb: str) -> list[Measurement]:
    """Get the list of measurements to be shown for a given summary plot.

    :param comb: Identifier for a given summary plot.
    :type comb: str
    """
    measures = [
        Measurement("Belle 2009", 1.1e-3, 6.1e-3, 5.2e-3, arxiv="0905.4185"),
        Measurement("BESIII 2015", -20.0e-3, 13.0e-3, 7.0e-3, arxiv="1501.01378"),
    ]
    measures.append(AVG_MEASURES[comb])
    return measures


def plot_average(comb: str, out_dir: Path, *, arxiv: bool, pub: bool) -> None:
    measures = get_measures(comb)
    units = 1e-3
    plot_measurements(
        measures,
        (-38.0, 47.0 if pub else 40.0),
        rf"$y_{{CP}}$ $[{get_units_label(units)}]$",
        out_dir / f"ycp-{comb}.pdf",
        units=units,
        arxiv=arxiv,
        pub=pub,
    )


def parse_args() -> argparse.Namespace:
    parser = blue_parser("ycp")
    parser.add_argument(
        "-c",
        "--comb",
        type=str,
        default=list(AVG_MEASURES.keys())[-1],
        help="Combination id",
        choices=list(AVG_MEASURES.keys()),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    setup_matplotlib(usetex=args.latex)
    plot_average(args.comb, args.outdir, arxiv=args.arxiv, pub=args.pub)
