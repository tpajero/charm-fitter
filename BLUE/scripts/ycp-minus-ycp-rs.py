#!/usr/bin/env python3
"""Plot the values of all existing measurements of yCP - yCP(K- pi+).

To check what options are available, run:

   python ycp-minus-ycp-rs.py -h

"""

import argparse
from pathlib import Path

from charm_fitter.blue import Measurement, blue_parser, get_units_label, plot_measurements
from charm_fitter.utils import setup_matplotlib

AVG_MEASURES = {
    "wa-2016": Measurement("World average 2016", 9.15e-3, 1.39e-3, 0.77e-3, is_average=True),
    "wa-2018": Measurement("World average 2018", 7.42e-3, 0.95e-3, 0.59e-3, is_average=True),
    "wa-2022": Measurement("World average 2022", 6.99e-3, 0.25e-3, 0.13e-3, is_average=True),
}


def get_measures(comb: str, pre_bf: bool = False) -> list[Measurement]:
    """Get the list of measurements to be shown for a given summary plot.

    :param comb: Identifier for a given summary plot.
    :type comb: str
    """
    measures = []
    if pre_bf:
        measures += [
            Measurement("E791 1999", 7.32e-3, 28.91e-3, 10.32e-3, arxiv="hep-ex/9903012"),
            Measurement("FOCUS 2000", 34.2e-3, 13.9e-3, 7.4e-3, arxiv="hep-ex/0004034"),
            Measurement("CLEO 2002", -12.0e-3, 25.0e-3, 14.0e-3, arxiv="hep-ex/0111024"),
        ]
    measures += [
        Measurement("BaBar 2012", 7.2e-3, 1.8e-3, 1.2e-3, arxiv="1209.3896"),
        Measurement("Belle 2016", 11.1e-3, 2.2e-3, 0.9e-3, arxiv="1509.08266"),
    ]
    if comb in ("wa-2018", "wa-2022"):
        measures.append(Measurement("LHCb 2018", 5.7e-3, 1.3e-3, 0.9e-3, arxiv="1810.06874"))
    if comb == "wa-2022":
        measures.append(Measurement("LHCb 2022", 6.96e-3, 0.26e-3, 0.13e-3, arxiv="2202.09106"))
    measures.append(AVG_MEASURES[comb])
    return measures


def get_xrange(comb: str, pre_bf: bool, pub: bool) -> tuple[float, float]:
    """Get the range of the x axis in 1e-3 units."""
    if pre_bf:
        return (-45.0, (105.0 if comb == "wa-2016" else 110.0) if pub else 95.0)
    xrange = {
        "wa-2016": (4.5, 18.5 if pub else 17.5),
        "wa-2018": (3.5, 20.0 if pub else 18.0),
        "wa-2022": (3.5, 20.0 if pub else 18.5),
    }
    return xrange[comb]


def plot_average(comb: str, out_dir: Path, *, pre_bf: bool, arxiv: bool, pub: bool) -> None:
    measures = get_measures(comb, pre_bf=pre_bf)
    units = 1e-3
    fig_name = f"ycp-minus-ycp-rs-{comb}" + ("-pre-bf" if pre_bf else "")
    plot_measurements(
        measures,
        get_xrange(comb, pre_bf, pub),
        rf"$y_{{CP}} - y_{{CP}}^{{K^- \pi^+}}$ $[{get_units_label(units)}]$",
        out_dir / f"{fig_name}.pdf",
        figsize=(8, len(measures)),
        units=units,
        arxiv=arxiv,
        pub=pub,
    )


def parse_args() -> argparse.Namespace:
    parser = blue_parser("ycp-minus-ycp-rs")
    parser.add_argument(
        "-c",
        "--comb",
        type=str,
        default=list(AVG_MEASURES.keys())[0],
        help="Combination id",
        choices=list(AVG_MEASURES.keys()),
    )
    parser.add_argument(
        "--pre-bf",
        action="store_true",
        help="Plot the (not so precise) measurements made before B factories",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    setup_matplotlib(usetex=args.latex)
    plot_average(args.comb, args.outdir, pre_bf=args.pre_bf, arxiv=args.arxiv, pub=args.pub)
