#!/usr/bin/env python3
"""Produce the plot with the world average of the DeltaY parameter.

To check what options are available, run:

   python dy.py -h

"""

import argparse
import sys
from pathlib import Path

from charm_fitter.blue import Measurement, blue_parser, plot_measurements
from charm_fitter.utils import setup_matplotlib

AVG_MEASURES = {
    # LHCb averages ----------------------------------------------------------------------------------------------------
    "lhcb-avg-2021": {
        fs: Measurement(r"LHCb average 2021", val, stat, syst, is_average=True)
        for fs, (val, stat, syst) in [
            ("all", (-1.08, 1.13, 0.33)),
            ("kk", (-0.35, 1.28, 0.32)),
            ("pipi", (-3.61, 2.38, 0.40)),
        ]
    },
    "lhcb-avg-2024": {"all": Measurement("LHCb average", -1.09, 1.11, 0.32, is_average=True)},
    # World averages ---------------------------------------------------------------------------------------------------
    "wa-2021": {
        fs: Measurement("World average", val, stat, syst, is_average=True)
        for fs, (val, stat, syst) in [
            ("all", (-0.96, 1.12, 0.32)),
            ("kk", (-0.20, 1.28, 0.32)),
            ("pipi", (-3.53, 2.36, 0.39)),
        ]
    },
    "wa-2024": {"all": Measurement("World average", -0.97, 1.11, 0.32, is_average=True)},
}

MEASURES = {
    "babar-2012": {"all": Measurement("BaBar 2012", -8.8, 25.5, 5.8, arxiv="1209.3896")},
    "cdf-2014": {
        fs: Measurement("CDF 2014", val, stat, syst, arxiv="1410.5435")
        for fs, (val, stat, syst) in [
            ("all", (12.0, 12.0, None)),
            ("kk", (19.0, 15.0, 4.0)),
            ("pipi", (1.0, 18.0, 3.0)),
        ]
    },
    "belle-2016": {"all": Measurement("Belle 2016", 3.0, 20.0, 7.0, arxiv="1509.08266")},
    "lhcb-2015": {
        fs: Measurement(r"LHCb 2015 $\mu^{-}$ tag (3 fb$^{-1}$)", val, stat, syst, arxiv="1501.06777")
        for fs, (val, stat, syst) in [
            ("all", (12.5, 7.3, None)),
            ("kk", (13.4, 7.7, 2.6)),
            ("pipi", (9.2, 14.5, 2.5)),
        ]
    },
    "lhcb-2017": {
        fs: Measurement("LHCb 2017 $D^{*+}$ tag (3 fb$^{-1}$)", val, stat, syst, arxiv="1702.06490")
        for fs, (val, stat, syst) in [
            ("all", (1.3, 2.8, 1.0)),
            ("kk", (3.0, 3.2, 1.0)),
            ("pipi", (-4.6, 5.8, 1.2)),
        ]
    },
    "lhcb-2019": {  # CONF note
        fs: Measurement("LHCb 2019 $D^{*+}$ tag (2 fb$^{-1}$)", val, stat, syst)
        for fs, (val, stat, syst) in [
            ("all", (-3.4, 3.1, 0.6)),
            ("kk", (-1.3, 3.5, 0.7)),
            ("pipi", (-11.3, 6.9, 0.8)),
        ]
    },
    "lhcb-2020": {
        fs: Measurement(r"LHCb 2020 $\mu^{-}$ tag (5.4 fb$^{-1}$)", val, stat, syst, arxiv="1911.01114")
        for fs, (val, stat, syst) in [
            ("all", (3.3, 3.6, 0.3)),
            ("kk", (4.8, 4.0, 0.3)),
            ("pipi", (-2.5, 7.8, 0.3)),
        ]
    },
    "lhcb-2021": {
        fs: Measurement("LHCb 2021 $D^{*+}$ tag (6 fb$^{-1}$)", val, stat, syst, arxiv="2105.09889")
        for fs, (val, stat, syst) in [
            ("all", (-2.70, 1.34, 0.30)),
            ("kk", (-2.32, 1.52, 0.32)),
            ("pipi", (-4.01, 2.81, 0.39)),
        ]
    },
    "lhcb-2024": {
        "all": Measurement(r"LHCb 2024 $D^0\to\pi^{+}\pi^{-}\pi^{0}$ (6 fb$^{-1}$)", -1.3, 6.3, 2.4, arxiv="2405.06556")
    },
} | AVG_MEASURES


def get_measures(comb: str, final_state: str) -> list[Measurement]:
    """Get the list of measurements to be shown for a given summary plot.

    :param comb: Identifier for a given summary plot.
    :type comb: str
    :type final_state: str
    """
    comb_measures = []
    if "lhcb" not in comb:
        if final_state == "all":
            comb_measures.append(MEASURES["babar-2012"]["all"])
        comb_measures.append(MEASURES["cdf-2014"][final_state])
    comb_measures.append(MEASURES["lhcb-2015"][final_state])
    if "lhcb" not in comb and final_state == "all":
        comb_measures.append(MEASURES["belle-2016"]["all"])
    comb_measures.extend(
        [
            MEASURES["lhcb-2017"][final_state],
            MEASURES["lhcb-2020"][final_state],
            MEASURES["lhcb-2021"][final_state],
        ]
    )
    if "2024" in comb and final_state == "all":
        comb_measures.append(MEASURES["lhcb-2024"]["all"])
    try:
        comb_measures.append(MEASURES[comb][final_state])  # World/LHCb average
    except KeyError:
        sys.exit(
            f"ERROR: Only the following final states are available for combination {comb}: {list(MEASURES[comb].keys())}. Exiting..."
        )
    return comb_measures


def get_xrange(comb: str, final_state: str, dy_notation: bool) -> tuple[float, float]:
    """Get the x ranges for a given summary plot, in 10^-4 units."""
    xrange = {
        "lhcb-avg-2021": {
            "all": (-7.0 if dy_notation else -23.0, 54.0 if dy_notation else 38.0),
            "kk": (-6.0 if dy_notation else -24.0, 56.0 if dy_notation else 40.0),
            "pipi": (-13.0 if dy_notation else -27.0, 70.0 if dy_notation else 58.0),
        },
        "lhcb-avg-2024": {
            "all": (-12.0 if dy_notation else -23.0, 80.0 if dy_notation else 67.0),
        },
        "wa-2021": {
            "all": (-40.0 if dy_notation else -27.0, 105.0 if dy_notation else 115.0),
            "kk": (-10.0 if dy_notation else -38.0, 90.0 if dy_notation else 57.0),
            "pipi": (-20.0 if dy_notation else -27.0, 84.0 if dy_notation else 72.0),
        },
        "wa-2024": {
            "all": (-40.0 if dy_notation else -30.0, 145.0 if dy_notation else 157.0),
        },
    }
    return xrange[comb][final_state]


def make_plot(
    comb: str,
    final_state: str,
    out_dir: Path,
    *,
    dy_notation: bool,
    arxiv: bool,
    pub: bool,
) -> None:
    measures = get_measures(comb, final_state)

    decay_label = {"all": "", "kk": r"(D^0 \to K^+ K^-)", "pipi": r"(D^0 \to \pi^+ \pi^-)"}
    observable = r"\Delta Y" if dy_notation else r"A_{\Gamma}"
    xlabel = f"${observable}{decay_label[final_state]}\\;[10^{{-4}}]$"

    obs = "dy" if dy_notation else "agamma"
    fig_name = f"{obs}-{comb}"
    if final_state != "all":
        fig_name += f"-{final_state}"

    plot_measurements(
        measures,
        get_xrange(comb, final_state, dy_notation),
        xlabel,
        out_dir / f"{fig_name}.pdf",
        val_transform=lambda v: v if dy_notation else -v,
        color_fn=lambda meas: meas.color if "lhcb" not in comb else "k",
        arxiv=arxiv,
        pub=pub,
    )


def parse_args() -> argparse.Namespace:
    parser = blue_parser("dy")
    parser.add_argument(
        "-c",
        "--comb",
        type=str,
        default=list(AVG_MEASURES.keys())[0],
        help="Combination id",
        choices=list(AVG_MEASURES.keys()),
    )
    parser.add_argument(
        "--fs",
        type=str,
        default="all",
        help="Final state for the combination",
        choices=["all", "kk", "pipi"],
    )
    parser.add_argument(
        "--agamma",
        action="store_false",
        dest="dy_notation",
        help="Use the A_Gamma notation rather than DeltaY",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    setup_matplotlib(usetex=args.latex)
    make_plot(args.comb, args.fs, args.outdir, dy_notation=args.dy_notation, arxiv=args.arxiv, pub=args.pub)
