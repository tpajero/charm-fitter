#!/usr/bin/env python3
"""Produce the plot with the world average of the DeltaY parameter.

To check what options are available, run:

   python dy.py -h

"""

import argparse
import itertools
import sys
from pathlib import Path

from charm_fitter.blue import Combination, Measurement, blue_parser, plot_measurements
from charm_fitter.utils import setup_matplotlib

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
}

# One entry per combination id in BLUE/main/dy.cpp. Organised as: {avg-id: {fs: combination}}
COMBINATIONS: dict[str, dict[str, Combination]] = {
    # LHCb averages ----------------------------------------------------------------------------------------------------
    "lhcb-2021": {
        fs: Combination(
            Measurement(r"LHCb average 2021", val, stat, syst, is_average=True),
            ["lhcb-2015", "lhcb-2017", "lhcb-2020", "lhcb-2021"],
        )
        for fs, (val, stat, syst) in [
            ("all", (-1.08, 1.13, 0.33)),
            ("kk", (-0.35, 1.28, 0.32)),
            ("pipi", (-3.61, 2.38, 0.40)),
        ]
    },
    "lhcb-2024": {
        "all": Combination(
            Measurement("LHCb average", -1.09, 1.11, 0.32, is_average=True),
            ["lhcb-2015", "lhcb-2017", "lhcb-2020", "lhcb-2021", "lhcb-2024"],
        )
    },
    # World averages ---------------------------------------------------------------------------------------------------
    "wa-2021": {
        fs: Combination(
            Measurement("World average", val, stat, syst, is_average=True),
            ["babar-2012", "cdf-2014", "belle-2016", "lhcb-2015", "lhcb-2017", "lhcb-2020", "lhcb-2021"],
        )
        for fs, (val, stat, syst) in [
            ("all", (-0.96, 1.12, 0.32)),
            ("kk", (-0.20, 1.28, 0.32)),
            ("pipi", (-3.53, 2.36, 0.39)),
        ]
    },
    "wa-2024": {
        "all": Combination(
            Measurement("World average", -0.97, 1.11, 0.32, is_average=True),
            ["babar-2012", "cdf-2014", "lhcb-2015", "belle-2016", "lhcb-2017", "lhcb-2020", "lhcb-2021", "lhcb-2024"],
        )
    },
}


def get_measures(comb: str, final_state: str) -> list[Measurement]:
    """Get the list of measurements to be shown for a given summary plot.

    :param comb: Identifier for a given summary plot.
    :type comb: str
    :type final_state: str
    """
    try:
        combination = COMBINATIONS[comb][final_state]
    except KeyError:
        err_str = (
            f"Available combinations: {list(COMBINATIONS.keys())}"
            if comb not in COMBINATIONS
            else f"Available final states for combination {comb}: {list(COMBINATIONS[comb].keys())}"
        )
        sys.exit(f"ERROR: Combination {comb!r} with final state {final_state!r} not available. {err_str}. Exiting...")
    return [MEASURES[name][final_state] for name in combination.members if final_state in MEASURES[name]] + [
        combination.average
    ]


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
        default=list(COMBINATIONS.keys())[-1],
        help="Combination id",
        choices=list(COMBINATIONS.keys()),
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
    parser.add_argument(
        "--run-all",
        action="store_true",
        help="Run all available combinations (ignores --comb, --fs, --agamma, --arxiv and --pub).",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    setup_matplotlib(usetex=args.latex)
    if args.run_all:
        for comb_id in COMBINATIONS:
            for fs in COMBINATIONS[comb_id]:
                for dy_notation, (arxiv, pub) in itertools.product(
                    [False, True], [(False, False), (True, False), (False, True)]
                ):
                    make_plot(comb_id, fs, args.outdir, dy_notation=dy_notation, arxiv=arxiv, pub=pub)
    else:
        make_plot(args.comb, args.fs, args.outdir, dy_notation=args.dy_notation, arxiv=args.arxiv, pub=args.pub)
