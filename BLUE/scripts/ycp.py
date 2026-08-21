#!/usr/bin/env python3
"""Plot the values of all existing measurements of yCP (using the sign convention based on CP-even final states) and of
the similar observables yCP - yCP(Kpi) and yCP - yCP(K- pi+), as combined in BLUE/main/ycp.cpp.

To check what options are available, run:

   python ycp.py -h

"""

import argparse
import sys
from pathlib import Path

from charm_fitter.blue import Combination, Measurement, blue_parser, get_units_label, plot_measurements
from charm_fitter.utils import setup_matplotlib

# Pre-B-factory measurements: much less precise than everything that came after, so by default they are omitted from
# the plots (pass --pre-bf to show them) to keep the rest of the figure readable.
PRE_BF_MEASURES = {"E791 1999", "FOCUS 2000", "CLEO 2001"}

# Individual measurements, matching the `estimates` map in BLUE/main/ycp.cpp.
MEASURES = {
    m.label: m
    for m in [
        Measurement("E791 1999", 7.32e-3, 28.91e-3, 10.32e-3, arxiv="hep-ex/9903012"),
        Measurement("FOCUS 2000", 34.2e-3, 13.9e-3, 7.4e-3, arxiv="hep-ex/0004034"),
        Measurement("CLEO 2001", -12.0e-3, 25.0e-3, 14.0e-3, arxiv="hep-ex/0111024"),
        Measurement("Belle 2009", 1.1e-3, 6.1e-3, 5.2e-3, arxiv="0905.4185"),
        Measurement("LHCb 2011", 5.5e-3, 6.3e-3, 4.1e-3, arxiv="1112.4698"),
        Measurement("BaBar 2012", 7.2e-3, 1.8e-3, 1.2e-3, arxiv="1209.3896"),
        Measurement("BESIII 2015", -20.0e-3, 13.0e-3, 7.0e-3, arxiv="1501.01378"),
        Measurement("Belle 2015", 11.1e-3, 2.2e-3, 0.9e-3, arxiv="1509.08266"),
        Measurement("LHCb 2018", 5.7e-3, 1.3e-3, 0.9e-3, arxiv="1810.06874"),
        Measurement("Belle 2019", 9.6e-3, 9.1e-3, 6.43e-3, arxiv="1912.10912"),
        Measurement("LHCb 2022", 6.96e-3, 0.26e-3, 0.13e-3, arxiv="2202.09106"),
    ]
}


# One entry per combination id in BLUE/main/ycp.cpp. Organised as: {observable: {id-for-out-file: combination}}
COMBINATIONS: dict[str, dict[str, Combination]] = {
    "ycp": {
        "wa-2015": Combination(
            Measurement("World average 2015", -3.70e-3, 5.56e-3, 4.32e-3, is_average=True),
            ["Belle 2009", "BESIII 2015"],
        ),
    },
    "ycp-minus-ycp-kp": {
        "wa-2015": Combination(
            Measurement("World average 2015", 7.70e-3, 1.78e-3, 1.18e-3, is_average=True),
            ["E791 1999", "FOCUS 2000", "BaBar 2012"],
        ),
    },
    "ycp-minus-ycp-rs": {
        "wa-2015": Combination(
            Measurement("World average 2015", 10.45e-3, 2.07e-3, 0.90e-3, is_average=True),
            ["CLEO 2001", "LHCb 2011", "Belle 2015"],
        ),
        "wa-2018": Combination(
            Measurement("World average 2018", 7.26e-3, 1.11e-3, 0.67e-3, is_average=True),
            ["CLEO 2001", "LHCb 2011", "Belle 2015", "LHCb 2018"],
        ),
        "wa-2022": Combination(
            Measurement("World average 2022", 6.97e-3, 0.25e-3, 0.13e-3, is_average=True),
            ["CLEO 2001", "LHCb 2011", "Belle 2015", "LHCb 2018", "LHCb 2022"],
        ),
        "wa-no-lhcb-2015": Combination(
            Measurement("World average (no LHCb) 2015", 10.94e-3, 2.19e-3, 0.90e-3, is_average=True),
            ["CLEO 2001", "Belle 2015"],
        ),
    },
    # "All" combinations mix yCP, yCP - yCP(Kpi) and yCP +/- yCP(K-pi+) measurements as if they were all measuring the
    # same observable, without correcting for the different control channels: the resulting averages are biased.
    "all": {
        "wa-biased-2015-01": Combination(
            Measurement("World average Jan 2015", 6.56e-3, 1.64e-3, 1.10e-3, is_average=True),
            ["E791 1999", "FOCUS 2000", "CLEO 2001", "Belle 2009", "LHCb 2011", "BaBar 2012", "BESIII 2015"],
        ),
        "wa-biased-2015-09": Combination(
            Measurement("World average Sep 2015", 8.41e-3, 1.32e-3, 0.75e-3, is_average=True),
            [
                "E791 1999",
                "FOCUS 2000",
                "CLEO 2001",
                "Belle 2009",
                "LHCb 2011",
                "BaBar 2012",
                "BESIII 2015",
                "Belle 2015",
            ],
        ),
        "wa-biased-2018": Combination(
            Measurement("World average 2018", 7.11e-3, 0.93e-3, 0.58e-3, is_average=True),
            [
                "E791 1999",
                "FOCUS 2000",
                "CLEO 2001",
                "Belle 2009",
                "LHCb 2011",
                "BaBar 2012",
                "BESIII 2015",
                "Belle 2015",
                "LHCb 2018",
            ],
        ),
        "wa-biased-2019": Combination(
            Measurement("World average 2019", 7.14e-3, 0.92e-3, 0.58e-3, is_average=True),
            [
                "E791 1999",
                "FOCUS 2000",
                "CLEO 2001",
                "Belle 2009",
                "LHCb 2011",
                "BaBar 2012",
                "BESIII 2015",
                "Belle 2015",
                "LHCb 2018",
                "Belle 2019",
            ],
        ),
        "wa-biased-2022": Combination(
            Measurement("World average 2022", 6.97e-3, 0.25e-3, 0.13e-3, is_average=True),
            [
                "E791 1999",
                "FOCUS 2000",
                "CLEO 2001",
                "Belle 2009",
                "LHCb 2011",
                "BaBar 2012",
                "BESIII 2015",
                "Belle 2015",
                "LHCb 2018",
                "Belle 2019",
                "LHCb 2022",
            ],
        ),
        "wa-biased-no-lhcb-2019": Combination(
            Measurement("World average (no LHCb) 2019", 8.55e-3, 1.34e-3, 0.75e-3, is_average=True),
            [
                "E791 1999",
                "FOCUS 2000",
                "CLEO 2001",
                "Belle 2009",
                "BaBar 2012",
                "BESIII 2015",
                "Belle 2015",
                "Belle 2019",
            ],
        ),
    },
}

XLABELS = {
    "ycp": r"y_{CP}",
    "ycp-minus-ycp-kp": r"y_{CP} - y_{CP}^{K\pi}",
    "ycp-minus-ycp-rs": r"y_{CP} - y_{CP}^{K^- \pi^+}",
    "all": r"y_{CP}",
}


def get_measures(obs: str, comb: str, *, pre_bf: bool) -> list[Measurement]:
    """Get the list of measurements to be shown for a given summary plot.

    :param obs: Observable (key of `COMBINATIONS`).
    :param comb: Identifier of the combination within `obs` (key of `COMBINATIONS[obs]`).
    :param pre_bf: Whether to also show the pre-B-factory measurements in `PRE_BF_MEASURES`, when they are members of
        the combination.
    """
    try:
        combination = COMBINATIONS[obs][comb]
    except KeyError:
        err_str = (
            f"Available observables: {list(COMBINATIONS.keys())}"
            if obs not in COMBINATIONS
            else f"Available combinations: {list(COMBINATIONS[obs].keys())}"
        )
        sys.exit(f"ERROR: Combination {comb!r} not available for observable {obs!r}.\n{err_str}. Exiting...")
    members = combination.members if pre_bf else [name for name in combination.members if name not in PRE_BF_MEASURES]
    return [MEASURES[name] for name in members] + [combination.average]


def plot_average(obs: str, comb: str, out_dir: Path, *, pre_bf: bool, arxiv: bool, pub: bool) -> None:
    measures = get_measures(obs, comb, pre_bf=pre_bf)
    units = 1e-3
    fig_prefix = "ycp-no-corrections" if obs == "all" else obs
    fig_name = f"{fig_prefix}-{comb}" + ("-pre-bf" if pre_bf else "")
    plot_measurements(
        measures,
        rf"${XLABELS[obs]}$ $[{get_units_label(units)}]$",
        out_dir / f"{fig_name}.pdf",
        units=units,
        arxiv=arxiv,
        pub=pub,
    )


def parse_args() -> argparse.Namespace:
    parser = blue_parser("ycp")
    default_obs = "ycp-minus-ycp-rs"
    parser.add_argument(
        "--obs",
        type=str,
        help="Observable",
        choices=list(COMBINATIONS.keys()),
        default=default_obs,
    )
    parser.add_argument(
        "-c",
        "--comb",
        type=str,
        help="\n".join(
            ["Combination id. Available choices depend on --obs:"]
            + [f"- {k}: {list(v.keys())}" for k, v in COMBINATIONS.items()]
        ),
        default=list(COMBINATIONS[default_obs].keys())[-1],
    )
    parser.add_argument(
        "--pre-bf",
        action="store_true",
        help=f"Also show the pre-B-factory measurements ({', '.join(sorted(PRE_BF_MEASURES))}), when they are "
        "members of the chosen combination. Hidden by default: their large uncertainties would otherwise force a "
        "much wider x-axis range, making the rest of the figure hard to read.",
    )
    parser.add_argument(
        "--run-all",
        action="store_true",
        help="Run all available combinations (ignores --obs, --comb, --pre-bf, --arxiv and --pub).",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    setup_matplotlib(usetex=args.latex)

    if args.run_all:
        for obs, combs in COMBINATIONS.items():
            for comb in combs:
                for pre_bf in [False, True]:
                    for arxiv, pub in [(False, False), (True, False), (False, True)]:
                        plot_average(obs, comb, args.outdir, pre_bf=pre_bf, arxiv=arxiv, pub=pub)
    else:
        plot_average(args.obs, args.comb, args.outdir, pre_bf=args.pre_bf, arxiv=args.arxiv, pub=args.pub)
