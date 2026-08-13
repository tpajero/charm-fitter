#!/usr/bin/env python3
"""Plot the values of all existing measurements of ACP and BF of D(s)+ -> eta(') h+ decays,
plus the PDG 2025 average.

To check what options are available, run:

   python d-to-etah.py -h

"""

import argparse
import itertools
from pathlib import Path

from charm_fitter.blue import Measurement, blue_parser, get_units_label, plot_measurements
from charm_fitter.utils import setup_matplotlib


def get_measures(meas_type: str, decay: str, date: str) -> list[Measurement]:
    """Get the list of measurements to be shown for a given summary plot.

    :param meas_type: 'BF' or 'ACP'
    :param decay: Decay identifier.
    :param date: Date identifier for a given summary plot.
    """

    measures = {
        "BF": {
            "dp-to-etapi": [
                Measurement(
                    "CLEO", 3.63e-3, 0.08e-3, 0.18e-3, 0.06e-3, arxiv="0906.3198"
                ),  # Using PDG 2025 for D+ -> K- pi+ pi+; original values (3.54e-3, 0.08e-3, 0.18e-3, 0.08e-3)
                Measurement("BES III", 3.790e-3, 0.070e-3, 0.068e-3, arxiv="1802.03119"),
                Measurement("PDG 2025", 3.77e-3, 0.09e-3, is_average=True),
            ],
            "dp-to-etappi": [
                Measurement(
                    "CLEO", 4.80e-3, 0.16e-3, 0.23e-3, 0.08e-3, arxiv="0906.3198"
                ),  # Using PDG 2025 for D+ -> K- pi+ pi+; original values (4.68e-3, 0.16e-3, 0.23e-3, 0.10e-3)
                Measurement(
                    "BES III", 5.12e-3, 0.14e-3, 0.02e-3, arxiv="1802.03119"
                ),  # TODO typo in the PDG result (one missing zero in sys. unc.)
                Measurement("PDG 2025", 4.97e-3, 0.19e-3, is_average=True),
            ],
            "dp-to-etak": [
                Measurement("Belle", 1.15e-4, 0.16e-4, 0.05e-4, 0.03e-4, arxiv="1107.0553"),
                Measurement("BES III 2018", 1.51e-4, 0.25e-4, 0.14e-4, arxiv="1802.03119"),
                Measurement("BES III 2025", 1.17e-4, 0.10e-4, 0.03e-4, arxiv="2506.15533"),
                Measurement("PDG 2025", 1.25e-4, 0.16e-4, is_average=True),
            ],
            "dp-to-etapk": [
                Measurement(
                    "Belle", 1.87e-4, 0.19e-4, 0.05e-4, 0.07e-4, arxiv="1107.0553"
                ),  # Using PDG 2025 for D+ -> eta' pi
                Measurement("BES III 2018", 1.64e-4, 0.51e-4, 0.24e-4, arxiv="1802.03119"),
                Measurement("BES III 2025", 1.88e-4, 0.15e-4, 0.11e-4, arxiv="2506.15533"),
                Measurement("PDG 2025", 1.85e-4, 0.20e-4, is_average=True),
            ],
            "ds-to-etapi": [
                Measurement(
                    "CLEO", 1.2e-2, 0.3e-2, 0.2e-2, 0.2e-2, arxiv="hep-ex/9705006"
                ),  # Using PDG 2025 for D+ -> omega pi+
                Measurement("CLEO", 1.67e-2, 0.08e-2, 0.06e-2, arxiv="1306.5363"),
                Measurement("Belle", 1.82e-2, 0.14e-2, 0.07e-2, arxiv="1307.6240"),
                Measurement(
                    "BES III", 1.715e-2, 0.018e-2, 0.026e-2, 0.032e-2, arxiv="2005.05072"
                ),  # Using PDG 2025 for D+ -> K- K+ pi+; original values (1.741e-2, 0.018e-2, 0.027e-2, 0.054e-2)
                Measurement(
                    "Belle", 1.874e-2, 0.010e-2, 0.058e-2, 0.051e-2, arxiv="2103.09969"
                ),  # Using PDG 2025 for D+ -> phi pi+; original values (1.900e-2, 0.010e-2, 0.059e-2, 0.068e-2)
                Measurement("PDG 2025", 1.67e-2, 0.09e-2, is_average=True),
            ],
            "ds-to-etappi": [
                Measurement("CLEO", 3.94e-2, 0.15e-2, 0.20e-2, arxiv="1306.5363"),
                Measurement("BES III", 37.8e-3, 0.4e-3, 2.1e-3, 1.2e-3, arxiv="2005.05072"),
                Measurement("PDG 2025", 3.94e-2, 0.25e-2, is_average=True),
            ],
            "ds-to-etak": [
                Measurement("CLEO", 1.76e-3, 0.33e-3, 0.09e-3, 0.10e-3, arxiv="0906.3198"),
                Measurement("BES III", 1.62e-3, 0.10e-3, 0.03e-3, 0.05e-3, arxiv="2005.05072"),
                Measurement(
                    "Belle", 1.73e-3, 0.05e-3, 0.05e-3, 0.05e-3, arxiv="2103.09969"
                ),  # Using PDG 2025 for D+ -> phi pi+; original values (1.75e-3, 0.05e-3, 0.05e-3, 0.06e-3)
                Measurement("PDG 2025", 1.73e-3, 0.08e-3, is_average=True),
            ],
            "ds-to-etapk": [
                Measurement("CLEO", 1.8e-3, 0.5e-3, 0.1e-3, 0.1e-3, arxiv="0906.3198"),
                Measurement("BES III", 2.68e-3, 0.17e-3, 0.17e-3, 0.08e-3, arxiv="2005.05072"),
                Measurement("PDG 2025", 2.64e-3, 0.24e-3, is_average=True),
            ],
        },
        "ACP": {
            "dp-to-etapi": [
                Measurement("CLEO", -2.0e-2, 2.3e-2, 0.3e-2, arxiv="0906.3198"),
                Measurement("Belle", 1.74e-2, 1.13e-2, 0.19e-2, arxiv="1107.0553"),
                Measurement(r"LHCb $\eta\to\gamma\gamma$, 6 fb$^{-1}$", -0.2e-2, 0.8e-2, 0.4e-2, arxiv="2103.11058"),
                Measurement(
                    r"LHCb $\eta\to\pi^+\pi^-\gamma$, 6 fb$^{-1}$",
                    3.4e-3,
                    6.6e-3,
                    1.6e-3,
                    0.5e-3,
                    arxiv="2204.12228",
                ),
                Measurement("PDG 2025", 3e-3, 5e-3, is_average=True),
            ],
            "dp-to-etappi": [
                Measurement("CLEO", -4.0e-2, 3.4e-2, 0.3e-2, arxiv="0906.3198"),
                Measurement("Belle", -0.12e-2, 1.12e-2, 0.17e-2, arxiv="1107.0553"),
                Measurement(
                    r"LHCb $\eta^{\prime}\to\pi^+\pi^-\gamma$, 3 fb$^{-1}$",
                    -6.1e-3,
                    7.2e-3,
                    5.3e-3,
                    1.2e-3,
                    arxiv="1701.01871",
                ),
                Measurement(
                    r"LHCb $\eta^{\prime}\to\pi^+\pi^-\gamma$, 6 fb$^{-1}$",
                    4.9e-3,
                    1.8e-3,
                    0.6e-3,
                    0.5e-3,
                    arxiv="2204.12228",
                ),
                Measurement("PDG 2025", 4.1e-3, 2.3e-3, is_average=True),
            ],
            "dp-to-etak": [
                Measurement(r"LHCb $\eta\to\gamma\gamma$, 6 fb$^{-1}$", -6e-2, 10e-2, 4e-2, arxiv="2103.11058"),
                Measurement("PDG 2025", -6e-2, 11e-2, is_average=True),
            ],
            "dp-to-etapk": [],
            "ds-to-etapi": [
                Measurement("CLEO", 1.1e-2, 3.0e-2, 0.8e-2, arxiv="1306.5363"),
                Measurement("Belle", 2e-3, 3e-3, 3e-3, arxiv="2103.09969"),
                Measurement(r"LHCb $\eta\to\gamma\gamma$, 6 fb$^{-1}$", 0.8e-2, 0.7e-2, 0.5e-2, arxiv="2103.11058"),
                Measurement(r"LHCb $\eta\to\pi^+\pi^-\gamma$, 6 fb$^{-1}$", 3.2e-3, 5.1e-3, 1.2e-3, arxiv="2204.12228"),
                Measurement("PDG 2025", 3.2e-3, 3.1e-3, is_average=True),
            ],
            "ds-to-etappi": [
                Measurement("CLEO", -2.2e-2, 2.2e-2, 0.6e-2, arxiv="1306.5363"),
                Measurement(
                    r"LHCb $\eta^{\prime}\to\pi^+\pi^-\gamma$, 3 fb$^{-1}$",
                    -8.2e-3,
                    3.6e-3,
                    2.2e-3,
                    2.7e-3,
                    arxiv="1701.01871",
                ),
                Measurement(
                    r"LHCb $\eta^{\prime}\to\pi^+\pi^-\gamma$, 6 fb$^{-1}$",
                    0.1e-3,
                    1.2e-3,
                    0.8e-3,
                    arxiv="2204.12228",
                ),
                Measurement("PDG 2025", -0.6e-3, 2.2e-3, is_average=True),
            ],
            "ds-to-etak": [
                Measurement("CLEO", 9.3e-2, 15.2e-2, 0.9e-2, arxiv="0906.3198"),
                Measurement("Belle", 2.1e-2, 2.1e-2, 0.4e-2, arxiv="2103.09969"),
                Measurement(r"LHCb $\eta\to\gamma\gamma$, 6 fb$^{-1}$", 0.9e-2, 3.7e-2, 1.1e-2, arxiv="2103.11058"),
                Measurement("PDG 2025", 1.8e-2, 1.9e-2, is_average=True),
            ],
            "ds-to-etapk": [
                Measurement("CLEO", 6.0e-2, 18.9e-2, 0.9e-2, arxiv="0906.3198"),
                Measurement("PDG 2025", 6e-2, 19e-2, is_average=True),
            ],
        },
    }

    if date == "2025":
        return measures[meas_type][decay]
    else:
        raise RuntimeError(f"The combination for {date} is not supported")


def get_decay_label(decay: str) -> str:
    dp = "D^+"
    ds = "D^+_s"
    eta = r"\eta"
    etap = r"\eta^{\prime}"
    pi = r"\pi^+"
    k = "K^+"

    lab = dp if "dp" in decay else ds
    lab += r" \to "
    lab += etap if any(e in decay for e in ["etapp", "etapk"]) else eta
    lab += " "
    lab += k if "k" in decay else pi
    return lab


def get_xrange_and_units(meas_type: str, decay: str, pub: bool) -> tuple[tuple[float, float], float]:
    """Get the x range and the units for a given summary plot."""
    xranges = {
        "BF": {
            "dp-to-etapi": (3.3e-3, 4.7e-3),
            "dp-to-etappi": (4.4e-3, 6.3e-3),
            "dp-to-etak": (0.8e-4, 3e-4),
            "dp-to-etapk": (1e-4, 3.8e-4),
            "ds-to-etapi": (0.7e-2, 3.6e-2),
            "ds-to-etappi": (3.4e-2, 5.1e-2),
            "ds-to-etak": (1.3e-3, 3.1e-3),
            "ds-to-etapk": (1e-3, 5.3e-3),
        },
        "ACP": {
            "dp-to-etapi": (-5e-2, 13e-2),
            "dp-to-etappi": (-8e-2, 12e-2),
            "dp-to-etak": (-20e-2, 30e-2),
            "dp-to-etapk": (-10e-2, 10e-2),
            "ds-to-etapi": (-2.5e-2, 12e-2),
            "ds-to-etappi": (-5e-2, 6e-2),
            "ds-to-etak": (-8e-2, 55e-2 if pub else 50e-2),
            "ds-to-etapk": (-15e-2, 60e-2),
        },
    }
    units = {
        "BF": {
            "dp-to-etapi": 1e-3,
            "dp-to-etappi": 1e-3,
            "dp-to-etak": 1e-4,
            "dp-to-etapk": 1e-4,
            "ds-to-etapi": 1e-2,
            "ds-to-etappi": 1e-2,
            "ds-to-etak": 1e-3,
            "ds-to-etapk": 1e-3,
        },
        "ACP": {
            "dp-to-etapi": 1e-2,
            "dp-to-etappi": 1e-2,
            "dp-to-etak": 1e-2,
            "dp-to-etapk": 1e-2,
            "ds-to-etapi": 1e-2,
            "ds-to-etappi": 1e-2,
            "ds-to-etak": 1e-2,
            "ds-to-etapk": 1e-2,
        },
    }
    return xranges[meas_type][decay], units[meas_type][decay]


def plot_average(meas_type: str, decay: str, date: str, out_dir: Path, *, arxiv: bool, pub: bool) -> None:
    measures = get_measures(meas_type, decay, date)
    if not measures:
        print(f"No measurements available for {meas_type} of {decay}. Will not produce the plot")
        return

    (x_min, x_max), units = get_xrange_and_units(meas_type, decay, pub)
    meas_label = "BF" if meas_type == "BF" else "A_{CP}"
    xlabel = rf"${meas_label}({get_decay_label(decay)})$ $[{get_units_label(units)}]$"

    plot_measurements(
        measures,
        (x_min / units, x_max / units),
        xlabel,
        out_dir / f"{meas_type.lower()}-{decay}-{date}.pdf",
        figsize=(6, 6),
        units=units,
        arxiv=arxiv,
        pub=pub,
    )


def parse_args() -> argparse.Namespace:
    parser = blue_parser("d-to-etah")
    parser.add_argument(
        "-c",
        "--date",
        type=str,
        default="2025",
        help="Date of the combination",
        choices=["2025"],
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    setup_matplotlib(usetex=args.latex)

    for m, d, e, h in itertools.product(["BF", "ACP"], ["dp", "ds"], ["eta", "etap"], ["pi", "k"]):
        plot_average(m, f"{d}-to-{e}{h}", args.date, args.outdir, arxiv=args.arxiv, pub=args.pub)
