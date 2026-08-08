import argparse
import logging
from dataclasses import dataclass
from math import floor, log, sqrt
from pathlib import Path

from charm_fitter.utils import repo_path


@dataclass(frozen=True)
class Measurement:
    """Class to define the objects storing the information relative to a single measurement for BLUE combinations."""

    label: str
    arxiv: str
    val: float
    stat: float
    sys: float | None = None
    sys2: float | None = None

    _EXPERIMENT_COLORS = {
        "BaBar": "g",
        "Belle": "r",
        "BES": "g",
        "CDF": "m",
        "CLEO": "m",
        "CMS": "g",
        "LHCb": "b",
        "average": "k",
        "PDG": "k",
    }

    def __post_init__(self):
        color = next((c for s, c in self._EXPERIMENT_COLORS.items() if s in self.label), None)
        if color is None:
            raise RuntimeError(f"The label {self.label} is not supported")
        object.__setattr__(self, "color", color)

    def result_str(self, units: float = 1.0) -> str:
        def _ndigits_to_print(stat, sys=None, sys2=None):
            """Given three uncertainties, return the number of digits after the comma to be printed according to PDG
            conventions. Assumes that the input measurements have the right number of digits according to PDG
            conventions (where the number of digits is set by the least precise measurement).
            """

            # Special case for publications not adhering to PDG conventions
            if (
                (
                    self.label == "CLEO"
                    and (self.arxiv == "hep-ex/9705006" or (self.arxiv == "0906.3198" and stat > 0.2))
                )
                or (self.label == "Belle" and self.arxiv == "2103.09969" and stat > 0.1)
                or self.arxiv in ["1911.01114", "2105.01565", "2405.11606", "2411.00306", "La Thuile"]
            ):
                return 1

            min_unc = stat
            if sys:
                min_unc = min(min_unc, sys)
            if sys2:
                min_unc = min(min_unc, sys2)
            main_exp = floor(log(min_unc) / log(10))
            if main_exp > 0:
                logging.warning(f"There may be too many figures printed for the uncertainties ({stat}, {sys}, {sys2})")
                return 0
            min_unc = floor(min_unc * 10 ** (2 - main_exp))
            if min_unc < 354:
                ndig = 2
            else:
                ndig = 1
            logging.debug(
                f"The number of digits to be printed for ({stat}, {sys}, {sys2}) is {ndig} ({ndig - main_exp - 1} after the comma; main_exp = {main_exp})"
            )
            return ndig - main_exp - 1

        val = self.val / units
        stat = self.stat / units
        sys = self.sys / units if self.sys else None
        sys2 = self.sys2 / units if self.sys2 else None

        ndig = _ndigits_to_print(stat, sys, sys2)
        if ndig > 2:  # TODO (fix for measurements non conformant to PDG
            ndig = 2

        s = f"{{:+.{ndig}f}} $\\pm$ {{:.{ndig}f}}".format(val, stat)
        if sys:
            s += f" $\\pm$ {{:.{ndig}f}}".format(sys)
        if sys2:
            s += f" $\\pm$ {{:.{ndig}f}}".format(sys2)
        s = s.replace("-", "\u2013")
        return s

    def err(self) -> float:
        if self.sys is None:
            return self.stat
        else:
            err2 = (self.stat) ** 2 + (self.sys) ** 2
            if self.sys2:
                err2 += (self.sys2) ** 2
            return sqrt(err2)


def get_units_label(units: float) -> str:
    """Get the LaTeX label (e.g. "\\%" or "10^{-2}") for a given power-of-ten unit scale."""
    raw_exp = log(units) / log(10)
    exp = round(raw_exp)
    if abs(exp - raw_exp) > 1e-2:
        raise RuntimeError(f"Units {units} not supported")
    return r"\%" if exp == -2 else f"10^{{{exp}}}"


def blue_parser(default_outdir: str) -> argparse.ArgumentParser:
    """Create a argument parser for the BLUE scripts with the arguments shared by all scripts."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        default=repo_path / "plots" / "BLUE" / default_outdir,
        help="Output directory for saving the plots",
    )
    parser.add_argument(
        "--arxiv",
        default=False,
        action="store_true",
        help="Print arXiv IDs next to the measurements on the plot.",
    )
    parser.add_argument(
        "--no-latex",
        dest="latex",
        default=True,
        action="store_false",
        help="Disable LaTeX in Matplotlib text processing (needed for, e.g., Gitlab CI).",
    )
    return parser
