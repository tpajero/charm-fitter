import logging
from dataclasses import dataclass
from math import floor, log, sqrt
from pathlib import Path

repo_path = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class Measurement:
    """Class to define the objects storing the information relative to a single measurement for BLUE combinations."""

    label: str
    arxiv: str
    val: float
    stat: float
    sys: float | None = None
    sys2: float | None = None

    def __post_init__(self):
        if "BaBar" in self.label:
            color = "g"
        elif "Belle" in self.label:
            color = "r"
        elif "BES" in self.label:
            color = "g"
        elif "CDF" in self.label:
            color = "m"
        elif "CLEO" in self.label:
            color = "m"
        elif "CMS" in self.label:
            color = "g"
        elif "LHCb" in self.label:
            color = "b"
        elif any(s in self.label for s in ["average", "PDG"]):
            color = "k"
        else:
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

    def err(self):
        if self.sys is None:
            return self.stat
        else:
            err2 = (self.stat) ** 2 + (self.sys) ** 2
            if self.sys2:
                err2 += (self.sys2) ** 2
            return sqrt(err2)


class Parameter:
    """Class to define the objects related to a scan parameter in a single place."""

    def __init__(self, identifier, title, range1d, parametrisation, range2d=None, degrees=False):
        assert parametrisation in ["theo", "pheno", "both"]
        self.id = identifier
        self.title = title
        self.parametrisation = parametrisation
        self.range1d = range1d
        self.range2d = range2d if range2d else range1d
        self.degrees = degrees
        if self.degrees:
            self.title += r"\,[^\circ]$"


class Plot2d:
    """Class to define the properties of 2D plots."""

    def __init__(self, pars, logo="r", legpos="l"):
        self.pars = pars
        self.logo = logo
        self.legpos = legpos


class Subcombination:
    """Class to define the PDFs and title of a subcombination of measurements."""

    def __init__(self, title, pdfs):
        self.title = title
        self.pdfs = pdfs


def setup_matplotlib(*, style="lhcb", usetex: bool = True) -> None:
    """Set the style for matplotlib plots."""
    from importlib.resources import files

    import matplotlib
    import matplotlib.pyplot as plt

    plt.style.use(files("charm_fitter") / "styles" / f"{style}.mplstyle")
    plt.rcParams["text.usetex"] = usetex
    if usetex:
        # Fix problems with rendering of minus sign in PDF
        matplotlib.use("pgf")
        plt.rcParams.update({"pgf.texsystem": "pdflatex"})
