import argparse
import logging
from collections.abc import Callable
from dataclasses import KW_ONLY, dataclass, field
from math import floor, log, sqrt
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt

from charm_fitter.utils import repo_path

# arXiv ID -> (publication, DOI)
PUBLICATIONS: dict[str, tuple[str, str]] = {
    "hep-ex/9705006": ("Phys. Rev. Lett. 79 (1997) 1436", "10.1103/PhysRevLett.79.1436"),
    "hep-ex/9903012": ("Phys. Rev. Lett. 83 (1999) 32", "10.1103/PhysRevLett.83.32"),
    "hep-ex/0004034": ("Phys. Lett. B 485 (2000) 62", "10.1016/S0370-2693(00)00694-8"),
    "hep-ex/0012054": ("Phys. Rev. D 63 (2001) 071101", "10.1103/PhysRevD.63.071101"),
    "hep-ex/0111024": ("Phys. Rev. D 65 (2002) 092001", "10.1103/PhysRevD.65.092001"),
    "0905.4185": ("Phys. Rev. D 80 (2009) 052006", "10.1103/PhysRevD.80.052006"),
    "0906.3198": ("Phys. Rev. D 81 (2010) 052013", "10.1103/PhysRevD.81.052013"),
    "1107.0553": ("Phys. Rev. Lett. 107 (2011) 221801", "10.1103/PhysRevLett.107.221801"),
    "1209.3896": ("Phys. Rev. D 87 (2013) 012004", "10.1103/PhysRevD.87.012004"),
    "1306.5363": ("Phys. Rev. D 88 (2013) 032009", "10.1103/PhysRevD.88.032009"),
    "1307.6240": ("JHEP 09 (2013) 139", "10.1007/JHEP09(2013)139"),
    "1410.5435": ("Phys. Rev. D 90 (2014) 111103", "10.1103/PhysRevD.90.111103"),
    "1501.01378": ("Phys. Lett. B 744 (2015) 339", "10.1016/j.physletb.2015.04.008"),
    "1501.06777": ("JHEP 04 (2015) 043", "10.1007/JHEP04(2015)043"),
    "1508.06087": ("JHEP 10 (2015) 055", "10.1007/JHEP10(2015)055"),
    "1509.08266": ("Phys. Lett. B 753 (2016) 412", "10.1016/j.physletb.2015.12.025"),
    "1701.01871": ("Phys. Lett. B 771 (2017) 21", "10.1016/j.physletb.2017.05.013"),
    "1702.06490": ("Phys. Rev. Lett. 118 (2017) 261803", "10.1103/PhysRevLett.118.261803"),
    "1802.03119": ("Phys. Rev. D 97 (2018) 072004", "10.1103/PhysRevD.97.072004"),
    "1810.06874": ("Phys. Rev. Lett. 122 (2019) 011802", "10.1103/PhysRevLett.122.011802"),
    "1911.01114": ("Phys. Rev. D 101 (2020) 012005", "10.1103/PhysRevD.101.012005"),
    "2005.05072": ("JHEP 08 (2020) 146", "10.1007/JHEP08(2020)146"),
    "2103.09969": ("Phys. Rev. D 103 (2021) 112005", "10.1103/PhysRevD.103.112005"),
    "2103.11058": ("JHEP 06 (2021) 019", "10.1007/JHEP06(2021)019"),
    "2105.01565": ("Phys. Rev. D 104 (2021) L031102", "10.1103/PhysRevD.104.L031102"),
    "2105.09889": ("Phys. Rev. D 104 (2021) 072010", "10.1103/PhysRevD.104.072010"),
    "2202.09106": ("Phys. Rev. D 105 (2022) 092013", "10.1103/PhysRevD.105.092013"),
    "2204.12228": ("JHEP 04 (2023) 081", "10.1007/JHEP04(2023)081"),
    "2405.06556": ("Phys. Rev. Lett. 133 (2024) 101803", "10.1103/PhysRevLett.133.101803"),
    "2405.11606": ("Eur. Phys. J. C 84 (2024) 1264", "10.1140/epjc/s10052-024-13244-0"),
    "2411.00306": ("Phys. Rev. D 111 (2025) 012015", "10.1103/PhysRevD.111.012015"),
    "2504.15881": ("Phys. Rev. D 112 (2025) 012017", "10.1103/8x1h-39dp"),
    "2506.15533": ("JHEP 10 (2025) 178", "10.1007/JHEP10(2025)178"),
    "2510.14732": ("JHEP 02 (2026) 253", "10.1007/JHEP02(2026)253"),
}


@dataclass(frozen=True)
class Measurement:
    """Class to define the objects storing the information relative to a single measurement for BLUE combinations."""

    label: str
    val: float
    stat: float
    sys: float | None = None
    sys2: float | None = None
    _: KW_ONLY
    is_average: bool = False
    arxiv: str = ""
    publication: str = field(init=False, default="")
    doi: str = field(init=False, default="")

    _EXPERIMENT_COLORS = {
        "BaBar": "g",
        "Belle": "r",
        "BES": "tab:orange",
        "CDF": "m",
        "CLEO": "m",
        "CMS": "g",
        "E791": "y",
        "FOCUS": "c",
        "LHCb": "b",
        "average": "k",
        "PDG": "k",
    }

    def __post_init__(self):
        color = next((c for s, c in self._EXPERIMENT_COLORS.items() if s in self.label), None)
        if color is None:
            raise RuntimeError(f"The label {self.label} is not supported")
        object.__setattr__(self, "color", color)
        if self.arxiv and self.arxiv not in PUBLICATIONS:
            raise RuntimeError(
                f"The arXiv ID {self.arxiv} is not supported. Please update the PUBLICATIONS dictionary in src/charm_fitter/blue.py"
            )
        publication, doi = PUBLICATIONS.get(self.arxiv, ("", ""))
        object.__setattr__(self, "publication", publication)
        object.__setattr__(self, "doi", doi)

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

    def reference_str(self, *, pub: bool = False) -> str:
        """Get a reference string for the measurement.

        If `pub` is True, print the publication reference linked to its DOI, reverting to the arXiv ID if the
        measurement has no publication set. Otherwise, print the arXiv ID. Links are available only if Matplotlib has
        latex enabled.
        """
        latex = plt.rcParams["text.usetex"]
        if pub and self.publication:
            if latex and self.doi:
                return rf"\href{{https://doi.org/{self.doi}}}{{{self.publication}}}"
            return self.publication
        elif self.arxiv:
            return (
                rf"\href{{https://arxiv.org/abs/{self.arxiv}}}{{\texttt{{arXiv:{self.arxiv}}}}}"
                if latex
                else f"arXiv:{self.arxiv}"
            )
        return ""

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


def text_size_pixels(text: str, fontsize: float) -> tuple[float, float]:
    """Get the size of the bounding box of a rendered text string in pixels."""
    if not text:
        return 0, 0
    fig = plt.figure()
    t = plt.text(0, 0, text, fontsize=fontsize)
    fig.canvas.draw()
    bbox = t.get_window_extent(fig.canvas.get_renderer())
    plt.close(fig)
    return bbox.width, bbox.height


def plot_measurements(
    measures: list[Measurement],
    xlabel: str,
    out_path: Path,
    *,
    figsize: tuple[float, float] | None = None,
    xrange: tuple[float, float] | None = None,
    xlabel_fontsize: float | None = None,
    units: float = 1.0,
    val_transform: Callable[[float], float] | None = None,
    color_fn: Callable[[Measurement], str] | None = None,
    arxiv: bool,
    pub: bool,
    add_subdir: bool = True,
) -> None:
    """Draw a summary plot of a set of measurements and save it to `out_path`.

    A horizontal divider line is drawn above every measurement flagged `is_average`, and if the last measurement is
    flagged `is_average` a dashed vertical line is drawn at its value.

    :param xrange: Plot x-axis range, already expressed in units of `units`. If None, the range is fit automatically.
    :param val_transform: Transformation applied to `meas.val` before converting to plotting units. Needed to pass from
        DeltaY to A_Gamma measurements.
    :param color_fn: Function returning the color to be used for a given measurement. If not provided, the color is set
        from the `color` member of the `Measurement` object.
    :param add_subdir: Whether to add a subdirectory to `out_path` depending on the `arxiv` and `pub` flags.
    """

    def _default_val_transform(v: float) -> float:
        return v / units

    val_transform = val_transform or _default_val_transform

    def _default_color(meas: Measurement) -> str:
        return meas.color

    color_fn = color_fn or _default_color

    n_meas = len(measures)
    fig, ax = plt.subplots(figsize=figsize or (6, n_meas))
    y_min, y_max = -0.5, n_meas - (0.5 if n_meas > 2 else 0.4)
    plt.ylim(y_min, y_max)
    plt.tick_params(axis="y", which="both", right=False, left=False, labelleft=False)
    plt.xlabel(xlabel, ha="center", **({"fontsize": xlabel_fontsize} if xlabel_fontsize is not None else {}))

    vals = [val_transform(meas.val) for meas in measures]
    errs = [meas.err() / units for meas in measures]
    data_lo = min(v - e for v, e in zip(vals, errs))
    data_hi = max(v + e for v, e in zip(vals, errs))

    fontsize = matplotlib.rcParams["font.size"]
    labels = [f"{meas.label}\n{meas.result_str(units)}" for meas in measures]
    ref_labels = [(meas.reference_str(pub=pub) if (arxiv or pub) else "") for meas in measures]
    ref_fontsize = fontsize * (0.8 if pub else 0.9)
    height_px = ax.bbox.height
    fontsize_yscale = fig.dpi / 72 * (y_max - y_min) / height_px

    # Fraction of x-axis range between the highest of all estimates intervals and text labels
    f_centre = 0.05

    if xrange is None:
        text_max_width_px = max(
            max(text_size_pixels(l, fontsize)[0] for l in labels),
            max(text_size_pixels(l, ref_fontsize)[0] for l in ref_labels),
        )
        f_text = text_max_width_px / (
            fig.bbox.width - plt.rcParams["figure.constrained_layout.w_pad"] * 2 * plt.rcParams["figure.dpi"]
        )
        f_left = 0.03
        f_right = 0.03
        data_span = data_hi - data_lo
        assert 1 - f_left - f_centre - f_text - f_right > 0.2, (
            "There is not enough space for showing the average itself. Please reduce the length of the labels"
        )
        conversion = data_span / (1 - f_left - f_centre - f_text - f_right)
        x_min = data_lo - f_left * conversion
        x_max = data_hi + (f_centre + f_text + f_right) * conversion
    else:
        x_min, x_max = xrange

    plt.xlim(x_min, x_max)
    x_text = data_hi + f_centre * (x_max - x_min)

    for i, meas in enumerate(measures):
        y = n_meas - 1 - i
        col = color_fn(meas)
        plt.errorbar(vals[i], y, xerr=errs[i], fmt=".", markersize=8, capsize=7, color=col)
        plt.errorbar(vals[i], y, xerr=meas.stat / units, capsize=5, color=col)

        y_text = y - 0.8 * fontsize_yscale * fontsize

        if ref_labels[i]:
            y_text += 0.5 * fontsize_yscale * ref_fontsize
            plt.text(
                x_text,
                y_text - fontsize_yscale * ref_fontsize,
                ref_labels[i],
                fontsize=ref_fontsize,
                color="deepskyblue",
            )
        plt.text(
            x_text,
            y_text,
            labels[i],
            fontsize=fontsize,
            color=col,
        )

    # Divider line(s) separating averaged/combined results from the individual measurements feeding into them
    for i, meas in enumerate(measures):
        if meas.is_average:
            y = n_meas - 1 - i
            plt.plot([x_min, x_max], [y + 0.5, y + 0.5], linestyle="-", linewidth=1, color="k")

    # Dashed vertical line at the combined/world-average value, if the last measurement shown is one
    if measures[-1].is_average:
        avg = vals[-1]
        plt.plot([avg, avg], [y_min, y_max], linestyle="--", linewidth=1, color="k")

    if add_subdir:
        subdir = "arxiv" if arxiv else "pub" if pub else "no-refs"
        out_path = out_path.parent / subdir / out_path.name

    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Need more than one iteration in constrained_layout mode to ensure correct rendering (see dy.py).
    # https://matplotlib.org/stable/users/explain/axes/constrainedlayout_guide.html -> other caveats
    fig.canvas.draw()
    fig.savefig(out_path)
    plt.close(fig)


def blue_parser(default_outdir: str) -> argparse.ArgumentParser:
    """Create a argument parser for the BLUE scripts with the arguments shared by all scripts."""
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    outdir = repo_path / "plots" / "BLUE" / default_outdir
    parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        default=outdir,
        help=f"Output directory for saving the plots (defaults to {str(outdir)!r})",
    )
    ref_group = parser.add_mutually_exclusive_group()
    ref_group.add_argument(
        "--arxiv",
        default=False,
        action="store_true",
        help="Print arXiv IDs (with web link) next to the measurements on the plot.",
    )
    ref_group.add_argument(
        "--pub",
        default=False,
        action="store_true",
        help="Print the publication reference (linked to its DOI) next to the measurements on the plot. "
        "Falls back to the arXiv ID for measurements not yet published",
    )
    parser.add_argument(
        "--no-latex",
        dest="latex",
        default=True,
        action="store_false",
        help="Disable LaTeX in Matplotlib text processing (needed for, e.g., Gitlab CI).",
    )
    return parser
