import argparse
import importlib
import os
import re
import subprocess
from collections.abc import Callable
from contextlib import contextmanager
from dataclasses import KW_ONLY, dataclass, field
from enum import Enum
from multiprocessing import Pool
from pathlib import Path
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from gc_core.mpl_tools import (
    lhcb_2d_cols,
    lhcb_cols,
    lhcb_ls,
    plot1d,
    plot2d,
    print_interval,
    read_gc_scan,
)
from matplotlib.ticker import MaxNLocator
from scipy.stats import chi2

# Hardcoded constants --------------------------------------------------------------------------------------------------

repo_path = Path(__file__).resolve().parents[2]


class DYFsc(Enum):
    """Enumeration for the final-state-correction hypotheses for DeltaY(h- h+).

    Values must match those defined in `str_repr` in src/CharmUtils.cpp.
    """

    NONE = "no"
    PARTIAL = "partial"
    FULL = "full"


DY_FSC_LABELS: dict[DYFsc, str] = {
    DYFsc.NONE: r"$\Delta Y_{f} = - x_{12}\sin\phi^M_2$",
    DYFsc.PARTIAL: r"$\Delta Y_{f} = - x_{12}\sin\phi^M_2 + y_{12} a^{\rm d}_f$",
    DYFsc.FULL: r"$\Delta Y_{f} = - x_{12}\sin\phi^M_2 + y_{12} a^{\rm d}_f (1 + \cot\delta_f \; x_{12} / y_{12})$",
}


class AcpParam(Enum):
    """Enumeration for the aCP(h- h+) parametrisations.

    Values must match those defined in `str_repr` in src/CharmUtils.cpp.
    """

    ACP_DY = "acp-dy"
    ACP_COT = "acp-cot"
    R_DELTA = "r-delta"


class MixParam(Enum):
    """Enumeration for the mixing parametrisations.

    Values must match those defined in `str_repr` in src/CharmUtils.cpp.
    """

    PHENO = "pheno"
    THEO = "theo"
    D0_TO_KPI = "d0-to-kpi"


# Utility functions and classes for 1D and 2D scans and plots ----------------------------------------------------------


def _bring_angle_to_pi_pi(x: float) -> float:
    """Shift the angle x (in radians) to the range [-pi, pi]."""
    return (x + np.pi) % (2 * np.pi) - np.pi


@dataclass(frozen=True)
class Range:
    r"""Representation of a range of values in arbitrary units.

    Provides access also to the range in absolute units, a latex-compatible label for the units, and transformations to
    convert values between the two units.

    Attribute:
        range: Range in `unit` units.
        unit: Units that the range is defined in (e.g. "deg", "mrad", "%", "1e-3").
        unit_label: LaTeX-compatible label for the units (e.g. r"^\circ", r"\%", r"10^{-3}").
        abs_range: Range in absolute units used by charm-fitter (e.g. rad for angles).
        transf: Transformation to convert numbers from native charm-fitter units to the given unit.
        inv_transf: Inverse transformation to convert numbers from the given unit to native charm-fitter units.
    """

    range: tuple[float, float]
    _: KW_ONLY
    unit: str = ""
    unit_label: str = ""
    abs_range: tuple[float, float] = field(init=False)
    transf: Callable[[float], float] = field(init=False, default=staticmethod(lambda x: x))
    inv_transf: Callable[[float], float] = field(init=False, default=staticmethod(lambda x: x))

    def __post_init__(self):
        unit = self.unit.strip().lower()
        if unit in ["deg", "degree", "degrees"]:
            object.__setattr__(self, "transf", lambda x: np.degrees(x))
            object.__setattr__(self, "inv_transf", lambda x: np.radians(x))
            if not self.unit_label:
                object.__setattr__(self, "unit_label", r"^\circ")
        elif unit == "%":
            object.__setattr__(self, "transf", lambda x: x * 100.0)
            object.__setattr__(self, "inv_transf", lambda x: x / 100.0)
            if not self.unit_label:
                object.__setattr__(self, "unit_label", r"\%")
        elif unit in ["0.001", "mrad"] or re.fullmatch(r"1e-(\d+)", unit) is not None:
            if unit in ["0.001", "mrad"]:
                exp = 3
            else:
                match = re.fullmatch(r"1e-(\d+)", unit)
                assert match is not None
                exp = int(match.group(1))
            scale = 10.0**-exp
            object.__setattr__(self, "transf", lambda x: x / scale)
            object.__setattr__(
                self, "inv_transf", lambda x: _bring_angle_to_pi_pi(x * scale) if unit == "mrad" else x * scale
            )
            if not self.unit_label:
                object.__setattr__(self, "unit_label", "mrad" if unit == "mrad" else rf"10^{{-{exp}}}")
        elif unit in ["", "rad"]:
            if not self.unit_label:
                object.__setattr__(self, "unit_label", unit)
        else:
            raise ValueError(f"Unsupported unit {unit!r}")
        object.__setattr__(self, "abs_range", tuple(self.inv_transf(v) for v in self.range))

    def transformed_range(self, transf: Callable[[float], float]) -> "Range":
        """Transform the range using a transformation that acts on the variable defined in absolute units.

        Args:
            transf: Transformation to be applied to the range.

        Returns:
            Transformed range.
        """
        x1, x2 = transf(self.abs_range[0]), transf(self.abs_range[1])
        x_min, x_max = min(x1, x2), max(x1, x2)
        return Range(range=(self.transf(x_min), self.transf(x_max)), unit=self.unit)


# Classes for parameter scans ------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class _Scan:
    """Base class for a 1D or 2D scan of charm-fitter parameters.

    Attribute:
        parfile: Start parameter files to be used to fix the parameters.
        extra_args: Extra arguments to be passed to GammaComboEngine.
        dy_fsc: Final-state correction for DeltaY(h- h+) to be used.
        acp_param: ACP parametrisation that should be used.
        mix_param: Mixing parametrisation that should be used.
        fixed_pars: Parameters to be fixed in the subcombination.
        fix_parfile: Parameter file from which the values of the fixed parameters should be read.
    """

    _: KW_ONLY
    parfile: str = ""
    extra_args: str = ""
    dy_fsc: DYFsc = DYFsc.NONE
    acp_param: AcpParam = AcpParam.ACP_DY
    mix_param: MixParam = MixParam.THEO
    fixed_pars: list[str] = field(default_factory=list)
    fix_parfile: str = ""

    def __post_init__(self):

        # List of parameters that exist only in specific parametrisations of aCP(h- h+).
        ACP_DY_PARAMS = ["DY_KK", "DY_PP"]
        ACP_COT_PARAMS = ["cot_delta_KK", "cot_delta_PP"]
        R_DELTA_PARAMS = ["r_KK", "r_PP", "delta_KK", "delta_PP"]

        # List of parameters that exist only in specific parametrisations of mixing.
        PHENO_PARAMS = ["x", "y", "qop", "phi"]
        THEO_PARAMS = ["x12", "y12", "phiM", "phiG"]
        D0_TO_KPI_PARAMS = ["yp", "dyp", "xp2", "dxp2"]

        # List of parameters defined in CharmParameters.cpp
        SCAN_PARAMS = (
            ACP_DY_PARAMS
            + ACP_COT_PARAMS
            + R_DELTA_PARAMS
            + PHENO_PARAMS
            + THEO_PARAMS
            + D0_TO_KPI_PARAMS
            + [
                "Acp_KK",
                "Acp_PP",
                "Acp_KP",
                "Delta_Kpi",
                "Delta_Kpipi0",
                "Delta_K3pi",
                "r_Kpi",
                "r_K3pi",
                "r_K3pipi0",
                "k_K3pi",
                "k_Kpipi0",
                "F_pipipi0",
            ]
        )

        pars = self._pars()

        if any(p not in SCAN_PARAMS for p in pars):
            raise ValueError(f"Invalid parameter(s) {', '.join(p for p in pars if p not in SCAN_PARAMS)}")

        if len(pars) > 1:
            assert not (any(p in PHENO_PARAMS for p in pars) and any(p in THEO_PARAMS for p in pars)), (
                f"Cannot mix phenomenological and theoretical parameters in the same scan ({' , '.join(pars)})"
            )
            assert (
                sum(
                    [
                        any(p in ACP_DY_PARAMS for p in pars),
                        any(p in ACP_COT_PARAMS for p in pars),
                        any(p in R_DELTA_PARAMS for p in pars),
                    ]
                )
                <= 1
            ), f"Parameters require incompatible aCP parametrisations ({' , '.join(pars)})"

        if any(p in ACP_DY_PARAMS for p in pars):
            object.__setattr__(self, "acp_param", AcpParam.ACP_DY)
        if any(p in ACP_COT_PARAMS for p in pars):
            object.__setattr__(self, "acp_param", AcpParam.ACP_COT)
        if any(p in R_DELTA_PARAMS for p in pars):
            object.__setattr__(self, "acp_param", AcpParam.R_DELTA)

        if any(p in PHENO_PARAMS for p in pars):
            object.__setattr__(self, "mix_param", MixParam.PHENO)
        if any(p in D0_TO_KPI_PARAMS for p in pars):
            object.__setattr__(self, "mix_param", MixParam.D0_TO_KPI)

        if self.fixed_pars and not self.fix_parfile:
            raise ValueError("`fix_parfile` must be specified whenever `fixed_pars` is used")

        extra_args = self.extra_args
        if self.dy_fsc != DYFsc.NONE:
            extra_args += f" --dy-fsc {self.dy_fsc.value}"
        if self.acp_param != AcpParam.ACP_DY:
            extra_args += f" --acp-param {self.acp_param.value}"
        if self.mix_param == MixParam.PHENO:
            extra_args += " --mix pheno"
        if self.parfile:
            extra_args += f" --parfile {self.parfile}"
        if self.fixed_pars:
            pars_string = ",".join(self.fixed_pars)
            extra_args += f" --fix-from-parfile {pars_string} --fix-parfile {self.fix_parfile}"
        object.__setattr__(self, "extra_args", extra_args.strip())


@dataclass(frozen=True)
class Scan1D(_Scan):
    """Collection of parameters to be used in a 1D parameter scan.

    Attribute:
        range: Range for the parameter to be scanned.
    """

    par: str
    range: Range

    def _pars(self) -> list[str]:
        return [self.par]


@dataclass(frozen=True)
class Scan2D(_Scan):
    """Collection of parameters to be used in a 2D parameter scan.

    Attribute:
        xrange: Range for the x-axis parameter.
        yrange: Range for the y-axis parameter.
    """

    xpar: str
    ypar: str
    xrange: Range
    yrange: Range

    def _pars(self) -> list[str]:
        return [self.xpar, self.ypar]


# Combiner class ----------------------------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class Combiner:
    """Representation of a combiner, along with the 1D and 2D scans that should be performed with it.

    Args:
        id: Numerical ID of the combiner, or list of the IDs of the PDFs that it should be built from.
        title: Title that will be used in plots legends.
        title_long: Longer version of the title for e.g. the subtitle under the logo "CharmFitter". Defaults to `title`.
        col: Color to be used in plots.
        ls: Line style to be used in plots.
        marker: Marker style to be used in plots.
        scans_1d: Dictionary of parameter names to 1D scan parameters.
        scans_2d: Dictionary of (xparname, yparname) to 2D scan parameters.
    """

    id: int | list[int]
    title: str
    _: KW_ONLY
    title_long: str = ""
    col: str = ""
    ls: str = ""
    marker: str = ""
    scans_1d: list[Scan1D] = field(default_factory=list)
    scans_2d: list[Scan2D] = field(default_factory=list)

    def __post_init__(self):
        if not self.title_long:
            object.__setattr__(self, "title_long", self.title)


# Plotting classes --------------------------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class PlotParameter:
    """Representation of a parameter to be plotted.

    Attribute:
        name: Name of the parameter (for the output file name).
        cf_par: Name of the charm-fitter parameter that the parameter is calculated from.
        title: LaTeX-formatted title for plotting (excluding the units).
        transf: Transformation applied to values to move from `cf_var` to the parameter.
    """

    name: str
    _: KW_ONLY
    cf_par: str = ""
    title: str = ""
    transf: Callable[[float], float] = lambda x: x

    def __post_init__(self):
        if not self.cf_par:
            object.__setattr__(self, "cf_par", self.name)


@dataclass(frozen=True)
class _Plot:
    """Collection of parameters for making a 1D or 2D plot.

    Attribute:
        combiners: List of combiner ids to be included in the plot.
        combiners_label: Label identifying the set of combiners for the output folder.
        legpos: Legend position.
        legfill: Whether to fill the legend box.
        legfontsize: Font size for the legend.
        logo: Position of the "CharmFitter" logo.
        dy_fsc: List of DeltaY(h- h+) final-state-correction hypotheses to be included in the plot.
    """

    combiners: list[str]
    _: KW_ONLY
    combiners_label: str = ""
    legpos: str = "l"
    legfill: bool = False
    legfontsize: int = -1
    logo: str = "r"
    dy_fsc: list[DYFsc] = field(default_factory=lambda: [DYFsc.NONE])

    def __post_init__(self):
        if not self.combiners_label:
            object.__setattr__(self, "combiners_label", "-".join(self.combiners))
        assert len(self.dy_fsc) == 1 or len(self.combiners) == 1, (
            "dy_fsc comparisons can be made only with a single combiner"
        )


@dataclass(frozen=True)
class Plot1D(_Plot):
    """Collection of parameters for making a 1D plot.

    Attribute:
        param: Name of the parameter to be plotted.
        range: Parameter range. The plot will use the same units used to define the range.
    """

    par: str
    range: Range


@dataclass(frozen=True)
class Plot2D(_Plot):
    """Collection of parameters for making a 2D plot.

    In all functions, all optional parameters are overwritten by their defaults, if not explicitly given
    (see Parameter for xrange and yrange, and the configuration file for the list of combiners).

    Attribute:
        params: Names of the two parameters to be plotted.
        xrange: Range for the x-axis parameter. The plot will use the same units used to define the range.
        yrange: Range for the y-axis parameter. The plot will use the same units used to define the range.
        levels: Number of contour levels to be drawn in sigma-equivalent units.
    """

    pars: tuple[str, str]
    xrange: Range
    yrange: Range
    _: KW_ONLY
    levels: int = 3


@dataclass
class Plotter:
    """Main class to make 1D and 2D Matplotlib plots."""

    def _default_extensions():
        return ["pdf"]

    _: KW_ONLY
    dim: int = 1
    save: Path | None = None
    xtitle: str | None = None
    ytitle: str | None = None
    xlabelpad: float = 4.0
    ylabelpad: float = 4.0
    xrange: tuple[float, float] | None = None
    yrange: tuple[float, float] | None = None
    xtransf: Callable[[float], float] | None = None
    ytransf: Callable[[float], float] | None = None
    levels: int = 2
    logo: str = "l"
    logo_subtitle: str = ""
    legpos: str = ""
    legfill: bool = False
    legfontsize: int = -1
    cls: str = ""
    axes_origin: tuple[int, int] = (-1, -1)
    prune_xlabel: bool = False
    prune_ylabel: bool = False

    # Mutable defaults must use default_factory
    scanpoints: list[list[float]] = field(default_factory=list)
    bfs: list[float | None] = field(default_factory=list)
    legtitles: list[str] = field(default_factory=list)
    lopts: list[Any] = field(default_factory=list)
    fopts: list[Any] = field(default_factory=list)
    mopts: list[Any] = field(default_factory=list)
    extensions: list[str] = field(default_factory=_default_extensions)

    def add_scan(
        self,
        scanname,
        pars,
        label,
        bf=False,
        col=None,
        hatch=None,
        lw=None,
        ls=None,
        marker=None,
    ):
        if scanname is not None:
            x, y, z, pt = get_scan_res(scanname, *pars)

            if self.xtransf is not None:
                x = self.xtransf(x)
                pt[0] = self.xtransf(pt[0])

            if self.ytransf is not None:
                y = self.ytransf(y)
                pt[1] = self.ytransf(pt[1])

        else:
            x = y = z = np.empty((2, 2)) * np.nan  # tricks a fake
            pt = None

        # best fit point
        if bf:
            self.bfs.append(pt)
        else:
            self.bfs.append(None)

        # scan points
        if self.dim == 1:
            self.scanpoints.append([x, y])
        elif self.dim == 2:
            self.scanpoints.append([x, y, z])

        # label
        self.legtitles.append(label)

        # 1d opts
        lw = lw or 1
        ls = ls or "-"
        if self.dim == 1:
            if col is not None:
                self.lopts.append(dict(c=lhcb_cols["d" + col], lw=lw, ls=lhcb_ls[ls]))
                if hatch is not None:
                    self.fopts.append(dict(ec=lhcb_cols[col], fc="none", hatch=hatch))
                else:
                    self.fopts.append(dict(ec="none", fc=lhcb_cols[col]))

        # 2d opts
        elif self.dim == 2:
            if col is not None:
                self.lopts.append(
                    dict(
                        colors=[lhcb_cols["d" + col]],
                        linewidths=[lw],
                        linestyles=[lhcb_ls[ls]],
                    )
                )
                self.fopts.append(dict(colors=lhcb_2d_cols[col][: self.levels][::-1], alpha=0.9))
                if bf:
                    self.mopts.append(dict(marker=marker or "o", color=lhcb_cols["d" + col]))
                else:
                    self.mopts.append(dict())

    def plot(self, interactive=False):
        legopts = {}

        leg_top_after_logo = 0.88
        xleg_right = 0.78
        if self.logo == "l" and self.legpos == "l":
            legopts = dict(
                bbox_to_anchor=(
                    0.01,
                    leg_top_after_logo - 0.07 * len(self.scanpoints),
                    0.2,
                    0.07 * len(self.scanpoints),
                ),
                loc="upper left",
            )
        elif self.logo == "r" and self.legpos == "r":
            legopts = dict(
                bbox_to_anchor=(
                    xleg_right,
                    leg_top_after_logo - 0.07 * len(self.scanpoints),
                    0.2,
                    0.07 * len(self.scanpoints),
                ),
                loc="upper right",
            )
        elif self.logo == "l" and self.legpos == "r":
            legopts = dict(
                bbox_to_anchor=(
                    xleg_right,
                    0.997 - 0.07 * len(self.scanpoints),
                    0.2,
                    0.07 * len(self.scanpoints),
                ),
                loc="upper right",
            )
        elif self.logo in ["r", "br"] and self.legpos == "l":
            legopts = dict(
                bbox_to_anchor=(
                    0.01,
                    0.997 - 0.07 * len(self.scanpoints),
                    0.2,
                    0.07 * len(self.scanpoints),
                ),
                loc="upper left",
            )
        elif self.legpos == "bl":
            legopts = dict(
                bbox_to_anchor=(0.01, 0.02, 0.2, 0.07 * len(self.scanpoints)),
                loc="lower left",
            )
        elif self.legpos == "br":
            legopts = dict(
                bbox_to_anchor=(0.8, 0.01, 0.2, 0.07 * len(self.scanpoints)),
                loc="lower right",
            )

        if self.legfill:
            legopts["frameon"] = True
            legopts["framealpha"] = 0.8
            legopts["facecolor"] = "w"
            legopts["edgecolor"] = "0.7"
            legopts["borderpad"] = 0.2

        if self.legfontsize != -1:
            legopts["fontsize"] = self.legfontsize

        fig, ax = plt.subplots()
        if self.axes_origin != (-1, -1):
            kwargs = {"axes_origin": self.axes_origin}
        else:
            kwargs = {}
        if self.dim == 1:
            plot1d(
                self.scanpoints,
                self.lopts,
                self.fopts,
                xtitle=self.xtitle,
                legtitles=self.legtitles,
                ax=ax,
                legopts=legopts,
                **kwargs,
            )
        elif self.dim == 2:
            NO_CPV_MS = "X"
            no_cpv_marker = any(
                var1 in str(self.save) and var2 in str(self.save)
                for var1, var2 in [("phiM", "phiG"), ("qop", "phi"), ("Acp_KK", "Acp_PP"), ("Acp_KK", "phiM")]
            )
            u_spin_line = "Acp_KK" in str(self.save) and "Acp_PP" in str(self.save)

            if no_cpv_marker:
                self.legtitles.append(r"No $C\!P$ violation")
                self.lopts.append({})
                self.fopts.append({})
                self.mopts.append({"marker": NO_CPV_MS, "color": "k"})

            if u_spin_line:
                self.legtitles.extend(["U-spin symmetry", "Improved U-spin prediction"])
                self.lopts.extend(
                    [
                        {"color": "g", "lw": 2, "ls": "--"},
                        {"color": "b", "lw": 2, "ls": ":"},
                    ]
                )
                self.fopts.extend([{}, {}])
                self.mopts.extend([{}, {}])

            if "mrad" in self.ytitle:
                self.ylabelpad = -2.0

            plot2d(
                self.scanpoints,
                self.lopts,
                self.fopts,
                self.mopts,
                title=[self.xtitle, self.ytitle],
                labelpad=(self.xlabelpad, self.ylabelpad),
                levels=self.levels,
                legtitles=self.legtitles,
                ax=ax,
                bf=self.bfs,
                cl2d=True,
                legopts=legopts,
                **kwargs,
            )

            if u_spin_line:
                # Exact U-spin
                xmin, xmax = ax.get_xlim()
                ymin, ymax = ax.get_ylim()
                xline_max = min(xmax, -ymin)
                xline_min = max(xmin, -ymax)
                ax.axline(
                    (xline_min, -xline_min),
                    (xline_max, -xline_max),
                    color="g",
                    linestyle="--",
                )

                # Improved U-spin (cf. https://indico.cern.ch/event/1440982/contributions/6530703/)
                # Amplitude ratio A(KK)/A(pipi), computed by scripts/u-spin.py from the BF(D0 -> KK/pipi) PDG values
                akk_over_app = 1.8091
                xmin, xmax = ax.get_xlim()
                ymin, ymax = ax.get_ylim()
                xline_max = min(xmax, -ymin / akk_over_app)
                xline_min = max(xmin, -ymax / akk_over_app)
                ax.axline(
                    (xline_min, -xline_min * akk_over_app),
                    (xline_max, -xline_max * akk_over_app),
                    color="b",
                    linestyle=":",
                )

            if no_cpv_marker:
                ax.plot([0], [0], NO_CPV_MS, color="k", ms=8)

        logo_xright = 0.97
        if self.logo == "l":
            logo_x, logo_y = 0.03, 0.92
        elif self.logo == "br":
            logo_x, logo_y = logo_xright, 0.06
            if self.logo_subtitle:
                logo_y += 0.04
        else:
            logo_x, logo_y = logo_xright, 0.92
        charm_fitter_logo(
            pos=(logo_x, logo_y), subtitle=self.logo_subtitle, ax=ax, ha="left" if self.logo == "l" else "right"
        )

        # limits
        if self.xrange is not None:
            ax.set_xlim(*self.xrange)
        if self.yrange is not None:
            ax.set_ylim(*self.yrange)

        # axis labels
        ax.set_xlabel(ax.get_xlabel(), loc="right")
        if self.dim == 1:
            ax.set_ylabel("$1-$CL", loc="top")
        elif self.dim == 2:
            ax.set_ylabel(ax.get_ylabel(), loc="top")
        if self.prune_xlabel:
            ax.xaxis.set_major_locator(MaxNLocator(prune="lower"))
        if self.prune_ylabel:
            ax.yaxis.set_major_locator(MaxNLocator(prune="lower"))

        # CL lines
        if self.cls:
            assert self.dim == 1, "CL lines only make sense for 1D plots"
            ax.axhline(chi2.sf(1, 1), c="k", ls=":", lw=1)
            ax.axhline(chi2.sf(4, 1), c="k", ls=":", lw=1)
            if self.cls == "l":
                ax.text(
                    0.02,
                    0.32,
                    r"$68.3\%$",
                    ha="left",
                    va="bottom",
                    transform=ax.transAxes,
                    fontsize=16,
                )
                ax.text(
                    0.02,
                    0.05,
                    r"$95.4\%$",
                    ha="left",
                    va="bottom",
                    transform=ax.transAxes,
                    fontsize=16,
                )
            elif self.cls == "r":
                ax.text(
                    0.98,
                    0.32,
                    r"$68.3\%$",
                    ha="right",
                    va="bottom",
                    transform=ax.transAxes,
                    fontsize=16,
                )
                ax.text(
                    0.98,
                    0.05,
                    r"$95.4\%$",
                    ha="right",
                    va="bottom",
                    transform=ax.transAxes,
                    fontsize=16,
                )

        self.save.parent.mkdir(parents=True, exist_ok=True)
        for ext in self.extensions:
            path = self.save.with_suffix(f".{ext}")
            fig.savefig(path)
            print(f"Saved plot: {str(path)!r}")

        if not interactive:
            fig.clf()

        plt.close()


# Configuration classes -----------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class PlottingConfig:
    """Configuration consumed by the `scripts/ws-combo.py` plotting script.

    A configuration file is expected to build one of these (or `CharmPlottingConfig`) and expose it as `CONFIG`.

    Attribute:
        parameters: Parameters made available for plotting; keyed by name after construction.
        baseline_combiners: Combiners used for the "1d"/"2d" actions. The keys should match the combiner names in the
            CPP executables if they are based on a core combiner.
        combiners: All combiners known to this configuration (baseline plus, for `CharmPlottingConfig`, breakdown
            combiners), used to resolve the combiner labels referenced by `plots_1d`/`plots_2d` entries.
        plots_1d: 1D plots to be produced.
        plots_2d: 2D plots to be produced.
    """

    parameters: list[PlotParameter]
    baseline_combiners: dict[str, Combiner]
    _: KW_ONLY
    plots_1d: list[Plot1D] = field(default_factory=list)
    plots_2d: list[Plot2D] = field(default_factory=list)
    combiners: dict[str, Combiner] = field(init=False, default_factory=dict)

    def __post_init__(self):
        object.__setattr__(self, "parameters", {p.name: p for p in self.parameters})
        object.__setattr__(self, "combiners", dict(self.baseline_combiners))


@dataclass(frozen=True)
class CharmPlottingConfig(PlottingConfig):
    """Configuration for `charm-combo.py`, additionally supporting breakdown and DY-FSC-comparison plots.

    Attribute:
        combiners_breakdown: Combiners used only for the "breakdown" action.
        plots_breakdown: 2D breakdown plots to be produced.
        plots_dy_fsc_1d: 1D plots comparing DeltaY(h- h+) final-state-correction hypotheses.
        plots_dy_fsc_2d: 2D plots comparing DeltaY(h- h+) final-state-correction hypotheses.
    """

    _: KW_ONLY
    combiners_breakdown: dict[str, Combiner] = field(default_factory=dict)
    plots_breakdown: list[Plot2D] = field(default_factory=list)
    plots_dy_fsc_1d: list[Plot1D] = field(default_factory=list)
    plots_dy_fsc_2d: list[Plot2D] = field(default_factory=list)

    def __post_init__(self):
        super().__post_init__()
        object.__setattr__(self, "combiners", {**self.combiners, **self.combiners_breakdown})


# Utility functions for producing scans, plots etc --------------------------------------------------------------------


def getfnames(prefix: str, xpar: str, ypar: str | None = None):
    yext = f"_{ypar}" if ypar is not None else ""
    fname = f"plots/scanner/{prefix}_{xpar}{yext}.root"

    bfname = fname.replace("scanner/", "par/").replace("_scanner", "").replace(".root", ".dat")

    if not Path(fname).exists():
        raise FileNotFoundError(f"Cannot find scan file {fname!r}")

    if not Path(bfname).exists():
        raise FileNotFoundError(f"Cannot find fit result file {bfname!r}")

    return fname, bfname


def print_cl(prefix: str, xpar: str, prob: bool = True):
    pref = prefix.split("scanner")[1]
    suff = "Prob" if prob else "Plugin"
    fname = f"plots/cl/clintervals{pref}_{xpar}_{suff}.py"
    if not Path(fname).exists():
        return

    try:
        print_interval(fname, nsigma=1)
    except Exception:
        return


def get_scan_res(prefix: str, xpar: str, ypar: str | None = None):
    pars = [xpar]
    if ypar is not None:
        pars.append(ypar)
        print_cl(prefix, xpar)

    fname, bfname = getfnames(prefix, xpar, ypar)

    header_str = f"{prefix} - {xpar}"
    if ypar is not None:
        header_str += f" , {ypar}"
    print(header_str)

    return read_gc_scan(fname, bfname, pars)


def _combiner_string(id: int | list[int]) -> str:
    """Get the cobiner string to be used by the combination executable."""
    if isinstance(id, int):
        return str(id)
    else:
        return "0:+" + ",+".join([str(i) for i in sorted(id)])


def _combiner_file_string(pdf_ids: list[int]):
    return "empty+" + "+".join([str(pdf_id) for pdf_id in sorted(pdf_ids)])


def _get_prefix(
    args_prefix: str,
    *,
    dy_fsc: DYFsc,
    acp_param: AcpParam,
    mix_param: MixParam,
    dcs_cpv: bool,
):
    """Get the prefix of the scanner ROOT file to read to produce Matplotlib plots."""
    if mix_param == MixParam.D0_TO_KPI:
        return args_prefix

    args_prefix = args_prefix.replace("_scanner", f"_{dy_fsc.value}-dyfsc_{acp_param.value}_{mix_param.value}_scanner")

    if dcs_cpv:
        args_prefix = args_prefix.replace("_scanner", "_dcs-cpv_scanner")

    return args_prefix


def _dcs_cpv_label(dcs_cpv: bool | None, args: argparse.Namespace, cfg: PlottingConfig):
    if dcs_cpv == args.dcs_cpv_default or dcs_cpv is None:
        return ""
    return ", {} {} 0$".format(cfg.parameters["Acp_KP"].title.strip()[:-1], r"\neq" if dcs_cpv else "=")


def _dcs_cpv_vals(args: argparse.Namespace, compare_dcs_hypos: bool) -> list[bool | None]:
    """Get the list of dcs_cpv values to be plotted/compared, for combos (e.g. ws-combo) that don't have this option."""
    if not hasattr(args, "dcs_cpv"):
        return [None]
    if compare_dcs_hypos:
        return [not args.dcs_cpv_default, args.dcs_cpv_default]
    return [args.dcs_cpv]


def _dcs_cpv_subdir(args: argparse.Namespace, compare_dcs_hypos: bool) -> str:
    """Get the output subdirectory encoding the dcs_cpv hypothesis of a plot."""
    if compare_dcs_hypos:
        return "comparison"
    if not hasattr(args, "dcs_cpv"):
        return ""
    return "with-dcs-cpv" if args.dcs_cpv else "no-dcs-cpv"


def scans_1d(combiners: list[Combiner], args: argparse.Namespace, dy_fsc_comparison: bool = False) -> None:
    """Run the 1D scans for all combiners (based on the scans_1d list of each combiner)."""

    if not (args.rescan or args.plugin):
        return

    cmds = []
    for combiner in combiners:
        for scan in combiner.scans_1d:
            if (scan.dy_fsc != DYFsc.NONE and not dy_fsc_comparison) or (
                scan.dy_fsc == DYFsc.NONE and dy_fsc_comparison
            ):
                continue
            par = scan.par
            if par == "Acp_KP" and not args.dcs_cpv:
                continue
            extra_args = " ".join([args.extra_args, scan.extra_args])
            if args.plugin:
                if args.submit:
                    extra_args += f" -a pluginbatch --ntoys 50 --nbatchjobs 200 {args.batchopts}"
                else:
                    extra_args += " -a plugin -j 1-200"
            scan_range = scan.range.abs_range
            cmd = (
                f"bin/{args.execfile} -c {_combiner_string(combiner.id):<31s}"
                f" --var {par:<12s}"
                f" --scanrange {scan_range[0]}:{scan_range[1]} {extra_args}"
            )
            cmds.append(cmd)

    if args.plugin and args.submit:
        for cmd in cmds:
            run_command(cmd)
        return
    else:
        run_commands(cmds)


def scans_2d(combiners: list[Combiner], args: argparse.Namespace, dy_fsc_comparison: bool = False) -> None:
    """Run the 2D scans for all combiners (based on the scans_2d list of each combiner)."""

    if not args.rescan:
        return

    cmds = []
    for combiner in combiners:
        for scan in combiner.scans_2d:
            if (scan.dy_fsc != DYFsc.NONE and not dy_fsc_comparison) or (
                scan.dy_fsc == DYFsc.NONE and dy_fsc_comparison
            ):
                continue
            (xpar, ypar) = scan.xpar, scan.ypar
            if "Acp_KP" in (xpar, ypar) and not args.dcs_cpv:
                continue
            extra_args = " ".join([args.extra_args, scan.extra_args])
            xrange = scan.xrange.abs_range
            yrange = scan.yrange.abs_range
            cmd = (
                f"bin/{args.execfile} -c {_combiner_string(combiner.id):<31s}"
                f" --var {xpar:<7s} --var {ypar:<7s}"
                f" --scanrange  {xrange[0]}:{xrange[1]}"
                f" --scanrangey {yrange[0]}:{yrange[1]}"
                f" {extra_args}"
            )
            cmds.append(cmd)
    run_commands(cmds)


def _add_unit_label(title: str, unit: str) -> str:
    if not unit:
        return title
    elif unit.isalpha():
        return title + f" [{unit}]"
    else:
        if title.endswith("$"):
            return title[:-1] + rf"\, [{unit}]$"
        else:
            return title + f" $[{unit}]$"


@dataclass
class _PlotterStyle:
    """Class to hold the colors, line styles and markers for plotting."""

    colors: list[str] = field(default_factory=list)
    line_styles: list[str] = field(default_factory=list)
    markers: list[str] = field(default_factory=list)


def _plot_style(user: _PlotterStyle, combiners: list[Combiner], n_plots: int, n_dim: int = 1) -> _PlotterStyle:
    """Get the list of colors, line styles or markers to be used for plotting.

    The first choice is given by user-provided options. If these are not provided, fall back to combiner properties
    (only if there are no duplicates). Otherwise, use the default options defined in the function body.

    Args:
        user: Plot options provided by the user.
        combiners: List of combiners that will be plotted.
        n_plots: Total number of plots that will be produced.
        n_dim: Dimension of the plot (1 or 2).
    """

    # Set default colors, line styles and markers.
    def_col = [
        ["b", "r", "y", "o", "g", "p", "lb"],
        ["lb", "r", "g", "y", "p"],
    ]
    def_ls = [
        ["longdash", "-", "shortdash", "dashdot", "dash3dots", "dash2dots", "dots"],
        ["-", "longdash", "shortdash", "dashdot", "dash3dots", "dash2dots", "dots"],
    ]
    def_mark = [
        ["o", "^", "s", "*", "D", "v", "<"],
        ["o", "x", "s", "*", "D", "v", "<"],
    ]
    if any(len(lst) < n_plots for lst in [def_col[n_dim - 1], def_ls[n_dim - 1], def_mark[n_dim - 1]]):
        raise ValueError(f"Not enough default plot options for {n_plots} plots in {n_dim}D")

    rt = _PlotterStyle(
        colors=[def_col[n_dim - 1][n_plots - 1 - i] for i in range(n_plots)],
        line_styles=[def_ls[n_dim - 1][n_plots - 1 - i] for i in range(n_plots)],
        markers=[def_mark[n_dim - 1][n_plots - 1 - i] for i in range(n_plots)],
    )

    if user.colors:
        if len(user.colors) < n_plots:
            raise ValueError(f"Not enough colors provided by the user ({len(user)} < {n_plots})")
        rt.colors = user.colors
    elif len(combiners) == n_plots and all(combiner.col for combiner in combiners):
        colors = [combiner.col for combiner in combiners]
        if len(set(colors)) == len(colors):
            rt.colors = colors

    if user.line_styles:
        if len(user.line_styles) < n_plots:
            raise ValueError(f"Not enough line styles provided by the user ({len(user.line_styles)} < {n_plots})")
        rt.line_styles = user.line_styles
    elif len(combiners) == n_plots and all(combiner.ls for combiner in combiners):
        line_styles = [combiner.ls for combiner in combiners]
        if len(set(line_styles)) == len(line_styles):
            rt.line_styles = line_styles

    if user.markers:
        if len(user.markers) < n_plots:
            raise ValueError(f"Not enough markers provided by the user ({len(user.markers)} < {n_plots})")
        rt.markers = user.markers
    elif len(combiners) == n_plots and all(combiner.marker for combiner in combiners):
        markers = [combiner.marker for combiner in combiners]
        if len(set(markers)) == len(markers):
            rt.markers = markers

    return rt


def plots_1d(
    plots: list[Plot1D],
    args: argparse.Namespace,
    cfg: PlottingConfig,
    *,
    plots_cat: str = "1d",
    compare_dcs_hypos: bool = False,
    colors: list[str] = [],
    line_styles: list[str] = [],
    dy_fsc_labels: dict[DYFsc, str] = DY_FSC_LABELS,
) -> None:
    """Produce a set of 1D CL plots.

    :param compare_dcs_hypos: Compare the results for different hypothesis of CPV in D0 -> K+ pi- decays
        (two distributions are plotted for each combiner).
    """

    for plot in plots:
        if plot.par == "Acp_KP" and (not args.dcs_cpv or compare_dcs_hypos):
            continue
        par = cfg.parameters[plot.par]

        plotter = Plotter(
            dim=1,
            save=args.savedir
            / plot.combiners_label
            / plots_cat
            / _dcs_cpv_subdir(args, compare_dcs_hypos)
            / f"{plot.par}.pdf",
            xtitle=_add_unit_label(par.title, plot.range.unit_label),
            xrange=plot.range.range,
            xtransf=lambda x: plot.range.transf(par.transf(x)),
            cls="r",
            logo=plot.logo,
            legpos=plot.legpos,
            legfill=plot.legfill,
            legfontsize=plot.legfontsize,
        )

        dcs_cpv_vals = _dcs_cpv_vals(args, compare_dcs_hypos)
        assert len(dcs_cpv_vals) == 1 or len(plot.dy_fsc) == 1, (
            "Comparison of different DeltaY(h- h+) final-state correction hypotheses and of allowing for CPV in "
            "DCS D0 -> K+ pi- decays are not allowed at the same time"
        )
        n_plots = len(plot.combiners) * len(dcs_cpv_vals) * len(plot.dy_fsc)

        opts = _plot_style(
            _PlotterStyle(colors=colors, line_styles=line_styles),
            [cfg.combiners[combiner_id] for combiner_id in plot.combiners],
            n_plots,
        )

        i = 0
        for cname in plot.combiners:
            combiner = cfg.combiners[cname]
            combiner_arg = cname if isinstance(combiner.id, int) else _combiner_file_string(combiner.id)
            for dcs_cpv in dcs_cpv_vals:
                for dy_fsc in plot.dy_fsc:
                    scan = next((x for x in combiner.scans_1d if x.par == par.cf_par and x.dy_fsc == dy_fsc), None)
                    if scan is None:
                        raise ValueError(f"Cannot find 1D scan for parameter {par.cf_par} in combiner {cname}")
                    prefix = _get_prefix(
                        args.prefix,
                        dy_fsc=scan.dy_fsc,
                        acp_param=scan.acp_param,
                        mix_param=scan.mix_param,
                        dcs_cpv=dcs_cpv,
                    )
                    label = (
                        (combiner.title + _dcs_cpv_label(dcs_cpv, args, cfg))
                        if len(plot.dy_fsc) == 1
                        else dy_fsc_labels[dy_fsc]
                    )
                    plotter.add_scan(
                        f"{prefix}_{combiner_arg}",
                        pars=[par.cf_par],
                        label=label,
                        col=opts.colors[i],
                        ls=opts.line_styles[i],
                    )
                    i += 1
        plotter.plot()


def plots_2d(
    plots: list[Plot2D],
    args: argparse.Namespace,
    cfg: PlottingConfig,
    *,
    plots_cat: str = "2d",
    compare_dcs_hypos: bool = False,
    colors: list[str] = [],
    line_styles: list[str] = [],
    markers: list[str] = [],
    dy_fsc_labels: dict[DYFsc, str] = DY_FSC_LABELS,
) -> None:
    """Produce a set of 2D CL plots.

    :param compare_dcs_hypos: Compare the results for different hypothesis of CPV in D0 -> K+ pi- decays
        (two distributions are plotted for each combiner).
    """
    for plot in plots:
        xname, yname = plot.pars
        if any(name == "Acp_KP" for name in [xname, yname]) and (not args.dcs_cpv or compare_dcs_hypos):
            continue

        xpar = cfg.parameters[xname]
        ypar = cfg.parameters[yname]

        plotter = Plotter(
            dim=2,
            save=args.savedir
            / plot.combiners_label
            / plots_cat
            / _dcs_cpv_subdir(args, compare_dcs_hypos)
            / f"{yname}-vs-{xname}.pdf",
            xtitle=_add_unit_label(xpar.title, plot.xrange.unit_label),
            ytitle=_add_unit_label(ypar.title, plot.yrange.unit_label),
            xtransf=lambda x: plot.xrange.transf(xpar.transf(x)),
            ytransf=lambda x: plot.yrange.transf(ypar.transf(x)),
            xrange=plot.xrange.range,
            yrange=plot.yrange.range,
            levels=plot.levels,
            logo=plot.logo,
            legpos=plot.legpos,
            legfill=plot.legfill,
            legfontsize=plot.legfontsize,
        )

        dcs_cpv_vals = _dcs_cpv_vals(args, compare_dcs_hypos)
        assert len(dcs_cpv_vals) == 1 or len(plot.dy_fsc) == 1, (
            "Comparison of different DeltaY(h- h+) final-state correction hypotheses and of allowing for CPV in DCS D0 -> K+ pi- decays are not allowed at the same time"
        )
        n_plots = len(plot.combiners) * len(dcs_cpv_vals) * len(plot.dy_fsc)

        opts = _plot_style(
            _PlotterStyle(colors=colors, line_styles=line_styles, markers=markers),
            [cfg.combiners[combiner_id] for combiner_id in plot.combiners],
            n_plots,
            n_dim=2,
        )

        i = 0
        for cname in plot.combiners:
            combiner = cfg.combiners[cname]
            combiner_arg = cname if isinstance(combiner.id, int) else _combiner_file_string(combiner.id)
            for dcs_cpv in dcs_cpv_vals:
                for dy_fsc in plot.dy_fsc:
                    scan = next(
                        (
                            x
                            for x in combiner.scans_2d
                            if x.xpar == xpar.name and x.ypar == ypar.name and x.dy_fsc == dy_fsc
                        ),
                        None,
                    )
                    if scan is None:
                        raise ValueError(
                            f"Cannot find 2D scan for parameters {xpar.name} and {ypar.name} in combiner {cname}"
                        )
                    prefix = _get_prefix(
                        args.prefix,
                        dy_fsc=scan.dy_fsc,
                        acp_param=scan.acp_param,
                        mix_param=scan.mix_param,
                        dcs_cpv=dcs_cpv,
                    )
                    plotter.add_scan(
                        f"{prefix}_{combiner_arg}",
                        pars=[xpar.name, ypar.name],
                        label=(combiner.title + _dcs_cpv_label(dcs_cpv, args, cfg))
                        if len(plot.dy_fsc) == 1
                        else dy_fsc_labels[dy_fsc],
                        col=opts.colors[i],
                        ls=opts.line_styles[i],
                        marker=opts.markers[i],
                        bf=(plots_cat != "breakdown" or i == n_plots - 1),
                    )
                    i += 1
        plotter.plot()


def parse_args(*, combo: str = "charm", dcs_cpv_default: bool = False) -> argparse.ArgumentParser:
    actions = {
        "all": "Run all actions",
        "1d": "Make the 1D plots",
        "2d": "Make the 2D plots ({} etc)".format(
            "x vs y, phi vs |q/p|" if combo == "charm" else "x'2 vs y', dx'2 vs dy'"
        ),
    }
    if combo == "charm":
        actions["dy-fsc"] = (
            "Make the 1D and 2D plots comparing the sensitivity for different hypotheses on the DeltaY(h- h+) "
            "final-state correction"
        )
        actions["breakdown"] = "Make the 2D breakdown plots, divided by D0 final-state (x vs y, phi vs |q/p| etc)"
    default_execname = f"{combo}-combo"

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

    N = max(len(k) for k in actions.keys())

    def _formatted_key(key: str) -> str:
        s = f'"{key}":'
        return f"{s:<{N + 3}s}"

    parser.add_argument(
        "-a",
        "--actions",
        default=[],
        action="append",
        choices=list(actions.keys()),
        help="\n".join(
            ["Run specific action(s). Can be parsed multiple times. Available actions:"]
            + [f"- {_formatted_key(k)} {v}" for k, v in actions.items()]
        ),
    )
    parser.add_argument(
        "-r", "--rescan", default=False, action="store_true", help="Rerun the scans (rather than just make plots)"
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=repo_path / f"config/{combo}-combo/plotting" / ("2024.py" if combo == "ws" else "WA-2025-06.py"),
        help="Path of the configuration file",
    )
    parser.add_argument("-e", "--execfile", type=str, default=default_execname, help="Executable name")
    parser.add_argument(
        "--extra-args",
        type=str,
        default="--ps 1",
        help="Additional arguments to be passed to the executable",
    )
    parser.add_argument(
        "-s",
        "--savedir",
        type=Path,
        default=repo_path / "plots/matplotlib" / default_execname,
        help="Output directory for Matplotlib plots",
    )
    parser.add_argument(
        "-p", "--prefix", type=str, default=f"{default_execname}_scanner", help="Prefix of scan file locations"
    )
    parser.add_argument(
        "-i", "--interactive", default=False, action="store_true", help="Show plots interactively at the end"
    )
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="Verbose output of scan commands")
    parser.add_argument("-P", "--plugin", default=False, action="store_true", help="Use plugin scans")
    parser.add_argument("-S", "--submit", default=False, action="store_true", help="Submit plugin batch jobs")
    parser.add_argument(
        "-B",
        "--batchopts",
        metavar="",
        default="",
        help="Additional batch submission options",
    )
    parser.add_argument(
        "--no-latex",
        dest="latex",
        default=True,
        action="store_false",
        help="Disable LaTeX in Matplotlib text processing (needed for, eg, Gitlab CI)",
    )
    if combo != "ws":
        parser.add_argument(
            "--dcs-cpv",
            action="store_true",
            dest="dcs_cpv",
            default=False,
            help="Allow for CP violation in doubly Cabibbo-suppressed D0 -> K+ pi- decays"
            + (" (default)" if dcs_cpv_default else ""),
        )
        parser.add_argument(
            "--no-dcs-cpv",
            action="store_false",
            dest="dcs_cpv",
            help="Do not allow for CP violation in doubly Cabibbo-suppressed D0 -> K+ pi- decays"
            + ("" if dcs_cpv_default else " (default)"),
        )
    args = parser.parse_args()
    args.dcs_cpv_default = dcs_cpv_default
    if combo != "ws" and args.dcs_cpv:
        args.extra_args += " --dcs-cpv"
    if not args.config.is_absolute():
        args.config = repo_path / args.config
    if not args.savedir.is_absolute():
        args.savedir = repo_path / args.savedir
    return args


# Generic utilities ----------------------------------------------------------------------------------------------------


@contextmanager
def cwd(path):
    oldpwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(oldpwd)


def charm_fitter_logo(
    pos: tuple[float, float] = (0.02, 0.88),
    subtitle: str = "",
    ax: matplotlib.axes.Axes | None = None,
    ha: str = "left",
) -> None:
    """Add the CharmFitter logo to the current Matplotlib plot (or ax, if given)."""
    ax = ax or plt.gca()
    props = dict(fc="none", ec="none", boxstyle="square,pad=0.1")
    font = {"family": "Times New Roman", "weight": 400}
    usetex = plt.rcParams["text.usetex"]
    ax.text(
        *pos,
        r"$\textsc{CharmFitter}$" if usetex else "CharmFitter",
        transform=ax.transAxes,
        size=16,
        ha=ha,
        bbox=props,
        fontdict=font,
        usetex=usetex,
    )
    if subtitle:
        ax.text(
            pos[0],
            pos[1] - 0.05,
            subtitle,
            transform=ax.transAxes,
            size=min(12.2, 12.2 * 20 / len(subtitle)),
            ha=ha,
            bbox=props,
            fontdict=font,
            usetex=usetex,
        )


def get_configuration(config_fpath: Path) -> PlottingConfig:
    """Load a configuration file as a Python module and return the `PlottingConfig` it exposes as `CONFIG`."""
    if not config_fpath.is_file():
        raise RuntimeError(f'Could not open configuration file "{config_fpath}"')
    spec = importlib.util.spec_from_file_location("module", str(config_fpath))
    cfg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cfg)
    return cfg.CONFIG


def run_command(cmd: str | list[str], verbose: bool = False, shell: bool = True) -> subprocess.CompletedProcess[bytes]:
    """Run a command in the shell via subprocess."""
    cmd_str = cmd if isinstance(cmd, str) else " ".join(cmd)
    print(f"Running: {cmd_str!r}")
    if shell and isinstance(cmd, list):
        cmd = " ".join(cmd)
    elif not shell and isinstance(cmd, str):
        cmd = cmd.split()
    return subprocess.run(cmd, shell=shell, capture_output=not verbose, check=True)


def run_commands(cmds: list[str]) -> list[subprocess.CompletedProcess[bytes]]:
    """Run a list of commands in parallel via subprocess."""
    pool = Pool()
    return pool.map(run_command, cmds)


def setup_matplotlib(*, style="lhcb", usetex: bool = True) -> None:
    """Set the style for matplotlib plots."""
    import matplotlib
    import matplotlib.pyplot as plt

    plt.style.use(importlib.resources.files("charm_fitter") / "styles" / f"{style}.mplstyle")
    plt.rcParams["text.usetex"] = usetex
    if usetex:
        # Fix problems with rendering of minus sign in PDF
        matplotlib.use("pgf")
        plt.rcParams.update({"pgf.texsystem": "pdflatex"})
