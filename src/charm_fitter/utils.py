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
from types import ModuleType
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

repo_path = Path(__file__).resolve().parents[2]


# CharmFitter scans and plotting ---------------------------------------------------------------------------------------


class MixParam(Enum):
    """Enumeration for the mixing parametrisations."""

    PHENO = "pheno"
    THEO = "theo"
    D0_TO_KPI = "d0-to-kpi"


class DYFsc(Enum):
    """Enumeration for the DeltaY final-state-correction hypotheses."""

    NONE = "no"
    PARTIAL = "partial"
    FULL = "full"


@dataclass(frozen=True)
class PlotOpts1D:
    _: KW_ONLY
    legpos: str = "l"
    legfill: bool = False


@dataclass(frozen=True)
class Parameter:
    """Representation of a parameter to be scanned and plotted.

    Args:
        name: Id of the parameter.
        scan_range: Tuple defining (min, max) scan range.
        plot_range: Tuple defining plotting range (defaults to scan_range).
        scan_range_2d: Range used for 2D scans (defaults to scan_range).
        plot_range_2d: Range used for 2D plots (defaults to scan_range_2d).
        title: LaTeX-formatted title for plotting (excluding the units).
        unit: Unit of the parameter (e.g. "deg", "%"). All ranges are defined in terms of these units.
        description: Human-readable description.
        transf: Optional transformation applied to values to move from scan units to plot units (and possibly variables).
        mix_params: List of mixing parametrisations that can be used when scanning the parameter likelihood.
        dy_fsc_param: List of DeltaY final-state-correction hypotheses that can be used when scanning the parameter likelihood.
        plot_opts_1d: Options for 1D plots.
        scan_1d: Sets whether a 1D scan for the parameter should be performed, or if the parameter depends on a
            different one that will be independently scanned.
    """

    def _default_mix():
        return [MixParam.PHENO, MixParam.THEO]

    def _default_dyfsc():
        return [DYFsc.NONE, DYFsc.PARTIAL]

    name: str
    scan_range: tuple[float, float]
    _: KW_ONLY
    plot_range: tuple[float, float] | None = None
    scan_range_2d: tuple[float, float] | None = None
    plot_range_2d: tuple[float, float] | None = None
    title: str = ""
    unit: str = ""
    description: str = ""
    transf: Callable[[float], float] | None = None
    mix_params: list[MixParam] = field(default_factory=_default_mix)
    dy_fsc_param: list[DYFsc] = field(default_factory=_default_dyfsc)
    plot_opts_1d: PlotOpts1D = PlotOpts1D()
    scan_1d: bool = True

    def __post_init__(self):
        scan_range = self.scan_range
        plot_range = self.plot_range if self.plot_range is not None else self.scan_range
        scan_range_2d = self.scan_range_2d if self.scan_range_2d is not None else self.scan_range
        plot_range_2d = self.plot_range_2d if self.plot_range_2d is not None else scan_range_2d
        unit = self.unit.strip()
        title = self.title.strip()
        object.__setattr__(self, "latex_name", title)
        default_transf = self.transf if self.transf is not None else lambda x: x
        set_transf = False

        if unit != "":
            set_transf = True
            if unit.lower() in ["deg", "degree", "degrees"]:
                object.__setattr__(self, "transf", lambda x: np.degrees(default_transf(x)))
                scan_range = tuple(np.radians(x) for x in scan_range)
                scan_range_2d = tuple(np.radians(x) for x in scan_range_2d)
                unit = r"^\circ"
            elif unit == "%":
                object.__setattr__(self, "transf", lambda x: default_transf(x) * 100.0)
                scan_range = tuple(x / 100.0 for x in scan_range)
                scan_range_2d = tuple(x / 100.0 for x in scan_range_2d)
                unit = r"\%"
            elif unit in ["0.001", "mrad"] or re.fullmatch(r"1e-(\d+)", unit) is not None:
                if unit in ["0.001", "mrad"]:
                    exp = 3
                else:
                    match = re.fullmatch(r"1e-(\d+)", unit)
                    assert match is not None
                    exp = int(match.group(1))
                scale = 10.0**exp
                object.__setattr__(
                    self,
                    "transf",
                    lambda x: (
                        (_bring_angle_to_pi_pi(default_transf(x)) if unit == "mrad" else default_transf(x)) * scale
                    ),
                )
                scan_range = tuple(x / scale for x in scan_range)
                scan_range_2d = tuple(x / scale for x in scan_range_2d)
                if unit != "mrad":
                    unit = rf"10^{{-{exp}}}"
            else:
                set_transf = False

            if unit.isalpha():
                title += f" [{unit}]"
            else:
                if title.endswith("$"):
                    title = title[:-1] + rf"\, [{unit}]$"
                else:
                    title += f" $[{unit}]$"

        object.__setattr__(self, "scan_range", scan_range)
        object.__setattr__(self, "plot_range", plot_range)
        object.__setattr__(self, "scan_range_2d", scan_range_2d)
        object.__setattr__(self, "plot_range_2d", plot_range_2d)
        object.__setattr__(self, "title", title)
        object.__setattr__(self, "unit", unit)
        if set_transf is False:
            object.__setattr__(self, "transf", default_transf)


@dataclass(frozen=True)
class ScanParams:
    """Collects a bunch of parameters that should be changed in scans of (sub)combinations.

    Args:
        fixed_pars: Parameters to be fixed in the subcombination.
        parfile: Start parameter files to be used to fix the parameters.
    """

    _: KW_ONLY
    fixed_pars: list[str] = field(default_factory=list)
    parfile: str | None = None
    extra_opts: str = ""


@dataclass(frozen=True)
class ScanParams1D(ScanParams):
    par: str
    range: tuple[float, float]


@dataclass(frozen=True)
class ScanParams2D(ScanParams):
    pars: tuple[str, str]
    xrange: tuple[float, float]
    yrange: tuple[float, float]
    _: KW_ONLY
    mix_param: MixParam | None = None


@dataclass(frozen=True)
class Plot2D:
    """Collects the parameters for a 2D plot.

    In all functions, all optional parameters are overwritten by their defaults, if not explicitly given
    (see Parameter for xrange and yrange, and the configuration file for the list of combiners).

    Args:
        params: Parameters to be plotted.
        combiners: List of combiners to be included in the plot.
        xrange: Range for the x-axis.
        yrange: Range for the y-axis.
        scan: Whether to perform a scan for the given parameters.
        label: Optional label identifying the set of combiners for the output folder.
        legpos: Position of the legend.
        legfill: Whether to fill the legend box.
        levels: Number of contour levels in sigma-equivalent units to be drawn.
    """

    params: tuple[str, str]
    _: KW_ONLY
    combiners: list[str] | None = None
    xrange: tuple[float, float] | None = None
    yrange: tuple[float, float] | None = None
    scan: bool = True
    label: str | None = None
    legpos: str = "l"
    legfill: bool = False
    legfontsize: int = -1
    logo: str = "r"
    levels: int = 3


@dataclass(frozen=True)
class Combiner:
    """Representation of a combiner.

    Args:
        id: Numerical ID of the combiner, or list of numerical ids of the PDFs that it should be built from.
        col: Color to be used in plots.
        ls: Linestyle to be used in plots.
    """

    id: int | list[int]
    title: str
    _: KW_ONLY
    title_long: str | None = None  # Longer version of the title for e.g. the subtitle under the logo "CharmFitter"
    scanrange_1d: tuple[float, float] | None = None  # 1D scan range in scan units (radians for angles); used by comp_1d
    scanparams_1d: list[ScanParams1D] = field(default_factory=list)
    scanparams_2d: list[ScanParams2D] = field(default_factory=list)
    external_file: str = ""
    col: str | None = None
    ls: str | None = None

    def __post_init__(self):
        if self.title_long is None:
            object.__setattr__(self, "title_long", self.title)


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
    logo_subtitle: str | None = None
    legpos: str = ""
    legfill: bool = False
    legfontsize: int = -1
    cls: str | None = None
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
            no_cpv_marker = any(
                var1 in str(self.save) and var2 in str(self.save)
                for var1, var2 in [("phiM", "phiG"), ("qop", "phi"), ("Acp_KK", "Acp_PP"), ("Acp_KK", "phiM")]
            )
            u_spin_line = "Acp_KK" in str(self.save) and "Acp_PP" in str(self.save)

            if no_cpv_marker:
                self.legtitles.append(r"No $C\!PV$")
                self.lopts.append({})
                self.fopts.append({})
                self.mopts.append({"marker": "X", "color": "k"})

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
                ax.plot([0], [0], "X", color="k", ms=8)

        logo_xright = 0.65
        if self.logo == "l":
            logo_x, logo_y = 0.03, 0.92
        elif self.logo == "br":
            logo_x, logo_y = logo_xright, 0.06
            if self.logo_subtitle is not None:
                logo_y += 0.04
        else:
            logo_x, logo_y = logo_xright, 0.92
        charm_fitter_logo(pos=(logo_x, logo_y), subtitle=self.logo_subtitle, ax=ax)

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
        if self.cls is not None:
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
    xpar: Parameter,
    ypar: Parameter | None = None,
    *,
    dcs_cpv: bool | None = None,
    mix: MixParam | None = None,
    dy_fsc: DYFsc = DYFsc.NONE,
):
    """Get the prefix of the scanner ROOT file that is read to produce Matplotlib plots.

    Args:
        mix: Mixing parametrisation to be used. If None, the default is used (PHENO if THEO is not available for the
            given parameter(s), otherwise THEO).
    """
    if MixParam.D0_TO_KPI in xpar.mix_params or (ypar is not None and MixParam.D0_TO_KPI in ypar.mix_params):
        return args_prefix

    param = (
        mix
        if mix is not None
        else MixParam.PHENO
        if MixParam.THEO not in xpar.mix_params or (ypar is not None and MixParam.THEO not in ypar.mix_params)
        else MixParam.THEO
    )
    args_prefix = args_prefix.replace("_scanner", f"_{param.value}_{dy_fsc.value}_dyfsc_scanner")
    if dcs_cpv is not None and not dcs_cpv:
        args_prefix = args_prefix.replace("_scanner", "_noDcsCpv_scanner")
    return args_prefix


def _dcs_cpv_label(dcs_cpv: bool | None, args: argparse.Namespace, cfg: ModuleType):
    if dcs_cpv == args.dcs_cpv_default or dcs_cpv is None:
        return ""
    return ", {} {} 0$".format(cfg.parameters["Acp_KP"].latex_name.strip()[:-1], r"\neq" if dcs_cpv else "=")


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


def scans_1d(args: argparse.Namespace, cfg: ModuleType, combiners_ids: list[str] | None = None) -> None:
    """Run the 1D scans for all parameters."""

    if not (args.rescan or args.plugin):
        return
    if combiners_ids is None:
        combiners_ids = cfg.baseline_combiners

    cmds = []
    for par in cfg.parameters.values():
        if not par.scan_1d or (par.name == "Acp_KP" and not args.dcs_cpv):
            continue
        extra_opts = args.extra_opts
        if MixParam.PHENO in par.mix_params and MixParam.THEO not in par.mix_params:
            extra_opts += " --mix pheno"
        if args.plugin:
            if args.submit:
                extra_opts += f" -a pluginbatch --ntoys 50 --nbatchjobs 200 {args.batchopts}"
            else:
                extra_opts += " -a plugin -j 1-200"
        for combiner_id in combiners_ids:
            scanparams = next((x for x in cfg.combiners[combiner_id].scanparams_1d if x.par == par.name), None)
            scan_range = scanparams.range if scanparams is not None else par.scan_range
            cmd = (
                f"bin/{args.execfile} -c {_combiner_string(cfg.combiners[combiner_id].id):<31s}"
                f" --var {par.name:<12s}"
                f" --scanrange {scan_range[0]}:{scan_range[1]} {extra_opts}"
            )
            cmds.append(cmd)

    if args.plugin and args.submit:
        for cmd in cmds:
            run_command(cmd)
        return
    else:
        run_commands(cmds)


def scans_2d(args: argparse.Namespace, cfg: ModuleType, plots_2d: list[Plot2D] | None = None) -> None:
    """Run the 2D scans for all parameters."""

    if not args.rescan:
        return
    if plots_2d is None:
        plots_2d = cfg.plots_2d

    cmds = []
    for plot_params in plots_2d:
        if plot_params.scan is False:
            continue
        xname, yname = plot_params.params
        xpar, ypar = cfg.parameters[xname], cfg.parameters[yname]
        combiners_ids = plot_params.combiners if plot_params.combiners is not None else cfg.baseline_combiners
        for combiner_id in combiners_ids:
            if combiners_ids != cfg.baseline_combiners and combiner_id in cfg.baseline_combiners:
                continue
            scanparams = next((x for x in cfg.combiners[combiner_id].scanparams_2d if x.pars == (xname, yname)), None)
            extra_opts = args.extra_opts
            if any(
                MixParam.PHENO in params and MixParam.THEO not in params
                for params in [xpar.mix_params, ypar.mix_params]
            ) or (scanparams is not None and scanparams.mix_param == MixParam.PHENO):
                extra_opts += " --mix pheno"
            xrange = (
                scanparams.xrange
                if scanparams is not None
                else plot_params.xrange
                if plot_params.xrange is not None
                else xpar.scan_range_2d
            )
            yrange = (
                scanparams.yrange
                if scanparams is not None
                else plot_params.yrange
                if plot_params.yrange is not None
                else ypar.scan_range_2d
            )
            if scanparams is not None:
                extra_opts += f" {scanparams.extra_opts}"
                if scanparams.fixed_pars:
                    pars_string = ",".join(scanparams.fixed_pars)
                    extra_opts += f" --fix-from-parfile {pars_string} --fix-parfile {cfg.baseline_parfile}"
                parfile = scanparams.parfile
                if parfile is None:
                    if hasattr(cfg, "baseline_parfile"):
                        parfile = cfg.baseline_parfile
                    elif hasattr(cfg, "baseline_combiner"):
                        prefix = _get_prefix(
                            args.prefix,
                            xpar,
                            dcs_cpv=args.dcs_cpv if hasattr(args, "dcs_cpv") else None,
                            dy_fsc=cfg.dy_fsc_baseline if hasattr(cfg, "dy_fsc_baseline") else None,
                        ).replace("_scanner", "")
                        parfile = f"plots/par/{prefix}_{cfg.baseline_combiner}_{xpar.name}.dat"
                if parfile is not None and Path(parfile).is_file():
                    extra_opts += f" --parfile {parfile}"
            cmd = (
                f"bin/{args.execfile} -c {_combiner_string(cfg.combiners[combiner_id].id):<31s}"
                f" --var {xpar.name:<7s} --var {ypar.name:<7s}"
                f" --scanrange  {xrange[0]}:{xrange[1]}"
                f" --scanrangey {yrange[0]}:{yrange[1]}"
                f" {extra_opts}"
            )
            cmds.append(cmd)
    run_commands(cmds)


def compare_dy_fsc_hypotheses_scans_2d(
    args: argparse.Namespace,
    cfg: ModuleType,
    combiners_ids: list[str] | None = None,
    dy_fscs: list[DYFsc] = [DYFsc.FULL, DYFsc.NONE],
) -> None:
    """Run the 2D scans for different hypotheses for the final-state dependence of DeltaY(h- h+)."""

    if not args.rescan:
        return
    if combiners_ids is None:
        combiners_ids = cfg.baseline_combiners

    cmds = []
    for plot_params in cfg.dy_fsc_comparison_plots_2d:
        if plot_params.scan is False:
            continue
        xname, yname = plot_params.params
        xpar, ypar = cfg.parameters[xname], cfg.parameters[yname]
        for dy_fsc in dy_fscs:
            if dy_fsc == DYFsc.NONE:
                continue
            extra_opts = args.extra_opts
            if any(MixParam.THEO not in params for params in [xpar.mix_params, ypar.mix_params]):
                extra_opts += " --mix pheno"
            extra_opts += f" --dy-fsc {dy_fsc.value}"
            for combiner_id in combiners_ids:
                xrange = plot_params.xrange if plot_params.xrange is not None else xpar.scan_range_2d
                yrange = plot_params.yrange if plot_params.yrange is not None else ypar.scan_range_2d
                cmd = (
                    f"bin/{args.execfile} -c {cfg.combiners[combiner_id].id}"
                    f" --var {xpar.name:<7s} --var {ypar.name:<7s}"
                    f" --scanrange  {xrange[0]}:{xrange[1]}"
                    f" --scanrangey {yrange[0]}:{yrange[1]}"
                    f" {extra_opts}"
                )
                cmds.append(cmd)
    run_commands(cmds)


def plots_1d(
    args: argparse.Namespace,
    cfg: ModuleType,
    *,
    combiners_ids: list[str] | None = None,
    compare_dcs_hypos: bool = False,
    colors: list[str] | None = None,
    line_styles: list[str] | None = None,
) -> None:
    """Compare the 1D CL plots of all parameters for a set of combiners.

    :param compare_dcs_hypos: Compare the results for different hypothesis of CPV in D0 -> K+ pi- decays
        (two distributions are plotted for each combiner).
    """

    if combiners_ids is None:
        combiners_ids = cfg.baseline_combiners
        combiners_label = "all" if len(combiners_ids) > 3 else "-".join(combiners_ids)
    else:
        combiners_label = "-".join(combiners_ids)

    for parname, par in cfg.parameters.items():
        if par.name == "Acp_KP" and (not args.dcs_cpv or compare_dcs_hypos):
            continue

        plot = Plotter(
            dim=1,
            save=args.savedir / combiners_label / "1d" / _dcs_cpv_subdir(args, compare_dcs_hypos) / f"{parname}.pdf",
            xtitle=par.title,
            xrange=par.plot_range,
            xtransf=par.transf,
            cls="r",
            logo="l",
            legpos=par.plot_opts_1d.legpos,
            legfill=par.plot_opts_1d.legfill,
        )

        dcs_cpv_vals = _dcs_cpv_vals(args, compare_dcs_hypos)
        prefixes = [_get_prefix(args.prefix, par, dcs_cpv=dcs_cpv) for dcs_cpv in dcs_cpv_vals] * len(combiners_ids)
        labels = [
            cfg.combiners[combiner_id].title + _dcs_cpv_label(dcs_cpv, args, cfg)
            for dcs_cpv in dcs_cpv_vals
            for combiner_id in combiners_ids
        ]

        for i, prefix in enumerate(prefixes):
            combiner_name = combiners_ids[i // (2 if compare_dcs_hypos else 1)]
            id = cfg.combiners[combiner_name].id
            combiner_arg = combiner_name if isinstance(id, int) else _combiner_file_string(id)
            plot.add_scan(
                f"{prefix}_{combiner_arg}",
                pars=[par.name],
                label=labels[i],
                col=colors[i] if colors is not None and i < len(colors) else cfg.colors[len(prefixes) - 1 - i],
                ls=line_styles[i]
                if line_styles is not None and i < len(line_styles)
                else cfg.ls[len(prefixes) - 1 - i],
            )
        plot.plot()


def _get_range(par_range, plot_params_range, par_transf):
    if plot_params_range is not None:
        return (par_transf(plot_params_range[0]), par_transf(plot_params_range[1]))
    else:
        return par_range


def plots_2d(
    args: argparse.Namespace,
    cfg: ModuleType,
    plots_2d: list[Plot2D] | None = None,
    *,
    compare_dcs_hypos: bool = False,
    colors: list[str] | None = None,
    line_styles: list[str] | None = None,
    markers: list[str] | None = None,
    breakdown: bool = False,
) -> None:
    """Compare the 2D CL plots of all parameters for a set of combiners.

    :param compare_dcs_hypos: Compare the results for different hypothesis of CPV in D0 -> K+ pi- decays
        (two distributions are plotted for each combiner).
    """
    if plots_2d is None:
        plots_2d = cfg.breakdowns if breakdown else cfg.plots_2d

    for plot_params in plots_2d:
        xname, yname = plot_params.params
        if any(name == "Acp_KP" for name in [xname, yname]) and (not args.dcs_cpv or compare_dcs_hypos):
            continue
        xpar, ypar = cfg.parameters[xname], cfg.parameters[yname]

        if plot_params.combiners is None:
            combiners_ids = cfg.baseline_combiners
            combiners_label = "all" if len(combiners_ids) > 3 else "-".join(combiners_ids)
        else:
            combiners_ids = plot_params.combiners
            combiners_label = "-".join(combiners_ids)
        if plot_params.label is not None:
            combiners_label = plot_params.label

        plot = Plotter(
            dim=2,
            save=args.savedir
            / (cfg.baseline_combiner if breakdown else combiners_label)
            / ("breakdown" if breakdown else "2d")
            / _dcs_cpv_subdir(args, compare_dcs_hypos)
            / f"{yname}-vs-{xname}.pdf",
            xtitle=xpar.title,
            ytitle=ypar.title,
            xtransf=xpar.transf,
            ytransf=ypar.transf,
            xrange=_get_range(xpar.plot_range_2d, plot_params.xrange, xpar.transf),
            yrange=_get_range(ypar.plot_range_2d, plot_params.yrange, ypar.transf),
            levels=2 if breakdown else plot_params.levels,
            logo=plot_params.logo,
            legpos=plot_params.legpos,
            legfill=plot_params.legfill,
            legfontsize=plot_params.legfontsize,
        )

        dcs_cpv_vals = _dcs_cpv_vals(args, compare_dcs_hypos)
        prefixes = []
        for combiner_id in combiners_ids:
            combiner = cfg.combiners[combiner_id]
            scanparams = next((x for x in combiner.scanparams_2d if x.pars == (xname, yname)), None)
            mix = scanparams.mix_param if scanparams is not None else None
            for dcs_cpv in dcs_cpv_vals:
                prefixes.append(
                    _get_prefix(
                        args.prefix,
                        xpar,
                        ypar,
                        dcs_cpv=dcs_cpv,
                        mix=mix,
                        dy_fsc=cfg.dy_fsc_baseline if hasattr(cfg, "dy_fsc_baseline") else None,
                    )
                )

        labels = [
            cfg.combiners[combiner_id].title + _dcs_cpv_label(dcs_cpv, args, cfg)
            for dcs_cpv in dcs_cpv_vals
            for combiner_id in combiners_ids
        ]

        for i, prefix in enumerate(prefixes):
            combiner_name = combiners_ids[i // (2 if compare_dcs_hypos else 1)]
            combiner = cfg.combiners[combiner_name]
            combiner_arg = combiner_name if isinstance(combiner.id, int) else _combiner_file_string(combiner.id)
            col = (
                colors[i]
                if colors is not None and i < len(colors)
                else combiner.col
                if combiner.col is not None and not compare_dcs_hypos
                else cfg.colors_2d[len(prefixes) - 1 - i]
            )
            ls = (
                line_styles[i]
                if line_styles is not None and i < len(line_styles)
                else combiner.ls
                if combiner.ls is not None and not compare_dcs_hypos
                else cfg.ls_2d[len(prefixes) - 1 - i]
            )
            marker = (
                "x"
                if breakdown
                else markers[i]
                if markers is not None and i < len(markers)
                else cfg.markers[len(prefixes) - 1 - i]
            )
            plot.add_scan(
                f"{prefix}_{combiner_arg}",
                pars=[xpar.name, ypar.name],
                label=labels[i],
                col=col,
                ls=ls,
                marker=marker,
                bf=(combiner_name == cfg.baseline_combiner) if breakdown else True,
            )
        plot.plot()


def compare_dy_fsc_hypotheses_plots_2d(
    args: argparse.Namespace,
    cfg: ModuleType,
    *,
    combiners_ids: list[str] | None = None,
    dy_fscs: list[DYFsc] = [DYFsc.FULL, DYFsc.NONE],
    colors: list[str] | None = None,
    line_styles: list[str] | None = None,
    markers: list[str] | None = None,
) -> None:
    """Compare the results for different hypotheses for the final-state dependence of DeltaY(h- h+)."""

    if combiners_ids is None:
        combiners_ids = cfg.baseline_combiners

    for plot_params in cfg.dy_fsc_comparison_plots_2d:
        xname, yname = plot_params.params
        if any(name == "Acp_KP" for name in [xname, yname]) and not args.dcs_cpv:
            continue
        xpar, ypar = cfg.parameters[xname], cfg.parameters[yname]

        for combiner_id in combiners_ids:
            plot = Plotter(
                dim=2,
                save=args.savedir
                / combiner_id
                / "dy-fsc-comparisons"
                / ("with-dcs-cpv" if args.dcs_cpv else "no-dcs-cpv")
                / f"{yname}-vs-{xname}.pdf",
                xtitle=xpar.title,
                ytitle=ypar.title,
                xtransf=xpar.transf,
                ytransf=ypar.transf,
                xrange=_get_range(xpar.plot_range_2d, plot_params.xrange, xpar.transf),
                yrange=_get_range(ypar.plot_range_2d, plot_params.yrange, ypar.transf),
                levels=3,
                logo="br",
                logo_subtitle=cfg.combiners[combiner_id].title_long,
                legpos=plot_params.legpos,
                legfill=plot_params.legfill,
                legfontsize=12,
            )

            for i, dy_fsc in enumerate(dy_fscs):
                prefix = _get_prefix(args.prefix, xpar, ypar, dcs_cpv=args.dcs_cpv, dy_fsc=dy_fsc)
                label = cfg.dy_fsc_labels[dy_fsc] + _dcs_cpv_label(args.dcs_cpv, args, cfg)
                col = colors[i] if colors is not None and i < len(colors) else cfg.colors_2d[len(dy_fscs) - 1 - i]
                ls = (
                    line_styles[i]
                    if line_styles is not None and i < len(line_styles)
                    else cfg.ls_2d[len(dy_fscs) - 1 - i]
                )
                marker = markers[i] if markers is not None and i < len(markers) else cfg.markers[len(dy_fscs) - 1 - i]
                plot.add_scan(
                    f"{prefix}_{combiner_id}",
                    pars=[xpar.name, ypar.name],
                    label=label,
                    col=col,
                    ls=ls,
                    marker=marker,
                    bf=True,
                )
            plot.plot()


def parse_args(*, combo: str = "charm", dcs_cpv_default: bool = False) -> argparse.ArgumentParser:
    actions = {
        "all": "Run all actions",
        "1d": "Make the 1D plots",
        "2d": "Make the 2D plots ({} etc)".format(
            "x vs y, phi vs |q/p|" if combo == "charm" else "x'2 vs y', dx'2 vs dy'"
        ),
    }
    if combo == "charm":
        actions["breakdown"] = "Make the 2D breakdown plots, divided by D0 final-state (x vs y, phi vs |q/p| etc)"
    default_execname = f"{combo}-combo"

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-a",
        "--actions",
        default=[],
        action="append",
        choices=list(actions.keys()),
        help="\n".join(
            ["Run specific action(s). Can be parsed multiple times. Available actions:"]
            + [f"{k:<17s} {v}" for k, v in actions.items()]
        ),
    )
    parser.add_argument(
        "-r", "--rescan", default=False, action="store_true", help="Rerun the scans (rather than just make plots)"
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=repo_path / f"config/plotting/{combo}-2025.py",
        help="Path of the configuration file",
    )
    parser.add_argument("-e", "--execfile", type=str, default=default_execname, help="Executable name")
    parser.add_argument(
        "--extra-opts",
        type=str,
        default="--ps 1",
        help="Additional line options to be passed to the executable",
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
    if combo != "ws" and not args.dcs_cpv:
        args.extra_opts += " --no-dcs-cpv"
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
    subtitle: str | None = None,
    ax: matplotlib.axes.Axes | None = None,
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
        ha="left",
        bbox=props,
        fontdict=font,
        usetex=usetex,
    )
    if subtitle is not None:
        ax.text(
            pos[0],
            pos[1] - 0.05,
            subtitle,
            transform=ax.transAxes,
            size=min(12.2, 12.2 * 20 / len(subtitle)),
            ha="left",
            bbox=props,
            fontdict=font,
            usetex=usetex,
        )


def get_configuration(config_fpath: Path) -> ModuleType:
    """Load a configuration file as a Python module."""
    if not config_fpath.is_file():
        raise RuntimeError(f'Could not open configuration file "{config_fpath}"')
    spec = importlib.util.spec_from_file_location("module", str(config_fpath))
    cfg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cfg)
    return cfg


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


def _bring_angle_to_pi_pi(x: float) -> float:
    """Shift the angle x (in radians) to the range [-pi, pi]."""
    return (x + np.pi) % (2 * np.pi) - np.pi
