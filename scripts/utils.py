"""Make and read scans, and plot them in matplotlib."""

import os
import subprocess
import sys
from multiprocessing import Pool

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator
from scipy.stats import chi2

charm_fitter_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), ".."))
sys.path.append(os.path.abspath(os.path.join(charm_fitter_dir, "..")))
from scripts.gc_mpl_tools import (
    lhcb_2d_cols,
    lhcb_cols,
    lhcb_ls,
    plot1d,
    plot2d,
    print_interval,
    read_gc_scan,
)


def run_command(cmd, verbose=False):
    print(f"Running: {cmd}")
    return subprocess.run(cmd, shell=True, capture_output=not verbose, check=True)


def run_commands(cmds):
    pool = Pool()
    [pool.map(run_command, cmds)]


def getfnames(prefix, xpar, ypar=None):
    yext = f"_{ypar}" if ypar is not None else ""
    fname = f"plots/scanner/{prefix}_{xpar}{yext}.root"

    bfname = fname.replace("scanner/", "par/").replace("_scanner", "").replace(".root", ".dat")

    if not os.path.exists(fname):
        raise FileNotFoundError(f"Cannot find scan file {fname}")

    if not os.path.exists(bfname):
        raise FileNotFoundError(f"Cannot find fit result file {bfname}")

    return fname, bfname


def print_cl(prefix, xpar, ypar=None, prob=True):
    if ypar is not None:
        return

    pref = prefix.split("scanner")[1]
    suff = "Prob" if prob else "Plugin"
    fname = f"plots/cl/clintervals{pref}_{xpar}_{suff}.py"
    if not os.path.exists(fname):
        return

    try:
        print_interval(fname, nsigma=1)
    except Exception:
        return


def get_scan_res(prefix, xpar, ypar=None):
    pars = [xpar]
    if ypar is not None:
        pars.append(ypar)

    fname, bfname = getfnames(prefix, xpar, ypar)

    header_str = prefix + " - " + xpar
    if ypar is not None:
        header_str += " , " + ypar

    print(header_str)
    print_cl(prefix, xpar, ypar)

    return read_gc_scan(fname, bfname, pars)


def charm_fitter_logo(pos=[0.02, 0.88], date=None, ax=None):
    """Add the CharmFitter logo to the current Matplotlib plot (or ax, if given)."""
    ax = ax or plt.gca()
    props = dict(fc="none", ec="none", boxstyle="square,pad=0.1")
    font = {"family": "Times New Roman", "weight": 400}
    ax.text(
        *pos,
        r"$\textsc{CharmFitter}$",
        transform=ax.transAxes,
        size=16,
        ha="left",
        bbox=props,
        fontdict=font,
        usetex=True,
    )
    if date is not None:
        ax.text(
            pos[0],
            pos[1] - 0.05,
            date,
            transform=ax.transAxes,
            size=12.2,
            ha="left",
            bbox=props,
            fontdict=font,
            usetex=False,
        )


class plotter:
    """Main class to make 1D and 2D Matplotlib plots."""

    def __init__(
        self,
        dim=1,
        save=None,
        xtitle=None,
        ytitle=None,
        xangle=False,
        yangle=False,
        xrange=None,
        yrange=None,
        xtransf=None,
        ytransf=None,
        levels=2,
        logo="l",
        legpos=None,
        legfill=False,
        legfontsize=-1,
        cls=None,
        axes_origin=(-1, -1),
        prune_xlabel=False,
        prune_ylabel=False,
    ):
        self.dim = dim
        self.save = save
        self.xtitle = xtitle
        self.ytitle = ytitle
        self.xangle = xangle
        self.yangle = yangle
        self.xrange = xrange
        self.yrange = yrange
        self.xtransf = xtransf
        self.ytransf = ytransf
        self.levels = levels
        self.logo = logo
        self.legpos = legpos
        self.legfill = legfill
        self.legfontsize = legfontsize
        self.cls = cls
        self.axes_origin = axes_origin
        self.prune_xlabel = prune_xlabel
        self.prune_ylabel = prune_ylabel

        self.scanpoints = []
        self.bfs = []
        self.legtitles = []
        self.lopts = []
        self.fopts = []
        self.mopts = []

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

        if self.logo == "l" and self.legpos == "l":
            legopts = dict(
                bbox_to_anchor=(
                    0.01,
                    0.83 - 0.07 * len(self.scanpoints),
                    0.2,
                    0.07 * len(self.scanpoints),
                ),
                loc="upper left",
            )
        elif self.logo == "r" and self.legpos == "r":
            legopts = dict(
                bbox_to_anchor=(
                    0.8,
                    0.83 - 0.07 * len(self.scanpoints),
                    0.2,
                    0.07 * len(self.scanpoints),
                ),
                loc="upper right",
                fontsize=12,
            )
        elif self.logo == "l" and self.legpos == "r":
            legopts = dict(
                bbox_to_anchor=(
                    0.8,
                    0.997 - 0.07 * len(self.scanpoints),
                    0.2,
                    0.07 * len(self.scanpoints),
                ),
                loc="upper right",
                fontsize=12,
            )
        elif self.logo == "r" and self.legpos == "l":
            legopts = dict(
                bbox_to_anchor=(
                    0.01,
                    0.997 - 0.07 * len(self.scanpoints),
                    0.2,
                    0.07 * len(self.scanpoints),
                ),
                loc="upper left",
                fontsize=12,
            )
        elif self.legpos == "bl":
            legopts = dict(
                bbox_to_anchor=(0.01, 0.02, 0.2, 0.07 * len(self.scanpoints)),
                loc="lower left",
                fontsize=12,
            )
        elif self.legpos == "br":
            legopts = dict(
                bbox_to_anchor=(0.8, 0.01, 0.2, 0.07 * len(self.scanpoints)),
                loc="lower right",
                fontsize=12,
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
                angle=self.xangle,
                ax=ax,
                legopts=legopts,
                **kwargs,
            )
        elif self.dim == 2:
            no_cpv_marker = any(var in self.save for var in ["phiM-phiG", "qop-phi", "Acp_KK-Acp_PP", "Acp_KK-phiM"])
            u_spin_line = "Acp_KK-Acp_PP" in self.save

            if no_cpv_marker:
                self.legtitles.append(r"No $C\!PV$")
                self.lopts.append({})
                self.fopts.append({})
                self.mopts.append({"marker": "X", "color": "k"})

            if u_spin_line:
                self.legtitles.extend(["$U$-spin symmetry", "Improved $U$-spin prediction"])
                self.lopts.extend(
                    [
                        {"color": "g", "lw": 2, "ls": "--"},
                        {"color": "b", "lw": 2, "ls": ":"},
                    ]
                )
                self.fopts.extend([{}, {}])
                self.mopts.extend([{}, {}])

            plot2d(
                self.scanpoints,
                self.lopts,
                self.fopts,
                self.mopts,
                title=[self.xtitle, self.ytitle],
                levels=self.levels,
                legtitles=self.legtitles,
                angle=[self.xangle, self.yangle],
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
                akk_over_app = 1.8091  # Amplitude ratio
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

        if self.logo == "l":
            logo_x, logo_y = 0.03, 0.92
        elif self.logo == "br":
            logo_x, logo_y = 0.69, 0.08
        else:
            logo_x, logo_y = 0.69, 0.92
        charm_fitter_logo(pos=[logo_x, logo_y], date=None, ax=ax)

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

        fig.savefig(self.save)

        if not interactive:
            fig.clf()

        plt.close()
