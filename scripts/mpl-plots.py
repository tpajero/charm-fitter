"""Make and read scans, and plot them in matplotlib."""

import importlib
import os
import sys
from argparse import ArgumentParser
from contextlib import contextmanager

import matplotlib.pyplot as plt
from utils import plotter, run_commands

charm_fitter_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), ".."))
sys.path.append(os.path.abspath(os.path.join(charm_fitter_dir, "..")))


@contextmanager
def cwd(path):
    oldpwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(oldpwd)


def parse_args():
    actions = ["all", "scans-1d", "scans-1d-no-dcs-cpv", "scans-2d", "scans-2d-no-dcs-cpv"]
    parser = ArgumentParser()
    parser.add_argument("-e", "--execfile", type=str, default="charmcombo", help="Executable name")
    parser.add_argument(
        "-c",
        "--combiner-id",
        type=int,
        default=55,
        help="Numerical ID of the combiner to be run",
    )
    parser.add_argument("--config", type=str, default="config/2025.py", help="Configuration file")
    parser.add_argument("-s", "--savedir", type=str, default="plots/matplotlib", help="Output directory")
    parser.add_argument(
        "-i",
        "--interactive",
        default=False,
        action="store_true",
        help="Show plots interactively at the end",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        default=False,
        action="store_true",
        help="Verbose output of scan commands",
    )
    parser.add_argument(
        "-r",
        "--rescan",
        default=False,
        action="store_true",
        help="Rerun the scans (rather than just make plots)",
    )
    parser.add_argument("-P", "--plugin", default=False, action="store_true", help="Use plugin scans")
    parser.add_argument(
        "-S",
        "--submit",
        default=False,
        action="store_true",
        help="Submit plugin batch jobs",
    )
    parser.add_argument(
        "-B",
        "--batchopts",
        metavar="",
        default="",
        help="Additional batch submission options",
    )
    parser.add_argument(
        "-a",
        "--actions",
        metavar="",
        default=[],
        action="append",
        choices=actions,
        help="Run specific action(s). Can be parsed multiple times.",
    )
    return parser.parse_args()


def run_scans_1d(args, cfg, no_dcs_cpv=False):
    combiner_id = args.combiner_id + 1000 if no_dcs_cpv else args.combiner_id

    cmds = []
    for par in cfg.params.values():
        if no_dcs_cpv and par.id == "Acp_KP":
            continue
        cmd = f"bin/{args.execfile} -c {combiner_id} --var {par.id} --scanrange {par.range1d[0]}:{par.range1d[1]} --pr --ps 1"
        cmd += " --param phenomenological" if par.parametrisation == "pheno" else " --param theoretical"
        if no_dcs_cpv:
            cmd += " --fix Acp_KP=0"
        if not args.rescan:
            cmd += " -a plot"
        cmds.append(cmd)
    run_commands(cmds)


def run_scans_2d(args, cfg, no_dcs_cpv=False):
    combiner_id = args.combiner_id + 1000 if no_dcs_cpv else args.combiner_id

    cmds = []
    for plot_2d in cfg.plots_2d:
        xpar, ypar = cfg.params[plot_2d.pars[0]], cfg.params[plot_2d.pars[1]]
        if no_dcs_cpv and "Acp_KP" in [xpar.id, ypar.id]:
            continue
        cmd = (
            f"bin/{args.execfile} -c {combiner_id} --var {xpar.id} --var {ypar.id} --pr --ps 1"
            f" --scanrange  {xpar.range2d[0]}:{xpar.range2d[1]}"
            f" --scanrangey {ypar.range2d[0]}:{ypar.range2d[1]}"
        )
        cmd += (
            " --param phenomenological"
            if "pheno" in [xpar.parametrisation, ypar.parametrisation]
            else " --param theoretical"
        )
        if no_dcs_cpv:
            cmd += " --fix Acp_KP=0"
        if not args.rescan:
            cmd += " -a plot"
        cmds.append(cmd)
    run_commands(cmds)


def dirname_from_combiners(combiners):
    """Set the name of the output directory, starting from a list of combiners."""

    if len(combiners) == 1:
        return f"{combiners[0] % 1000}/{'baseline' if combiners[0] < 1000 else 'no-dcs-cpv'}"

    if len(combiners) == 2 and combiners[0] % 1000 == combiners[1] % 1000:
        assert combiners[0] != combiners[1], "Comparing the same combination twice does not make sense"
        return f"{combiners[0] % 1000}/baseline-vs-no-dcs-cpv"

    return "-".join(combiners)


def make_plots_1d(combiners, cfg, savedir):
    """Plot the 1D comparison of a set of combiners, for all parameters in the combination."""

    out_dir = os.path.join(savedir, dirname_from_combiners(combiners), "1d")
    os.makedirs(out_dir, exist_ok=True)

    for par in cfg.params.values():
        if any(combiner > 1000 for combiner in combiners) and par.id == "Acp_KP":
            continue
        param = "phenomenological" if par.parametrisation == "pheno" else "theoretical"
        plot = plotter(
            dim=1,
            save=f"{out_dir}/{par.id}.pdf",
            xtitle=par.title,
            xangle=par.degrees,
            cls="r",
            logo="l",
        )
        for i, combiner in enumerate(combiners):
            scan_name = f"{param}_scanner_{cfg.combiners_info[combiner][0]}"
            plot.add_scan(
                scan_name,
                pars=[par.id],
                lw=1,
                col=cfg.colors[i],
                ls=cfg.ls[i],
                label=cfg.combiners_info[combiner][1],
            )
        plot.plot()


def make_plots_2d(combiners, cfg, savedir):
    """Plot the 1D comparison of a set of combiners, for all pairs parameters in the configuration file."""

    out_dir = os.path.join(savedir, dirname_from_combiners(combiners), "2d")
    os.makedirs(out_dir, exist_ok=True)

    for plot_2d in cfg.plots_2d:
        xpar, ypar = cfg.params[plot_2d.pars[0]], cfg.params[plot_2d.pars[1]]
        if any(combiner > 1000 for combiner in combiners) and "Acp_KP" in [xpar.id, ypar.id]:
            continue
        param = "phenomenological" if "pheno" in [xpar.parametrisation, ypar.parametrisation] else "theoretical"
        plot = plotter(
            dim=2,
            save=f"{out_dir}/{xpar.id}-{ypar.id}.pdf",
            xtitle=xpar.title,
            ytitle=ypar.title,
            xrange=xpar.range2d,
            yrange=ypar.range2d,
            levels=3,
            logo=plot_2d.logo,
            legpos=plot_2d.legpos,
        )
        for i, combiner in enumerate(combiners):
            scan_name = f"{param}_scanner_{cfg.combiners_info[combiner][0]}"
            plot.add_scan(
                scan_name,
                pars=[xpar.id, ypar.id],
                label=cfg.combiners_info[combiner][1],
                bf=True,
                col=cfg.colors_2d[i],
                ls=cfg.ls_2d[i],
                marker=cfg.markers[i],
            )
        plot.plot()


if __name__ == "__main__":
    plt.style.use(os.path.join(charm_fitter_dir, "scripts/lhcb.mplstyle"))

    # Parse the arguments
    args = parse_args()
    if not os.path.isabs(args.config):
        args.config = os.path.join(charm_fitter_dir, args.config)
    if not os.path.isabs(args.savedir):
        args.savedir = os.path.join(charm_fitter_dir, args.savedir)
    os.makedirs(args.savedir, exist_ok=True)
    os.makedirs(os.path.join(charm_fitter_dir, "plots/corr"), exist_ok=True)

    # Import plots configuration
    spec = importlib.util.spec_from_file_location("module", args.config)
    cfg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cfg)

    with cwd(charm_fitter_dir):
        # Run 1D scans and make plots
        if any(x in ["all", "scans-1d"] for x in args.actions):
            run_scans_1d(args, cfg)
            make_plots_1d([args.combiner_id], cfg, args.savedir)
        if any(x in ["all", "scans-1d-no-dcs-cpv"] for x in args.actions):
            run_scans_1d(args, cfg, args.savedir)
            make_plots_1d([args.combiner_id + 1000], cfg, args.savedir)
            make_plots_1d([args.combiner_id, args.combiner_id + 1000], cfg, args.savedir)

        # Run 2D scans and make plots
        if any(x in ["all", "scans-2d"] for x in args.actions):
            run_scans_2d(args, cfg)
            make_plots_2d([args.combiner_id], cfg, args.savedir)
        if any(x in ["all", "scans-2d-no-dcs-cpv"] for x in args.actions):
            run_scans_2d(args, cfg, no_dcs_cpv=True)
            make_plots_2d([args.combiner_id + 1000], cfg, args.savedir)
            make_plots_2d([args.combiner_id, args.combiner_id + 1000], cfg, args.savedir)

    if args.interactive:
        plt.show()
