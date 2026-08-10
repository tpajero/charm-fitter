#!/usr/bin/env python3

"""Main script to run the likelihood scans and produce all the plots in Matplotlib."""

import matplotlib.pyplot as plt
from gc_core.utils import load_gc_core_lib

from charm_fitter.utils import (
    cwd,
    get_configuration,
    parse_args,
    plots_1d,
    plots_2d,
    repo_path,
    scans_1d,
    scans_2d,
    setup_matplotlib,
)

if __name__ == "__main__":
    args = parse_args()

    # Setup environment and Matplotlib, import plotting configuration
    load_gc_core_lib()
    setup_matplotlib(usetex=args.latex)
    cfg = get_configuration(args.config)

    # Create the output directory and run the requested actions
    args.savedir.mkdir(parents=True, exist_ok=True)

    with cwd(repo_path):  # Need to call the executable from the repository root due to GammaCombo core limitations
        if any(x in ["all", "1d"] for x in args.actions):
            scans_1d(list(cfg.baseline_combiners.values()), args)
            plots_1d(cfg.plots_1d, args, cfg, compare_dcs_hypos=args.dcs_cpv != args.dcs_cpv_default)
        if any(x in ["all", "2d"] for x in args.actions):
            scans_2d(list(cfg.baseline_combiners.values()), args)
            plots_2d(cfg.plots_2d, args, cfg, compare_dcs_hypos=args.dcs_cpv != args.dcs_cpv_default)
        if any(x in ["all", "dy-fsc"] for x in args.actions):
            scans_1d(list(cfg.baseline_combiners.values()), args, dy_fsc_comparison=True)
            scans_2d(list(cfg.baseline_combiners.values()), args, dy_fsc_comparison=True)
            plots_1d(cfg.plots_dy_fsc_1d, args, cfg)
            plots_2d(cfg.plots_dy_fsc_2d, args, cfg, plots_cat="compare-dy-fsc")
        if any(x in ["all", "breakdown"] for x in args.actions):
            scans_2d(list(cfg.combiners_breakdown.values()), args)
            plots_2d(cfg.plots_breakdown, args, cfg, plots_cat="breakdown")

    if args.interactive:
        plt.show()
