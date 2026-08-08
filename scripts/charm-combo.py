#!/usr/bin/env python3

"""Main script to run the likelihood scans and produce all the plots in Matplotlib."""

import matplotlib.pyplot as plt
from gc_core.utils import load_gc_core_lib

from charm_fitter.utils import (
    compare_dy_fsc_hypotheses_plots_2d,
    compare_dy_fsc_hypotheses_scans_2d,
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
            scans_1d(args, cfg)
            plots_1d(args, cfg, compare_dcs_hypos=args.dcs_cpv != args.dcs_cpv_default)
        if any(x in ["all", "2d"] for x in args.actions):
            scans_2d(args, cfg)
            plots_2d(args, cfg, compare_dcs_hypos=args.dcs_cpv != args.dcs_cpv_default)
            if args.dcs_cpv == args.dcs_cpv_default:
                compare_dy_fsc_hypotheses_scans_2d(args, cfg)
                compare_dy_fsc_hypotheses_plots_2d(args, cfg)
        if any(x in ["all", "breakdown"] for x in args.actions):
            scans_2d(args, cfg, cfg.breakdowns)
            plots_2d(args, cfg, breakdown=True)

    if args.interactive:
        plt.show()
