#!/usr/bin/env python3

"""Main script to run the likelihood scans and produce plots for the WS/RS D0 -> K+ pi- measurements."""

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
    args = parse_args(combo="ws")

    # Setup environment and Matplotlib, import plotting configuration
    load_gc_core_lib()
    setup_matplotlib(usetex=args.latex)
    cfg = get_configuration(args.config)

    # Create the output directory and run the requested actions
    args.savedir.mkdir(parents=True, exist_ok=True)

    with cwd(repo_path):  # Need to call the executable from the repository root due to GammaCombo core limitations
        if any(x in ["all", "1d"] for x in args.actions):
            scans_1d(args, cfg)
            plots_1d(args, cfg)
        if any(x in ["all", "2d"] for x in args.actions):
            scans_2d(args, cfg)
            plots_2d(args, cfg)

    if args.interactive:
        plt.show()
