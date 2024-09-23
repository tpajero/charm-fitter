#!/usr/bin/env python3
"""Produce the plot with the world average of the DeltaY parameter."""

from argparse import ArgumentParser
import logging
import os
import matplotlib.pyplot as plt
from utils import Measurement


def get_measures(comb):
    """Get the list of measurements to be shown for a given summary plot.

    :param comb: Identifier for a given summary plot.
    :type comb: str
    """
    measures = []
    if "lhcb" not in comb:
        measures.append(Measurement("BaBar 2012", "1209.3896", -8.8, 25.5, 5.8))
        measures.append(Measurement("CDF 2014", "1410.5435", 12.0, 12.0))
    measures.append(Measurement(r"LHCb 2015 $\mu^{-}$ tag (3 fb$^{-1}$)", "1501.06777", 12.5, 7.3))
    if "lhcb" not in comb:
        measures.append(Measurement("Belle 2016", "1509.08266", 3.0, 20.0, 7.0))
    measures.append(Measurement("LHCb 2017 $D^{*+}$ tag (3 fb$^{-1}$)", "1702.06490", 1.3, 2.8, 1.0))
    measures.append(Measurement(r"LHCb 2020 $\mu^{-}$ tag (5.4 fb$^{-1}$)", "1911.01114", 2.9, 3.2, 0.5))
    if comb == "wa-summer-2020":
        measures.insert(Measurement("LHCb 2019 $D^{*+}$ tag (2 fb$^{-1}$)", "", -3.4, 3.1, 0.6), len(measures) - 2)
        measures.append(Measurement("World average", None, 1.9, 1.6, 0.5))
    else:
        measures.append(Measurement("LHCb 2021 $D^{*+}$ tag (6 fb$^{-1}$)", "2105.09889", -2.70, 1.34, 0.30))
        if comb == "wa-2021":
            measures.append(Measurement("World average", None, -0.9, 1.1, 0.3))
        elif comb == "lhcb-2021":
            measures.append(Measurement("LHCb average", None, -1.0, 1.1, 0.3))
        elif "2024" in comb:
            measures.append(
                Measurement(r"LHCb 2024 $D^0\to\pi^{+}\pi^{-}\pi^{0}$ (6 fb$^{-1}$)", "2405.06556", -1.3, 6.3, 2.4))
            if comb == "wa-2024":
                measures.append(Measurement("World average", None, -0.96, 1.10, 0.32))
            elif comb == "lhcb-2024":
                measures.append(Measurement("LHCb average", None, -1.09, 1.11, 0.32))
    return measures


def get_xy_ranges(comb, dy_notation):
    """Get the x and y ranges for a given summary plot."""
    n_meas = len(get_measures(comb))
    if 'lhcb-only' in comb:
        x_min = -10. if dy_notation else -27.
        x_max = 70.
    else:
        x_min = -40. if dy_notation else -27.
        if '2024' in comb:
            x_max = 120
        else:
            x_max = 90. if dy_notation else 70.
    y_min, y_max = -0.5, n_meas - 0.5
    return (x_min, x_max), (y_min, y_max)


def make_plot(comb, dy_notation, out_dir):

    if comb == 'lhcb-only':
        fig, ax = plt.subplots(figsize=(6, 6))
    else:
        fig, ax = plt.subplots(figsize=(8 if dy_notation else 6, 8))
    plt.xlabel(r'$\Delta Y\;[10^{-4}]$' if dy_notation else r'$A_{\Gamma}\;[10^{-4}]$', fontsize=24, ha='center')
    (x_min, x_max), (y_min, y_max) = get_xy_ranges(comb, dy_notation)
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    plt.tick_params(axis='y', which='both', right=False, left=False, labelleft=False)

    measures = get_measures(comb)
    n_meas = len(measures)
    for i in range(n_meas):
        meas = measures[i]
        if not dy_notation:
            meas.m = -meas.m
        plt.errorbar(meas.m,
                     n_meas - 1 - i,
                     xerr=meas.err(),
                     fmt='.',
                     markersize=8,
                     capsize=7,
                     color=meas.color if comb != 'lhcb-only' else 'k')
        plt.errorbar(meas.m,
                     n_meas - 1 - i,
                     xerr=meas.stat,
                     capsize=7,
                     color=meas.color if comb != 'lhcb-only' else 'k')

        if dy_notation:
            x_text = 26 if 'lhcb-only' in comb else 30
            y_text = n_meas - i - 1.25
        else:
            x_text = 22 if i == 3 else 9
            y_text = n_meas - 0.8 if i == 0 else n_meas - i - 1.25
        plt.text(x_text, y_text, f'{meas.label}\n{meas.result_str()}', color=meas.color if comb != 'lhcb-only' else 'k')

    # Vertical line for world average
    plt.plot([measures[-1].m, measures[-1].m], [y_min, y_max], linestyle='--', linewidth=1, color='k')

    # Horizontal line to divide world average from other measurements
    plt.plot([x_min, x_max], [0.5, 0.5], linestyle='-', linewidth=1, color='k')

    # Save figure
    os.makedirs(out_dir, exist_ok=True)
    obs = 'dy' if dy_notation else 'agamma'
    fig_name = f'{obs}-{comb}'
    for ext in ['pdf']:
        plt.savefig(f'{out_dir}/{fig_name}.{ext}')


def parse_args():
    cwd = os.path.dirname(os.path.realpath(__file__))
    parser = ArgumentParser()
    parser.add_argument('-c',
                        '--comb',
                        type=str,
                        default='wa-2024',
                        help='Combination id',
                        choices=['wa-summer-2020', 'wa-2021', 'lhcb-2021', 'wa-2024', 'lhcb-2024'])
    parser.add_argument('--agamma',
                        action='store_false',
                        dest='dy_notation',
                        help='Use the A_Gamma notation rather than DeltaY')
    parser.add_argument('-o',
                        '--outdir',
                        type=str,
                        default=os.path.join(cwd, 'figs'),
                        help='Output directory for saving the plots')
    return parser.parse_args()


if __name__ == '__main__':

    logging.basicConfig(encoding='utf-8', level=logging.INFO)
    cwd = os.path.dirname(os.path.realpath(__file__))
    plt.style.use(os.path.join(cwd, 'lhcb.mplstyle'))

    args = parse_args()

    make_plot(args.comb, args.dy_notation, args.outdir)
