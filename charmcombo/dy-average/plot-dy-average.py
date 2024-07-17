#!/usr/bin/env python3
"""Produce the plot with the world average of the DeltaY parameter."""

from argparse import ArgumentParser
import os
from math import sqrt
import matplotlib.pyplot as plt

plt.style.use('lhcb.mplstyle')


class Measurement:
    """Class to define the objects storing the information relative to a single measurement of DY/AGamma."""

    def __init__(self, m, stat, sys, label, colour):
        self.m = m
        self.stat = stat
        self.sys = sys
        self.label = label
        self.colour = colour

    def value(self):
        sys_label = '' if (self.sys is None) else f' $\\pm$ {self.sys:1.1f}'
        if self.stat / 10 > 1:
            value = f'{self.m:+2.1f} $\\pm$ {self.stat:2.1f}{sys_label}'
        else:
            value = f'{self.m:+2.1f} $\\pm$ {self.stat:1.1f}{sys_label}'
        value = value.replace('-', '\u2013')
        return value

    def err(self):
        if self.sys is None:
            return self.stat
        else:
            return sqrt((self.stat)**2 + (self.sys)**2)


def get_measures(comb):
    """Get the list of measurements to be shown for a given summary plot.

    :param comb: Identifier for a given summary plot.
    :type comb: str
    """
    measures = []
    if 'lhcb' not in comb:
        measures.append(Measurement(-8.8, 25.5, 5.8, 'BaBar 2012', 'g'))
        measures.append(Measurement(12.0, 12.0, None, 'CDF 2014', 'r'))
    measures.append(Measurement(12.5, 7.3, None, r'LHCb 2015 $\mu^{-}$ tag (3 fb$^{-1}$)', 'b'))
    if 'lhcb' not in comb:
        measures.append(Measurement(3.0, 20.0, 7.0, 'Belle 2016', 'm'))
    measures.append(Measurement(1.3, 2.8, 1.0, 'LHCb 2017 $D^{*+}$ tag (3 fb$^{-1}$)', 'b'))
    measures.append(Measurement(2.9, 3.2, 0.5, r'LHCb 2020 $\mu^{-}$ tag (5.4 fb$^{-1}$)', 'b'))
    if comb == 'wa-summer-2020':
        measures.insert(Measurement(-3.4, 3.1, 0.6, 'LHCb 2019 $D^{*+}$ tag (2 fb$^{-1}$)', 'b'), len(measures) - 2)
        measures.append(Measurement(1.9, 1.6, 0.5, 'World average', 'k'))
    else:
        measures.append(Measurement(-2.7, 1.3, 0.3, 'LHCb 2021 $D^{*+}$ tag (6 fb$^{-1}$)', 'b'))
        if comb == 'wa-2021':
            measures.append(Measurement(-0.9, 1.1, 0.3, 'World average', 'k'))
        elif comb == 'lhcb-2021':
            measures.append(Measurement(-1.0, 1.1, 0.3, 'LHCb average', 'k'))
        elif '2024' in comb:
            measures.append(Measurement(-1.3, 6.3, 2.4, r'LHCb 2024 $D^0\to\pi^{+}\pi^{-}\pi^{0}$ (6 fb$^{-1}$)', 'b'))
            if comb == 'wa-2024':
                measures.append(Measurement(-1.0, 1.1, 0.3, 'World average', 'k'))
            elif comb == 'lhcb-2024':
                measures.append(Measurement(-1.1, 1.1, 0.3, 'LHCb average', 'k'))
    return measures


def get_xy_ranges(comb, dy_notation):
    """Get the x and y ranges for a given summary plot."""
    n_meas = len(get_measures(comb))
    if 'lhcb-only' in comb:
        x_min = -10. if dy_notation else -27.
        x_max = 70.
    else:
        x_min = -40. if dy_notation else -27.
        if comb == '2024':
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
                     color=meas.colour if comb != 'lhcb-only' else 'k')
        plt.errorbar(meas.m,
                     n_meas - 1 - i,
                     xerr=meas.stat,
                     capsize=7,
                     color=meas.colour if comb != 'lhcb-only' else 'k')

        if dy_notation:
            x_text = 26 if 'lhcb-only' in comb else 30
            y_text = n_meas - i - 1.25
        else:
            x_text = 22 if i == 3 else 9
            y_text = n_meas - 0.8 if i == 0 else n_meas - i - 1.25
        plt.text(x_text, y_text, f'{meas.label}\n{meas.value()}', color=meas.colour if comb != 'lhcb-only' else 'k')

    # Vertical line for world average
    plt.plot([measures[-1].m, measures[-1].m], [y_min, y_max], linestyle='--', linewidth=1, color='k')

    # Horizontal line to divide world average from other measurements
    plt.plot([x_min, x_max], [0.5, 0.5], linestyle='-', linewidth=1, color='k')

    # Save figure
    os.makedirs(out_dir, exist_ok=True)
    obs = 'deltay' if dy_notation else 'agamma'
    fig_name = f'{obs}-{comb}'
    for ext in ['pdf']:
        plt.savefig(f'{out_dir}/{fig_name}.{ext}')


if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument('-c',
                        '--comb',
                        type=str,
                        default='2024',
                        help='Combination id',
                        choices=['wa-summer-2020', 'wa-2021', 'lhcb-2021', 'wa-2024', 'lhcb-2024'])
    parser.add_argument('--agamma',
                        action='store_false',
                        dest='dy_notation',
                        help='Use the A_Gamma notation rather than DeltaY')
    parser.add_argument('-o', '--outdir', type=str, default='./figs', help='Output directory for saving the plots')
    args = parser.parse_args()

    make_plot(args.comb, args.dy_notation, args.outdir)
