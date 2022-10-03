#!/usr/bin/env python3
'''
Produce 1D and 2D plots for a single combination, or comparing different
combinations. The prob method is employed.

Example (plot the current averages):
    python3 make_plots.py --1d --2d -c 40 --param theoretical
'''

import argparse
from copy import copy
import subprocess


not_phen_params = ['x12', 'y12', 'phiG', 'phiM']
not_supw_params = ['x', 'y', 'qop', 'phi', 'phiG', 'Acp_KK', 'Acp_PP', 'Acp_KP', 'cot_delta_KK', 'cot_delta_PP']
not_theo_params = ['x', 'y', 'qop', 'phi']
full_fsc_params = ['cot_delta_KK', 'cot_delta_PP']
cpv_decay_params = ['Acp_KK', 'Acp_PP', 'cot_delta_KK', 'cot_delta_PP']
full_fsc_combinations = [42]


class Scan1d:
    def __init__(self, var, x_min=None, x_max=None, npoints=200):
        self.var = var
        self.x_min = x_min
        self.x_max = x_max
        self.npoints = npoints
        self.phenomenological = var not in not_phen_params
        self.superweak = var not in not_supw_params
        self.theoretical = var not in not_theo_params
        self.full_fsc_only = var in full_fsc_params


class Scan2d:
    def __init__(self, varx, vary, x_min=None, x_max=None, y_min=None, y_max=None, npoints_x=80, npoints_y=80):
        self.id = f'{varx}_{vary}'
        self.var_x = varx
        self.var_y = vary
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.npoints_x = npoints_x
        self.npoints_y = npoints_y
        self.phenomenological = varx not in not_phen_params and vary not in not_phen_params
        self.superweak = varx not in not_supw_params and vary not in not_supw_params
        self.theoretical = varx not in not_theo_params and vary not in not_theo_params
        self.full_fsc_only = varx in full_fsc_params or vary in full_fsc_params


def get_scans_1d(scan_dy_rs=False):
    scans = [
        Scan1d('x12',           0.,   0.8 ),
        Scan1d('y12',           0.4,  0.8 ),
        Scan1d('phiM',         -0.2,  0.2 ),
        Scan1d('phiG',         -0.2,  0.2 ),
        Scan1d('qop',          -0.15, 0.15),
        Scan1d('phi',          -0.2,  0.2 ),
        Scan1d('x',             0.,   0.8 ),
        Scan1d('y',             0.4,  0.8 ),
        Scan1d('Delta_Kpi',    -0.7,  0.7 ),
        Scan1d('Delta_Kpipi',  -2.,   2.  ),
        Scan1d('R_Kpi'                    ),
        Scan1d('Acp_KP',                  ),
        Scan1d('Acp_KK',       -0.5,  0.5 ),
        Scan1d('Acp_PP',       -0.4,  0.7 ),
        Scan1d('cot_delta_KK', -4e2,  4e2 ),
        Scan1d('cot_delta_PP', -4e2,  4e2 ),
        ]
    if scan_dy_rs:
        scans.append(Scan1d('DY_RS', 0., 0.4))
    return scans


def get_scans_2d():
    scans = [
        Scan2d('x12', 'y12', 0., 0.8, 0.4, 0.8),
        Scan2d('phiM', 'phiG', -0.2, 0.2, -0.2, 0.2),
        Scan2d('x', 'y', 0., 0.8, 0.4, 0.8),
        Scan2d('qop', 'phi', -0.2, 0.2, -0.2, 0.2),
        Scan2d('Acp_KK', 'Acp_PP', -0.3, 0.3, -0.3, 0.5),
        ]
    return scans


def get_base_options(param, combinations, titles, parfile, already_plotted):
    options = [f'--param {param}', '--ps 1', '--pulls', '--pr', '--qh 29', '--group off']
    for i in range(len(combinations)):
        options.append('-c {}{}'.format(
                    combinations[i],
                    '' if titles is None else ' --title "{}"'.format(titles[i])))
    if len(combinations) == 1:
        options.append('--leg off')
    if parfile is not None:
        options.append(f'--parfile {parfile}')
    if already_plotted:
        options.append('-a plot')
    return options


def get_options_1d(param, combinations, titles, parfile, already_plotted):
    return get_base_options(param, combinations, titles, parfile, already_plotted)


def get_options_2d(param, combinations, titles, parfile, already_plotted, n_contours=2):
    options = get_base_options(param, combinations, titles, parfile, already_plotted)
    options.extend(['--magnetic', '--2dcl 1', f'--ncontours {n_contours}', '--ps 2'])
    if len(combinations) > 1:
        options.append('--leg 0.18:0.25')
    return options


def execute_cmd(command):
    print(f'Execute `{command}`')
    subprocess.run(command, shell=True)


def plot_1d(param, combinations, titles, parfile, already_plotted, scan_dy_rs):
    for scan in get_scans_1d():
        cc = copy(combinations)
        tt = copy(titles)
        i = 0
        while i < len(cc):
            if ((scan.full_fsc_only and cc[i] not in full_fsc_combinations) or
                    (scan.var in cpv_decay_params and ((not isinstance(cc[i], int)) or cc[i] < 40))):  # TODO
                cc.pop(i)
                if tt is not None:
                    tt.pop(i)
            else:
                i += 1
        if len(cc) == 0 or not getattr(scan, param):
            continue
        options = get_options_1d(param, cc, tt, parfile, already_plotted)
        options.append(f'--var {scan.var}')
        options.append(f'--npoints {scan.npoints}')
        if scan.x_min is not None and scan.x_max is not None:
            options.append(f'--scanrange {scan.x_min}:{scan.x_max}')
        if scan.var == 'DY_RS':
            if scan_dy_rs:
                options.append('--teststat 2')  # one-sided CL
            else:
                continue
        execute_cmd(f'bin/charmcombo {" ".join(options)}')


def plot_2d(param, combinations, titles, parfile, already_plotted):
    for scan in get_scans_2d():
        cc = copy(combinations)
        tt = copy(titles)
        i = 0
        while i < len(cc):
            if ((scan.full_fsc_only and cc[i] not in full_fsc_combinations) or
                    (any(v in cpv_decay_params for v in [scan.var_x, scan.var_y]) and
                        ((not isinstance(cc[i], int)) or cc[i] < 40))):  # TODO
                cc.pop(i)
                if tt is not None:
                    tt.pop(i)
            else:
                i += 1
        if len(cc) == 0 or not getattr(scan, param):
            continue
        options = get_options_2d(param, cc, tt, parfile, already_plotted)
        options.append(f'--var {scan.var_x}')
        options.append(f'--var {scan.var_y}')
        options.append(f'--npoints2dx {scan.npoints_x}')
        options.append(f'--npoints2dy {scan.npoints_y}')
        if scan.x_min is not None and scan.x_max is not None:
            options.append(f'--scanrange {scan.x_min}:{scan.x_max}')
        if scan.y_min is not None and scan.y_max is not None:
            options.append(f'--scanrangey {scan.y_min}:{scan.y_max}')
        execute_cmd(f'bin/charmcombo {" ".join(options)}')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--param', choices=['phenomenological', 'theoretical', 'superweak'], default='theoretical', help='Parametrisation to be used in the combination')
    parser.add_argument('-c', '--combinations', nargs='+', default=[41], help='Combinations to plot')
    parser.add_argument('-t', '--titles', nargs='+', default=None, help='Titles for the combinations to plot. If none, the default values are used.')
    parser.add_argument('--1d', action='store_true', dest='plot_1d', help='Make the 1d plots')
    parser.add_argument('--2d', action='store_true', dest='plot_2d', help='Make the 2d plots')
    parser.add_argument('--parfile', default=None, help='File with initial value of the parameters')
    parser.add_argument('-a', '--already-plotted', action='store_true', help='Do not rerun the combination and just plot the results retrieving them from a previous execution.')
    parser.add_argument('--scan-dy-rs', action='store_true', help='Provide predictions for the value of DY(K- pi+).')
    args = parser.parse_args()
    if args.titles is not None:
        assert(len(args.titles) == len(args.combinations))
    fn_args = copy(vars(args))
    del fn_args['plot_1d'], fn_args['plot_2d']
    if args.plot_1d:
        plot_1d(**fn_args)
    del fn_args['scan_dy_rs']
    if args.plot_2d:
        plot_2d(**fn_args)


if __name__ == '__main__':
    parse_args()
