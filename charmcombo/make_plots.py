#!/usr/bin/env python3
'''
Produce the 1D and 2D plots comparing the world averages in 2020 and 2021.
The prob method is employed.
'''

import subprocess
from copy import copy

# --- configuration -----------------------------------------------------------
th_cfg = 'theoretical'  # 'theoretical', 'phenomenological', 'superweak'
combo_to_plot = [3, 1]
titles = [
    'World average (Dec 2020)',
    'World average (Feb 2022)',
]  # or `None` for default combination names
plot_1d = True
plot_2d = False
n_contours = 2
already_plotted = False  # save time if True
use_start_file = False
scan_dy_rs = False
# -----------------------------------------------------------------------------


assert(th_cfg in ['phenomenological', 'theoretical', 'superweak'])
assert(len(combo_to_plot) > 0)


class Parameter:
    def __init__(self, var, x_min=None, x_max=None, npoints=200):
        self.var = var
        self.x_min = x_min
        self.x_max = x_max
        self.npoints = npoints


class Scan1d:
    def __init__(self, par, pheno=True, theo=True, superweak=True, options=None):
        self.var = par.var
        self.x_min = par.x_min
        self.x_max = par.x_max
        self.npoints = par.npoints
        self.phenomenological = pheno
        self.theoretical = theo
        self.superweak = superweak
        self.options = options

        
class Scan2d:
    def __init__(self, par1, par2, pheno=True, theo=True, superweak=True, options=None):
        self.par1 = par1
        self.par2 = par2
        self.phenomenological = pheno
        self.theoretical = theo
        self.superweak = superweak
        self.options = options


params = {
    'x12':         Parameter('x12',          0.,   0.8 ),
    'y12':         Parameter('y12',          0.4,  0.8 ),
    'x12_for2d':   Parameter('x12',          0.,   0.8 ,npoints=100),
    'y12_for2d':   Parameter('y12',          0.4,  0.8 ,npoints=100),
    'y12_for_x12': Parameter('y12',          0.2,  1   ),
    'phiM':        Parameter('phiM',        -0.2,  0.2 ),
    'phiG':        Parameter('phiG',        -0.2,  0.2 ),
    'qop':         Parameter('qop',         -0.2,  0.2 ),
    'phi':         Parameter('phi',         -0.5,  0.5 ),
    'x':           Parameter('x',            0.,   1.  ),
    'y':           Parameter('y',            0.,   1.  ),
    'Delta_Kpi':   Parameter('Delta_Kpi',   -0.7,  0.7 ),
    'Delta_Kpipi': Parameter('Delta_Kpipi', -2.,   2.  ),
    'R_Kpi':       Parameter('R_Kpi'),
    'A_Kpi':       Parameter('A_Kpi'),
    'DY_RS':       Parameter('DY_RS', 0., 0.4, npoints=400),
}


scan_1d = {
    'x12':       Scan1d(params['x12'], pheno=False),
    'y12':       Scan1d(params['y12'], pheno=False),
    'phiM':      Scan1d(params['phiM'], pheno=False),
    'phiG':      Scan1d(params['phiG'], pheno=False, superweak=False),
    'qop':       Scan1d(params['qop'], theo=False, superweak=False),
    'phi':       Scan1d(params['phi'], theo=False, superweak=False),
    'x':         Scan1d(params['x'], theo=False, superweak=False),
    'y':         Scan1d(params['y'], theo=False, superweak=False),
    'Delta_Kpi': Scan1d(params['Delta_Kpi']),
    'Delta_Kpipi': Scan1d(params['Delta_Kpipi']),
    'R_Kpi':     Scan1d(params['R_Kpi']),
    'A_Kpi':     Scan1d(params['A_Kpi'], superweak=False),
}
if scan_dy_rs:
    scan_1d['DY_RS'] = Scan1d(params['DY_RS'])


scan_2d = {
    'x11_y12':   Scan2d(params['x12'], params['y12'], pheno=False),
    'phiM_phiG': Scan2d(params['phiM'], params['phiG'], pheno=False, superweak=False),
    'x12_phiM':  Scan2d(params['x12_for2d'], params['phiM'], pheno=False),
    'x12_phiG':  Scan2d(params['x12_for2d'], params['phiG'], pheno=False, superweak=False),
    'y12_phiM':  Scan2d(params['y12_for2d'], params['phiM'], pheno=False),
    'y12_phiG':  Scan2d(params['y12_for2d'], params['phiG'], pheno=False, superweak=False),
    'qop_phi':   Scan2d(params['qop'], params['phi'], theo=False, superweak=False),
    'x_y':       Scan2d(params['x'], params['y'], theo=False, superweak=False),
    'x_phi':     Scan2d(params['x'], params['phi'], theo=False, superweak=False),
    'y_phi':     Scan2d(params['y'], params['phi'], theo=False, superweak=False),
    'x_qop':     Scan2d(params['x'], params['qop'], theo=False, superweak=False),
    'y_qop':     Scan2d(params['y'], params['qop'], theo=False, superweak=False),
}


base_options = ['--ps 1', '--pulls', '--pr', '--qh 29', '--group off']
if use_start_file:
    base_options.append('--parfile start_files/{}_{}start.dat'.format(
                th_cfg, 'DY_RS_' if scan_dy_rs else ''))
if already_plotted:
    base_options.append('-a plot')
for i in range(len(combo_to_plot)):
    base_options.append('-c {}{}'.format(
                combo_to_plot[i], '' if titles is None else ' --title "{}"'.format(titles[i])))
if len(combo_to_plot) == 1:
    base_options.append('--leg off')
base_options_2d = copy(base_options)
base_options_2d.extend(
        ['--magnetic',
         '--2dcl 1',
         '--ncontours {}'.format(n_contours),
         '--leg 0.18:0.25',
         '--ps 2'])


if plot_1d:
    for k, par in scan_1d.items():
        do_plot = getattr(par, th_cfg)
        if not do_plot:
            continue
        options = copy(base_options)
        options.append('--var {}'.format(par.var))
        options.append('--npoints {}'.format(par.npoints))
        if par.x_min is not None and par.x_max is not None:
            options.append('--scanrange {}:{}'.format(
                        par.x_min, par.x_max))
        if par.var == 'DY_RS':
            if not scan_dy_rs:
                continue
            options.append('--teststat 2')  # one-sided CL
        command = 'bin/charmcombo_{th_cfg} {options}'.format(
                th_cfg=th_cfg,
                options=' '.join(options))
        print('Execute {}'.format(command))
        subprocess.run(command, shell=True)

    
if plot_2d:
    for k, par in scan_2d.items():
        do_plot = getattr(par, th_cfg)
        if not do_plot:
            continue
        options = copy(base_options_2d)
        options.append('--var {}'.format(par.par1.var))
        options.append('--var {}'.format(par.par2.var))
        options.append('--npoints2dx {}'.format(par.par1.npoints))
        options.append('--npoints2dy {}'.format(par.par2.npoints))
        if par.par1.x_min is not None and par.par1.x_max is not None:
            options.append('--scanrange {}:{}'.format(
                        par.par1.x_min, par.par1.x_max))
        if par.par2.x_min is not None and par.par2.x_max is not None:
            options.append('--scanrangey {}:{}'.format(
                        par.par2.x_min, par.par2.x_max))
        if par.options is not None:
            options.extend(par.options)
        command = 'bin/charmcombo_{th_cfg} {options}'.format(
                th_cfg=th_cfg,
                options=' '.join(options))
        print('Execute {}'.format(command))
        subprocess.run(command, shell=True)
