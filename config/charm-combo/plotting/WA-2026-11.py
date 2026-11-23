"""Configuration file for summary plots of the World Average of November 2026."""

from charm_fitter.config import BREAKDOWN_COMBINERS_KPI_RD, BREAKDOWN_STYLES, PARAMETERS, RANGES_1D, RANGES_2D
from charm_fitter.utils import CharmPlottingConfig, Combiner, DYFsc, Plot1D, Plot2D, Range, Scan1D, Scan2D

# Parameters and combiners ---------------------------------------------------------------------------------------------

_PARAMETERS_BY_NAME = {p.name: p for p in PARAMETERS}

_BASELINE_COMBINER = "WA-2026-11"
_DY_FSCS = [DYFsc.FULL, DYFsc.NONE]
_DY_FSC_1D_VARS = ["phiM", "phiG", "qop", "phi", "Acp_KK", "Acp_PP"]
_DY_FSC_2D_VARS = [("phiM", "phiG"), ("qop", "phi"), ("Acp_KK", "Acp_PP"), ("Acp_KK", "phiM")]

_FIX_PARFILE = "config/charm-combo/start/WA-2026-11.dat"

# Ranges of the LHCb Run 1+2 average, adapted to the larger value of x (x12) and the improved precision on CPV in mixing
_RANGES_1D = RANGES_1D["LHCb-R12"] | {
    "x": Range((2.0, 7.5), unit="1e-3"),
    "x12": Range((2.0, 7.5), unit="1e-3"),
    "phi": Range((-90.0, 70.0), unit="mrad"),
    "phiM": Range((-60.0, 60.0), unit="mrad"),
    "phiG": Range((-120.0, 150.0), unit="mrad"),
}
_RANGES_2D = RANGES_2D["LHCb-R12"] | {
    ("x12", "y12"): (Range((2.5, 7.0), unit="1e-3"), Range((5.4, 7.2), unit="1e-3")),
    ("phiM", "phiG"): (Range((-50.0, 55.0), unit="mrad"), Range((-80.0, 130.0), unit="mrad")),
    ("x", "y"): (Range((2.5, 7.0), unit="1e-3"), Range((5.4, 7.2), unit="1e-3")),
    ("qop", "phi"): (Range((0.93, 1.05)), Range((-70.0, 50.0), unit="mrad")),
    ("Acp_KK", "phiM"): (Range((-2.0, 3.5), unit="1e-3"), Range((-50.0, 55.0), unit="mrad")),
}

_BASELINE_COMBINERS = {
    _BASELINE_COMBINER: Combiner(
        13,
        "WA November 2026",
        title_long="World Average November 2026",
        scans_1d=[Scan1D(par, r, parfile=_FIX_PARFILE) for par, r in _RANGES_1D.items()]
        + [
            Scan1D(par, _RANGES_1D[par], dy_fsc=dy_fsc, parfile=_FIX_PARFILE)
            for dy_fsc in _DY_FSCS
            if dy_fsc != DYFsc.NONE
            for par in _DY_FSC_1D_VARS
        ],
        scans_2d=[
            Scan2D(xpar, ypar, xrange, yrange, parfile=_FIX_PARFILE)
            for (xpar, ypar), (xrange, yrange) in _RANGES_2D.items()
        ]
        + [
            Scan2D(
                xpar,
                ypar,
                _RANGES_2D[(xpar, ypar)][0],
                _RANGES_2D[(xpar, ypar)][1],
                dy_fsc=dy_fsc,
                parfile=_FIX_PARFILE,
            )
            for dy_fsc in _DY_FSCS
            if dy_fsc != DYFsc.NONE
            for (xpar, ypar) in _DY_FSC_2D_VARS
        ],
        col="r",
        ls="-",
    ),
}

_COMBINERS_BREAKDOWN = {
    **BREAKDOWN_COMBINERS_KPI_RD,
    "kpi-mix": Combiner(
        [1, 21, 40, 41, 250, 252, 260, 286, 288],
        **BREAKDOWN_STYLES["kpi"],
        scans_2d=[
            Scan2D(
                "x",
                "y",
                xrange=Range((-8e-3, 8e-3)),
                yrange=Range((3e-3, 8e-3)),
                fixed_pars=["qop", "phi"],
                fix_parfile=_FIX_PARFILE,
                parfile="config/charm-combo/start/WA-2025-06-breakdown-kpi-x-y.dat",
                extra_args="--npoints2dx 100 --npoints2dy 100",
            ),
            Scan2D(
                "x12",
                "y12",
                xrange=Range((0.0, 8e-3)),
                yrange=Range((3e-3, 8e-3)),
                fixed_pars=["phiM", "phiG"],
                fix_parfile=_FIX_PARFILE,
                parfile="config/charm-combo/start/WA-2025-06-breakdown-kpi-x12-y12.dat",
                extra_args="--npoints2dx 100 --npoints2dy 100",
            ),
        ],
    ),
    "kpi-cpv": Combiner(
        [1, 21, 40, 41, 260, 261, 286, 288],
        **BREAKDOWN_STYLES["kpi"],
        scans_2d=[
            Scan2D(
                "qop",
                "phi",
                xrange=Range((0.7, 1.4)),
                yrange=Range((-0.4, 0.35)),
                fixed_pars=["x", "y"],
                fix_parfile=_FIX_PARFILE,
                extra_args="--npoints2dx 100 --npoints2dy 100",
            ),
            Scan2D(
                "phiM",
                "phiG",
                xrange=Range((-0.1, 0.1)),
                yrange=Range((-0.41, 0.3)),
                fixed_pars=["x12", "y12"],
                fix_parfile=_FIX_PARFILE,
            ),
        ],
    ),
    "k3pi": Combiner(
        [60, 371, 372],
        **BREAKDOWN_STYLES["k3pi"],
        scans_2d=[
            Scan2D(
                "x",
                "y",
                xrange=Range((-2e-3, 10e-3)),
                yrange=Range((0.0, 10e-3)),
                fixed_pars=["qop", "phi"],
                fix_parfile=_FIX_PARFILE,
            ),
            Scan2D(
                "x12",
                "y12",
                xrange=Range((0.0, 10e-3)),
                yrange=Range((0.0, 10e-3)),
                fixed_pars=["phiM", "phiG"],
                fix_parfile=_FIX_PARFILE,
            ),
            Scan2D(
                "qop",
                "phi",
                xrange=Range((0.9, 1.1)),
                yrange=Range((-0.2, 0.2)),
                fixed_pars=["x", "y"],
                fix_parfile=_FIX_PARFILE,
            ),
            Scan2D(
                "phiM",
                "phiG",
                xrange=Range((-0.2, 0.2)),
                yrange=Range((-0.2, 0.2)),
                fixed_pars=["x12", "y12"],
                fix_parfile=_FIX_PARFILE,
            ),
        ],
    ),
    "kspipi": Combiner(
        [320, 322, 327, 328, 329],
        **BREAKDOWN_STYLES["kspipi"],
        scans_2d=[
            Scan2D(
                "x",
                "y",
                xrange=Range((2e-3, 7e-3)),
                yrange=Range((2e-3, 9e-3)),
                fixed_pars=["qop", "phi"],
                fix_parfile=_FIX_PARFILE,
            ),
            Scan2D(
                "x12",
                "y12",
                xrange=Range((2e-3, 7e-3)),
                yrange=Range((2e-3, 9e-3)),
                fixed_pars=["phiM", "phiG"],
                fix_parfile=_FIX_PARFILE,
            ),
            Scan2D(
                "qop",
                "phi",
                xrange=Range((0.85, 1.2)),
                yrange=Range((-0.25, 0.1)),
                fixed_pars=["x", "y"],
                fix_parfile=_FIX_PARFILE,
            ),
            Scan2D(
                "phiM",
                "phiG",
                xrange=Range((-0.2, 0.3)),
                yrange=Range((-0.1, 0.3)),
                fixed_pars=["x12", "y12"],
                fix_parfile=_FIX_PARFILE,
            ),
        ],
    ),
    "hh-mix": Combiner(
        [110, 130, 162],  # yCP only
        **BREAKDOWN_STYLES["hh-mix"],
        scans_2d=[
            Scan2D(
                "x",
                "y",
                xrange=Range((-8e-3, 12.5e-3)),
                yrange=Range((5e-3, 8e-3)),
                fixed_pars=["r_Kpi", "Delta_Kpi", "qop", "phi"],
                fix_parfile=_FIX_PARFILE,
                extra_args="--npoints2dx 100",
            ),
            Scan2D(
                "x12",
                "y12",
                xrange=Range((0.0, 10.5e-3)),
                yrange=Range((5e-3, 8e-3)),
                fixed_pars=["r_Kpi", "Delta_Kpi", "phiM", "phiG"],
                fix_parfile=_FIX_PARFILE,
                extra_args="--npoints2dx 100",
            ),
        ],
    ),
    "hh": Combiner(
        [110, 130, 162, 202],  # yCP + DeltaY
        **BREAKDOWN_STYLES["hh"],
        scans_2d=[
            Scan2D(
                "qop",
                "phi",
                xrange=Range((0.78, 1.3)),
                yrange=Range((-0.4, 0.4)),
                fixed_pars=["x", "y", "r_Kpi", "Delta_Kpi"],
                fix_parfile=_FIX_PARFILE,
                extra_args="--npoints2dx 80 --npoints2dy 80",
            ),
            Scan2D(
                "phiM",
                "phiG",
                xrange=Range((-0.1, 0.15)),
                yrange=Range((-0.5, 0.5)),
                fixed_pars=["x12", "y12", "r_Kpi", "Delta_Kpi"],
                fix_parfile=_FIX_PARFILE,
            ),
        ],
    ),
}

# 1D plots -------------------------------------------------------------------------------------------------------------

_1D_LEG_RIGHT = ["x", "y", "qop", "x12", "y12", "Delta_Kpi", "r_Kpi", "r_Kpi", "r_K3pi", "r_Kpipi0", "k_K3pi", "R_Kpi_"]
_1D_LEGFILL = ["phi", "phiM", "Acp_KK", "Acp_PP", "r_K3pi"]
_PLOTS_1D = [
    Plot1D(
        [_BASELINE_COMBINER],
        name,
        _RANGES_1D[par.cf_par].transformed_range(_PARAMETERS_BY_NAME[name].transf),
        legpos="r" if name in _1D_LEG_RIGHT else "l",
        legfill=name in _1D_LEGFILL,
        legfontsize=14,
    )
    for name, par in _PARAMETERS_BY_NAME.items()
    if par.cf_par in _RANGES_1D
]

# 2D plots -------------------------------------------------------------------------------------------------------------

_PLOTS_2D = [
    Plot2D(
        [_BASELINE_COMBINER],
        (xpar, ypar),
        xrange.transformed_range(_PARAMETERS_BY_NAME[xpar].transf),
        yrange.transformed_range(_PARAMETERS_BY_NAME[ypar].transf)
        if ypar != "Acp_PP"
        else Range((-0.5, 5.5), unit="1e-3"),
        legfill=xpar in ["Acp_KK", "phiM"],
    )
    for (xpar, ypar), (xrange, yrange) in _RANGES_2D.items()
]


# Breakdown 2D plots ---------------------------------------------------------------------------------------------------

# N.B.: This is tuned only for the hypothesis of no CPV in D0 -> K+ pi- decays.
_BREAKDOWN_CPV_COMBOS = ["kspipi", "k3pi", "hh", "kpi-cpv", _BASELINE_COMBINER]
_BREAKDOWN_MIX_COMBOS = ["kpi-mix", "k3pi", "kspipi", "hh-mix", _BASELINE_COMBINER]
_PLOTS_BREAKDOWN = [
    Plot2D(
        _BREAKDOWN_MIX_COMBOS,
        ("x", "y"),
        xrange=Range((-8.0, 10.0), unit="1e-3"),
        yrange=Range((2.0, 10.0), unit="1e-3"),
        combiners_label=_BASELINE_COMBINER,
        levels=2,
        logo="bl",
        legncol=2,
    ),
    Plot2D(
        _BREAKDOWN_MIX_COMBOS,
        ("x12", "y12"),
        xrange=Range((0.0, 10.0), unit="1e-3"),
        yrange=Range((2.0, 10.5), unit="1e-3"),
        combiners_label=_BASELINE_COMBINER,
        levels=2,
        logo="bl",
        legncol=2,
    ),
    Plot2D(
        _BREAKDOWN_CPV_COMBOS,
        ("qop", "phi"),
        xrange=Range((-0.25, 0.3)),
        yrange=Range((-350.0, 400.0), unit="mrad"),
        logo="br",
        combiners_label=_BASELINE_COMBINER,
        levels=2,
    ),
    Plot2D(
        _BREAKDOWN_CPV_COMBOS,
        ("phiM", "phiG"),
        xrange=Range((-320.0, 250.0), unit="mrad"),
        yrange=Range((-400.0, 480.0), unit="mrad"),
        combiners_label=_BASELINE_COMBINER,
        levels=2,
    ),
    Plot2D(
        ["cf-kpi", "d0-kpi", "hhprime", _BASELINE_COMBINER],
        ("Delta_Kpi", "r_Kpi"),
        xrange=Range((-1.3, 0.7), unit="rad"),
        yrange=Range((5.6, 6.1), unit="1e-2"),
        combiners_label=_BASELINE_COMBINER,
        levels=2,
        logo="l",
        legpos="bl",
        legfill=True,
    ),
]

# Plots for comparing the impact of different DeltaY(h- h+) final-state correction hypotheses --------------------------

_PLOTS_DY_FSC_1D = [
    Plot1D(
        [_BASELINE_COMBINER],
        par,
        _RANGES_1D[par].transformed_range(_PARAMETERS_BY_NAME[par].transf),
        legpos="r",
        legfontsize=12,
        dy_fsc=_DY_FSCS,
    )
    for par in _DY_FSC_1D_VARS
]

# N.B.: This is tuned only for the hypothesis of no CPV in D0 -> K+ pi- decays.
_RANGES_DY_FSC = {
    ("phiM", "phiG"): (Range((-60.0, 80.0), unit="mrad"), Range((-80.0, 180.0), unit="mrad")),
    ("qop", "phi"): (Range((-0.08, 0.05)), Range((-90.0, 60.0), unit="mrad")),
    ("Acp_KK", "Acp_PP"): (Range((-1.5, 3.0), unit="1e-3"), Range((-0.5, 6.5), unit="1e-3")),
    ("Acp_KK", "phiM"): (Range((-2.0, 3.5), unit="1e-3"), Range((-60.0, 70.0), unit="mrad")),
}

_PLOTS_DY_FSC_2D = [
    Plot2D(
        [_BASELINE_COMBINER],
        (xpar, ypar),
        xrange,
        yrange,
        logo="br",
        legpos="l",
        legfontsize=12,
        dy_fsc=_DY_FSCS,
    )
    for (xpar, ypar), (xrange, yrange) in _RANGES_DY_FSC.items()
]

CONFIG = CharmPlottingConfig(
    PARAMETERS,
    _BASELINE_COMBINERS,
    combiners_breakdown=_COMBINERS_BREAKDOWN,
    plots_1d=_PLOTS_1D,
    plots_2d=_PLOTS_2D,
    plots_breakdown=_PLOTS_BREAKDOWN,
    plots_dy_fsc_1d=_PLOTS_DY_FSC_1D,
    plots_dy_fsc_2d=_PLOTS_DY_FSC_2D,
)
