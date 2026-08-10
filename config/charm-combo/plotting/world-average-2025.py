"""Configuration file for summary plots of the World Average of 2025."""

from charm_fitter.config import PARAMETERS, RANGES_1D, RANGES_2D
from charm_fitter.utils import CharmPlottingConfig, Combiner, DYFsc, MixParam, Plot1D, Plot2D, Range, Scan1D, Scan2D

# Parameters and combiners ---------------------------------------------------------------------------------------------

_PARAMETERS_BY_NAME = {p.name: p for p in PARAMETERS}

_BASELINE_COMBINER = "WA-2025-10"
_DY_FSCS = [DYFsc.FULL, DYFsc.NONE]
_DY_FSC_1D_VARS = ["phiM", "phiG", "qop", "phi", "Acp_KK", "Acp_PP"]
_DY_FSC_2D_VARS = [("phiM", "phiG"), ("qop", "phi"), ("Acp_KK", "Acp_PP"), ("Acp_KK", "phiM")]

_BASELINE_COMBINERS = {
    _BASELINE_COMBINER: Combiner(
        55,
        "WA Oct 2025",
        title_long="World Average October 2025",
        scans_1d=[Scan1D(par, r, parfile="config/start/wa-2025.dat") for par, r in RANGES_1D["LHCb-R12"].items()]
        + [
            Scan1D(par, RANGES_1D["LHCb-R12"][par], dy_fsc=dy_fsc, parfile="config/start/wa-2025.dat")
            for dy_fsc in _DY_FSCS
            if dy_fsc != DYFsc.NONE
            for par in _DY_FSC_1D_VARS
        ],
        scans_2d=[
            Scan2D(xpar, ypar, xrange, yrange, parfile="config/start/wa-2025.dat")
            for (xpar, ypar), (xrange, yrange) in RANGES_2D["LHCb-R12"].items()
        ]
        + [
            Scan2D(
                xpar,
                ypar,
                RANGES_2D["LHCb-R12"][(xpar, ypar)][0],
                RANGES_2D["LHCb-R12"][(xpar, ypar)][1],
                dy_fsc=dy_fsc,
                parfile="config/start/wa-2025.dat",
            )
            for dy_fsc in _DY_FSCS
            if dy_fsc != DYFsc.NONE
            for (xpar, ypar) in _DY_FSC_2D_VARS
        ],
        col="r",
        ls="-",
    ),
}

_FIX_PARFILE = "config/start/wa-2025.dat"

_COMBINERS_BREAKDOWN = {
    "ext-kpi": Combiner(
        [50, 53, 55, 56],
        r"CLEO/BES $D^0 \to K^\pm \pi^\mp$",
        col="y",
        ls="dash3dots",
        scans_2d=[
            Scan2D(
                "Delta_Kpi",
                "r_Kpi",
                xrange=Range((-0.6, 0.3)),
                yrange=Range((5.55e-2, 6.15e-2)),
                mix_param=MixParam.PHENO,
            )
        ],
    ),
    "lhcb-kpi": Combiner(
        [39, 41],
        r"LHCb $D^0 \to K^\pm \pi^\mp$",
        col="g",
        ls="dash2dots",
        scans_2d=[
            Scan2D(
                "Delta_Kpi",
                "r_Kpi",
                xrange=Range((-1.35, 0.55)),
                yrange=Range((5.81e-2, 5.92e-2)),
                parfile="config/start/charm-2025-breakdown-kpi-rd.dat",
                extra_args="--npoints2dx 200",
            ),
        ],
    ),
    "kpi": Combiner(
        [39, 41, 50, 53, 55, 56],
        r"$D^0 \to K^\pm \pi^\mp$",
        col="g",
        ls="dash2dots",
        scans_2d=[
            Scan2D(
                "x",
                "y",
                xrange=Range((-8e-3, 8e-3)),
                yrange=Range((3e-3, 8e-3)),
                fixed_pars=["qop", "phi"],
                fix_parfile=_FIX_PARFILE,
                parfile="config/start/charm-2025-breakdown-kpi-x-y.dat",
                extra_args="--npoints2dx 100 --npoints2dy 100",
            ),
            Scan2D(
                "x12",
                "y12",
                xrange=Range((0.0, 8e-3)),
                yrange=Range((3e-3, 8e-3)),
                fixed_pars=["phiM", "phiG"],
                fix_parfile=_FIX_PARFILE,
                parfile="config/start/charm-2025-breakdown-kpi-x12-y12.dat",
                extra_args="--npoints2dx 100 --npoints2dy 100",
            ),
            Scan2D(
                "qop",
                "phi",
                xrange=Range((0.7, 1.4)),
                yrange=Range((-0.4, 0.35)),
                fixed_pars=["x", "y"],
                fix_parfile=_FIX_PARFILE,
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
        [5, 54],
        r"$D^0 \to K^\pm \pi^\mp \pi^+ \pi^-$",
        col="p",
        ls="dash3dots",
        scans_2d=[
            Scan2D(
                "x",
                "y",
                xrange=Range((1e-3, 6e-3)),
                yrange=Range((2e-3, 10e-3)),
                fixed_pars=["qop", "phi"],
                fix_parfile=_FIX_PARFILE,
            ),
            Scan2D(
                "x12",
                "y12",
                xrange=Range((1e-3, 6e-3)),
                yrange=Range((2e-3, 10e-3)),
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
                xrange=Range((-0.2, 0.35)),
                yrange=Range((-0.1, 0.25)),
                fixed_pars=["x12", "y12"],
                fix_parfile=_FIX_PARFILE,
            ),
        ],
    ),
    "kspipi": Combiner(
        [1, 3, 6, 21, 24],
        r"$D^0 \to K_S^0 \pi^+ \pi^-$",
        col="lb",
        ls="shortdash",
        scans_2d=[
            Scan2D(
                "x",
                "y",
                xrange=Range((2e-3, 6e-3)),
                yrange=Range((2e-3, 7e-3)),
                fixed_pars=["qop", "phi"],
                fix_parfile=_FIX_PARFILE,
            ),
            Scan2D(
                "x12",
                "y12",
                xrange=Range((2e-3, 6e-3)),
                yrange=Range((2e-3, 10e-3)),
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
                xrange=Range((-0.2, 0.35)),
                yrange=Range((-0.05, 0.15)),
                fixed_pars=["x12", "y12"],
                fix_parfile=_FIX_PARFILE,
            ),
        ],
    ),
    "hh-mix": Combiner(
        [61, 62, 65],  # yCP only
        r"$D^0 \to h^+ h^-$",
        col="y",
        ls="dashdot",
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
        [61, 62, 65, 72],  # yCP + DeltaY
        r"$D^0 \to h^+ h^-$",
        col="y",
        ls="dashdot",
        scans_2d=[
            Scan2D(
                "qop",
                "phi",
                xrange=Range((0.8, 1.3)),
                yrange=Range((-0.4, 0.4)),
                fixed_pars=["x", "y", "r_Kpi", "Delta_Kpi"],
                fix_parfile=_FIX_PARFILE,
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
        RANGES_1D["LHCb-R12"][par.cf_par].transformed_range(_PARAMETERS_BY_NAME[name].transf),
        legpos="r" if name in _1D_LEG_RIGHT else "l",
        legfill=name in _1D_LEGFILL,
        legfontsize=14,
    )
    for name, par in _PARAMETERS_BY_NAME.items()
    if par.cf_par in RANGES_1D["LHCb-R12"]
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
    for (xpar, ypar), (xrange, yrange) in RANGES_2D["LHCb-R12"].items()
]


# Breakdown 2D plots ---------------------------------------------------------------------------------------------------

# N.B.: This is tuned only for the hypothesis of no CPV in D0 -> K+ pi- decays.
_BREAKDOWN_CPV_COMBOS = ["kspipi", "hh", "kpi", _BASELINE_COMBINER]
_BREAKDOWN_MIX_COMBOS = ["kpi", "kspipi", "hh-mix", _BASELINE_COMBINER]
_PLOTS_BREAKDOWN = [
    Plot2D(
        _BREAKDOWN_MIX_COMBOS,
        ("x", "y"),
        xrange=Range((-7.0, 12.0), unit="1e-3"),
        yrange=Range((2.0, 10.0), unit="1e-3"),
        combiners_label=_BASELINE_COMBINER,
        levels=2,
    ),
    Plot2D(
        _BREAKDOWN_MIX_COMBOS,
        ("x12", "y12"),
        xrange=Range((0.0, 10.0), unit="1e-3"),
        yrange=Range((2.0, 10.0), unit="1e-3"),
        combiners_label=_BASELINE_COMBINER,
        levels=2,
    ),
    Plot2D(
        _BREAKDOWN_CPV_COMBOS,
        ("qop", "phi"),
        xrange=Range((-0.2, 0.3)),
        yrange=Range((-400.0, 350.0), unit="mrad"),
        logo="br",
        combiners_label=_BASELINE_COMBINER,
        levels=2,
    ),
    Plot2D(
        _BREAKDOWN_CPV_COMBOS,
        ("phiM", "phiG"),
        xrange=Range((-400.0, 350.0), unit="mrad"),
        yrange=Range((-400.0, 400.0), unit="mrad"),
        combiners_label=_BASELINE_COMBINER,
        levels=2,
    ),
    Plot2D(
        ["ext-kpi", "lhcb-kpi", _BASELINE_COMBINER],
        ("Delta_Kpi", "r_Kpi"),
        xrange=Range((-1.3, 0.5), unit="rad"),
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
        RANGES_1D["LHCb-R12"][par].transformed_range(_PARAMETERS_BY_NAME[par].transf),
        legpos="r",
        legfontsize=12,
        dy_fsc=_DY_FSCS,
    )
    for par in _DY_FSC_1D_VARS
]

# N.B.: This is tuned only for the hypothesis of no CPV in D0 -> K+ pi- decays.
_RANGES_DY_FSC = {
    ("phiM", "phiG"): (Range((-75.0, 110.0), unit="mrad"), Range((-80.0, 200.0), unit="mrad")),
    ("qop", "phi"): (Range((-0.1, 0.05)), Range((-120.0, 60.0), unit="mrad")),
    ("Acp_KK", "Acp_PP"): (Range((-1.5, 3.0), unit="1e-3"), Range((-0.5, 6.0), unit="1e-3")),
    ("Acp_KK", "phiM"): (Range((-2.0, 3.5), unit="1e-3"), Range((-90.0, 80.0), unit="mrad")),
}

_PLOTS_DY_FSC_2D = [
    Plot2D(
        [_BASELINE_COMBINER],
        (xpar, ypar),
        xrange,
        yrange,
        logo="br",
        legfill="Acp_PP" == ypar,
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
