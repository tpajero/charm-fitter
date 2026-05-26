"""Configuration file for LHCb Upgrade I and II sensitivity projections.

All scans are made allowing for final-state corrections to DeltaY(h- h+).
"""

from charm_fitter.config import PARAMETERS
from charm_fitter.utils import (
    ACP_COT_PARAMS,
    R_DELTA_PARAMS,
    CharmPlottingConfig,
    Combiner,
    DYFsc,
    Plot1D,
    Plot2D,
    Range,
    Scan1D,
    Scan2D,
)

_FULL_FSC_PARAMS = ACP_COT_PARAMS + R_DELTA_PARAMS

_EXTRAPOLATIONS = {
    "R12": 1.0,
    "UI": 0.27,
    "UII": 0.12,
}


def _extrapolated_range(period, val, unc, unit=None, is_2d=False):
    factor = _EXTRAPOLATIONS[period] * (5 if is_2d else 3.5)
    return Range((val - factor * unc, val + factor * unc), unit=unit)


_RUN2_VALS = {
    "x12": (3.9, 0.44, "1e-3"),
    "y12": (6.31, 0.224, "1e-3"),
    "phiM": (-6, 15.2, "mrad"),
    "phiG": (40, 27.2, "mrad"),
    "Delta_Kpi": (-0.25, 0.056, "rad"),
    "r_Kpi": (5.862, 0.01, "%"),
    "Acp_KK": (0.9, 0.52, "1e-3"),
    "Acp_PP": (2.8, 0.68, "1e-3"),
    # "cot_delta_KK": Range((-15.0, 15.0)),
    # "cot_delta_PP": Range((-15.0, 15.0)),
    # "r_KK": Range((0.0, 100.0), unit="%"),
    # "r_PP": Range((0.0, 100.0), unit="%"),
    # "delta_KK": Range((-math.pi, 0.0), unit="rad"),
    # "delta_PP": Range((-math.pi, 0.0), unit="rad"),
}

_RANGES_1D = {
    period: {p: _extrapolated_range(period, val, unc, unit) for p, (val, unc, unit) in _RUN2_VALS.items()}
    for period in _EXTRAPOLATIONS
}

_RANGES_2D = {
    period: {
        (xpar, ypar): (
            _extrapolated_range(period, _RUN2_VALS[xpar][0], _RUN2_VALS[xpar][1], unit=_RUN2_VALS[xpar][2], is_2d=True),
            _extrapolated_range(period, _RUN2_VALS[ypar][0], _RUN2_VALS[ypar][1], unit=_RUN2_VALS[ypar][2], is_2d=True),
        )
        for xpar, ypar in [
            ("x12", "y12"),
            ("phiM", "phiG"),
            ("Delta_Kpi", "r_Kpi"),
            ("Acp_KK", "Acp_PP"),
            # ("cot_delta_KK", "cot_delta_PP"),
            # ("r_KK", "r_PP"),
            # ("delta_KK", "delta_PP"),
        ]
    }
    for period in _RANGES_1D
}

# Parameters and combiner ------------------------------------------------------------------------------------------

_PARAMETERS = [p for p in PARAMETERS if p.name in _RANGES_1D["R12"]]
_PARAMETERS_BY_NAME = {p.name: p for p in _PARAMETERS}


_EXTRA_ARGS_UPGRADES = (
    "--pr --truth-parfile config/charm-combo/start/LHCb-UII.dat --parfile config/charm-combo/start/LHCb-UII.dat"
)
_BASELINE_COMBINERS = {
    f"LHCb-{period}": Combiner(
        combiner_id,
        title,
        title_long=title_long,
        col=col,
        scans_1d=[
            Scan1D(
                par,
                r,
                dy_fsc=DYFsc.FULL,
                extra_args=extra_args,
                parfile=f"config/charm-combo/start/LHCb-{'UII' if period == 'UI' else period}.dat",
            )
            for par, r in _RANGES_1D[period].items()
        ],
        scans_2d=[
            Scan2D(
                xpar,
                ypar,
                xrange,
                yrange,
                dy_fsc=DYFsc.FULL,
                parfile=f"config/charm-combo/start/LHCb-{'UII' if period == 'UI' else period}.dat",
                extra_args=extra_args,
            )
            for (xpar, ypar), (xrange, yrange) in _RANGES_2D[period].items()
        ],
    )
    for period, combiner_id, title, title_long, col, extra_args in [
        ("R12", 301, "LHCb Run 1+2", "LHCb Run 1+2", "y", "--pr"),
        ("UI", 400, "LHCb UI", "LHCb Upgrade I", "r", _EXTRA_ARGS_UPGRADES),
        ("UII", 450, "LHCb UII", "LHCb Upgrade II", "lb", _EXTRA_ARGS_UPGRADES),
    ]
}

# 1D plots ----------------------------------------------------------------------------------------------------------

_PLOTS_1D = [
    Plot1D(
        # LHCb-R12 does not have enough precision for determining the strong parameters relevant for measuring CP violation in D0 -> hh decays
        [
            combiner
            for combiner in _BASELINE_COMBINERS
            if not (combiner == "LHCb-R12" and _PARAMETERS_BY_NAME[name].cf_par in _FULL_FSC_PARAMS)
        ],
        name,
        _RANGES_1D["R12"][_PARAMETERS_BY_NAME[name].cf_par].transformed_range(_PARAMETERS_BY_NAME[name].transf),
        dy_fsc=[DYFsc.FULL],
    )
    for name in _PARAMETERS_BY_NAME
] + [
    Plot1D(
        ["LHCb-UI", "LHCb-UII"],
        name,
        _RANGES_1D["UI"][_PARAMETERS_BY_NAME[name].cf_par].transformed_range(_PARAMETERS_BY_NAME[name].transf),
        dy_fsc=[DYFsc.FULL],
    )
    for name in _PARAMETERS_BY_NAME
]

# 2D plots ------------------------------------------------------------------------------------------------------------

_PLOTS_2D = [
    Plot2D(
        # LHCb-R12 does not have enough precision for determining the strong parameters relevant for measuring CP violation in D0 -> hh decays
        [
            combiner
            for combiner in _BASELINE_COMBINERS
            if not (
                combiner == "LHCb-R12"
                and (
                    _PARAMETERS_BY_NAME[xpar].cf_par in _FULL_FSC_PARAMS
                    or _PARAMETERS_BY_NAME[ypar].cf_par in _FULL_FSC_PARAMS
                )
            )
        ],
        (xpar, ypar),
        xrange.transformed_range(_PARAMETERS_BY_NAME[xpar].transf),
        yrange.transformed_range(_PARAMETERS_BY_NAME[ypar].transf),
        dy_fsc=[DYFsc.FULL],
        legfill=True,
    )
    for (xpar, ypar), (xrange, yrange) in _RANGES_2D["R12"].items()
]

CONFIG = CharmPlottingConfig(_PARAMETERS, _BASELINE_COMBINERS, plots_1d=_PLOTS_1D, plots_2d=_PLOTS_2D)
