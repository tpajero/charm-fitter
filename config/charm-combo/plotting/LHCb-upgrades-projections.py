"""Configuration file for LHCb Upgrade I and II sensitivity projections."""

import math

from charm_fitter.config import PARAMETERS
from charm_fitter.utils import CharmPlottingConfig, Combiner, Plot1D, Plot2D, Range, Scan1D, Scan2D

# Parameters and combiner ------------------------------------------------------------------------------------------

_PARAMETERS = PARAMETERS
_PARAMETERS_BY_NAME = {p.name: p for p in _PARAMETERS}

# Projected ranges, keyed by the underlying charm-fitter parameter name (i.e. Parameter.cf_par). Parameters that
# are also defined in charm_fitter.config use the same units as there.
_RANGES = {
    "x12": Range((3.5, 4.5), unit="1e-3"),
    "y12": Range((6.0, 6.5), unit="1e-3"),
    "phiM": Range((-10.0, 10.0), unit="mrad"),
    "phiG": Range((-10.0, 10.0), unit="mrad"),
    "Delta_Kpi": Range((-0.35, -0.1), unit="rad"),
    "r_Kpi": Range((5.84, 5.88), unit="%"),  # R_Kpi_ = r_Kpi^2 ~ (0.341, 0.346)%, matching the old R_Kpi range
    "F_pipipi0": Range((0.9, 1.0)),
    "Acp_KP": Range((-0.4, 0.4), unit="%"),
    "Acp_KK": Range((0.1, 1.1), unit="1e-3"),
    "Acp_PP": Range((1.7, 2.7), unit="1e-3"),
    "cot_delta_KK": Range((-15.0, 15.0)),
    "cot_delta_PP": Range((-15.0, 15.0)),
    "r_KK": Range((0.0, 100.0), unit="%"),
    "r_PP": Range((0.0, 100.0), unit="%"),
    "delta_KK": Range((-math.pi, 0.0), unit="rad"),
    "delta_PP": Range((-math.pi, 0.0), unit="rad"),
}

_BASELINE_COMBINER = "LHCb-UII"
_BASELINE_COMBINERS = {
    _BASELINE_COMBINER: Combiner(
        600,
        "LHCb UII",
        title_long="LHCb Upgrade II",
        scans_1d=[Scan1D(par, r) for par, r in _RANGES.items()],
        scans_2d=[
            Scan2D(xpar, ypar, _RANGES[xpar], _RANGES[ypar])
            for xpar, ypar in [
                ("x12", "y12"),
                ("phiM", "phiG"),
                ("Delta_Kpi", "r_Kpi"),
                ("Acp_KK", "Acp_PP"),
                ("cot_delta_KK", "cot_delta_PP"),
                ("r_KK", "r_PP"),
                ("delta_KK", "delta_PP"),
            ]
        ],
    ),
}

# 1D plots ----------------------------------------------------------------------------------------------------------

_PLOT_NAMES_1D = [
    "x12",
    "y12",
    "phiM",
    "phiG",
    "Delta_Kpi",
    "R_Kpi_",
    "F_pipipi0",
    "Acp_KP",
    "Acp_KK",
    "Acp_PP",
    "cot_delta_KK",
    "cot_delta_PP",
    "r_KK",
    "r_PP",
    "delta_KK",
    "delta_PP",
]
_PLOTS_1D = [
    Plot1D(
        [_BASELINE_COMBINER],
        name,
        _RANGES[_PARAMETERS_BY_NAME[name].cf_par].transformed_range(_PARAMETERS_BY_NAME[name].transf),
    )
    for name in _PLOT_NAMES_1D
]

# 2D plots ------------------------------------------------------------------------------------------------------------

_PLOTS_2D = [
    Plot2D([_BASELINE_COMBINER], ("x12", "y12"), _RANGES["x12"], _RANGES["y12"]),
    Plot2D([_BASELINE_COMBINER], ("phiM", "phiG"), _RANGES["phiM"], _RANGES["phiG"], logo="br", legpos="bl"),
    Plot2D([_BASELINE_COMBINER], ("Delta_Kpi", "r_Kpi"), _RANGES["Delta_Kpi"], _RANGES["r_Kpi"]),
    Plot2D([_BASELINE_COMBINER], ("Acp_KK", "Acp_PP"), _RANGES["Acp_KK"], _RANGES["Acp_PP"]),
    Plot2D([_BASELINE_COMBINER], ("cot_delta_KK", "cot_delta_PP"), _RANGES["cot_delta_KK"], _RANGES["cot_delta_PP"]),
    Plot2D([_BASELINE_COMBINER], ("r_KK", "r_PP"), _RANGES["r_KK"], _RANGES["r_PP"]),
    Plot2D([_BASELINE_COMBINER], ("delta_KK", "delta_PP"), _RANGES["delta_KK"], _RANGES["delta_PP"]),
]

CONFIG = CharmPlottingConfig(_PARAMETERS, _BASELINE_COMBINERS, plots_1d=_PLOTS_1D, plots_2d=_PLOTS_2D)
