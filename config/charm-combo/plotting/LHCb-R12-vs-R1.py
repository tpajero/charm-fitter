"""Configuration file for LHCb Run 1 vs. LHCb Run 1+2 comparison plots.

Combiners also include the latest determinations of hadronic parameters from charm factories.
"""

from charm_fitter.config import COMBINERS, PARAMETERS, RANGES_1D, RANGES_2D
from charm_fitter.utils import CharmPlottingConfig, Plot1D, Plot2D, Range

# Parameters and combiners ---------------------------------------------------------------------------------------------

_BASELINE_COMBINERS = {
    "LHCb-R1": COMBINERS["LHCb-R1"],
    "LHCb-R12": COMBINERS["LHCb-R12"],
}

_PARAMETERS_BY_NAME = {p.name: p for p in PARAMETERS}

# 1D plots -------------------------------------------------------------------------------------------------------------

_RIGHT_1D = ["x12", "phiM", "qop", "r_Kpi", "R_Kpi_"]
_PLOTS_1D = [
    Plot1D(
        list(_BASELINE_COMBINERS.keys()),
        name,
        RANGES_1D["LHCb-R1"][par.cf_par].transformed_range(_PARAMETERS_BY_NAME[name].transf),
        logo="r" if name in _RIGHT_1D else "l",
        legpos="r" if name in _RIGHT_1D else "l",
        legfill=True,
    )
    for name, par in _PARAMETERS_BY_NAME.items()
    if par.cf_par in RANGES_1D["LHCb-R1"]
]

# 2D plots -------------------------------------------------------------------------------------------------------------

_PLOTS_2D = [
    Plot2D(
        list(_BASELINE_COMBINERS.keys()),
        (xpar, ypar),
        xrange.transformed_range(_PARAMETERS_BY_NAME[xpar].transf),
        yrange if ypar != "Acp_PP" else Range((-5.5, 9.5), unit="1e-3"),
        legfill=xpar in ["Acp_KK", "phiM"],
    )
    for (xpar, ypar), (xrange, yrange) in RANGES_2D["LHCb-R1"].items()
]

CONFIG = CharmPlottingConfig(
    PARAMETERS, _BASELINE_COMBINERS, plots_1d=_PLOTS_1D, plots_2d=_PLOTS_2D, compare_dcs_hypos=False
)
