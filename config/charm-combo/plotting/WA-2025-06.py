"""Configuration file for summary plots of the World Average of June 2025."""

from charm_fitter.config import (
    BREAKDOWN_COMBINERS_KPI_RD,
    CPV_PAIRS,
    MIX_PAIRS,
    PARAMETERS,
    RANGES_1D,
    RANGES_2D,
    breakdown_combiner,
    breakdown_plot_kpi_rd,
    dy_fsc_plots_1d,
    dy_fsc_plots_2d,
    wa_combiner,
    wa_plots_1d,
    wa_plots_2d,
)
from charm_fitter.utils import CharmPlottingConfig, Plot2D, Range

# Parameters and combiners ---------------------------------------------------------------------------------------------

_BASELINE_COMBINER = "WA-2025-06"
_FIX_PARFILE = "config/charm-combo/start/WA-2025-06.dat"

_RANGES_1D = RANGES_1D["LHCb-R12"]
_RANGES_2D = RANGES_2D["LHCb-R12"]

_BASELINE_COMBINERS = {
    _BASELINE_COMBINER: wa_combiner(12, "World average", _FIX_PARFILE, _RANGES_1D, _RANGES_2D),
}

_COMBINERS_BREAKDOWN = {
    **BREAKDOWN_COMBINERS_KPI_RD,
    "kpi-mix": breakdown_combiner("kpi", [1, 21, 40, 41, 250, 252, 260, 286, 288], _FIX_PARFILE, MIX_PAIRS),
    "kpi-cpv": breakdown_combiner("kpi", [1, 21, 40, 41, 260, 261, 286, 288], _FIX_PARFILE, CPV_PAIRS),
    "kspipi": breakdown_combiner("kspipi", [320, 322, 323, 326, 327], _FIX_PARFILE),
    "hh-mix": breakdown_combiner("hh", [110, 130, 162], _FIX_PARFILE, MIX_PAIRS),  # yCP only
    "hh-cpv": breakdown_combiner("hh", [110, 130, 162, 202], _FIX_PARFILE, CPV_PAIRS),  # yCP + DeltaY
}

# 1D and 2D plots ------------------------------------------------------------------------------------------------------

_PLOTS_1D = wa_plots_1d(_BASELINE_COMBINER, _RANGES_1D)
_PLOTS_2D = wa_plots_2d(
    _BASELINE_COMBINER,
    _RANGES_2D | {("Acp_KK", "Acp_PP"): (_RANGES_2D[("Acp_KK", "Acp_PP")][0], Range((-0.5, 5.5), unit="1e-3"))},
)

# Breakdown 2D plots ---------------------------------------------------------------------------------------------------

# N.B.: This is tuned only for the hypothesis of no CPV in D0 -> K+ pi- decays.
_BREAKDOWN_CPV_COMBOS = ["kspipi", "hh-cpv", "kpi-cpv", _BASELINE_COMBINER]
_BREAKDOWN_MIX_COMBOS = ["kpi-mix", "kspipi", "hh-mix", _BASELINE_COMBINER]
_PLOTS_BREAKDOWN = [
    Plot2D(
        _BREAKDOWN_MIX_COMBOS,
        ("x", "y"),
        xrange=Range((-8.0, 10.0), unit="1e-3"),
        yrange=Range((2.0, 9.0), unit="1e-3"),
        combiners_label=_BASELINE_COMBINER,
        levels=2,
        logo="bl",
        legncol=2,
    ),
    Plot2D(
        _BREAKDOWN_MIX_COMBOS,
        ("x12", "y12"),
        xrange=Range((0.0, 10.0), unit="1e-3"),
        yrange=Range((2.0, 9.0), unit="1e-3"),
        combiners_label=_BASELINE_COMBINER,
        levels=2,
        logo="br",
        legncol=2,
    ),
    Plot2D(
        _BREAKDOWN_CPV_COMBOS,
        ("qop", "phi"),
        xrange=Range((-0.25, 0.3)),
        yrange=Range((-350.0, 350.0), unit="mrad"),
        logo="br",
        combiners_label=_BASELINE_COMBINER,
        levels=2,
    ),
    Plot2D(
        _BREAKDOWN_CPV_COMBOS,
        ("phiM", "phiG"),
        xrange=Range((-320.0, 350.0), unit="mrad"),
        yrange=Range((-400.0, 450.0), unit="mrad"),
        combiners_label=_BASELINE_COMBINER,
        levels=2,
    ),
    breakdown_plot_kpi_rd(_BASELINE_COMBINER),
]

# Plots for comparing the impact of different DeltaY(h- h+) final-state correction hypotheses --------------------------

_PLOTS_DY_FSC_1D = dy_fsc_plots_1d(_BASELINE_COMBINER, _RANGES_1D)

# N.B.: This is tuned only for the hypothesis of no CPV in D0 -> K+ pi- decays.
_PLOTS_DY_FSC_2D = dy_fsc_plots_2d(
    _BASELINE_COMBINER,
    {
        ("phiM", "phiG"): (Range((-80.0, 80.0), unit="mrad"), Range((-100.0, 220.0), unit="mrad")),
        ("qop", "phi"): (Range((-0.1, 0.05)), Range((-120.0, 60.0), unit="mrad")),
        ("Acp_KK", "Acp_PP"): (Range((-1.5, 3.0), unit="1e-3"), Range((-0.5, 6.5), unit="1e-3")),
        ("Acp_KK", "phiM"): (Range((-2.0, 3.5), unit="1e-3"), Range((-80.0, 70.0), unit="mrad")),
    },
)

CONFIG = CharmPlottingConfig(
    PARAMETERS,
    _BASELINE_COMBINERS,
    combiners_breakdown=_COMBINERS_BREAKDOWN,
    plots_1d=_PLOTS_1D,
    plots_2d=_PLOTS_2D,
    plots_breakdown=_PLOTS_BREAKDOWN,
    plots_dy_fsc_1d=_PLOTS_DY_FSC_1D,
    plots_dy_fsc_2d=_PLOTS_DY_FSC_2D,
    logo_date="June 2025",
)
