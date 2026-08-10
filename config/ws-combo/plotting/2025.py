from charm_fitter.utils import Combiner, Plot1D, Plot2D, PlotParameter, PlottingConfig, Range, Scan1D, Scan2D

# Parameters definitions -----------------------------------------------------------------------------------------------

_PARAMETERS = [
    PlotParameter("yp", title=r"$y^{\prime}$"),
    PlotParameter("dyp", title=r"$\Delta y^{\prime}$"),
    PlotParameter("xp2", title=r"$x^{\prime 2}$"),
    PlotParameter("dxp2", title=r"$\Delta x^{\prime 2}$"),
]

# Combiners definitions ------------------------------------------------------------------------------------------------

_RANGES = {
    "yp": Range((1.5, 11.0), unit="1e-3"),
    "dyp": Range((-10.0, 8.0), unit="1e-3"),
    "xp2": Range((-30.0, 32.0), unit="1e-5"),
    "dxp2": Range((-40.0, 57.0), unit="1e-5"),
}
_RANGES_WA = {
    "yp": Range((3.5, 7.0), unit="1e-3"),
    "dyp": Range((-2.0, 2.0), unit="1e-3"),
    "xp2": Range((-4.0, 9.0), unit="1e-5"),
    "dxp2": Range((-5.0, 5.0), unit="1e-5"),
}
_RANGES_LHCB = {
    "yp": Range((1.5, 9.5), unit="1e-3"),
    "xp2": Range((-28.0, 32.0), unit="1e-5"),
    "dyp": Range((-5.2, 2.5), unit="1e-3"),
    "dxp2": Range((-26.0, 36.0), unit="1e-5"),
}

_2D_COMBINATIONS = [("yp", "xp2"), ("dyp", "dxp2")]

_SCANS_1D_DEFAULT = [Scan1D(par, r) for par, r in _RANGES.items()]
_SCANS_2D_DEFAULT = [Scan2D(xpar, ypar, _RANGES[xpar], _RANGES[ypar]) for xpar, ypar in _2D_COMBINATIONS]

_SCANS_1D_WA = [Scan1D(par, r) for par, r in _RANGES_WA.items()]
_SCANS_2D_WA = [Scan2D(xpar, ypar, _RANGES_WA[xpar], _RANGES_WA[ypar]) for xpar, ypar in _2D_COMBINATIONS]

_BASELINE_COMBINERS = {
    "NonLHCb2025": Combiner(1, "CDF + BaBar + Belle", scans_1d=_SCANS_1D_DEFAULT, scans_2d=_SCANS_2D_DEFAULT),
    "LHCbPrompt2025": Combiner([23], "LHCb Run 1+2 prompt", scans_1d=_SCANS_1D_WA, scans_2d=_SCANS_2D_WA),
    "LHCbDT2025": Combiner(
        [25],
        "LHCb Run 1+2 double tag",
        scans_1d=_SCANS_1D_DEFAULT,
        scans_2d=[Scan2D(xpar, ypar, _RANGES_LHCB[xpar], _RANGES_LHCB[ypar]) for xpar, ypar in _2D_COMBINATIONS],
    ),
    "LHCb2025": Combiner(11, "LHCb Run 1+2", scans_1d=_SCANS_1D_WA, scans_2d=_SCANS_2D_WA),
    "WA2025": Combiner(20, "WA 2025", scans_1d=_SCANS_1D_WA, scans_2d=_SCANS_2D_WA),
}

# 1D plots -------------------------------------------------------------------------------------------------------------

_PLOTS_1D = [
    Plot1D(list(_BASELINE_COMBINERS.keys()), par.name, _RANGES[par.name], legpos="r", legfill=True)
    for par in _PARAMETERS
]

# 2D plots -------------------------------------------------------------------------------------------------------------

_PLOTS_2D = [
    Plot2D(list(_BASELINE_COMBINERS.keys()), (xpar, ypar), _RANGES[xpar], _RANGES[ypar], levels=2, legpos="bl")
    for xpar, ypar in _2D_COMBINATIONS
] + [
    Plot2D(
        ["LHCbDT2025", "LHCbPrompt2025", "LHCb2025"],
        (xpar, ypar),
        _RANGES_LHCB[xpar],
        _RANGES_LHCB[ypar],
        combiners_label="lhcb",
        levels=2,
        legpos="bl",
    )
    for xpar, ypar in _2D_COMBINATIONS
]

CONFIG = PlottingConfig(_PARAMETERS, _BASELINE_COMBINERS, plots_1d=_PLOTS_1D, plots_2d=_PLOTS_2D)
