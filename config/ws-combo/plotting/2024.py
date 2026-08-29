from charm_fitter.utils import Combiner, Plot1D, Plot2D, PlotParameter, PlottingConfig, Range, Scan1D, Scan2D

# Parameters definitions -----------------------------------------------------------------------------------------------

_PARAMETERS = [
    PlotParameter("yp", title=r"$y^{\prime}$"),
    PlotParameter("dyp", title=r"$\Delta y^{\prime}$"),
    PlotParameter("xp2", title=r"$x^{\prime 2}$"),
    PlotParameter("dxp2", title=r"$\Delta x^{\prime 2}$"),
]

# Combiners definitions ------------------------------------------------------------------------------------------------

_RANGES = {  # For CDF + BaBar + Belle
    "yp": Range((1.0, 10.5), unit="1e-3"),
    "dyp": Range((-10.0, 7.0), unit="1e-3"),
    "xp2": Range((-34.0, 34.0), unit="1e-5"),
    "dxp2": Range((-40.0, 60.0), unit="1e-5"),
}
_RANGES_WA_1D = {  # For the world average and LHCb Run 1+2
    "yp": Range((3.8, 6.6), unit="1e-3"),
    "dyp": Range((-1.5, 1.5), unit="1e-3"),
    "xp2": Range((-5.0, 9.0), unit="1e-5"),
    "dxp2": Range((-8.0, 7.0), unit="1e-5"),
}
_RANGES_WA_2D = {  # For the world average and LHCb Run 1+2
    "yp": Range((4.2, 6.5), unit="1e-3"),
    "dyp": Range((-1.5, 1.5), unit="1e-3"),
    "xp2": Range((-3.0, 7.0), unit="1e-5"),
    "dxp2": Range((-8.0, 7.0), unit="1e-5"),
}
_RANGES_LHCB_R1 = {  # For LHCb Run 1 2D plots
    "yp": Range((2.5, 8.0), unit="1e-3"),
    "dyp": Range((-2.5, 3.0), unit="1e-3"),
    "xp2": Range((-14.0, 15.0), unit="1e-5"),
    "dxp2": Range((-15.0, 15.0), unit="1e-5"),
}
_RANGES_LHCB_DT = {  # For LHCb Run 1+2 double tag 2D plots
    "yp": Range((1.0, 9.5), unit="1e-3"),
    "dyp": Range((-5.2, 2.5), unit="1e-3"),
    "xp2": Range((-28.0, 32.0), unit="1e-5"),
    "dxp2": Range((-26.0, 36.0), unit="1e-5"),
}

_2D_COMBINATIONS = [("yp", "xp2"), ("dyp", "dxp2")]

_SCANS_1D_DEFAULT = [Scan1D(par, r) for par, r in _RANGES.items()]
_SCANS_2D_DEFAULT = [Scan2D(xpar, ypar, _RANGES[xpar], _RANGES[ypar]) for xpar, ypar in _2D_COMBINATIONS]

_SCANS_1D_WA = [Scan1D(par, r) for par, r in _RANGES_WA_1D.items()]
_SCANS_2D_WA = [Scan2D(xpar, ypar, _RANGES_WA_2D[xpar], _RANGES_WA_2D[ypar]) for xpar, ypar in _2D_COMBINATIONS]

_BASELINE_COMBINERS = {
    "NonLHCb-2013": Combiner(1, "CDF + BaBar + Belle", scans_1d=_SCANS_1D_DEFAULT, scans_2d=_SCANS_2D_DEFAULT),
    "NonLHCb-2014-ForCPV": Combiner(2, "CDF + BaBar + Belle", scans_1d=_SCANS_1D_DEFAULT, scans_2d=_SCANS_2D_DEFAULT),
    "LHCb-R12-prompt": Combiner([35], "LHCb Run 1+2 prompt", scans_1d=_SCANS_1D_WA, scans_2d=_SCANS_2D_WA),
    "LHCb-R12-DT": Combiner(
        [37],
        "LHCb Run 1+2 double tag",
        scans_1d=_SCANS_1D_DEFAULT,
        scans_2d=[Scan2D(xpar, ypar, _RANGES_LHCB_DT[xpar], _RANGES_LHCB_DT[ypar]) for xpar, ypar in _2D_COMBINATIONS],
    ),
    "LHCb-R1": Combiner(
        200,
        "LHCb Run 1",
        scans_2d=[Scan2D(xpar, ypar, _RANGES_LHCB_R1[xpar], _RANGES_LHCB_R1[ypar]) for xpar, ypar in _2D_COMBINATIONS],
    ),
    "LHCb-R12": Combiner(201, "LHCb Run 1+2", scans_2d=_SCANS_2D_WA),
    "WA-2024": Combiner(300, "WA 2024", scans_1d=_SCANS_1D_WA, scans_2d=_SCANS_2D_WA),
}

_MIX_COMBINERS = ["NonLHCb-2013", "LHCb-R12-DT", "LHCb-R12-prompt", "WA-2024"]
_CPV_COMBINERS = ["NonLHCb-2014-ForCPV", "LHCb-R12-DT", "LHCb-R12-prompt", "WA-2024"]

_BASELINE_LABEL = "WA-2024"

# 1D plots -------------------------------------------------------------------------------------------------------------

_PLOTS_1D = [
    Plot1D(
        _MIX_COMBINERS if par.name in ["yp", "xp2"] else _CPV_COMBINERS,
        par.name,
        _RANGES[par.name],
        combiners_label=_BASELINE_LABEL,
        legpos="l",
        legfontsize=13,
        legfill=True,
    )
    for par in _PARAMETERS
]

# 2D plots -------------------------------------------------------------------------------------------------------------

_PLOTS_2D = (
    [
        Plot2D(
            combiners,
            (xpar, ypar),
            _RANGES[xpar],
            _RANGES[ypar],
            levels=2,
            combiners_label=_BASELINE_LABEL,
            legpos="bl",
        )
        for (xpar, ypar), combiners in [(("yp", "xp2"), _MIX_COMBINERS), (("dyp", "dxp2"), _CPV_COMBINERS)]
    ]
    + [
        Plot2D(
            ["LHCb-R12-DT", "LHCb-R12-prompt", "LHCb-R12"],
            (xpar, ypar),
            _RANGES_LHCB_DT[xpar],
            _RANGES_LHCB_DT[ypar],
            combiners_label="LHCb-prompt-vs-DT",
            levels=2,
            legpos="bl",
        )
        for xpar, ypar in _2D_COMBINATIONS
    ]
    + [
        Plot2D(
            ["LHCb-R1", "LHCb-R12"],
            (xpar, ypar),
            _RANGES_LHCB_R1[xpar],
            _RANGES_LHCB_R1[ypar],
            combiners_label="LHCb-R12-vs-R1",
            levels=2,
            legpos="bl",
        )
        for xpar, ypar in _2D_COMBINATIONS
    ]
)

CONFIG = PlottingConfig(_PARAMETERS, _BASELINE_COMBINERS, plots_1d=_PLOTS_1D, plots_2d=_PLOTS_2D)
