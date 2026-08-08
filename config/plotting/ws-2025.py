from charm_fitter.utils import (
    Combiner,
    MixParam,
    Parameter,
    Plot2D,
    PlotOpts1D,
    ScanParams1D,
    ScanParams2D,
)

# Style for 1D and 2D plots --------------------------------------------------------------------------------------------

colors = ["b", "r", "g", "y", "p", "lb"]
ls = ["-", "longdash", "shortdash", "dashdot", "dash3dots", "dash2dots", "dots"]

colors_2d = ["lb", "r", "g", "y", "p"]
ls_2d = ["-", "longdash", "shortdash", "dashdot", "dash3dots", "dash2dots", "dots"]
markers = ["o", "x", "s", "^", "D", "v", "P"]

# Parameters definitions -----------------------------------------------------------------------------------------------

parameters = {
    "yp": Parameter(
        "yp",
        (1.5, 11.0),
        title=r"$y^{\prime}$",
        unit="1e-3",
        mix_params=[MixParam.D0_TO_KPI],
        plot_opts_1d=PlotOpts1D(legpos="r", legfill=True),
    ),
    "dyp": Parameter(
        "dyp",
        (-10.0, 8.0),
        title=r"$\Delta y^{\prime}$",
        unit="1e-3",
        mix_params=[MixParam.D0_TO_KPI],
        plot_opts_1d=PlotOpts1D(legpos="r", legfill=True),
    ),
    "xp2": Parameter(
        "xp2",
        (-30.0, 32.0),
        title=r"$x^{\prime 2}$",
        unit="1e-5",
        mix_params=[MixParam.D0_TO_KPI],
        plot_opts_1d=PlotOpts1D(legpos="r", legfill=True),
    ),
    "dxp2": Parameter(
        "dxp2",
        (-40.0, 57.0),
        title=r"$\Delta x^{\prime 2}$",
        unit="1e-5",
        mix_params=[MixParam.D0_TO_KPI],
        plot_opts_1d=PlotOpts1D(legpos="r", legfill=True),
    ),
}

# Combiners definitions ------------------------------------------------------------------------------------------------

# Reduced ranges for some of the scans
_scanparams_1d_wa = [
    ScanParams1D(par="yp", range=(3.5e-3, 7e-3)),
    ScanParams1D(par="dyp", range=(-2e-3, 2e-3)),
    ScanParams1D(par="xp2", range=(-4e-5, 9e-5)),
    ScanParams1D(par="dxp2", range=(-5e-5, 5e-5)),
]
_scanparams_2d_wa = [
    ScanParams2D(("yp", "xp2"), xrange=(4e-3, 7e-3), yrange=(-4e-5, 9e-5)),
    ScanParams2D(("dyp", "dxp2"), xrange=(-1e-3, 1.5e-3), yrange=(-7e-5, 6e-5)),
]
_lhcb_yp_range = (1.5e-3, 9.5e-3)
_lhcb_xp2_range = (-28e-5, 32e-5)
_lhcb_dyp_range = (-5.2e-3, 2.5e-3)
_lhcb_dxp2_range = (-26e-5, 36e-5)

combiners = {
    "NonLHCb2025": Combiner(1, "CDF + BaBar + Belle"),
    "LHCbPrompt2025": Combiner(
        [23], "LHCb Run 1+2 prompt", scanparams_1d=_scanparams_1d_wa, scanparams_2d=_scanparams_2d_wa
    ),
    "LHCbDT2025": Combiner(
        [25],
        "LHCb Run 1+2 double tag",
        scanparams_2d=[
            ScanParams2D(("yp", "xp2"), xrange=_lhcb_yp_range, yrange=_lhcb_xp2_range),
            ScanParams2D(("dyp", "dxp2"), xrange=_lhcb_dyp_range, yrange=_lhcb_dxp2_range),
        ],
    ),
    "LHCb2025": Combiner(11, "LHCb Run 1+2", scanparams_1d=_scanparams_1d_wa, scanparams_2d=_scanparams_2d_wa),
    "WA2025": Combiner(20, "WA 2025", scanparams_1d=_scanparams_1d_wa, scanparams_2d=_scanparams_2d_wa),
}

baseline_combiners = list(combiners.keys())

# Parameters and options for 2D plots ----------------------------------------------------------------------------------

_lhcb_combiners = ["LHCbDT2025", "LHCbPrompt2025", "LHCb2025"]

# List all 2D plots to be produced.
plots_2d = [
    Plot2D(("yp", "xp2"), levels=2, legpos="bl"),
    Plot2D(("dyp", "dxp2"), levels=2, legpos="bl"),
    Plot2D(
        ("yp", "xp2"),
        label="lhcb",
        combiners=_lhcb_combiners,
        xrange=_lhcb_yp_range,
        yrange=_lhcb_xp2_range,
        levels=2,
        legpos="bl",
        scan=False,
    ),
    Plot2D(
        ("dyp", "dxp2"),
        label="lhcb",
        combiners=_lhcb_combiners,
        xrange=_lhcb_dyp_range,
        yrange=_lhcb_dxp2_range,
        levels=2,
        legpos="bl",
        scan=False,
    ),
]
