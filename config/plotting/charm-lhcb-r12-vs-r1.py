from charm_fitter.utils import (
    AcpParam,
    Combiner,
    DYFsc,
    MixParam,
    Parameter,
    Plot2D,
    PlotOpts1D,
)

# Parameters definitions -----------------------------------------------------------------------------------------------

_delta_kpi_range = (-0.7, 0.3)

parameters = {
    # --- Mixing in the phenomenological parametrisation ---
    "x": Parameter(
        "x",
        (-2.0, 8.0),
        title=r"$x$",
        unit="1e-3",
        mix_params=[MixParam.PHENO],
        plot_opts_1d=PlotOpts1D(legpos="r"),
    ),
    "y": Parameter(
        "y",
        (4.0, 9.0),
        title=r"$y$",
        unit="1e-3",
        mix_params=[MixParam.PHENO],
        plot_opts_1d=PlotOpts1D(legpos="r"),
    ),
    "qop": Parameter(
        "qop",
        (0.8, 1.2),
        plot_range=(-0.2, 0.2),
        scan_range_2d=(0.90, 1.07),
        plot_range_2d=(-0.1, 0.07),
        title=r"$|q/p|-1$",
        transf=lambda x: x - 1.0,
        mix_params=[MixParam.PHENO],
        plot_opts_1d=PlotOpts1D(legpos="r"),
    ),
    "phi": Parameter(
        "phi",
        (-500.0, 200.0),
        title=r"$\phi_2$",
        unit="mrad",
        mix_params=[MixParam.PHENO],
        plot_opts_1d=PlotOpts1D(legpos="r", legfill=True),
    ),
    # --- Mixing in the theoretical parametrisation ---
    "x12": Parameter(
        "x12",
        (0.0, 8.0),
        title=r"$x_{12}$",
        unit="1e-3",
        mix_params=[MixParam.THEO],
        plot_opts_1d=PlotOpts1D(legpos="r"),
    ),
    "y12": Parameter(
        "y12",
        (4.0, 9.0),
        title=r"$y_{12}$",
        unit="1e-3",
        mix_params=[MixParam.THEO],
        plot_opts_1d=PlotOpts1D(legpos="r"),
    ),
    "phiM": Parameter(
        "phiM",
        (-800.0, 800.0),
        scan_range_2d=(-80.0, 110.0),
        title=r"$\phi_{2}^{M}$",
        unit="mrad",
        mix_params=[MixParam.THEO],
        plot_opts_1d=PlotOpts1D(legfill=True),
    ),
    "phiG": Parameter(
        "phiG",
        (-300.0, 500.0),
        scan_range_2d=(-100.0, 200.0),
        title=r"$\phi_{2}^{\Gamma}$",
        unit="mrad",
        mix_params=[MixParam.THEO],
    ),
    # --- CPV in the decay ---
    "Acp_KK": Parameter(
        "Acp_KK",
        (-5.0, 5.0),
        plot_range_2d=(-3.0, 3.0),
        title=r"$a^{\rm d}_{K^+K^-}$",
        unit="1e-3",
        plot_opts_1d=PlotOpts1D(legfill=True),
        acp_params=[AcpParam.ACP_DY, AcpParam.ACP_COT],
    ),
    "Acp_PP": Parameter(
        "Acp_PP",
        (-5.0, 6.0),
        plot_range_2d=(-0.5, 6.0),
        title=r"$a^{\rm d}_{\pi^+\pi^-}$",
        unit="1e-3",
        plot_opts_1d=PlotOpts1D(legpos="r"),
        acp_params=[AcpParam.ACP_DY, AcpParam.ACP_COT],
    ),
    "Acp_KP": Parameter("Acp_KP", (-2.0, 2.0), title=r"$a^{\rm d}_{K^+\pi^-}$", unit="%"),
    # --- Hadronic nuisance parameters ---
    "Delta_Kpi": Parameter(
        "Delta_Kpi",
        _delta_kpi_range,
        title=r"$\Delta_{K\pi}$",
        unit="rad",
        plot_opts_1d=PlotOpts1D(legpos="r"),
    ),
    # "Delta_K3pi": Parameter(
    #     "Delta_K3pi", (-0.5, 1.5), title=r"$\Delta_{K\pi\pi\pi}$", unit="rad", plot_opts_1d=PlotOpts1D(legpos="r")
    # ),
    "r_Kpi": Parameter(
        "r_Kpi",
        (5.5, 6.7),
        scan_range_2d=(5.81, 5.92),
        title=r"$r_{K\pi}$",
        unit="%",
        plot_opts_1d=PlotOpts1D(legpos="r"),
    ),
    "r_K3pi": Parameter("r_K3pi", (5.1, 6.0), title=r"$r_{K\pi\pi\pi}$", unit="%", plot_opts_1d=PlotOpts1D(legpos="r")),
    "k_K3pi": Parameter("k_K3pi", (0.2, 0.75), title=r"$k_{K\pi\pi\pi}$", plot_opts_1d=PlotOpts1D(legpos="r")),
    "F_pipipi0": Parameter("F_pipipi0", (0.91, 0.96), title=r"$F_{\pi\pi\pi^0}$"),
    # --- Reparametrisations ---
    # NB: Avoid redefining existing parameters with different case (e.g. R_Kpi vs r_Kpi), since this will overwrite
    #     output files in case-insensitive file systems.
    "RKpi": Parameter(
        "r_Kpi",
        (3.32, 3.5),
        plot_range_2d=(3.38, 3.5),
        title=r"$R_{K\pi}$",
        unit="1e-3",
        transf=lambda x: x * x,
        scan_1d=False,
    ),
    "deltaKpi": Parameter(
        "Delta_Kpi",
        (-_delta_kpi_range[1], -_delta_kpi_range[0]),
        title=r"$\delta_{K\pi}$",
        unit="rad",
        transf=lambda x: -x,
        scan_1d=False,
    ),
}

# Combiners definitions ------------------------------------------------------------------------------------------------

baseline_combiners = ["LHCb-R1", "LHCb-R12"]
baseline_combiner = baseline_combiners[0]
baseline_parfile = "config/start/charm-2025.dat"

combiners = {
    "LHCb-R1": Combiner(300, "LHCb Run 1", parfile="config/start/charm-lhcb-r1.dat"),
    "LHCb-R12": Combiner(301, "LHCb Run 1+2", parfile="config/start/charm-lhcb-r12.dat"),
}

# Style for 1D and 2D plots --------------------------------------------------------------------------------------------

colors = ["b", "r", "y", "o", "g", "p", "lb"]
ls = ["longdash", "-", "shortdash", "dashdot", "dash3dots", "dash2dots", "dots"]

colors_2d = ["lb", "r", "y", "o", "g", "p", "b"]
ls_2d = ["-", "longdash", "shortdash", "dashdot", "dash3dots", "dash2dots", "dots"]
markers = ["o", "x", "s"]

# Parameters and options for 2D plots ----------------------------------------------------------------------------------

# List all pairs of parameters for which 2D scans should be prepared.
plots_2d = [
    Plot2D(("x12", "y12")),
    Plot2D(("phiM", "phiG")),
    Plot2D(("x", "y")),
    Plot2D(("qop", "phi")),
    Plot2D(("Delta_Kpi", "r_Kpi")),
    Plot2D(("Delta_Kpi", "RKpi"), scan=False),
    Plot2D(("deltaKpi", "r_Kpi"), scan=False),
    Plot2D(("Acp_KK", "Acp_PP"), legfill=True),
    Plot2D(("Acp_KK", "phiM")),
]

# Parameters and options for comparisons of different DeltaY(h- h+) final-state correction hypotheses
dy_fsc_baseline = DYFsc.NONE
dy_fsc_labels = {
    DYFsc.NONE: r"$\Delta Y_{f} = - x_{12}\sin\phi^M_2$",
    DYFsc.PARTIAL: r"$\Delta Y_{f} = - x_{12}\sin\phi^M_2 + y_{12} a^{\rm d}_f$",
    DYFsc.FULL: r"$\Delta Y_{f} = - x_{12}\sin\phi^M_2 + y_{12} a^{\rm d}_f (1 + \cot\delta_f \; x_{12} / y_{12})$",
}
dy_fsc_comparison_plots_2d = [
    Plot2D(("phiM", "phiG")),
    Plot2D(("qop", "phi")),
    Plot2D(("Acp_KK", "Acp_PP"), legfill=True),
    Plot2D(("Acp_KK", "phiM"), xrange=(-2e-3, 4e-3), yrange=(-0.1, 0.1)),
]

# Breakdown plots
breakdowns = []
