from charm_fitter.utils import (
    AcpParam,
    Combiner,
    DYFsc,
    MixParam,
    Parameter,
    Plot2D,
    PlotOpts1D,
    ScanParams2D,
)

# Parameters definitions -----------------------------------------------------------------------------------------------

_delta_kpi_range = (-0.5, 0.1)

parameters = {
    # --- Mixing in the phenomenological parametrisation ---
    "x": Parameter(
        "x",
        (1.5, 7.0),
        title=r"$x$",
        unit="1e-3",
        mix_params=[MixParam.PHENO],
        plot_opts_1d=PlotOpts1D(legpos="r"),
    ),
    "y": Parameter(
        "y",
        (5.0, 8.0),
        title=r"$y$",
        unit="1e-3",
        mix_params=[MixParam.PHENO],
        plot_opts_1d=PlotOpts1D(legpos="r"),
    ),
    "qop": Parameter(
        "qop",
        (0.92, 1.08),
        plot_range=(-0.08, 0.08),
        scan_range_2d=(0.90, 1.07),
        plot_range_2d=(-0.1, 0.07),
        title=r"$|q/p|-1$",
        transf=lambda x: x - 1.0,
        mix_params=[MixParam.PHENO],
        plot_opts_1d=PlotOpts1D(legpos="r"),
    ),
    "phi": Parameter(
        "phi",
        (-120.0, 70.0),
        title=r"$\phi_2$",
        unit="mrad",
        mix_params=[MixParam.PHENO],
        plot_opts_1d=PlotOpts1D(legpos="r", legfill=True),
    ),
    # --- Mixing in the theoretical parametrisation ---
    "x12": Parameter(
        "x12",
        (1.5, 7.0),
        title=r"$x_{12}$",
        unit="1e-3",
        mix_params=[MixParam.THEO],
        plot_opts_1d=PlotOpts1D(legpos="r"),
    ),
    "y12": Parameter(
        "y12",
        (5.0, 8.0),
        title=r"$y_{12}$",
        unit="1e-3",
        mix_params=[MixParam.THEO],
        plot_opts_1d=PlotOpts1D(legpos="r"),
    ),
    "phiM": Parameter(
        "phiM",
        (-120.0, 120.0),
        scan_range_2d=(-80.0, 110.0),
        title=r"$\phi_{2}^{M}$",
        unit="mrad",
        mix_params=[MixParam.THEO],
        plot_opts_1d=PlotOpts1D(legfill=True),
    ),
    "phiG": Parameter(
        "phiG",
        (-150.0, 150.0),
        scan_range_2d=(-100.0, 200.0),
        title=r"$\phi_{2}^{\Gamma}$",
        unit="mrad",
        mix_params=[MixParam.THEO],
    ),
    # --- CPV in the decay ---
    "Acp_KK": Parameter(
        "Acp_KK",
        (-1.5, 3.0),
        plot_range_2d=(-3.0, 3.0),
        title=r"$a^{\rm d}_{K^+K^-}$",
        unit="1e-3",
        plot_opts_1d=PlotOpts1D(legfill=True),
        acp_params=[AcpParam.ACP_DY, AcpParam.ACP_COT],
    ),
    "Acp_PP": Parameter(
        "Acp_PP",
        (0.0, 5.5),
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
    "Delta_Kpipi0": Parameter("Delta_Kpipi0", (-1.5, 0.5), title=r"$\Delta_{K\pi\pi^0}$", unit="rad"),
    "Delta_K3pi": Parameter(
        "Delta_K3pi", (-0.5, 1.5), title=r"$\Delta_{K\pi\pi\pi}$", unit="rad", plot_opts_1d=PlotOpts1D(legpos="r")
    ),
    "r_Kpi": Parameter(
        "r_Kpi",
        (5.79, 5.95),
        scan_range_2d=(5.81, 5.92),
        title=r"$r_{K\pi}$",
        unit="%",
        plot_opts_1d=PlotOpts1D(legpos="r"),
    ),
    "r_K3pi": Parameter("r_K3pi", (5.1, 6.0), title=r"$r_{K\pi\pi\pi}$", unit="%", plot_opts_1d=PlotOpts1D(legpos="r")),
    "r_K3pipi0": Parameter(
        "r_Kpipi0", (3.7, 5.5), title=r"$r_{K\pi\pi^0}$", unit="%", plot_opts_1d=PlotOpts1D(legpos="r")
    ),
    "k_K3pi": Parameter("k_K3pi", (0.2, 0.75), title=r"$k_{K\pi\pi\pi}$", plot_opts_1d=PlotOpts1D(legpos="r")),
    "k_Kpipi0": Parameter("k_Kpipi0", (0.4, 1.0), title=r"$k_{K\pi\pi^0}$"),
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

baseline_combiners = ["WA-2025-10"]
baseline_combiner = baseline_combiners[0]
baseline_parfile = "config/start/charm-2025.dat"

combiners = {
    "WA-2025-10": Combiner(55, "WA Oct 2025", col="r", ls="-", title_long="World Average October 2025"),
    # Breakdowns by D0 final state and experiment
    "ext-kpi": Combiner(
        [50, 53, 55, 56],
        r"CLEO/BES $D^0 \to K^\pm \pi^\mp$",
        col="y",
        ls="dash3dots",
        scanparams_2d=[
            ScanParams2D(
                ("Delta_Kpi", "r_Kpi"), xrange=(-0.6, 0.3), yrange=(5.55e-2, 6.15e-2), mix_param=MixParam.PHENO
            )
        ],
    ),
    "lhcb-kpi": Combiner(
        [39, 41],
        r"LHCb $D^0 \to K^\pm \pi^\mp$",
        col="g",
        ls="dash2dots",
        scanparams_2d=[
            ScanParams2D(
                ("Delta_Kpi", "r_Kpi"),
                xrange=(-1.35, 0.55),
                yrange=(5.81e-2, 5.92e-2),
                parfile="config/start/charm-2025-breakdown-kpi-rd.dat",
                extra_opts="--npoints2dx 200",
            ),
        ],
    ),
    "kpi": Combiner(
        [39, 41, 50, 53, 55, 56],
        r"$D^0 \to K^\pm \pi^\mp$",
        col="g",
        ls="dash2dots",
        scanparams_2d=[
            ScanParams2D(
                ("x", "y"),
                xrange=(-8e-3, 8e-3),
                yrange=(3e-3, 8e-3),
                fixed_pars=["qop", "phi"],
                parfile="config/start/charm-2025-breakdown-kpi-x-y.dat",
                extra_opts="--npoints2dx 100 --npoints2dy 100",
            ),
            ScanParams2D(
                ("x12", "y12"),
                xrange=(0.0, 8e-3),
                yrange=(3e-3, 8e-3),
                fixed_pars=["phiM", "phiG"],
                parfile="config/start/charm-2025-breakdown-kpi-x12-y12.dat",
                extra_opts="--npoints2dx 100 --npoints2dy 100",
            ),
            ScanParams2D(("qop", "phi"), xrange=(0.7, 1.4), yrange=(-0.4, 0.35), fixed_pars=["x", "y"]),
            ScanParams2D(("phiM", "phiG"), xrange=(-0.1, 0.1), yrange=(-0.41, 0.3), fixed_pars=["x12", "y12"]),
        ],
    ),
    "k3pi": Combiner(
        [5, 54],
        r"$D^0 \to K^\pm \pi^\mp \pi^+ \pi^-$",
        col="p",
        ls="dash3dots",
        scanparams_2d=[
            ScanParams2D(("x", "y"), xrange=(1e-3, 6e-3), yrange=(2e-3, 10e-3), fixed_pars=["qop", "phi"]),
            ScanParams2D(("x12", "y12"), xrange=(1e-3, 6e-3), yrange=(2e-3, 10e-3), fixed_pars=["phiM", "phiG"]),
            ScanParams2D(("qop", "phi"), xrange=(0.85, 1.2), yrange=(-0.25, 0.1), fixed_pars=["x", "y"]),
            ScanParams2D(("phiM", "phiG"), xrange=(-0.2, 0.35), yrange=(-0.1, 0.25), fixed_pars=["x12", "y12"]),
        ],
    ),
    "kspipi": Combiner(
        [1, 3, 6, 21, 24],
        r"$D^0 \to K_S^0 \pi^+ \pi^-$",
        col="lb",
        ls="shortdash",
        scanparams_2d=[
            ScanParams2D(("x", "y"), xrange=(2e-3, 6e-3), yrange=(2e-3, 7e-3), fixed_pars=["qop", "phi"]),
            ScanParams2D(("x12", "y12"), xrange=(2e-3, 6e-3), yrange=(2e-3, 10e-3), fixed_pars=["phiM", "phiG"]),
            ScanParams2D(("qop", "phi"), xrange=(0.85, 1.2), yrange=(-0.25, 0.1), fixed_pars=["x", "y"]),
            ScanParams2D(("phiM", "phiG"), xrange=(-0.2, 0.35), yrange=(-0.05, 0.15), fixed_pars=["x12", "y12"]),
        ],
    ),
    "hh-mix": Combiner(
        [61, 62, 64],  # yCP only
        r"$D^0 \to h^+ h^-$",
        col="y",
        ls="dashdot",
        scanparams_2d=[
            ScanParams2D(
                ("x", "y"),
                xrange=(-8e-3, 12.5e-3),
                yrange=(5e-3, 8e-3),
                fixed_pars=["r_Kpi", "Delta_Kpi", "qop", "phi"],
                extra_opts="--npoints2dx 100",
            ),
            ScanParams2D(
                ("x12", "y12"),
                xrange=(0.0, 10.5e-3),
                yrange=(5e-3, 8e-3),
                fixed_pars=["r_Kpi", "Delta_Kpi", "phiM", "phiG"],
                extra_opts="--npoints2dx 100",
            ),
        ],
    ),
    "hh": Combiner(
        [61, 62, 64, 72],  # yCP + DeltaY
        r"$D^0 \to h^+ h^-$",
        col="y",
        ls="dashdot",
        scanparams_2d=[
            ScanParams2D(
                ("qop", "phi"), xrange=(0.8, 1.3), yrange=(-0.4, 0.4), fixed_pars=["x", "y", "r_Kpi", "Delta_Kpi"]
            ),
            ScanParams2D(
                ("phiM", "phiG"),
                xrange=(-0.1, 0.15),
                yrange=(-0.5, 0.5),
                fixed_pars=["x12", "y12", "r_Kpi", "Delta_Kpi"],
            ),
        ],
    ),
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
_breakdown_cpv_combos = ["kspipi", "hh", "kpi", baseline_combiner]
_breakdown_mix_combos = ["kpi", "kspipi", "hh-mix", baseline_combiner]
breakdowns = [
    Plot2D(("x", "y"), combiners=_breakdown_mix_combos, xrange=(-7e-3, 12e-3), yrange=(2e-3, 10e-3)),
    Plot2D(("x12", "y12"), combiners=_breakdown_mix_combos, xrange=(0.0, 10e-3), yrange=(2e-3, 10e-3)),
    Plot2D(("qop", "phi"), combiners=_breakdown_cpv_combos, xrange=(0.8, 1.3), yrange=(-0.4, 0.35), logo="br"),
    Plot2D(("phiM", "phiG"), combiners=_breakdown_cpv_combos, xrange=(-0.4, 0.35), yrange=(-0.4, 0.4)),
    Plot2D(
        ("Delta_Kpi", "r_Kpi"),
        combiners=["ext-kpi", "lhcb-kpi", baseline_combiner],
        xrange=(-1.3, 0.5),
        yrange=(5.6e-2, 6.1e-2),
        logo="l",
        legpos="bl",
        legfill=True,
    ),
]
