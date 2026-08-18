"""Define default titles for parameters, default combiners and relative ranges for 1D and 2D scans."""

from charm_fitter.utils import PlotParameter, Range

PARAMETERS = [
    # --- Mixing in the phenomenological parametrisation ---
    PlotParameter("x", title=r"$x$"),
    PlotParameter("y", title=r"$y$"),
    PlotParameter("qop", title=r"$|q/p|-1$", transf=lambda x: x - 1.0),
    PlotParameter("phi", title=r"$\phi_2$"),
    # --- Mixing in the theoretical parametrisation ---
    PlotParameter("x12", title=r"$x_{12}$"),
    PlotParameter("y12", title=r"$y_{12}$"),
    PlotParameter("phiM", title=r"$\phi_{2}^{M}$"),
    PlotParameter("phiG", title=r"$\phi_{2}^{\Gamma}$"),
    # --- CPV in the decay ---
    PlotParameter("Acp_KP", title=r"$a^{\rm d}_{K^+\pi^-}$"),
    PlotParameter("Acp_KK", title=r"$a^{\rm d}_{K^+K^-}$"),
    PlotParameter("Acp_PP", title=r"$a^{\rm d}_{\pi^+\pi^-}$"),
    PlotParameter("cot_delta_KK", title=r"$\cot\delta_{KK}$"),
    PlotParameter("cot_delta_PP", title=r"$\cot\delta_{\pi\pi}$"),
    PlotParameter("r_KK", title=r"$(P/T)_{KK}$"),
    PlotParameter("r_PP", title=r"$(P/T)_{\pi\pi}$"),
    PlotParameter("delta_KK", title=r"$\delta_{KK}$"),
    PlotParameter("delta_PP", title=r"$\delta_{\pi\pi}$"),
    # --- Hadronic nuisance parameters ---
    PlotParameter("Delta_Kpi", title=r"$\Delta_{K\pi}$"),
    PlotParameter("Delta_Kpipi0", title=r"$\Delta_{K\pi\pi^0}$"),
    PlotParameter("Delta_K3pi", title=r"$\Delta_{K\pi\pi\pi}$"),
    PlotParameter("r_Kpi", title=r"$r_{K\pi}$"),
    PlotParameter("r_K3pi", title=r"$r_{K\pi\pi\pi}$"),
    PlotParameter("r_Kpipi0", title=r"$r_{K\pi\pi^0}$"),
    PlotParameter("k_K3pi", title=r"$k_{K\pi\pi\pi}$"),
    PlotParameter("k_Kpipi0", title=r"$k_{K\pi\pi^0}$"),
    PlotParameter("F_pipipi0", title=r"$F_{\pi\pi\pi^0}$"),
    # --- Reparametrisations ---
    # NB: Avoid redefining existing parameters with different case (e.g. R_Kpi vs r_Kpi), since this will overwrite
    #     output files in case-insensitive file systems.
    PlotParameter("R_Kpi_", cf_par="r_Kpi", title=r"$R_{K\pi}$", transf=lambda x: x * x),
    PlotParameter("delta_Kpi_", cf_par="Delta_Kpi", title=r"$\delta_{K\pi}$", transf=lambda x: -x),
]

# Default scan parameters for 1D scans, tuned for different combiners.
RANGES_1D = {
    "LHCb-R1": {
        "x": Range((-2.0, 8.0), unit="1e-3"),
        "y": Range((4.0, 9.0), unit="1e-3"),
        "qop": Range((0.8, 1.2)),
        "phi": Range((-500.0, 250.0), unit="mrad"),
        "x12": Range((0.0, 8.0), unit="1e-3"),
        "y12": Range((4.0, 9.0), unit="1e-3"),
        "phiM": Range((-300.0, 500.0), unit="mrad"),
        "phiG": Range((-300.0, 500.0), unit="mrad"),
        "Acp_KK": Range((-5.0, 5.0), unit="1e-3"),
        "Acp_PP": Range((-5.0, 6.0), unit="1e-3"),
        "Acp_KP": Range((-0.2, 0.2), unit="%"),
        "Delta_Kpi": Range((-0.7, 0.3), unit="rad"),
        "r_Kpi": Range((5.8, 6.05), unit="%"),
        "r_K3pi": Range((5.2, 5.7), unit="%"),
        "k_K3pi": Range((0.2, 0.8)),
        "F_pipipi0": Range((0.91, 0.96)),
    },
    "LHCb-R12": {
        "x": Range((1.5, 7.0), unit="1e-3"),
        "y": Range((5.0, 8.0), unit="1e-3"),
        "qop": Range((0.92, 1.08)),
        "phi": Range((-120.0, 70.0), unit="mrad"),
        "x12": Range((1.5, 7.0), unit="1e-3"),
        "y12": Range((5.0, 8.0), unit="1e-3"),
        "phiM": Range((-120.0, 120.0), unit="mrad"),
        "phiG": Range((-150.0, 150.0), unit="mrad"),
        "Acp_KK": Range((-1.5, 3.0), unit="1e-3"),
        "Acp_PP": Range((0.0, 5.0), unit="1e-3"),
        "Acp_KP": Range((-2.0, 2.0), unit="%"),
        "Delta_Kpi": Range((-0.5, 0.1), unit="rad"),
        "r_Kpi": Range((5.82, 5.92), unit="%"),
        "r_K3pi": Range((5.25, 5.7), unit="%"),
        "k_K3pi": Range((0.2, 0.8)),
        "F_pipipi0": Range((0.92, 0.96)),
    },
}

RANGES_2D = {
    "LHCb-R12": {
        ("x12", "y12"): (Range((2.0, 6.0), unit="1e-3"), Range((5.4, 7.4), unit="1e-3")),
        ("phiM", "phiG"): (Range((-80.0, 100.0), unit="mrad"), Range((-80.0, 160.0), unit="mrad")),
        ("x", "y"): (Range((2.0, 6.0), unit="1e-3"), Range((5.4, 7.4), unit="1e-3")),
        ("qop", "phi"): (Range((0.9, 1.05)), Range((-120.0, 60.0), unit="mrad")),
        ("Delta_Kpi", "r_Kpi"): (Range((-0.5, 0.0), unit="rad"), Range((5.82, 5.92), unit="%")),
        ("Acp_KK", "Acp_PP"): (Range((-1.5, 3.0), unit="1e-3"), Range((0.0, 5.5), unit="1e-3")),
        ("Acp_KK", "phiM"): (Range((-2.0, 3.5), unit="1e-3"), Range((-75.0, 95.0), unit="mrad")),
    },
}
