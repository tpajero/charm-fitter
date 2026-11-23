"""Define default titles for parameters, default combiners and relative ranges for 1D and 2D scans."""

import numpy as np

from charm_fitter.utils import Combiner, MixParam, PlotParameter, Range, Scan1D, Scan2D

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
        "x": Range((-3.0, 8.0), unit="1e-3"),
        "y": Range((4.0, 8.0), unit="1e-3"),
        "qop": Range((0.75, 1.2)),
        "phi": Range((-500.0, 250.0), unit="mrad"),
        "x12": Range((0.0, 8.0), unit="1e-3"),
        "y12": Range((4.0, 8.0), unit="1e-3"),
        "phiM": Range((-np.pi, np.pi), unit="rad"),
        "phiG": Range((-300.0, 500.0), unit="mrad"),
        "Acp_KK": Range((-5.0, 5.0), unit="1e-3"),
        "Acp_PP": Range((-5.0, 6.0), unit="1e-3"),
        "Acp_KP": Range((-8.0, 3.0), unit="%"),
        "Delta_Kpi": Range((-0.7, 0.3), unit="rad"),
        "r_Kpi": Range((5.74, 6.03), unit="%"),
        "r_K3pi": Range((5.2, 5.7), unit="%"),
        "k_K3pi": Range((0.2, 0.8)),
        "F_pipipi0": Range((0.92, 0.97)),
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
        "r_K3pi": Range((5.2, 5.7), unit="%"),
        "k_K3pi": Range((0.2, 0.8)),
        "F_pipipi0": Range((0.92, 0.97)),
    },
}

RANGES_2D = {
    "LHCb-R1": {
        ("Delta_Kpi", "r_Kpi"): (Range((-0.7, 0.4), unit="rad"), Range((5.7, 6.1), unit="%")),
        ("x12", "y12"): (Range((0.0, 8.0), unit="1e-3"), Range((3.5, 8.5), unit="1e-3")),
        ("phiM", "phiG"): (Range((-np.pi, np.pi), unit="rad"), Range((-400.0, 600.0), unit="mrad")),
        ("x", "y"): (Range((-4.0, 8.0), unit="1e-3"), Range((3.5, 8.5), unit="1e-3")),
        ("qop", "phi"): (Range((0.7, 1.2)), Range((-550.0, 400.0), unit="mrad")),
        ("Acp_KK", "Acp_PP"): (Range((-7.0, 7.0), unit="1e-3"), Range((-7.0, 8.0), unit="1e-3")),
        ("Acp_KK", "phiM"): (Range((-6.0, 6.0), unit="1e-3"), Range((-np.pi, np.pi), unit="rad")),
    },
    "LHCb-R12": {
        ("Delta_Kpi", "r_Kpi"): (Range((-0.5, 0.0), unit="rad"), Range((5.82, 5.92), unit="%")),
        ("x12", "y12"): (Range((2.0, 6.0), unit="1e-3"), Range((5.4, 7.4), unit="1e-3")),
        ("phiM", "phiG"): (Range((-80.0, 100.0), unit="mrad"), Range((-80.0, 160.0), unit="mrad")),
        ("x", "y"): (Range((2.0, 6.0), unit="1e-3"), Range((5.4, 7.4), unit="1e-3")),
        ("qop", "phi"): (Range((0.9, 1.05)), Range((-120.0, 60.0), unit="mrad")),
        ("Acp_KK", "Acp_PP"): (Range((-1.5, 3.0), unit="1e-3"), Range((0.0, 5.5), unit="1e-3")),
        ("Acp_KK", "phiM"): (Range((-2.0, 3.5), unit="1e-3"), Range((-75.0, 95.0), unit="mrad")),
    },
}

COMBINERS = {
    "LHCb-R1": Combiner(
        300,
        "LHCb Run 1",
        scans_1d=[
            Scan1D(
                par,
                r,
                parfile="config/charm-combo/start/LHCb-R1.dat",
                extra_args="--pr " + ("--npoints 300" if par == "phiM" else ""),
            )
            for par, r in RANGES_1D["LHCb-R1"].items()
        ],
        scans_2d=[
            Scan2D(
                xpar,
                ypar,
                xrange,
                yrange,
                parfile="config/charm-combo/start/"
                + (
                    "LHCb-R1-AcpKK-phiM.dat"
                    if (xpar, ypar) == ("Acp_KK", "phiM")
                    else "LHCb-R1-phiM-phiG.dat"
                    if (xpar, ypar) == ("phiM", "phiG")
                    else "LHCb-R1.dat"
                ),
                extra_args="--pr "
                + ("--npoints2dx 150" if xpar == "phiM" else "--npoints2dy 150" if ypar == "phiM" else ""),
            )
            for (xpar, ypar), (xrange, yrange) in RANGES_2D["LHCb-R1"].items()
        ],
    ),
    "LHCb-R12": Combiner(
        301,
        "LHCb Run 1+2",
        scans_1d=[
            Scan1D(par, r, parfile="config/charm-combo/start/LHCb-R12.dat") for par, r in RANGES_1D["LHCb-R12"].items()
        ],
        scans_2d=[
            Scan2D(xpar, ypar, xrange, yrange, parfile="config/charm-combo/start/LHCb-R12.dat")
            for (xpar, ypar), (xrange, yrange) in RANGES_2D["LHCb-R12"].items()
        ],
    ),
}

# Legend labels and plotting styles of the subsets of measurements entering breakdown plots. Use them as
# `Combiner(pdf_ids, **BREAKDOWN_STYLES[key], scans_2d=...)`.
BREAKDOWN_STYLES = {
    "ext-kpi": {"title": r"CLEO/BESIII $D \to K^\pm \pi^\mp$", "col": "y", "ls": "dash3dots"},
    "lhcb-kpi": {"title": r"LHCb $D^0 \to K^\pm \pi^\mp$", "col": "g", "ls": "dash2dots"},
    "lhcb-hhprime": {"title": r"LHCb $D^0 \to h^+ h^{'-}$", "col": "b", "ls": "shortdash"},
    "kpi": {"title": r"$D^0 \to K^\pm \pi^\mp$", "col": "g", "ls": "dash2dots"},
    "k3pi": {"title": r"$D^0 \to K^\pm \pi^\mp \pi^+ \pi^-$", "col": "p", "ls": "dash3dots"},
    "kspipi": {"title": r"$D^0 \to K_{\rm S}^0 \pi^+ \pi^-$", "col": "lb", "ls": "shortdash"},
    "hh-mix": {"title": r"$D^0 \to h^+ h^-$", "col": "y", "ls": "dashdot"},
    "hh": {"title": r"$D^0 \to h^+ h^-$", "col": "y", "ls": "dashdot"},
}

# Breakdown of the constraints on the D0 -> K+ pi- hadronic parameters, independent of the specific average.
_KPI_RD_SCAN_LHCB = Scan2D(
    "Delta_Kpi",
    "r_Kpi",
    xrange=Range((-1.35, 0.75)),
    yrange=Range((5.81e-2, 5.92e-2)),
    parfile="config/charm-combo/start/WA-2025-06-breakdown-kpi-rd.dat",
    extra_args="--npoints2dx 200",
)
BREAKDOWN_COMBINERS_KPI_RD = {
    "ext-kpi": Combiner(
        [1, 21, 40, 41],
        **BREAKDOWN_STYLES["ext-kpi"],
        scans_2d=[
            Scan2D(
                "Delta_Kpi",
                "r_Kpi",
                xrange=Range((-0.6, 0.3)),
                yrange=Range((5.55e-2, 6.15e-2)),
                mix_param=MixParam.PHENO,
            )
        ],
    ),
    "lhcb-kpi": Combiner([286, 288], **BREAKDOWN_STYLES["lhcb-kpi"], scans_2d=[_KPI_RD_SCAN_LHCB]),
    # LHCb D0 -> K pi and yCP(hh) - yCP(RS)
    "lhcb-hhprime": Combiner([140, 141, 286, 288], **BREAKDOWN_STYLES["lhcb-hhprime"], scans_2d=[_KPI_RD_SCAN_LHCB]),
}
