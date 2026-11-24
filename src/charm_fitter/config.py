"""Define default titles for parameters, default combiners and relative ranges for 1D and 2D scans."""

from typing import Any

import numpy as np

from charm_fitter.utils import Combiner, DYFsc, MixParam, Plot1D, Plot2D, PlotParameter, Range, Scan1D, Scan2D

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
PARAMETERS_BY_NAME = {p.name: p for p in PARAMETERS}

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

# Ranges tuned to work both with --dcs-cpv and without it
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
        ("x12", "y12"): (Range((2.0, 7.0), unit="1e-3"), Range((5.3, 7.4), unit="1e-3")),
        ("phiM", "phiG"): (Range((-50.0, 80.0), unit="mrad"), Range((-80.0, 140.0), unit="mrad")),
        ("x", "y"): (Range((2.0, 7.0), unit="1e-3"), Range((5.3, 7.4), unit="1e-3")),
        ("qop", "phi"): (Range((0.93, 1.05)), Range((-100.0, 70.0), unit="mrad")),
        ("Acp_KK", "Acp_PP"): (Range((-1.5, 2.5), unit="1e-3"), Range((-0.5, 5.0), unit="1e-3")),
        ("Acp_KK", "phiM"): (Range((-2.0, 3.0), unit="1e-3"), Range((-50.0, 80.0), unit="mrad")),
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

# World averages: baseline combiner and plots --------------------------------------------------------------------------

# DeltaY(h- h+) final-state-correction hypotheses compared in the plots, and parameters for which they are compared.
DY_FSCS = [DYFsc.FULL, DYFsc.NONE]
DY_FSC_1D_VARS = ["phiM", "phiG", "qop", "phi", "Acp_KK", "Acp_PP"]
DY_FSC_2D_VARS = [("phiM", "phiG"), ("qop", "phi"), ("Acp_KK", "Acp_PP"), ("Acp_KK", "phiM")]

_WA_1D_LEG_RIGHT = ["x", "y", "qop", "x12", "y12", "Delta_Kpi", "r_Kpi", "r_K3pi", "r_Kpipi0", "k_K3pi", "R_Kpi_"]
_WA_1D_LEGFILL = ["phi", "phiM", "Acp_KK", "Acp_PP", "r_K3pi"]


def wa_combiner(
    id: int,
    title: str,
    parfile: str,
    ranges_1d: dict[str, Range],
    ranges_2d: dict[tuple[str, str], tuple[Range, Range]],
) -> Combiner:
    """Build the combiner of a world average, including the 1D and 2D scans for all DeltaY(h- h+) FSC hypotheses.

    Args:
        id: Numerical ID of the combiner.
        title: Title used in plot legends.
        parfile: Start-parameter file used for all scans.
        ranges_1d: Ranges of the 1D scans, by parameter name.
        ranges_2d: Ranges of the 2D scans, by (xpar, ypar).

    Returns:
        The world-average combiner.
    """
    dy_fscs = [dy_fsc for dy_fsc in DY_FSCS if dy_fsc != DYFsc.NONE]
    return Combiner(
        id,
        title,
        scans_1d=[Scan1D(par, r, parfile=parfile) for par, r in ranges_1d.items()]
        + [Scan1D(par, ranges_1d[par], dy_fsc=dy_fsc, parfile=parfile) for dy_fsc in dy_fscs for par in DY_FSC_1D_VARS],
        scans_2d=[
            Scan2D(xpar, ypar, xrange, yrange, parfile=parfile) for (xpar, ypar), (xrange, yrange) in ranges_2d.items()
        ]
        + [
            Scan2D(xpar, ypar, *ranges_2d[(xpar, ypar)], dy_fsc=dy_fsc, parfile=parfile)
            for dy_fsc in dy_fscs
            for (xpar, ypar) in DY_FSC_2D_VARS
        ],
        col="r",
        ls="-",
    )


def wa_plots_1d(combiner: str, ranges_1d: dict[str, Range]) -> list[Plot1D]:
    """1D plots of all parameters (including reparametrisations) whose charm-fitter parameter has a range defined."""
    return [
        Plot1D(
            [combiner],
            name,
            ranges_1d[par.cf_par].transformed_range(par.transf),
            legpos="r" if name in _WA_1D_LEG_RIGHT else "l",
            legfill=name in _WA_1D_LEGFILL,
            legfontsize=14,
        )
        for name, par in PARAMETERS_BY_NAME.items()
        if par.cf_par in ranges_1d
    ]


def wa_plots_2d(combiner: str, ranges_2d: dict[tuple[str, str], tuple[Range, Range]]) -> list[Plot2D]:
    """2D plots for all the given pairs of parameters."""
    return [
        Plot2D(
            [combiner],
            (xpar, ypar),
            xrange.transformed_range(PARAMETERS_BY_NAME[xpar].transf),
            yrange.transformed_range(PARAMETERS_BY_NAME[ypar].transf),
        )
        for (xpar, ypar), (xrange, yrange) in ranges_2d.items()
    ]


def dy_fsc_plots_1d(combiner: str, ranges_1d: dict[str, Range]) -> list[Plot1D]:
    """1D plots comparing the DeltaY(h- h+) final-state-correction hypotheses."""
    return [
        Plot1D(
            [combiner],
            par,
            ranges_1d[par].transformed_range(PARAMETERS_BY_NAME[par].transf),
            legpos="r",
            legfontsize=12,
            dy_fsc=DY_FSCS,
        )
        for par in DY_FSC_1D_VARS
    ]


def dy_fsc_plots_2d(
    combiner: str,
    ranges_2d: dict[tuple[str, str], tuple[Range, Range]],
    legfill: list[tuple[str, str]] | None = None,
) -> list[Plot2D]:
    """2D plots comparing the DeltaY(h- h+) final-state-correction hypotheses.

    Args:
        combiner: Name of the combiner.
        ranges_2d: Plot ranges by (xpar, ypar), used as they are (no transformation applied).
        legfill: Pairs of parameters for which the legend box should be filled.
    """
    return [
        Plot2D(
            [combiner],
            (xpar, ypar),
            xrange,
            yrange,
            logo="br",
            legfill=(xpar, ypar) in (legfill or []),
            legpos="l",
            legfontsize=14 if (xpar, ypar) != ("Acp_KK", "Acp_PP") else 13.5,
            dy_fsc=DY_FSCS,
        )
        for (xpar, ypar), (xrange, yrange) in ranges_2d.items()
    ]


# Breakdown plots ------------------------------------------------------------------------------------------------------

# Legend labels and plotting styles of the subsets of measurements entering breakdown plots.
BREAKDOWN_STYLES = {
    "cf-kpi": {"title": r"CLEO/BESIII $D \to K^\pm \pi^\mp$", "col": "y", "ls": "dash3dots"},
    "d0-kpi": {"title": r"$D^0 \to K^\pm \pi^\mp$", "col": "g", "ls": "dash2dots"},
    "hhprime": {"title": r"$D^0 \to h^+ h^{'-}$", "col": "b", "ls": "shortdash"},
    "kpi": {"title": r"$D^0 \to K^\pm \pi^\mp$", "col": "g", "ls": "dash2dots"},
    "k3pi": {"title": r"$D^0 \to K^\pm \pi^\mp \pi^+ \pi^-$", "col": "p", "ls": "dash3dots"},
    "kspipi": {"title": r"$D^0 \to K_{\rm S}^0 \pi^+ \pi^-$", "col": "lb", "ls": "shortdash"},
    "hh": {"title": r"$D^0 \to h^+ h^-$", "col": "y", "ls": "dashdot"},
}

# Pairs of mixing parameters scanned in the breakdown plots.
MIX_PAIRS = [("x", "y"), ("x12", "y12")]
CPV_PAIRS = [("qop", "phi"), ("phiM", "phiG")]

# Mixing parameters (complementary to the scanned ones) that are fixed in breakdown scans.
_BREAKDOWN_FIXED_MIX_PARS = {
    ("x", "y"): ["qop", "phi"],
    ("x12", "y12"): ["phiM", "phiG"],
    ("qop", "phi"): ["x", "y"],
    ("phiM", "phiG"): ["x12", "y12"],
}
# Additional parameters that are fixed in breakdown scans, as they are not constrained by the subset of measurements.
_BREAKDOWN_FIXED_PARS = {"hh": ["r_Kpi", "Delta_Kpi"]}

# Default settings (`Scan2D` keyword arguments) of the 2D breakdown scans, tuned on the latest world average. Ranges may
# need retuning when the subset of measurements changes: override them via `breakdown_combiner(..., overrides=...)`.
BREAKDOWN_SCANS_2D: dict[str, dict[tuple[str, str], dict[str, Any]]] = {
    "kpi": {
        ("x", "y"): {
            "xrange": Range((-8e-3, 8e-3)),
            "yrange": Range((3e-3, 8e-3)),
            "parfile": "config/charm-combo/start/WA-2025-06-breakdown-kpi-x-y.dat",
            "extra_args": "--npoints2dx 100 --npoints2dy 100",
        },
        ("x12", "y12"): {
            "xrange": Range((0.0, 8e-3)),
            "yrange": Range((3e-3, 8e-3)),
            "parfile": "config/charm-combo/start/WA-2025-06-breakdown-kpi-x12-y12.dat",
            "extra_args": "--npoints2dx 100 --npoints2dy 100",
        },
        ("qop", "phi"): {
            "xrange": Range((0.7, 1.4)),
            "yrange": Range((-0.4, 0.35)),
            "extra_args": "--npoints2dx 100 --npoints2dy 100",
        },
        ("phiM", "phiG"): {"xrange": Range((-0.1, 0.1)), "yrange": Range((-0.41, 0.3))},
    },
    "k3pi": {
        ("x", "y"): {"xrange": Range((-2e-3, 10e-3)), "yrange": Range((0.0, 10e-3))},
        ("x12", "y12"): {"xrange": Range((0.0, 10e-3)), "yrange": Range((0.0, 10e-3))},
        ("qop", "phi"): {"xrange": Range((0.9, 1.1)), "yrange": Range((-0.2, 0.2))},
        ("phiM", "phiG"): {"xrange": Range((-0.2, 0.2)), "yrange": Range((-0.2, 0.2))},
    },
    "kspipi": {
        ("x", "y"): {"xrange": Range((2e-3, 6.5e-3)), "yrange": Range((2e-3, 8e-3))},
        ("x12", "y12"): {"xrange": Range((2e-3, 6.5e-3)), "yrange": Range((2e-3, 8e-3))},
        ("qop", "phi"): {"xrange": Range((0.85, 1.15)), "yrange": Range((-0.2, 0.1))},
        ("phiM", "phiG"): {"xrange": Range((-0.2, 0.35)), "yrange": Range((-0.1, 0.2))},
    },
    "hh": {
        ("x", "y"): {
            "xrange": Range((-8e-3, 12.5e-3)),
            "yrange": Range((5e-3, 8e-3)),
            "extra_args": "--npoints2dx 100",
        },
        ("x12", "y12"): {
            "xrange": Range((0.0, 10.5e-3)),
            "yrange": Range((5e-3, 8e-3)),
            "extra_args": "--npoints2dx 100",
        },
        ("qop", "phi"): {
            "xrange": Range((0.78, 1.3)),
            "yrange": Range((-0.4, 0.4)),
            "extra_args": "--npoints2dx 80 --npoints2dy 80",
        },
        ("phiM", "phiG"): {"xrange": Range((-0.1, 0.15)), "yrange": Range((-0.5, 0.5))},
    },
}


def breakdown_combiner(
    key: str,
    pdf_ids: list[int],
    fix_parfile: str,
    pairs: list[tuple[str, str]] | None = None,
    overrides: dict[tuple[str, str], dict[str, Any]] | None = None,
) -> Combiner:
    """Build the combiner of a subset of measurements entering the breakdown plots.

    The mixing parameters complementary to the scanned ones (e.g. qop and phi when scanning x and y) are fixed to their
    values in `fix_parfile`, together with any parameter not constrained by the subset (see `_BREAKDOWN_FIXED_PARS`).

    Args:
        key: Subset of measurements, defining the plotting style and the default scan settings (`BREAKDOWN_SCANS_2D`).
        pdf_ids: IDs of the PDFs of the measurements in the subset.
        fix_parfile: Parameter file (typically, that of the world average) used to fix the parameters.
        pairs: Pairs of parameters to be scanned. Defaults to all pairs in `BREAKDOWN_SCANS_2D[key]`.
        overrides: `Scan2D` keyword arguments (e.g. `xrange`, `yrange`, `extra_args`) overriding the defaults, by pair.

    Returns:
        The breakdown combiner.
    """
    defaults = BREAKDOWN_SCANS_2D[key]
    overrides = overrides or {}
    return Combiner(
        pdf_ids,
        **BREAKDOWN_STYLES[key],
        scans_2d=[
            Scan2D(
                *pair,
                fixed_pars=_BREAKDOWN_FIXED_MIX_PARS[pair] + _BREAKDOWN_FIXED_PARS.get(key, []),
                fix_parfile=fix_parfile,
                **(defaults.get(pair, {}) | overrides.get(pair, {})),
            )
            for pair in (pairs if pairs is not None else list(defaults))
        ],
    )


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
    "cf-kpi": Combiner(
        [1, 21, 40, 41],
        **BREAKDOWN_STYLES["cf-kpi"],
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
    # The theoretical parametrisation gets stuck in local minima of (x12, y12, phiM, phiG) for this subset, which
    # makes the contours jagged: use the phenomenological one (the plotted parameters are independent of it)
    "d0-kpi": Combiner(
        [250, 252, 260, 286, 288],
        **BREAKDOWN_STYLES["d0-kpi"],
        scans_2d=[
            Scan2D(
                "Delta_Kpi",
                "r_Kpi",
                xrange=_KPI_RD_SCAN_LHCB.xrange,
                yrange=_KPI_RD_SCAN_LHCB.yrange,
                parfile="config/charm-combo/start/WA-2026-11-breakdown-kpi-rd-pheno.dat",
                extra_args="--npoints2dx 200",
                mix_param=MixParam.PHENO,
            )
        ],
    ),
    # LHCb D0 -> K pi and yCP(hh) - yCP(RS)
    "hhprime": Combiner([140, 141, 286, 288], **BREAKDOWN_STYLES["hhprime"], scans_2d=[_KPI_RD_SCAN_LHCB]),
}


def breakdown_plot_kpi_rd(combiner: str) -> Plot2D:
    """Breakdown plot of the constraints on the D0 -> K+ pi- hadronic parameters (`BREAKDOWN_COMBINERS_KPI_RD`)."""
    return Plot2D(
        [*BREAKDOWN_COMBINERS_KPI_RD, combiner],
        ("Delta_Kpi", "r_Kpi"),
        xrange=Range((-1.3, 0.7), unit="rad"),
        yrange=Range((5.6, 6.1), unit="1e-2"),
        combiners_label=combiner,
        levels=2,
        logo="l",
        legpos="bl",
        legfill=True,
    )
