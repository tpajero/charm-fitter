import os
import sys

config_dir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.abspath(config_dir))

from classes import Parameter, Plot2d, Subcombination

# List all CharmCombo parameters with the properties to be used for plotting.
params = {
    "x12": Parameter("x12", r"$x_{12}$ [\%]", (0.1, 0.7), "theo"),
    "y12": Parameter("y12", r"$y_{12}$ [\%]", (0.5, 0.8), "theo"),
    "phiM": Parameter("phiM", r"$\phi_2^M$ [rad]", (-0.125, 0.125), "theo"),
    "phiG": Parameter("phiG", r"$\phi_2^{\Gamma}$ [rad]", (-0.15, 0.15), "theo"),
    "x": Parameter("x", r"$x$ [\%]", (0.1, 0.7), "pheno"),
    "y": Parameter("y", r"$y$ [\%]", (0.5, 0.8), "pheno"),
    "qop": Parameter("qop", r"$|q/p|-1$", (-0.08, 0.08), "pheno"),
    "phi": Parameter("phi", r"$\phi_2$ [rad]", (-0.125, 0.125), "pheno"),
    "Acp_KK": Parameter("Acp_KK", r"$a^{\rm d}_{K^+K^-}\, [\%]$", (-0.2, 0.4), "both"),
    "Acp_PP": Parameter("Acp_PP", r"$a^{\rm d}_{\pi^+\pi^-}\, [\%]$", (-0.1, 0.6), "both"),
    "Acp_KP": Parameter("Acp_KP", r"$a^{\rm d}_{K^+\pi^-}\, [\%]$", (-2.0, 2.0), "both"),
    "Delta_Kpi": Parameter("Delta_Kpi", r"$\Delta_{K\pi}$ [rad]", (-0.5, 0.2), "both"),
    "Delta_Kpipi0": Parameter("Delta_Kpipi0", r"$\Delta_{K\pi\pi^0}$ [rad]", (-2.0, 2.0), "both"),
    "Delta_K3pi": Parameter("Delta_K3pi", r"$\Delta_{K\pi\pi\pi}$ [rad]", (-0.5, 1.5), "both"),
    "R_Kpi": Parameter("R_Kpi", r"$R_{K\pi}$ [\%]", (0.335, 0.355), "both"),
    "r_K3pi": Parameter("r_K3pi", r"$r_{K\pi\pi\pi}$ [\%]", (4.0, 6.0), "both"),
    "r_K3pipi0": Parameter("r_Kpipi0", r"$r_{K\pi\pi^0}$ [\%]", (4.0, 6.0), "both"),
    "k_K3pi": Parameter("k_K3pi", r"$k_{K\pi\pi\pi}$", (0.0, 1.0), "both"),
    "k_Kpipi0": Parameter("k_Kpipi0", r"$k_{K\pi\pi^0}$", (0.0, 1.0), "both"),
    "F_pipipi0": Parameter("F_pipipi0", r"$F_{\pi\pi\pi^0}$", (0.9, 1.0), "both"),
}


# List all pairs of parameters for which 2D scans should be prepared.
plots_2d = [
    Plot2d(("x12", "y12")),
    Plot2d(("phiM", "phiG"), "br", "bl"),
    Plot2d(("x", "y")),
    Plot2d(("qop", "phi")),
    Plot2d(("Delta_Kpi", "R_Kpi")),
    Plot2d(("Acp_KK", "Acp_PP")),
    Plot2d(("Acp_KK", "phiM")),
]


# List the IDs and titles for all combinations that may be passed to the main plotting script.
combiners_info = {
    55: ("WAOct2025NoFSC", "WA 2025"),
    56: ("WAOct2025PartialFSC", ""),
    57: ("WAOct2025FullFSC", ""),
}
for k in list(combiners_info.keys()):
    combiners_info[k + 1000] = (
        combiners_info[k][0] + "NoDcsCPV",
        combiners_info[k][1] + r", $a^{\rm d}_{K^+\pi^-} = 0$",
    )


def combiners_acp(no_dcs_cpv=False):
    combiners = {
        57: r"$\Delta Y_{f} = - x_{12}\sin\phi^M_2 + y_{12} a^{\rm d}_f (1 + \cot\delta x_{12} / y_{12})$",
        56: r"$\Delta Y_{f} = - x_{12}\sin\phi^M_2 + y_{12} a^{\rm d}_f$",
        55: r"$\Delta Y_{f} = - x_{12}\sin\phi^M_2$",
    }
    if no_dcs_cpv:
        tmp = combiners.copy()
        combiners = {}
        for k, v in tmp.items():
            combiners[k + 1000] = v
    return combiners


# Dictionary of all subcombinations, of their titles and PDFs.
sub_comb_dict = {
    "cc_kpi": Subcombination(r"CLEO/BESIII $D^0 \to K^\pm \pi^\mp$", [50, 55, 56]),
    "lhcb_kpi": Subcombination(r"LHCb $D^0 \to K^\pm \pi^\mp$", [39, 41]),
    "all_kpi": Subcombination(r"$D^0 \to K^\pm \pi^\mp$", [39, 41, 50, 55, 56]),
    "k3pi": Subcombination(r"$D^0 \to K^\pm \pi^\mp \pi^+ \pi^-$", [5, 54]),
    "kspipi": Subcombination(r"$D^0 \to K_S^0 \pi^+\pi^-$", [1, 3, 21, 24]),
    "hh": Subcombination(r"$D^0 \to h^+ h^-$", [61, 62, 72, 64, 72]),
}

# Style for 1D and 2D plots.
colors = ["r", "lb", "y", "o", "g", "p", "b"]
ls = ["longdash", "-", "shortdash", "dashdot", "dash3dots", "dash2dots", "dots"]

colors_2d = ["r", "lb", "y", "o", "g", "p", "b"]
ls_2d = ["longdash", "-", "shortdash", "dashdot", "dash3dots", "dash2dots", "dots"]
markers = ["o", "P", "s"]

colors_acp = ["y", "r", "lb"]
ls_acp = ["shortdash", "longdash", "-"]
markers_acp = ["o", "s", "P"]
