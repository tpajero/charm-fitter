class Parameter:
    """Class to define the objects related to a scan parameter in a single place."""

    def __init__(self, identifier, title, range1d, parametrisation, range2d=None, degrees=False):
        assert parametrisation in ["theo", "pheno", "both"]
        self.id = identifier
        self.title = title
        self.parametrisation = parametrisation
        self.range1d = range1d
        self.range2d = range2d if range2d else range1d
        self.degrees = degrees
        if self.degrees:
            self.title += r"\,[^\circ]$"


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


class Plot2d:
    """Class to define the properties of 2D plots."""

    def __init__(self, pars, logo="r", legpos="l"):
        self.pars = pars
        self.logo = logo
        self.legpos = legpos


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
}
for k in list(combiners_info.keys()):
    combiners_info[k + 1000] = (
        combiners_info[k][0] + "NoDcsCPV",
        combiners_info[k][1] + r", $a^{\rm d}_{K^+\pi^-} = 0$",
    )


class Subcombination:
    """Class to define the PDFs and title of a subcombination of measurements."""

    def __init__(self, title, pdfs):
        self.title = title
        self.pdfs = pdfs


# dict which lists which pdfs are included for each subcombination
sub_comb_dict = {
    "cc_kpi": Subcombination(r"CLEO/BES $D^0 \to K^\pm \pi^\mp$", [254, 256]),
    "lhcb_kpi": Subcombination(r"LHCb $D^0 \to K^\pm \pi^\mp$", [220, 216]),
    "all_kpi": Subcombination(r"$D^0 \to K^\pm \pi^\mp$", [220, 216, 254, 256]),
    "k3pi": Subcombination(r"$D^0 \to K_S^0 \pi^+\pi^-$", [206, 188]),
    "kspipi": Subcombination(r"$D^0 \to K_S^0 \pi^+\pi^-$", [201, 202, 222, 188]),
    "hh": Subcombination(r"$D^0 \to K^\pm \pi^\mp$", [219, 213, 217, 214]),
}

colors = ["r", "lb", "y", "o", "g", "p", "b"]
ls = ["longdash", "-", "shortdash", "dashdot", "dash3dots", "dash2dots", "dots"]

colors_2d = ["r", "lb", "y", "o", "g", "p", "b"]
ls_2d = ["longdash", "-", "shortdash", "dashdot", "dash3dots", "dash2dots", "dots"]
markers = ["o", "P"]
