/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include "PDF_scan_DY_RS.h"


PDF_scan_DY_RS::PDF_scan_DY_RS(const theory_config& th_cfg)
: PDF_Abs(1)
{
    name = "scan_DY_RS";
    initParameters(th_cfg);
    initRelations(th_cfg);
    initObservables();
    setObservables();
    setUncertainties();
    setCorrelations();
    buildCov();
    buildPdf();
}


PDF_scan_DY_RS::~PDF_scan_DY_RS() {}


void PDF_scan_DY_RS::initParameters(const theory_config& th_cfg)
{
    ParametersCharmCombo p(th_cfg, kTRUE);
    parameters = new RooArgList("parameters");
    parameters->add(*(p.get("R_Kpi")));
    if (th_cfg != superweak)
        parameters->add(*(p.get("A_Kpi")));
    parameters->add(*(p.get("Delta_Kpi")));
    parameters->add(*(p.get("DY_RS")));

    switch (th_cfg) {
        case phenomenological:
            parameters->add(*(p.get("x")));
            parameters->add(*(p.get("y")));
            parameters->add(*(p.get("qop")));
            parameters->add(*(p.get("phi")));
            break;
        case theoretical:
            parameters->add(*(p.get("phiG")));
        case superweak:
            parameters->add(*(p.get("x12")));
            parameters->add(*(p.get("y12")));
            parameters->add(*(p.get("phiM")));
            break;
        default:
            cout << "PDF_scan_DY_RS::initParameters : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
}


void PDF_scan_DY_RS::initRelations(const theory_config& th_cfg)
{
    theory = new RooArgList("theory");
    switch (th_cfg) {
        case phenomenological:
            theory->add(
                    *(new RooFormulaVar(
                            "DY_RS_scan_th","DY_RS_scan_th",
                            "DY_RS - abs(100. * 0.5 * pow(R_Kpi/100, 0.5) * "
                            "(  (y*cos(Delta_Kpi) - x*sin(Delta_Kpi))*((qop+1) - 1/(qop+1) - A_Kpi/100)*cos(phi)"
                            " - (x*cos(Delta_Kpi) + y*sin(Delta_Kpi))*((qop+1) + 1/(qop+1)            )*sin(phi)))",
                            *(RooArgSet*)parameters)));
            break;
        case theoretical:
            theory->add(
                    *(new RooFormulaVar(
                            "DY_RS_scan_th","DY_RS_scan_th",
                            "DY_RS - abs(100 * pow(R_Kpi/100, 0.5) * "
                            "(  (-y12*cos(Delta_Kpi)*cos(phiG) + x12*sin(Delta_Kpi)*cos(phiM))*(A_Kpi/100)*0.5"
                            " + (+y12*sin(Delta_Kpi)*sin(phiG) + x12*cos(Delta_Kpi)*sin(phiM))                 ))",
                            *(RooArgSet*)parameters)));
            break;
        case superweak:
            theory->add(
                    *(new RooFormulaVar(
                            "DY_RS_scan_th","DY_RS_scan_th",
                            "DY_RS - abs(100 * pow(R_Kpi/100, 0.5) * "
                            "            x12 * cos(Delta_Kpi) * sin(phiM))",
                            *(RooArgSet*)parameters)));
            break;
        default:
            cout << "PDF_scan_DY_RS::initRelations : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
}


void PDF_scan_DY_RS::initObservables()
{
    observables = new RooArgList("observables");
    observables->add(
            *(new RooRealVar(
                    "DY_RS_scan_obs",
                    "scan   #Delta#it{Y}_{#it{K^{#minus}#pi^{+}}}",
                    0, -1e4, 1e4)));
}

void PDF_scan_DY_RS::setObservables()
{
    setObservable("DY_RS_scan_obs", 0.);
}


void PDF_scan_DY_RS::setUncertainties()
{
    StatErr[0] = 0.005;
    SystErr[0] = 0.;
}


void PDF_scan_DY_RS::setCorrelations()
{
    resetCorrelations();
    corSource = "No correlations for one observable";
}


void PDF_scan_DY_RS::buildPdf()
{
    pdf = new RooMultiVarGaussian(
            "pdf_" + name, "pdf_" + name, *(RooArgSet*)observables,
            *(RooArgSet*)theory, covMatrix);
}
