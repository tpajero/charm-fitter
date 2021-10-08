/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include "PDF_yCP_minus_yCP_RS.h"


PDF_yCP_minus_yCP_RS::PDF_yCP_minus_yCP_RS(
        TString cObs, TString cErr, TString cCor, const theory_config& th_cfg)
: PDF_Abs(1)
{
    name = "yCP_minus_yCP_RS_" + cObs;
    initParameters(th_cfg);
    initRelations(th_cfg);
    initObservables(cObs);
    setObservables(cObs);
    setUncertainties(cErr);
    setCorrelations(cCor);
    buildCov();
    buildPdf();
}


PDF_yCP_minus_yCP_RS::~PDF_yCP_minus_yCP_RS() {}


void PDF_yCP_minus_yCP_RS::initParameters(const theory_config& th_cfg)
{
    ParametersCharmCombo p(th_cfg);
    parameters = new RooArgList("parameters");

    parameters->add(*(p.get("R_Kpi")));
    parameters->add(*(p.get("Delta_Kpi")));
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
            parameters->add(*(p.get("phiM")));
            parameters->add(*(p.get("x12")));
            parameters->add(*(p.get("y12")));
            break;
        default:
            cout << "PDF_yCP::initParameters : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
}


void PDF_yCP_minus_yCP_RS::initRelations(const theory_config& th_cfg)
{
    theory = new RooArgList("theory");
    switch (th_cfg) {
        case phenomenological:
            theory->add(
                    *(new RooFormulaVar(
                            "yCP_minus_yCP_RS_th", "yCP_minus_yCP_RS_th",
                            "0.5*( "
                            "      y*(qop+1 + 1/(qop+1))*cos(phi)"
                            "    - x*(qop+1 - 1/(qop+1))*sin(phi)"
                            " + sqrt(R_Kpi/100) * ("
                            "      (y * cos(Delta_Kpi) - x * sin(Delta_Kpi)) * (qop+1 + 1/(qop+1)) * cos(phi)"
                            "    - (x * cos(Delta_Kpi) + y * sin(Delta_Kpi)) * (qop+1 - 1/(qop+1)) * sin(phi)))",
                            *(RooArgSet*)parameters)));
            break;
        case theoretical:
            theory->add(
                    *(new RooFormulaVar(
                            "yCP_minus_yCP_RS_th", "yCP_minus_yCP_RS_th",
                            " y12 * cos(phiG)"
                            " + sqrt(R_Kpi/100) * ("
                            "       y12 * cos(Delta_Kpi) * cos(phiG)"
                            "     - x12 * sin(Delta_Kpi) * cos(phiM))",
                            *(RooArgSet*)parameters)));
            break;
        case superweak:
            theory->add(
                    *(new RooFormulaVar(
                            "yCP_minus_yCP_RS_th", "yCP_minus_yCP_RS_th",
                            " y12"
                            " + sqrt(R_Kpi/100) * ("
                            "       y12 * cos(Delta_Kpi)"
                            "     - x12 * sin(Delta_Kpi) * cos(phiM))",
                            *(RooArgSet*)parameters)));
            break;
        default:
            cout << "PDF_yCP::initRelations : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
}


void PDF_yCP_minus_yCP_RS::initObservables(const TString& setName)
{
    observables = new RooArgList("observables");
    observables->add(*(new RooRealVar(
            "yCP_minus_yCP_RS_obs",
            setName + "   #it{y_{CP}}#minus#it{y_{CP}^{K^{#minus}#pi^{+}}}",
            0, -1e4, 1e4)));
}


void PDF_yCP_minus_yCP_RS::setObservables(TString c)
{
    if (c.EqualTo("truth")) {
        setObservablesTruth();
    }
    else if (c.EqualTo("toy")) {
        setObservablesToy();
    }
    else if (c.EqualTo("WA2020")) {
        obsValSource = "https://cds.cern.ch/record/2747731";
        setObservable("yCP_minus_yCP_RS_obs", 0.738);
    }
    else {
        cout << "PDF_yCP_minus_yCP_RS::setObservables() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_yCP_minus_yCP_RS::setUncertainties(TString c)
{
    if (c.EqualTo("WA2020")) {
        obsErrSource = "https://cds.cern.ch/record/2747731";
        StatErr[0] = 0.111;
        SystErr[0] = 0;
    }
    else {
        cout << "PDF_yCP_minus_yCP_RS::setUncertainties() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_yCP_minus_yCP_RS::setCorrelations(TString c)
{
    resetCorrelations();
    corSource = "No correlations for one observable";
}


void PDF_yCP_minus_yCP_RS::buildPdf()
{
    pdf = new RooMultiVarGaussian(
            "pdf_" + name, "pdf_" + name, *(RooArgSet*)observables,
            *(RooArgSet*)theory, covMatrix);
}
