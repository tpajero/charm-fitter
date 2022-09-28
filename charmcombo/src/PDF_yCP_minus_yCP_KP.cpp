/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include "PDF_yCP_minus_yCP_KP.h"


PDF_yCP_minus_yCP_KP::PDF_yCP_minus_yCP_KP(TString measurement_id, const theory_config& th_cfg)
    : PDF_Abs{1}, th_cfg{th_cfg}
{
    name = "yCP_minus_yCP_KP_" + measurement_id;
    initParameters();
    initRelations();
    initObservables(measurement_id);
    setObservables(measurement_id);
    setUncertainties(measurement_id);
    setCorrelations(measurement_id);
    buildCov();
    buildPdf();
}


PDF_yCP_minus_yCP_KP::~PDF_yCP_minus_yCP_KP() {}


void PDF_yCP_minus_yCP_KP::initParameters() {
    ParametersCharmCombo p;
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
            parameters->add(*(p.get("x12")));
            parameters->add(*(p.get("y12")));
            break;
        default:
            cout << "PDF_yCP::initParameters : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
}


void PDF_yCP_minus_yCP_KP::initRelations() {
    theory = new RooArgList("theory");
    switch (th_cfg) {
        case phenomenological:
            theory->add(
                    *(Utils::makeTheoryVar(
                            "yCP_minus_yCP_KP_th", "yCP_minus_yCP_KP_th",
                            " 0.5*( "
                            "       y*(qop+1 + 1/(qop+1))*cos(phi)"
                            "     - x*(qop+1 - 1/(qop+1))*sin(phi))"
                            " + sqrt(R_Kpi/100) * cos(Delta_Kpi) * ("
                            "      y * (qop+1 + 1/(qop+1)) * cos(phi)"
                            "    - x * (qop+1 - 1/(qop+1)) * sin(phi))",
                            parameters)));
            break;
        case theoretical:
            theory->add(
                    *(Utils::makeTheoryVar(
                            "yCP_minus_yCP_KP_th", "yCP_minus_yCP_KP_th",
                            "y12*cos(phiG)"
                            "+ 2 * sqrt(R_Kpi/100) * y12 * cos(Delta_Kpi) * cos(phiG)",
                            parameters)));
            break;
        case superweak:
            theory->add(
                    *(Utils::makeTheoryVar(
                            "yCP_minus_yCP_KP_th", "yCP_minus_yCP_KP_th",
                            "y12 + 2 * sqrt(R_Kpi/100) * y12 * cos(Delta_Kpi)",
                            parameters)));
            break;
        default:
            cout << "PDF_yCP::initRelations : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
}


void PDF_yCP_minus_yCP_KP::initObservables(const TString& setName) {
    observables = new RooArgList("observables");
    observables->add(*(new RooRealVar(
                    "yCP_minus_yCP_KP_obs",
                    setName + "   #it{y_{CP}}#minus#it{y_{CP}^{K#pi}}",
                    0, -1e4, 1e4)));
}


void PDF_yCP_minus_yCP_KP::setObservables(TString c) {
    if (c.EqualTo("truth")) setObservablesTruth();
    else if (c.EqualTo("toy")) setObservablesToy();
    else if (c.EqualTo("WA2020")) {
        obsValSource = "https://cds.cern.ch/record/2747731";
        setObservable("yCP_minus_yCP_KP_obs", 0.732);
    } else {
        cout << "PDF_yCP_minus_yCP_KP::setObservables() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_yCP_minus_yCP_KP::setUncertainties(TString c) {
    if (c.EqualTo("WA2020")) {
        obsErrSource = "https://cds.cern.ch/record/2747731";
        StatErr[0] = 3.068;
        SystErr[0] = 0;
    } else {
        cout << "PDF_yCP_minus_yCP_KP::setUncertainties() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_yCP_minus_yCP_KP::setCorrelations(TString c) {
    resetCorrelations();
    corSource = "No correlations for one observable";
}


void PDF_yCP_minus_yCP_KP::buildPdf() {
    pdf = new RooMultiVarGaussian(
            "pdf_" + name, "pdf_" + name, *(RooArgSet*)observables,
            *(RooArgSet*)theory, covMatrix);
}
