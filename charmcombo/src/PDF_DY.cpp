/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include "PDF_DY.h"


PDF_DY::PDF_DY(TString measurement_id, const theory_config& th_cfg)
    : PDF_Abs{1}, th_cfg{th_cfg}
{
    name = "DY_" + measurement_id;
    initParameters();
    initRelations();
    initObservables(measurement_id);
    setObservables(measurement_id);
    setUncertainties(measurement_id);
    setCorrelations(measurement_id);
    buildCov();
    buildPdf();
}


PDF_DY::~PDF_DY() {}


void PDF_DY::initParameters() {
    ParametersCharmCombo p(th_cfg);
    parameters = new RooArgList("parameters");
    switch (th_cfg) {
        case phenomenological:
            parameters->add(*(p.get("x")));
            parameters->add(*(p.get("y")));
            parameters->add(*(p.get("qop")));
            parameters->add(*(p.get("phi")));
            break;
        case theoretical:
        case superweak:
            parameters->add(*(p.get("x12")));
            parameters->add(*(p.get("y12")));
            parameters->add(*(p.get("phiM")));
            break;
        default:
            cout << "PDF_DY::initParameters : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
}


void PDF_DY::initRelations() {
    theory = new RooArgList("theory");
    switch (th_cfg) {
        case phenomenological:
            theory->add(
                    *(new RooFormulaVar(
                            "DY_th", "DY_th",
                            "0.5 * (- y*(qop+1 - 1/(qop+1))*cos(phi)  "
                            "       + x*(qop+1 + 1/(qop+1))*sin(phi) )",
                            *(RooArgSet*)parameters)));
            break;
        case theoretical:
        case superweak:
            theory->add(
                    *(new RooFormulaVar(
                            "DY_th", "DY_th", "-x12 * sin(phiM)",
                            *(RooArgSet*)parameters)));
            break;
        default:
            cout << "PDF_DY::initRelations : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
}


void PDF_DY::initObservables(const TString& setName) {
    observables = new RooArgList("observables");
    observables->add(*(new RooRealVar("DY_obs", setName + "   #Delta#it{Y}",  0, -1e4, 1e4)));
}


void PDF_DY::setObservables(TString c) {
    if (c.EqualTo("truth")) setObservablesTruth();
    else if (c.EqualTo("toy")) setObservablesToy();
    else if (c.EqualTo("WA2020")) {
        obsValSource = "https://cds.cern.ch/record/2747731";
        setObservable("DY_obs", 3.1e-2);
    } else if (c.EqualTo("WA2021")) {
        obsValSource = "https://cds.cern.ch/record/2747731";
        setObservable("DY_obs", -0.92e-2);
    } else {
        cout << "PDF_DY::setObservables() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_DY::setUncertainties(TString c) {
    if (c.EqualTo("WA2020")) {
        obsErrSource = "https://cds.cern.ch/record/2747731";
        StatErr[0] = 2.0e-2;
        SystErr[0] = 0.;
    } else if (c.EqualTo("WA2021")) {
        obsErrSource = "https://cds.cern.ch/record/2747731";
        StatErr[0] = 1.11e-2;
        SystErr[0] = 0.33e-2;
    } else {
        cout << "PDF_DY::setUncertainties() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_DY::setCorrelations(TString c) {
    resetCorrelations();
    corSource = "No correlations for one observable";
}


void PDF_DY::buildPdf() {
    pdf = new RooMultiVarGaussian(
            "pdf_" + name, "pdf_" + name, *(RooArgSet*)observables,
            *(RooArgSet*)theory, covMatrix);
}
