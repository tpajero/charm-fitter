/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include "PDF_yCP.h"


PDF_yCP::PDF_yCP(TString measurement_id, const theory_config& th_cfg)
    : PDF_Abs{1}, th_cfg{th_cfg}
{
    name = "yCP_" + measurement_id;
    initParameters();
    initRelations();
    initObservables(measurement_id);
    setObservables(measurement_id);
    setUncertainties(measurement_id);
    setCorrelations(measurement_id);
    buildCov();
    buildPdf();
}


PDF_yCP::~PDF_yCP() {}


void PDF_yCP::initParameters() {
    ParametersCharmCombo p;
    parameters = new RooArgList("parameters");

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
            parameters->add(*(p.get("y12")));
            break;
        default:
            cout << "PDF_yCP::initParameters : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
}


void PDF_yCP::initRelations() {
    theory = new RooArgList("theory");
    switch (th_cfg) {
        case phenomenological:
            theory->add(
                    *(Utils::makeTheoryVar(
                            "yCP_th", "yCP_th",
                            "0.5*(  y * (qop+1 + 1/(qop+1)) * cos(phi)"
                            "     - x * (qop+1 - 1/(qop+1)) * sin(phi))",
                            parameters)));
            break;
        case theoretical:
            theory->add(
                    *(Utils::makeTheoryVar(
                            "yCP_th", "yCP_th", "y12*cos(phiG)",
                            parameters)));
            break;
        case superweak:
            theory->add(
                    *(Utils::makeTheoryVar(
                            "yCP_th", "yCP_th", "y12",
                            parameters)));
            break;
        default:
            cout << "PDF_yCP::initRelations : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
}


void PDF_yCP::initObservables(const TString& setName) {
    observables = new RooArgList("observables");
    observables->add(*(new RooRealVar("yCP_obs", setName + "   #it{y_{CP}}",  0, -1e4, 1e4)));
}


void PDF_yCP::setObservables(TString c) {
    if (c.EqualTo("truth")) setObservablesTruth();
    else if (c.EqualTo("toy")) setObservablesToy();
    else if (c.EqualTo("WA2020")) {
        obsValSource = "https://cds.cern.ch/record/2747731";
        setObservable("yCP_obs", -0.370);
    } else if (c.EqualTo("WA2020_biased")) {
        obsValSource = "HFLAV";
        setObservable("yCP_obs", 0.719);
    } else if (c.EqualTo("LHCb2022_biased")) {
        obsValSource = "https://inspirehep.net/literature/2035063";
        setObservable("yCP_obs", 0.696);
    } else {
        cout << "PDF_yCP::setObservables() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_yCP::setUncertainties(TString c) {
    if (c.EqualTo("WA2020")) {
        obsErrSource = "https://cds.cern.ch/record/2747731";
        StatErr[0] = 0.704;
        SystErr[0] = 0;
    } else if (c.EqualTo("WA2020_biased")) {
        obsValSource = "HFLAV";
        StatErr[0] = 0.113;
        SystErr[0] = 0;
    } else if (c.EqualTo("LHCb2022_biased")) {
        obsValSource = "https://inspirehep.net/literature/2035063";
        StatErr[0] = 0.026;
        SystErr[0] = 0.013;
    } else {
        cout << "PDF_yCP::setUncertainties() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_yCP::setCorrelations(TString c) {
    resetCorrelations();
    corSource = "No correlations for one observable";
}


void PDF_yCP::buildPdf() {
    pdf = new RooMultiVarGaussian(
            "pdf_" + name, "pdf_" + name, *(RooArgSet*)observables,
            *(RooArgSet*)theory, covMatrix);
}
