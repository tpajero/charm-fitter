/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include "PDF_RM.h"


PDF_RM::PDF_RM(TString cObs, TString cErr, TString cCor,
               const theory_config& th_cfg)
: PDF_Abs(1)
{
    name = "RM_" + cObs;
    initParameters(th_cfg);
    initRelations(th_cfg);
    initObservables(cObs);
    setObservables(cObs);
    setUncertainties(cErr);
    setCorrelations(cCor);
    buildCov();
    buildPdf();
}


PDF_RM::~PDF_RM() {}


void PDF_RM::initParameters(const theory_config& th_cfg)
{
    ParametersCharmCombo p(th_cfg);
    parameters = new RooArgList("parameters");

    switch (th_cfg) {
        case phenomenological:
            parameters->add(*(p.get("x")));
            parameters->add(*(p.get("y")));
            break;
        case theoretical:
            parameters->add(*(p.get("x12")));
            parameters->add(*(p.get("y12")));
            parameters->add(*(p.get("phiM")));
            parameters->add(*(p.get("phiG")));
            break;
        case superweak:
            parameters->add(*(p.get("x12")));
            parameters->add(*(p.get("y12")));
            parameters->add(*(p.get("phiM")));
            break;
        default:
            cout << "PDF_RM::initParameters : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
}


void PDF_RM::initRelations(const theory_config& th_cfg)
{
    theory = new RooArgList("theory");
    switch (th_cfg) {
        case phenomenological:
            theory->add(
                    *(new RooFormulaVar(
                            "RM_th", "RM_th",
                            "(pow(x,2) + pow(y,2))/2",
                            *(RooArgSet*)parameters)));
            break;
        case theoretical:
            theory->add(
                    *(new RooFormulaVar(
                            "RM_th", "RM_th",
                            "0.5 * pow( "
                            "    + pow(pow(x12,2) + pow(y12,2),2)"
                            "    - pow(2 * x12 * y12 * sin(phiM - phiG),2), 0.5)",
                            *(RooArgSet*)parameters)));
            break;
        case superweak:
            theory->add(
                    *(new RooFormulaVar(
                            "RM_th", "RM_th",
                            "0.5 * pow( "
                            "    + pow(pow(x12,2) + pow(y12,2),2)"
                            "    - pow(2 * x12 * y12 * sin(phiM),2), 0.5)",
                            *(RooArgSet*)parameters)));
            break;
        default:
            cout << "PDF_RM::initRelations : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
}


void PDF_RM::initObservables(const TString& setName)
{
    observables = new RooArgList("observables");
    observables->add(*(new RooRealVar("RM_obs", setName + "   #it{R_{M}}",  0, 0., 1e4)));
}


void PDF_RM::setObservables(TString c)
{
    if (c.EqualTo("truth")) {
        setObservablesTruth();
    }
    else if (c.EqualTo("toy")) {
        setObservablesToy();
    }
    else if (c.EqualTo("HFLAV2016")) {
        obsValSource = "https://hflav-eos.web.cern.ch/hflav-eos/charm/CHARM21/results_mixing.html";
        setObservable("RM_obs", 1.30);
    }
    else if (c.EqualTo("LHCb_K3pi_Run1")) {
        obsValSource = "https://inspirehep.net/literature/1423070";
        setObservable("RM_obs", 2 * 0.48);
    }
    else {
        cout << "PDF_RM::setObservables() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_RM::setUncertainties(TString c)
{
    if (c.EqualTo("HFLAV2016")) {
        obsErrSource = "https://hflav-eos.web.cern.ch/hflav-eos/charm/CHARM21/results_mixing.html";
        StatErr[0] = 2.69;
        SystErr[0] = 0;
    }
    else if (c.EqualTo("LHCb_K3pi_Run1")) {
        obsErrSource = "https://inspirehep.net/literature/1423070";
        StatErr[0] = 2 * 0.18;
        SystErr[0] = 0;
    }
    else {
        cout << "PDF_RM::setUncertainties() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_RM::setCorrelations(TString c)
{
    resetCorrelations();
    corSource = "No correlations for one observable";
}


void PDF_RM::buildPdf()
{
    pdf = new RooMultiVarGaussian(
            "pdf_" + name, "pdf_" + name, *(RooArgSet*)observables,
            *(RooArgSet*)theory, covMatrix);
}
