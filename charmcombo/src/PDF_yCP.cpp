/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include "PDF_yCP.h"


PDF_yCP::PDF_yCP(TString cObs, TString cErr, TString cCor,
                 const theory_config& th_cfg)
: PDF_Abs(1)
{
    name = "yCP_" + cObs;
    initParameters(th_cfg);
    initRelations(th_cfg);
    initObservables(cObs);
    setObservables(cObs);
    setUncertainties(cErr);
    setCorrelations(cCor);
    buildCov();
    buildPdf();
}


PDF_yCP::~PDF_yCP() {}


void PDF_yCP::initParameters(const theory_config& th_cfg)
{
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


void PDF_yCP::initRelations(const theory_config& th_cfg)
{
    theory = new RooArgList("theory");
    switch (th_cfg) {
        case phenomenological:
            theory->add(
                    *(new RooFormulaVar(
                            "yCP_th", "yCP_th",
                            "0.5*(  y * (qop+1 + 1/(qop+1)) * cos(phi)"
                            "     - x * (qop+1 - 1/(qop+1)) * sin(phi))",
                            *(RooArgSet*)parameters)));
            break;
        case theoretical:
            theory->add(
                    *(new RooFormulaVar(
                            "yCP_th", "yCP_th", "y12*cos(phiG)",
                            *(RooArgSet*)parameters)));
            break;
        case superweak:
            theory->add(
                    *(new RooFormulaVar(
                            "yCP_th", "yCP_th", "y12",
                            *(RooArgSet*)parameters)));
            break;
        default:
            cout << "PDF_yCP::initRelations : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
}


void PDF_yCP::initObservables(const TString& setName)
{
    observables = new RooArgList("observables");
    observables->add(*(new RooRealVar("yCP_obs", setName + "   #it{y_{CP}}",  0, -1e4, 1e4)));
}


void PDF_yCP::setObservables(TString c)
{
    if (c.EqualTo("truth")) {
        setObservablesTruth();
    }
    else if (c.EqualTo("toy")) {
        setObservablesToy();
    }
    else if (c.EqualTo("WA2020")) {
        obsValSource = "https://cds.cern.ch/record/2747731";
        setObservable("yCP_obs", -0.370);
    }
    else if (c.EqualTo("WA2020_biased")) {
        obsValSource = "HFLAV";
        setObservable("yCP_obs", 0.719);
    }
    else if (c.EqualTo("LHCb2022_biased")) {
        obsValSource = "https://inspirehep.net/literature/2035063";
        setObservable("yCP_obs", 0.696);
    }
    else {
        cout << "PDF_yCP::setObservables() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_yCP::setUncertainties(TString c)
{
    if (c.EqualTo("WA2020")) {
        obsErrSource = "https://cds.cern.ch/record/2747731";
        StatErr[0] = 0.704;
        SystErr[0] = 0;
    }
    else if (c.EqualTo("WA2020_biased")) {
        obsValSource = "HFLAV";
        StatErr[0] = 0.113;
        SystErr[0] = 0;
    }
    else if (c.EqualTo("LHCb2022_biased")) {
        obsValSource = "https://inspirehep.net/literature/2035063";
        StatErr[0] = 0.026;
        SystErr[0] = 0.013;
    }
    else {
        cout << "PDF_yCP::setUncertainties() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_yCP::setCorrelations(TString c)
{
    resetCorrelations();
    corSource = "No correlations for one observable";
}


void PDF_yCP::buildPdf()
{
    pdf = new RooMultiVarGaussian(
            "pdf_" + name, "pdf_" + name, *(RooArgSet*)observables,
            *(RooArgSet*)theory, covMatrix);
}
