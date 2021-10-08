/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include "PDF_BES.h"


PDF_BES::PDF_BES(TString cObs, TString cErr, TString cCor,
                 const theory_config& th_cfg)
: PDF_Abs(1)
{
    name = "BES";
    initParameters(th_cfg);
    initRelations(th_cfg);
    initObservables(name);
    setObservables(cObs);
    setUncertainties(cErr);
    setCorrelations (cCor);
    buildCov();
    buildPdf();
}


PDF_BES::~PDF_BES() {}


void PDF_BES::initParameters(const theory_config& th_cfg)
{
    ParametersCharmCombo p(th_cfg);
    parameters = new RooArgList("parameters");
    parameters->add(*(p.get("Delta_Kpi")));
    parameters->add(*(p.get("R_Kpi")));

    switch (th_cfg) {
        case phenomenological:
            parameters->add(*(p.get("x")));
            parameters->add(*(p.get("y")));
            break;
        case theoretical:
            parameters->add(*(p.get("phiG")));
        case superweak:
            parameters->add(*(p.get("x12")));
            parameters->add(*(p.get("y12")));
            parameters->add(*(p.get("phiM")));
            break;
        default:
            cout << "PDF_BES::initParameters : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
}


void PDF_BES::initRelations(const theory_config& th_cfg)
{
    RooArgSet *p = (RooArgSet*)parameters;
    theory = new RooArgList("theory");
    switch (th_cfg) {
        case phenomenological:
            theory->add(
                    *(new RooFormulaVar(
                            "asymm_th", "asymm_th",
                            "(2 * 10 * sqrt(R_Kpi) * cos(Delta_Kpi) + y) / "
                            "(1 + R_Kpi/100 "
                            "   + sqrt(R_Kpi/100)*(y*cos(Delta_Kpi) + x*sin(Delta_Kpi))/100"
                            "   + (pow(x/100,2) + pow(y/100,2))/2)", *p)));
            break;
        case theoretical:
            theory->add(
                    *(new RooFormulaVar(
                            "asymm_th", "asymm_th",
                            "(2 * 10 * sqrt(R_Kpi) * cos(Delta_Kpi)"
                            " + pow(2,-0.5) * pow( "
                            "      pow(y12,2) - pow(x12,2) + pow( "
                            "         + pow(pow(x12,2) + pow(y12,2),2)"
                            "         - pow(2*x12*y12*sin(phiM - phiG),2),"
                            "      0.5)"
                            "   ,0.5)) / "
                            "(1 + R_Kpi/100 "
                            "   + sqrt(R_Kpi/100)*(  y12*cos(Delta_Kpi)*TMath::Sign(1.,cos(phiG))"
                            "                      + x12*sin(Delta_Kpi)*TMath::Sign(1.,cos(phiM)))/100"
                            "   + sqrt(  pow(pow(x12,2) + pow(y12,2),2)"
                            "          - pow(2*x12*y12*sin(phiM - phiG),2))/2e4)", *p)));
            break;
        case superweak:
            theory->add(
                    *(new RooFormulaVar(
                            "asymm_th", "asymm_th",
                            "(2 * 10 * sqrt(R_Kpi) * cos(Delta_Kpi)"
                            " + pow(2,-0.5) * pow( "
                            "      pow(y12,2) - pow(x12,2) + pow( "
                            "         + pow(pow(x12,2) + pow(y12,2),2)"
                            "         - pow(2*x12*y12*sin(phiM),2),"
                            "      0.5)"
                            "   ,0.5)) / "
                            "(1 + R_Kpi/100 "
                            "   + sqrt(R_Kpi/100)*(  y12*cos(Delta_Kpi) "
                            "                      + x12*sin(Delta_Kpi)*TMath::Sign(1.,cos(phiM)))/100"
                            "   + sqrt(  pow(pow(x12,2) + pow(y12,2),2)"
                            "          - pow(2*x12*y12*sin(phiM),2))/2e4)", *p)));
            break;
        default:
            cout << "PDF_BES::initRelations : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
}


void PDF_BES::initObservables(const TString& setName)
{
    observables = new RooArgList("observables");
    observables->add(*(new RooRealVar(
            "asymm_obs", setName + "   #it{A_{K#pi}^{CP}}", 0., -1e4, 1e4)));
}


void PDF_BES::setObservables(TString c)
{
    if (c.EqualTo("truth")) {
        setObservablesTruth();
    }
    else if (c.EqualTo("toy")) {
        setObservablesToy();
    }
    else if (c.EqualTo("BES")) {
        obsValSource = "http://inspirehep.net/record/1291279";
        setObservable("asymm_obs", 12.7);
    }
    else {
        cout << "PDF_BES::setObservables() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_BES::setUncertainties(TString c)
{
    if (c.EqualTo("BES")) {
        obsErrSource = "http://inspirehep.net/record/1291279";
        StatErr[0] = 1.3;
        SystErr[0] = 0.7;
    }
    else {
        cout << "PDF_BES::setUncertainties() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_BES::setCorrelations(TString c)
{
    resetCorrelations();
    corSource = "No correlations for one observable";
}


void PDF_BES::buildPdf()
{
  pdf = new RooMultiVarGaussian(
          "pdf_" + name, "pdf_" + name, *(RooArgSet*)observables,
          *(RooArgSet*)theory, covMatrix);
}
