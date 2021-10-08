/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include "PDF_XY.h"
#include "TMath.h"


PDF_XY::PDF_XY(TString cObs, TString cErr, TString cCor,
               const theory_config& th_cfg)
: PDF_Abs(2)
{
    name = "XY_"+cObs;

    TString label = cObs;
    if (cObs.EqualTo("BaBar_Kshh"))
        label = "BaBar #it{K}_{S}^{0}#it{h}^{+}#it{h}^{#minus}";
    else if (cObs.EqualTo("BaBar_pipipi0"))
        label = "BaBar #it{#pi}^{+}#it{#pi}^{#minus}#it{#pi}^{0}";
    else if (cObs.EqualTo("LHCb_KSpipi"))
        label = "LHCb #it{K}^{0}_{s}#it{#pi}^{+}#pi^{#minus}";
    initParameters(th_cfg);
    initRelations(th_cfg);
    initObservables(label);
    setObservables(cObs);
    setUncertainties(cErr);
    setCorrelations(cCor);
    buildCov();
    buildPdf();
}


PDF_XY::~PDF_XY() {}


void PDF_XY::initParameters(const theory_config& th_cfg)
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
            cout << "PDF_XY::initRelations : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
}


void PDF_XY::initRelations(const theory_config& th_cfg)
{
    RooArgSet *p = (RooArgSet*)parameters;
    theory = new RooArgList("theory"); ///< the order of this list must match that of the COR matrix!

    switch (th_cfg) {
        case phenomenological:
            theory->add(*(new RooFormulaVar("x_th", "x_th", "x", *p)));
            theory->add(*(new RooFormulaVar("y_th", "y_th", "y", *p)));
            break;
        case theoretical:
            theory->add(*(new RooFormulaVar(
                    "x_th", "x_th",
                    "pow(2,-0.5) * pow( "
                    "    pow(x12,2) - pow(y12,2) + pow( "
                    "       + pow(pow(x12,2) + pow(y12,2),2)"
                    "       - pow(2 * x12 * y12 * sin(phiM - phiG),2),"
                    "    0.5)"
                    ",0.5) * TMath::Sign(1., cos(phiM - phiG))",
                    *p)));
            theory->add(*(new RooFormulaVar(
                    "y_th", "y_th",
                    "pow(2,-0.5) * pow( "
                    "    pow(y12,2) - pow(x12,2) + pow( "
                    "       + pow(pow(x12,2) + pow(y12,2),2)"
                    "       - pow(2 * x12 * y12 * sin(phiM - phiG),2),"
                    "    0.5)"
                    ",0.5)",
                    *p)));
            break;
        case superweak:
            theory->add(*(new RooFormulaVar(
                    "x_th", "x_th",
                    "pow(2,-0.5) * pow( "
                    "    pow(x12,2) - pow(y12,2) + pow( "
                    "       + pow(pow(x12,2) + pow(y12,2),2)"
                    "       - pow(2 * x12 * y12 * sin(phiM),2),"
                    "    0.5)"
                    ",0.5) * TMath::Sign(1., cos(phiM))",
                    *p)));
            theory->add(*(new RooFormulaVar(
                    "y_th", "y_th",
                    "pow(2,-0.5) * pow( "
                    "    pow(y12,2) - pow(x12,2) + pow( "
                    "       + pow(pow(x12,2) + pow(y12,2),2)"
                    "       - pow(2 * x12 * y12 * sin(phiM),2),"
                    "    0.5)"
                    ",0.5)",
                    *p)));
            break;
        default:
            cout << "PDF_XY::initRelations : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
}


void PDF_XY::initObservables(const TString& setName)
{
    observables = new RooArgList("observables"); ///< the order of this list must match that of the COR matrix!
    observables->add(*(new RooRealVar("x_obs", setName + "   #it{x}",  0, -1e4, 1e4)));
    observables->add(*(new RooRealVar("y_obs", setName + "   #it{y}",  0, -1e4, 1e4)));
}


void PDF_XY::setObservables(TString c)
{
    if (c.EqualTo("truth")) {
        setObservablesTruth();
    }
    else if (c.EqualTo("toy")) {
        setObservablesToy();
    }
    else if (c.EqualTo("BaBar_Kshh")) {
        obsValSource = "https://inspirehep.net/literature/853279";
        setObservable("x_obs", 0.16);
        setObservable("y_obs", 0.57);
    }
    else if (c.EqualTo("BaBar_pipipi0")) {
        obsValSource = "https://inspirehep.net/literature/1441203";
        setObservable("x_obs", 1.5);
        setObservable("y_obs", 0.2);
    }
    else if (c.EqualTo("LHCb_KSpipi")) {
        obsValSource = "https://inspirehep.net/literature/1396327";
        setObservable("x_obs", -0.86);
        setObservable("y_obs", 0.03);
    }
    else {
        cout << "PDF_XY::setObservables() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_XY::setUncertainties(TString c)
{
    if (c.EqualTo("BaBar_Kshh")) {
        obsErrSource = "https://inspirehep.net/literature/853279";
        StatErr[0] = sqrt(pow(0.23,2) + pow(0.12,2) +pow(0.08,2)); // x
        StatErr[1] = sqrt(pow(0.20,2) + pow(0.13,2) +pow(0.07,2)); // y
        SystErr[0] = 0;
        SystErr[1] = 0;
    }
    else if (c.EqualTo("BaBar_pipipi0")) {
        obsErrSource = "https://inspirehep.net/literature/1441203";
        StatErr[0] = sqrt(pow(1.2,2) + pow(0.6,2)); // x
        StatErr[1] = sqrt(pow(0.9,2) + pow(0.5,2)); // y
        SystErr[0] = 0;
        SystErr[1] = 0;
    }
    else if (c.EqualTo("LHCb_KSpipi")) {
        obsErrSource = "https://inspirehep.net/literature/1396327";
        StatErr[0] = sqrt(pow(0.53,2) + pow(0.17,2)); // x
        StatErr[1] = sqrt(pow(0.46,2) + pow(0.13,2)); // y
        SystErr[0] = 0;
        SystErr[1] = 0;
    }
    else {
        cout << "PDF_XY::setUncertainties() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_XY::setCorrelations(TString c)
{
    resetCorrelations();
    if (c.EqualTo("BaBar_Kshh")) {
        corSource = "https://inspirehep.net/literature/853279";
        corStatMatrix[1][0] = 0.0586;
        corSystMatrix[1][0] = 0.0;
    }
    else if (c.EqualTo("BaBar_pipipi0")) {
        corSource = "https://inspirehep.net/literature/1441203";
        corStatMatrix[1][0] = -0.006;
        corSystMatrix[1][0] = 0.0;
    }
    else if (c.EqualTo("LHCb_KSpipi")) {
        corSource = "https://inspirehep.net/literature/1396327";
        corStatMatrix[1][0] = 0.37;
        corSystMatrix[1][0] = 0.0;
    }
    else {
        cout << "PDF_XY::setCorrelations() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_XY::buildPdf()
{
  pdf = new RooMultiVarGaussian("pdf_" + name, "pdf_" + name, *(RooArgSet*)observables, *(RooArgSet*)theory, covMatrix);
}
