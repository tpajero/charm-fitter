/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include "PDF_Kpipi0.h"


PDF_Kpipi0::PDF_Kpipi0(TString cObs, TString cErr, TString cCor,
                       const theory_config& th_cfg)
: PDF_Abs(2)
{
    name = cObs + "_Kpipi0";
    TString label;
    if ( cObs.EqualTo("BaBar") )
        label = "BaBar #it{K}^{+}#pi^{#minus}#pi^{0}";
    else
        exit(1);
    initParameters(th_cfg);
    initRelations(th_cfg);
    initObservables(label);
    setObservables(cObs);
    setUncertainties(cErr);
    setCorrelations(cCor);
    buildCov();
    buildPdf();
}


PDF_Kpipi0::~PDF_Kpipi0() {}


void PDF_Kpipi0::initParameters(const theory_config& th_cfg)
{
    ParametersCharmCombo p(th_cfg);
    parameters = new RooArgList("parameters");
    parameters->add(*(p.get("Delta_Kpipi")));

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
            cout << "PDF_Kpipi0::initParameters : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
}


void PDF_Kpipi0::initRelations(const theory_config& th_cfg)
{
    RooArgSet *p = (RooArgSet*)parameters;
    theory = new RooArgList("theory");
    switch (th_cfg) {
        case phenomenological:  // CP violation in the decay is missing!
                                // check the CP violation in x'
            theory->add(
                    *(new RooFormulaVar(
                            "xpp_th", "xpp_th",
                            "x*cos(Delta_Kpipi) - y*sin(Delta_Kpipi)",
                            *p)));
            theory->add(
                    *(new RooFormulaVar(
                            "ypp_th", "ypp_th",
                            "y*cos(Delta_Kpipi) + x*sin(Delta_Kpipi)", *p)));
            break;
        case theoretical:
            theory->add(
                    *(new RooFormulaVar(
                            "xpp_th", "xpp_th",
                            "  x12*cos(Delta_Kpipi) * TMath::Sign(1.,cos(phiM))"
                            "- y12*sin(Delta_Kpipi) * TMath::Sign(1.,cos(phiG))", *p)));
            theory->add(
                    *(new RooFormulaVar(
                            "ypp_th", "ypp_th",
                            "  y12 * cos(Delta_Kpipi) * TMath::Sign(1.,cos(phiG))"
                            "+ x12 * sin(Delta_Kpipi) * TMath::Sign(1.,cos(phiM))", *p)));
            break;
        case superweak:
            theory->add(
                    *(new RooFormulaVar(
                            "xpp_th", "xpp_th",
                            "  x12*cos(Delta_Kpipi) * TMath::Sign(1.,cos(phiM))"
                            "- y12*sin(Delta_Kpipi)                            ", *p)));
            theory->add(
                    *(new RooFormulaVar(
                            "ypp_th", "ypp_th",
                            "  y12 * cos(Delta_Kpipi)                            "
                            "+ x12 * sin(Delta_Kpipi) * TMath::Sign(1.,cos(phiM))", *p)));
            break;
        default:
            cout << "PDF_Kpipi0::initRelations : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
}


void PDF_Kpipi0::initObservables(const TString& setName)
{
    observables = new RooArgList("observables"); ///< the order of this list must match that of the COR matrix!
    observables->add(*(new RooRealVar("xpp_obs", setName + "   #it{x''}",  0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("ypp_obs", setName + "   #it{y''}",  0., -1e4, 1e4)));
}


void PDF_Kpipi0::setObservables(TString c)
{
    if (c.EqualTo("truth")) {
        setObservablesTruth();
    }
    else if (c.EqualTo("toy")) {
        setObservablesToy();
    }
    else if (c.EqualTo("BaBar")) {
        obsValSource = "https://inspirehep.net/literature/791715";
        setObservable("xpp_obs", 2.61);
        setObservable("ypp_obs", -0.06);
    }
    else {
        cout << "PDF_Kpipi0::setObservables() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_Kpipi0::setUncertainties(TString c)
{
    if (c.EqualTo("BaBar")) {
        obsErrSource = "https://inspirehep.net/literature/791715";
        StatErr[0] = pow(pow(0.625,2) + pow(0.39,2), 0.5);  // x''
        StatErr[1] = pow(pow(0.595,2) + pow(0.34,2), 0.5);  // y''
    }
    else {
        cout << "PDF_Kpipi0::setUncertainties() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_Kpipi0::setCorrelations(TString c)
{
    resetCorrelations();
    if (c.EqualTo("BaBar")) {
        obsErrSource = "https://inspirehep.net/literature/791715";
        double dataStat[] = {
             1.00, -0.75,
            -0.75,  1.00
        };
        corStatMatrix = TMatrixDSym(nObs,dataStat);
        double dataSyst[] = {
            1.00,  0.00,
            0.00,  1.00
        };
        corSystMatrix = TMatrixDSym(nObs,dataSyst);
    }
    else {
        cout << "PDF_Kpipi0::setCorrelations() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_Kpipi0::buildPdf()
{
    pdf = new RooMultiVarGaussian(
            "pdf_" + name, "pdf_" + name, *(RooArgSet*)observables,
            *(RooArgSet*)theory, covMatrix);
}
