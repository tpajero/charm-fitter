/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include "PDF_Kshh.h"


PDF_Kshh::PDF_Kshh(TString cObs, TString cErr, TString cCor,
                   const theory_config& th_cfg)
: PDF_Abs(4)
{
    name = "Kshh_"+cObs;
    initParameters(th_cfg);
    initRelations(th_cfg);

    TString label = cObs;
    if (cObs == "Belle")
        label = "Belle #it{K}_{S}^{0}#it{#pi^{+}#pi^{#minus}}";

    initObservables(label);
    setObservables(cObs);
    setUncertainties(cErr);
    setCorrelations (cCor);
    buildCov();
    buildPdf();
}


PDF_Kshh::~PDF_Kshh() {}


void PDF_Kshh::initParameters(const theory_config& th_cfg)
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
            parameters->add(*(p.get("x12")));
            parameters->add(*(p.get("y12")));
            parameters->add(*(p.get("phiM")));
            break;
        default:
            cout << "PDF_Kshh::initParameters : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
}


void PDF_Kshh::initRelations(const theory_config& th_cfg)
{
    RooArgSet *p = (RooArgSet*)parameters;
    theory = new RooArgList("theory"); ///< the order of this list must match that of the COR matrix!
    switch (th_cfg) {
        case phenomenological:
            theory->add(*(new RooFormulaVar("x_th",   "x_th",   "x",     *p)));
            theory->add(*(new RooFormulaVar("y_th",   "y_th",   "y",     *p)));
            theory->add(*(new RooFormulaVar("qop_th", "qop_th", "qop+1", *p)));
            theory->add(*(new RooFormulaVar("phi_th", "phi_th", "phi",   *p)));
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
            theory->add(*(new RooFormulaVar(
                            "qop_th", "qop_th",
                            "pow(  (pow(x12,2) + pow(y12,2) + 2 * x12 * y12 * sin(phiM - phiG))"
                            "    /  pow(  pow(pow(x12,2) + pow(y12,2),2)                       "
                            "           - pow(2 * x12 * y12 * sin(phiM - phiG),2), 0.5), 0.5)  ", *p)));
            theory->add(
                    *(new RooFormulaVar(
                            "phi_th", "phi_th",
                            "-0.5 * TMath::ATan("
                            "      (pow(x12,2) * sin(2*phiM) + pow(y12,2) * sin(2*phiG))"
                            "    / (pow(x12,2) * cos(2*phiM) + pow(y12,2) * cos(2*phiG)))", *p)));
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
            theory->add(
                    *(new RooFormulaVar(
                            "qop_th", "qop_th",
                            "pow(  (pow(x12,2) + pow(y12,2) + 2 * x12 * y12 * sin(phiM))"
                            "    / pow(                                                      "
                            "        pow(pow(x12,2) + pow(y12,2),2)                      "
                            "      - pow(2 * x12 * y12 * sin(phiM),2), 0.5), 0.5)", *p)));
            theory->add(
                    *(new RooFormulaVar(
                            "phi_th", "phi_th",
                            "-0.5 * TMath::ATan("
                            "       pow(x12,2) * sin(2*phiM)"
                            "    / (pow(x12,2) * cos(2*phiM) + pow(y12,2)))", *p)));
            break;
        default:
            cout << "PDF_Kshh::initRelations : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
}


void PDF_Kshh::initObservables(const TString& setName)
{
    observables = new RooArgList("observables"); ///< the order of this list must match that of the COR matrix!
    observables->add(*(new RooRealVar("x_obs"  , setName + "   #it{x}"    ,  0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("y_obs"  , setName + "   #it{y}"    ,  0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("qop_obs", setName + "   |#it{q}/#it{p}|",  0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("phi_obs", setName + "   #it{#phi}_{2}" ,  0., -1e4, 1e4)));
}


void PDF_Kshh::setObservables(TString c)
{
    if (c.EqualTo("truth")) {
        setObservablesTruth();
    }
    else if (c.EqualTo("toy")) {
        setObservablesToy();
    }
    else if (c.EqualTo("Belle")) {
        obsValSource = "https://inspirehep.net/literature/1289224";
        setObservable("x_obs", 0.56);
        setObservable("y_obs", 0.30);
        setObservable("qop_obs", 0.90);
        setObservable("phi_obs", DegToRad(-6.));
    }
    else {
        cout << "PDF_Kshh::setObservables() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_Kshh::setUncertainties(TString c)
{
    if (c.EqualTo("Belle")) {
        obsErrSource = "https://inspirehep.net/literature/1289224";
        StatErr[0] =          sqrt(pow(0.19,  2) + pow(0.093, 2));  // x
        StatErr[1] =          sqrt(pow(0.15,  2) + pow(0.068, 2));  // y
        StatErr[2] =          sqrt(pow(0.155, 2) + pow(0.071, 2));  // qop
        StatErr[3] = DegToRad(sqrt(pow(11,    2) + pow(4.6,  2)));  // phi
        SystErr[0] = 0;  // x
        SystErr[1] = 0;  // y
        SystErr[2] = 0;  // qop
        SystErr[3] = 0;  // phi
    }
    else {
        cout << "PDF_Kshh::setUncertainties() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_Kshh::setCorrelations(TString c)
{
    resetCorrelations();
    if (c.EqualTo("Belle")) {
        corSource = "hflav";
        double dataStat[]  = {
        //   x       y       qop     phi
             1.000,  0.054, -0.074, -0.031,
             0.054,  1.000,  0.034, -0.019,
            -0.074,  0.034,  1.000,  0.044,
            -0.031, -0.019,  0.044,  1.000
        };
        corStatMatrix = TMatrixDSym(nObs,dataStat);
        double dataSyst[]  = {
        //  x    y    qop  phi
            1.,  0.,  0.,  0.,
            0.,  1.,  0.,  0.,
            0.,  0.,  1.,  0.,
            0.,  0.,  0.,  1.
        };
        corSystMatrix = TMatrixDSym(nObs,dataSyst);
    }
    else {
        cout << "PDF_Kshh::setCorrelations() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_Kshh::buildPdf()
{
  pdf = new RooMultiVarGaussian(
          "pdf_" + name, "pdf_" + name, *(RooArgSet*)observables,
          *(RooArgSet*)theory, covMatrix);
}
