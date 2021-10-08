/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include "PDF_BinFlip.h"


PDF_BinFlip::PDF_BinFlip(TString c, TString err, TString corr,
                         const theory_config& th_cfg)
: PDF_Abs(4)
{
    if (err.IsNull()) err = c;
    if (corr.IsNull()) corr = c;
    name = "BinFlip_" + c;
    TString label;
    if (c.EqualTo("LHCb_Run1"))
        label = "LHCb Binflip Run 1";
    else if (c.EqualTo("LHCb_Run2"))
        label = "LHCb Binflip Run 2";
    initParameters(th_cfg);
    initRelations(th_cfg);
    initObservables(label);
    setObservables(c);
    setUncertainties(err);
    setCorrelations(corr);
    buildCov();
    buildPdf();
}


PDF_BinFlip::~PDF_BinFlip() {}


void PDF_BinFlip::initParameters(const theory_config& th_cfg)
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
            cout << "PDF_BinFlip::initParameters : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
}


void PDF_BinFlip::initRelations(const theory_config& th_cfg)
{
    RooArgSet *p = (RooArgSet*)parameters;
    theory = new RooArgList("theory"); ///< the order of this list must match that of the COR matrix!
    switch (th_cfg) {
        case phenomenological:
            theory->add(
                    *(new RooFormulaVar(
                            "x_th" , "x_th" ,
                            "0.5*(  x*cos(phi)*(qop+1 + 1/(qop+1))"
                            "     + y*sin(phi)*(qop+1 - 1/(qop+1)))", *p)));
            theory->add(
                    *(new RooFormulaVar(
                            "y_th" , "y_th" ,
                            "0.5*(  y*cos(phi)*(qop+1 + 1./(qop+1))"
                            "     - x*sin(phi)*(qop+1 - 1./(qop+1)))", *p)));
            theory->add(
                    *(new RooFormulaVar(
                            "dx_th", "dx_th",
                            "0.5*(  x*cos(phi)*(qop+1 - 1./(qop+1))"
                            "     + y*sin(phi)*(qop+1 + 1./(qop+1)))", *p)));
            theory->add(
                    *(new RooFormulaVar(
                            "dy_th", "dy_th",
                            "0.5*(  y*cos(phi)*(qop+1 - 1./(qop+1))"
                            "     - x*sin(phi)*(qop+1 + 1./(qop+1)))", *p)));
            break;
        case theoretical:
            theory->add(
                    *(new RooFormulaVar("x_th", "x_th", "x12*cos(phiM)", *p)));
            theory->add(
                    *(new RooFormulaVar("y_th", "y_th", "y12*cos(phiG)", *p)));
            theory->add(
                    *(new RooFormulaVar("dx_th", "dx_th", "-y12*sin(phiG)", *p)));
            theory->add(
                    *(new RooFormulaVar("dy_th", "dy_th", "x12*sin(phiM)", *p)));
            break;
        case superweak:
            theory->add(
                    *(new RooFormulaVar("x_th", "x_th", "x12*cos(phiM)", *p)));
            theory->add(
                    *(new RooFormulaVar("y_th", "y_th", "y12", *p)));
            theory->add(
                    *(new RooFormulaVar("dx_th", "dx_th", "0", *p)));
            theory->add(
                    *(new RooFormulaVar("dy_th", "dy_th", "x12*sin(phiM)", *p)));
            break;
        default:
            cout << "PDF_BinFlip::initRelations : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
}


void PDF_BinFlip::initObservables(const TString& setName)
{
    observables = new RooArgList("observables"); ///< the order of this list must match that of the COR matrix!
    observables->add(*(new RooRealVar("x_obs", setName + "   #it{x_{CP}}", 0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("y_obs", setName + "   #it{y_{CP}}", 0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("dx_obs", setName + "   #it{#Deltax}", 0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("dy_obs", setName + "   #it{#Deltay}", 0., -1e4, 1e4)));
}


void PDF_BinFlip::setObservables(TString c)
{
    if (c.EqualTo("truth")) {
        setObservablesTruth();
    }
    else if (c.EqualTo("toy")) {
        setObservablesToy();
    }
    else if (c.EqualTo("LHCb_Run1")) {
        obsValSource = "https://inspirehep.net/literature/1724179";
        setObservable("x_obs", 0.27);
        setObservable("y_obs", 0.74);
        setObservable("dx_obs", -0.053);
        setObservable("dy_obs", 0.06);
    }
    else if (c.EqualTo("LHCb_Run2")) {
        obsValSource = "https://inspirehep.net/literature/1867376";
        setObservable("x_obs", 0.3973);
        setObservable("y_obs", 0.4589);
        setObservable("dx_obs", -0.0271);
        setObservable("dy_obs", 0.0203);
    }
    else {
        cout << "PDF_BinFlip::setObservables() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_BinFlip::setUncertainties(TString c)
{
    if (c.EqualTo("LHCb_Run1")) {
        obsErrSource = "https://inspirehep.net/literature/1724179";
        StatErr[0] = 0.16;  // x
        StatErr[1] = 0.36;  // y
        StatErr[2] = 0.07;  // dx
        StatErr[3] = 0.16;  // dy
        SystErr[0] = 0.04;  // x
        SystErr[1] = 0.11;  // y
        SystErr[2] = 0.022; // dx
        SystErr[3] = 0.03;  // dy
    }
    else if (c.EqualTo("LHCb_Run2")) {
        obsErrSource = "https://inspirehep.net/literature/1867376";
        StatErr[0] = pow(pow(0.0459,2) + pow(0.029,2),0.5); // x
        StatErr[1] = pow(pow(0.1198,2) + pow(0.085,2),0.5); // y
        StatErr[2] = pow(pow(0.0182,2) + pow(0.001,2),0.5); // dx
        StatErr[3] = pow(pow(0.0365,2) + pow(0.011,2),0.5); // dy
        SystErr[0] = 0.; // x
        SystErr[1] = 0.; // y
        SystErr[2] = 0.; // dx
        SystErr[3] = 0.; // dy
    }
    else {
        cout << "PDF_BinFlip::setUncertainties() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_BinFlip::setCorrelations(TString c)
{
    resetCorrelations();
    if (c.EqualTo("LHCb_Run1")) {
        corSource = "https://inspirehep.net/literature/1724179";
        double dataStat[]  = {
             1.00, -0.17,  0.04, -0.02,  // x
            -0.17,  1.00, -0.03,  0.01,  // y
             0.04, -0.03,  1.00, -0.13,  // dx
            -0.02,  0.01, -0.13,  1.00   // dy
        };
        corStatMatrix = TMatrixDSym(nObs,dataStat);
        double dataSyst[]  = {
             1.00,  0.15,  0.01, -0.02,  // x
             0.15,  1.00, -0.05, -0.03,  // y
             0.01, -0.05,  1.00,  0.14,  // dx
            -0.02, -0.03,  0.14,  1.00   // dy
        };
        corSystMatrix = TMatrixDSym(nObs,dataSyst);
    }
    else if (c.EqualTo("LHCb_Run2")) {
        corSource = "https://inspirehep.net/literature/1867376";
        double dataStat[]  = {
             1.   ,  0.111,  -0.017, -0.010,  // x
             0.111,  1.   ,  -0.011, -0.051,  // y
            -0.017, -0.011,   1.   ,  0.077,  // dx
            -0.010, -0.051,  0.077,   1.      // dy
        };
        corStatMatrix = TMatrixDSym(nObs,dataStat);
        double dataSyst[]  = {
             1.,  0.,  0.,  0.,  // x
             0.,  1.,  0.,  0.,  // y
             0.,  0.,  1.,  0.,  // dx
             0.,  0.,  0.,  1.   // dy
        };
        corSystMatrix = TMatrixDSym(nObs,dataSyst);
    }
    else {
        cout << "PDF_BinFlip::setCorrelations() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_BinFlip::buildPdf()
{
  pdf = new RooMultiVarGaussian(
          "pdf_" + name, "pdf_" + name, *(RooArgSet*)observables,
          *(RooArgSet*)theory, covMatrix);
}
