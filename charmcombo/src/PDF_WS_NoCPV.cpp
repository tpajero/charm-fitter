/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include "PDF_WS_NoCPV.h"


PDF_WS_NoCPV::PDF_WS_NoCPV(TString cObs, TString cErr, TString cCor,
                           const theory_config& th_cfg)
: PDF_Abs(3)
{
    name = cObs + "_WS_NoCPV";
    initParameters(th_cfg);
    initRelations(th_cfg);
    initObservables(cObs);
    setObservables(cObs);
    setUncertainties(cErr);
    setCorrelations(cCor);
    buildCov();
    buildPdf();
}


PDF_WS_NoCPV::~PDF_WS_NoCPV() {}


void PDF_WS_NoCPV::initParameters(const theory_config& th_cfg)
{
    ParametersCharmCombo p(th_cfg);
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
            parameters->add(*(p.get("phiM")));
            break;
        default:
            cout << "PDF_WS_NoCPV::initParameters : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
}


void PDF_WS_NoCPV::initRelations(const theory_config& th_cfg)
{
    RooArgSet *p = (RooArgSet*)parameters;
    theory = new RooArgList("theory");
    theory->add(*(new RooFormulaVar("RD_th", "RD_th", "R_Kpi", *p)));
    switch (th_cfg) {
        case phenomenological:
            theory->add(
                    *(new RooFormulaVar(
                            "yp_th", "yp_th",
                            "y*cos(Delta_Kpi) + x*sin(Delta_Kpi)", *p)));
            theory->add(
                    *(new RooFormulaVar(
                            "xp2_th", "xp2_th",
                            "pow(x*cos(Delta_Kpi) - y*sin(Delta_Kpi),2)",
                            *p)));
            break;
        case theoretical:
            theory->add(
                    *(new RooFormulaVar(
                            "yp_th", "yp_th",
                            "  y12 * cos(Delta_Kpi) * TMath::Sign(1.,cos(phiG)) "
                            "+ x12 * sin(Delta_Kpi) * TMath::Sign(1.,cos(phiM))", *p)));
            theory->add(
                    *(new RooFormulaVar(
                            "xp2_th", "xp2_th",
                            "pow(- y12*sin(Delta_Kpi) * TMath::Sign(1.,cos(phiG))"
                            "    + x12*cos(Delta_Kpi) * TMath::Sign(1.,cos(phiM)), 2)", *p)));
            break;
        case superweak:
            theory->add(
                    *(new RooFormulaVar(
                            "yp_th", "yp_th",
                            "  y12 * cos(Delta_Kpi)                            "
                            "+ x12 * sin(Delta_Kpi) * TMath::Sign(1.,cos(phiM))", *p)));
            theory->add(
                    *(new RooFormulaVar(
                            "xp2_th", "xp2_th",
                            "pow(- y12*sin(Delta_Kpi)                            "
                            "    + x12*cos(Delta_Kpi) * TMath::Sign(1.,cos(phiM)), 2)", *p)));
            break;
        default:
            cout << "PDF_WS_NoCPV::initRelations : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
}


void PDF_WS_NoCPV::initObservables(const TString& setName)
{
    observables = new RooArgList("observables"); ///< the order of this list must match that of the COR matrix!
    observables->add(*(new RooRealVar("RD_obs" , setName + "   #it{R_{K#pi}}",  0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("yp_obs" , setName + "   #it{y'}",  0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("xp2_obs", setName + "   #it{x'}^{2}",  0., -1e4, 1e4)));
}


void PDF_WS_NoCPV::setObservables(TString c)
{
    if (c.EqualTo("truth")) {
        setObservablesTruth();
    }
    else if (c.EqualTo("toy")) {
        setObservablesToy();
    }
    else if (c.EqualTo("CDF")) {
        obsValSource = "https://inspirehep.net/literature/1254229";
        setObservable("RD_obs" , 0.351);
        setObservable("yp_obs" , 0.43);
        setObservable("xp2_obs", 0.8);
    }
    else if (c.EqualTo("BaBar")) {
        obsValSource = "https://inspirehep.net/literature/746245";
        setObservable("RD_obs" , 0.303);
        setObservable("yp_obs" , 0.97);
        setObservable("xp2_obs", -2.2);
    }
    else if (c.EqualTo("Belle")) {
        obsValSource = "https://inspirehep.net/literature/1277238";
        setObservable("RD_obs" , 0.353);
        setObservable("yp_obs" , 0.46);
        setObservable("xp2_obs", 0.9);
    }
    else {
        cout << "PDF_WS_NoCPV::setObservables() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_WS_NoCPV::setUncertainties(TString c)
{
    if (c.EqualTo("CDF")) {
        obsErrSource = "https://inspirehep.net/literature/1254229";
        StatErr[0] = 0.035; // RD
        StatErr[1] = 0.43;  // y'
        StatErr[2] = 1.8;   // x'2
        SystErr[0] = 0;  // RD
        SystErr[1] = 0;  // y'
        SystErr[2] = 0;  // x'2
    } else if (c.EqualTo("BaBar")) {
        obsErrSource = "https://inspirehep.net/literature/746245";
        StatErr[0] = pow(pow(0.016,2) + pow(0.010,2),0.5);  // RD
        StatErr[1] = pow(pow(0.44, 2) + pow(0.31, 2),0.5);  // y'+
        StatErr[2] = pow(pow(3.0,  2) + pow(2.1,  2),0.5);  // x'2+
        SystErr[0] = 0;  // RD
        SystErr[1] = 0;  // y'+
        SystErr[2] = 0;  // x'2+
    } else if (c.EqualTo("Belle")) {
        obsErrSource = "https://inspirehep.net/literature/1277238";
        StatErr[0] = 0.013; // RD
        StatErr[1] = 0.34;  // y'
        StatErr[2] = 2.2;   // x'2
        SystErr[0] = 0;  // RD
        SystErr[1] = 0;  // y'
        SystErr[2] = 0;  // x'2
    }
    else {
        cout << "PDF_WS_NoCPV::setUncertainties() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_WS_NoCPV::setCorrelations(TString c)
{
    resetCorrelations();
    if (c.EqualTo("CDF")) {
        corSource = "https://inspirehep.net/literature/1254229";
        cout << "The correlation matrix of paper\n"
                "https://inspirehep.net/literature/1254229\n"
                "(WS/RS CDF) is not positive definite;\n"
                "thus, it has been modified to avoid non convergence.\n";
        double dataStat[] = {
            //  RD      y'     x'2
             1.000, -0.967,  0.900,
            -0.967,  1.000, -0.975,
             0.900, -0.975,  1.000
        };
        corStatMatrix = TMatrixDSym(nObs,dataStat);
        double dataSyst[] = {
            1.,  0.,  0.,
            0.,  1.,  0.,
            0.,  0.,  1.
        };
        corSystMatrix = TMatrixDSym(nObs,dataSyst);
    }
    else if (c.EqualTo("BaBar")) {
        corSource = "https://inspirehep.net/literature/746245";
        double dataStat[] = {
            //  RD      y'     x'2
             1.000, -0.87 ,  0.77 ,
            -0.87 ,  1.000, -0.94 ,
             0.77 , -0.94 ,  1.000
        };
        corStatMatrix = TMatrixDSym(nObs,dataStat);
        double dataSyst[] = {
            1.,  0.,  0.,
            0.,  1.,  0.,
            0.,  0.,  1.
        };
        corSystMatrix = TMatrixDSym(nObs,dataSyst);
    }
    else if (c.EqualTo("Belle")) {
        corSource = "https://inspirehep.net/literature/1277238";
        double dataStat[] = {
            //  RD      y'     x'2
             1.000, -0.865,  0.737,
            -0.865,  1.000, -0.948,
             0.737, -0.948,  1.000
        };
        corStatMatrix = TMatrixDSym(nObs,dataStat);
        double dataSyst[] = {
            1.,  0.,  0.,
            0.,  1.,  0.,
            0.,  0.,  1.
        };
        corSystMatrix = TMatrixDSym(nObs,dataSyst);
    }
    else {
        cout << "PDF_WS_NoCPV::setCorrelations() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_WS_NoCPV::buildPdf()
{
    pdf = new RooMultiVarGaussian(
            "pdf_" + name, "pdf_" + name, *(RooArgSet*)observables,
            *(RooArgSet*)theory, covMatrix);
}
