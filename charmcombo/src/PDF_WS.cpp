/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include "PDF_WS.h"


PDF_WS::PDF_WS(TString cObs, TString cErr, TString cCor,
               const theory_config& th_cfg)
: PDF_Abs(6)
{
    TString label;
    if (cObs.EqualTo("LHCb_Prompt_2011_2016"))
        label = "WS/RS LHCb prompt";
    else if (cObs.EqualTo("LHCb_DT_Run1"))
        label = "WS/RS LHCb dt";
    else if (cObs.EqualTo("Belle"))
        label = "WS/RS Belle CPV";
    else
        exit(1);
    name = "WS_"+cObs;
    initParameters(th_cfg);
    initRelations(th_cfg);
    initObservables(label);
    setObservables(cObs);
    setUncertainties(cErr);
    setCorrelations(cCor);
    buildCov();
    buildPdf();
}


PDF_WS::~PDF_WS() {}


void PDF_WS::initParameters(const theory_config& th_cfg)
{
    ParametersCharmCombo p(th_cfg);
    parameters = new RooArgList("parameters");
    parameters->add(*(p.get("R_Kpi")));
    if (th_cfg != superweak)
        parameters->add(*(p.get("A_Kpi")));
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
            cout << "PDF_WS::initParameters : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
}


void PDF_WS::initRelations(const theory_config& th_cfg)
{
    RooArgSet *p = (RooArgSet*)parameters;
    theory = new RooArgList("theory");
    switch (th_cfg) {
        case phenomenological:
        case theoretical:
            theory->add(*(new RooFormulaVar("RD_p_th", "RD_p_th",
                                            "R_Kpi*(1+A_Kpi/100)", *p)));
            break;
        case superweak:
            theory->add(*(new RooFormulaVar("RD_p_th", "RD_p_th",
                                            "R_Kpi", *p)));
            break;
        default:
            cout << "PDF_WS::initRelations : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
    switch (th_cfg) {
        case phenomenological:
            theory->add(
                    *(new RooFormulaVar(
                            "yp_p_th", "yp_p_th",
                            "(qop+1)*(  y*cos(Delta_Kpi - phi)"
                            "         + x*sin(Delta_Kpi - phi))", *p)));
            theory->add(
                    *(new RooFormulaVar(
                            "xp2_p_th", "xp2_p_th",
                            "pow((qop+1)*(  x*cos(Delta_Kpi - phi)"
                            "             - y*sin(Delta_Kpi - phi)),2)", *p)));
            break;
        case theoretical:
            theory->add(
                    *(new RooFormulaVar(
                            "yp_p_th", "yp_p_th",
                            "  y12*cos(Delta_Kpi+phiG)"
                            "+ x12*sin(Delta_Kpi+phiM)",
                            *p)));
            theory->add(
                    *(new RooFormulaVar(
                            "xp2_p_th", "xp2_p_th",
                            "pow(- y12*sin(Delta_Kpi+phiG)"
                            "    + x12*cos(Delta_Kpi+phiM)"
                            "    ,2)", *p)));
            break;
        case superweak:
            theory->add(
                    *(new RooFormulaVar(
                            "yp_p_th", "yp_p_th",
                            "(  y12*cos(Delta_Kpi)"
                            " + x12*sin(Delta_Kpi+phiM))",
                            *p)));
            theory->add(
                    *(new RooFormulaVar(
                            "xp2_p_th", "xp2_p_th",
                            "pow((- y12*sin(Delta_Kpi)"
                            "     + x12*cos(Delta_Kpi+phiM))"
                            "    ,2)", *p)));
            break;
        default:
            cout << "PDF_WS::initRelations : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
    switch (th_cfg) {
        case phenomenological:
        case theoretical:
            theory->add(*(new RooFormulaVar("RD_m_th", "RD_m_th",
                                            "R_Kpi*(1-A_Kpi/100)", *p)));
            break;
        case superweak:
            theory->add(*(new RooFormulaVar("RD_m_th", "RD_m_th",
                                            "R_Kpi", *p)));
            break;
        default:
            cout << "PDF_WS::initRelations : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
    switch (th_cfg) {
        case phenomenological:
            theory->add(
                    *(new RooFormulaVar(
                            "yp_m_th", "yp_m_th",
                            "1/(qop+1)*(  y*cos(Delta_Kpi + phi)"
                            "           + x*sin(Delta_Kpi + phi))", *p)));
            theory->add(
                    *(new RooFormulaVar(
                            "xp2_m_th", "xp2_m_th",
                            "pow(1/(qop+1)*(  x*cos(Delta_Kpi + phi)"
                            "               - y*sin(Delta_Kpi + phi)),2)", *p)));
            break;
        case theoretical:
            theory->add(
                    *(new RooFormulaVar(
                            "yp_m_th", "yp_m_th",
                            "  y12*cos(Delta_Kpi-phiG)"
                            "+ x12*sin(Delta_Kpi-phiM)",
                            *p)));
            theory->add(
                    *(new RooFormulaVar(
                            "xp2_m_th", "xp2_m_th",
                            "pow(- y12*sin(Delta_Kpi-phiG)"
                            "    + x12*cos(Delta_Kpi-phiM)"
                            "    ,2)", *p)));
            break;
        case superweak:
            theory->add(
                    *(new RooFormulaVar(
                            "yp_m_th", "yp_m_th",
                            "(  y12*cos(Delta_Kpi)"
                            " + x12*sin(Delta_Kpi-phiM))",
                            *p)));
            theory->add(
                    *(new RooFormulaVar(
                            "xp2_m_th", "xp2_m_th",
                            "pow((- y12*sin(Delta_Kpi)"
                            "     + x12*cos(Delta_Kpi-phiM))"
                            "    ,2)", *p)));
            break;
        default:
            cout << "PDF_WS::initRelations : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
}


void PDF_WS::initObservables(const TString& setName)
{
    observables = new RooArgList("observables"); ///< the order of this list must match that of the COR matrix!
    observables->add(*(new RooRealVar("RD_p_obs" , setName + "   #it{R_{K#pi}^{+}}", 0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("yp_p_obs" , setName + "   #it{y'^{+}}", 0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("xp2_p_obs", setName + "   #it{x'^{+2}}", 0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("RD_m_obs" , setName + "   #it{R_{K#pi}^{#minus}}", 0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("yp_m_obs" , setName + "   #it{y'}^{#minus}", 0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("xp2_m_obs", setName + "   #it{x'}^{#minus2}", 0., -1e4, 1e4)));
}


void PDF_WS::setObservables(TString c)
{
    if (c.EqualTo("truth")) {
        setObservablesTruth();
    }
    else if (c.EqualTo("toy")) {
        setObservablesToy();
    }
    else if (c.EqualTo("LHCb_Prompt_2011_2016")) {
        obsValSource = "https://inspirehep.net/literature/1642234";
        setObservable("RD_p_obs" , 0.3454);
        setObservable("yp_p_obs" , 0.501);
        setObservable("xp2_p_obs", 0.61);
        setObservable("RD_m_obs" , 0.3454);
        setObservable("yp_m_obs" , 0.554);
        setObservable("xp2_m_obs", 0.16);
    }
    else if (c.EqualTo("LHCb_DT_Run1")) {
        obsValSource = "https://inspirehep.net/literature/1499047";
        setObservable("RD_p_obs" , 0.338);
        setObservable("yp_p_obs" , 0.581);
        setObservable("xp2_p_obs", -0.19);
        setObservable("RD_m_obs" , 0.360);
        setObservable("yp_m_obs" , 0.332);
        setObservable("xp2_m_obs", 0.79);
    }
    else if (c.EqualTo("Belle")) {
        obsValSource = "http://belle.kek.jp/belle/theses/doctor/lmzhang06/phd-mix-400.ps.gz";
        setObservable("RD_p_obs" ,  0.373);
        setObservable("yp_p_obs" , -0.12);
        setObservable("xp2_p_obs",  3.2);
        setObservable("RD_m_obs" ,  0.356);
        setObservable("yp_m_obs" ,  0.2);
        setObservable("xp2_m_obs",  0.6);
    }
    else {
        cout << "PDF_WS::setObservables() : ERROR : config " + c + " not found."
             << endl;
        exit(1);
    }
}


void PDF_WS::setUncertainties(TString c)
{
    if (c.EqualTo("LHCb_Prompt_2011_2016")) {
        obsErrSource = "https://inspirehep.net/literature/1642234";
        StatErr[0] = 0.0045;  // RD+
        StatErr[1] = 0.074;   // y'+
        StatErr[2] = 0.37;    // x'2+
        StatErr[3] = 0.0045;  // RD-
        StatErr[4] = 0.074;   // y'-
        StatErr[5] = 0.39;    // x'2-
        SystErr[0] = 0;  // RD+
        SystErr[1] = 0;  // y'+
        SystErr[2] = 0;  // x'2+
        SystErr[3] = 0;  // RD-
        SystErr[4] = 0;  // y'-
        SystErr[5] = 0;  // x'2-
    }
    else if (c.EqualTo("LHCb_DT_Run1")) {
        obsErrSource = "https://inspirehep.net/literature/1499047";
        StatErr[0] = pow(pow(0.15, 2) + pow(0.06, 2),0.5);  // RD+
        StatErr[1] = pow(pow(0.525,2) + pow(0.032,2),0.5);  // y'+
        StatErr[2] = pow(pow(4.46, 2) + pow(0.31, 2),0.5);  // x'2+
        StatErr[3] = pow(pow(0.15, 2) + pow(0.07, 2),0.5);  // RD-
        StatErr[4] = pow(pow(0.521,2) + pow(0.040,2),0.5);  // y'-
        StatErr[5] = pow(pow(4.31, 2) + pow(0.38, 2),0.5);  // x'2-
        SystErr[0] = 0;  // RD+
        SystErr[1] = 0;  // y'+
        SystErr[2] = 0;  // x'2+
        SystErr[3] = 0;  // RD-
        SystErr[4] = 0;  // y'-
        SystErr[5] = 0;  // x'2-
    }
    else if (c.EqualTo("Belle")) {
        obsErrSource = "http://belle.kek.jp/belle/theses/doctor/lmzhang06/phd-mix-400.ps.gz";
        StatErr[0] = 0.024;  // RD+
        StatErr[1] = 0.57;   // y'+
        StatErr[2] = 3.1;    // x'2+
        StatErr[3] = 0.024;  // RD-
        StatErr[4] = 0.54;   // y'-
        StatErr[5] = 2.9;    // x'2-
        SystErr[0] = 0;  // RD+
        SystErr[1] = 0;  // y'+
        SystErr[2] = 0;  // x'2+
        SystErr[3] = 0;  // RD-
        SystErr[4] = 0;  // y'-
        SystErr[5] = 0;  // x'2-
    }
    else {
        cout << "PDF_WS::setUncertainties() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}


void PDF_WS::setCorrelations(TString c)
{
    resetCorrelations();
    if (c.EqualTo("LHCb_Prompt_2011_2016")) {
        corSource = "https://inspirehep.net/literature/1642234";
        double dataStat[] = {
            //  RD+     y'+    x'2+     RD-     y'-    x'2-
             1.000, -0.935,  0.843, -0.012, -0.003, -0.002,  // RD+
            -0.935,  1.000, -0.963, -0.003,  0.004, -0.003,  // y'+
             0.843, -0.963,  1.000,  0.002, -0.003,  0.003,  // x'2+
            -0.012, -0.003,  0.002,  1.000, -0.935,  0.846,  // RD-
            -0.003,  0.004, -0.003, -0.935,  1.000, -0.964,  // y'-
            -0.002, -0.003,  0.003,  0.846, -0.964,  1.000   // x'2-
        };
        corStatMatrix = TMatrixDSym(nObs,dataStat);
        double dataSyst[] = {
        //   RD+  y'+  x'2+ RD-  y'-  x'2-
             1.,  0.,  0.,  0.,  0.,  0.,  // RD+
             0.,  1.,  0.,  0.,  0.,  0.,  // y'+
             0.,  0.,  1.,  0.,  0.,  0.,  // x'2+
             0.,  0.,  0.,  1.,  0.,  0.,  // RD-
             0.,  0.,  0.,  0.,  1.,  0.,  // y'-
             0.,  0.,  0.,  0.,  0.,  1.   // x'2-
        };
        corSystMatrix = TMatrixDSym(nObs,dataSyst);
    }
    else if (c.EqualTo("LHCb_DT_Run1")) {
        corSource = "https://inspirehep.net/literature/1499047";
        double dataStat[] = {
            //  RD+     y'+    x'2+     RD-     y'-    x'2-
             1.000, -0.658,  0.043, -0.005,  0.000,  0.000,  // RD+
            -0.658,  1.000,  0.438, -0.001,  0.000, -0.001,  // y'+
             0.043,  0.438,  1.000,  0.000,  0.000, -0.002,  // x'2+
            -0.005, -0.001,  0.000,  1.000, -0.621,  0.074,  // RD-
             0.000,  0.000,  0.000, -0.621,  1.000,  0.050,  // y'-
             0.000, -0.001, -0.002,  0.074,  0.050,  1.000   // x'2-
        };
        corStatMatrix = TMatrixDSym(nObs,dataStat);
        double dataSyst[] = {
        //   RD+  y'+  x'2+ RD-  y'-  x'2-
             1.,  0.,  0.,  0.,  0.,  0.,  // RD+
             0.,  1.,  0.,  0.,  0.,  0.,  // y'+
             0.,  0.,  1.,  0.,  0.,  0.,  // x'2+
             0.,  0.,  0.,  1.,  0.,  0.,  // RD-
             0.,  0.,  0.,  0.,  1.,  0.,  // y'-
             0.,  0.,  0.,  0.,  0.,  1.   // x'2-
        };
        corSystMatrix = TMatrixDSym(nObs,dataSyst);
    }
    else if (c.EqualTo("Belle")) {
        corSource = "http://belle.kek.jp/belle/theses/doctor/lmzhang06/phd-mix-400.ps.gz";
        double dataStat[] = {
            //  RD+     y'+    x'2+     RD-     y'-    x'2-
             1.000, -0.834,  0.655,  0.000,  0.000,  0.000,  // RD+
            -0.834,  1.000, -0.909,  0.000,  0.000,  0.000,  // y'+
             0.655, -0.909,  1.000,  0.000,  0.000,  0.000,  // x'2+
             0.000,  0.000,  0.000,  1.000, -0.834,  0.655,  // RD-
             0.000,  0.000,  0.000, -0.834,  1.000, -0.909,  // y'-
             0.000,  0.000,  0.000,  0.655, -0.909,  1.000   // x'2-
        };
        corStatMatrix = TMatrixDSym(nObs,dataStat);
        double dataSyst[] = {
        //   RD+  y'+  x'2+ RD-  y'-  x'2-
             1.,  0.,  0.,  0.,  0.,  0.,  // RD+
             0.,  1.,  0.,  0.,  0.,  0.,  // y'+
             0.,  0.,  1.,  0.,  0.,  0.,  // x'2+
             0.,  0.,  0.,  1.,  0.,  0.,  // RD-
             0.,  0.,  0.,  0.,  1.,  0.,  // y'-
             0.,  0.,  0.,  0.,  0.,  1.   // x'2-
        };
        corSystMatrix = TMatrixDSym(nObs,dataSyst);
    }
    else {
        cout << "PDF_WS::setCorrelations() : ERROR : config " + c
                + " not found." << endl;
        exit(1);
    }
}


void PDF_WS::buildPdf()
{
    pdf = new RooMultiVarGaussian(
            "pdf_" + name, "pdf_" + name, *(RooArgSet*)observables,
            *(RooArgSet*)theory, covMatrix);
}
