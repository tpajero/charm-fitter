/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include "CharmUtils.h"
#include "ParametersCharmCombo.h"
#include "PDF_Fp_pipipi0.h"

PDF_Fp_pipipi0::PDF_Fp_pipipi0() : PDF_Abs{1} {
    name = "Cleo-Fp-pipipi0";
    initParameters();
    initRelations();
    initObservables();
    setObservables("Cleo-c");
    setUncertainties("Cleo-c");
    setCorrelations("Cleo-c");
    buildCov();
    buildPdf();
}

PDF_Fp_pipipi0::~PDF_Fp_pipipi0(){}

void PDF_Fp_pipipi0::initParameters() {
    ParametersCharmCombo p;
    parameters = new RooArgList("parameters");
    parameters->add(*(p.get("F_pipipi0")));
}

void PDF_Fp_pipipi0::initRelations() {
    RooArgSet *p = (RooArgSet*)parameters;
    theory = new RooArgList("theory");
    theory->add(*(new RooFormulaVar("F_pipipi0_th", "F_pipipi0_th", "F_pipipi0", *p)));
}

void PDF_Fp_pipipi0::initObservables() {
    observables = new RooArgList("observables");
    observables->add(*(new RooRealVar("F_pipipi0_obs", "F_pipipi0",  0, -1e4, 1e4)));
}

void PDF_Fp_pipipi0::setObservables(TString c) {
    if (c.EqualTo("Cleo-c")) {
        obsValSource = "https://inspirehep.net/literature/2139827";
        setObservable("F_pipipi0_obs", 0.973);
    } else {
        cout << "PDF_Fp_pipipi0::setObservables() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}

void PDF_Fp_pipipi0::setUncertainties(TString c) {
    if (c.EqualTo("Cleo-c")) {
        obsErrSource = "https://inspirehep.net/literature/2139827";
        StatErr[0] = 0.017;
        SystErr[0] = 0.;
    } else {
        cout << "PDF_Fp_pipipi0::setObservables() : ERROR : config " + c + " not found." << endl;
        exit(1);
    }
}

void PDF_Fp_pipipi0::setCorrelations(TString c) {
    resetCorrelations();
    corSource = "No correlations for one observable";
}

void PDF_Fp_pipipi0::buildPdf()
{
    pdf = new RooMultiVarGaussian("pdf_" + name, "pdf_" + name, *(RooArgSet*)observables, *(RooArgSet*)theory, covMatrix);
}
