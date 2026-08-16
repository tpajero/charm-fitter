/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include <PDF_Fp_pipipi0.h>

#include <CharmUtils.h>

#include <Utils.h>

#include <RooFormulaVar.h>
#include <RooRealVar.h>

#include <TString.h>

#include <algorithm>

PDF_Fp_pipipi0::PDF_Fp_pipipi0(const TString measurement_id) : PDF_Charm{1}, measurement_id{measurement_id} {
  name = "Fp-pipipi0" + measurement_id;
  initialise(measurement_id, measurement_id, measurement_id);
}

std::set<std::string> PDF_Fp_pipipi0::getParameterNames() const { return {"F_pipipi0"}; }

void PDF_Fp_pipipi0::initRelations() {
  theory = new RooArgList("theory");
  theory->add(*(Utils::makeTheoryVar("F_pipipi0_th", "F_pipipi0", parameters)));
}

void PDF_Fp_pipipi0::initObservables() {
  observables = new RooArgList("observables");
  observables->add(*(new RooRealVar("F_pipipi0_obs", "F_{#pi#pi#pi^{0}} " + measurement_id, 0, -1e4, 1e4)));
}

void PDF_Fp_pipipi0::setObservables(const TString c) {
  if (c.EqualTo("truth"))
    setObservablesTruth();
  else if (c.EqualTo("toy"))
    setObservablesToy();
  else if (c.EqualTo("Cleo-c")) {
    obsValSource = "https://inspirehep.net/literature/2139827";
    setObservable("F_pipipi0_obs", 0.973);
  } else if (c.EqualTo("BESIII")) {
    obsValSource = "https://inspirehep.net/literature/2827201";
    setObservable("F_pipipi0_obs", 0.9406);
  } else {
    std::cout << "PDF_Fp_pipipi0::setObservables() : ERROR : config " + c + " not found." << std::endl;
    exit(1);
  }
}

void PDF_Fp_pipipi0::setUncertainties(const TString c) {
  if (c.EqualTo("Cleo-c")) {
    obsErrSource = "https://inspirehep.net/literature/2139827";
    StatErr = {0.017};
    SystErr = {0.0};
  } else if (c.EqualTo("BESIII")) {
    obsErrSource = "https://inspirehep.net/literature/2827201";
    StatErr = {0.0036};
    SystErr = {0.0021};
  } else {
    std::cout << "PDF_Fp_pipipi0::setObservables() : ERROR : config " + c + " not found." << std::endl;
    exit(1);
  }
}

void PDF_Fp_pipipi0::setCorrelations(const TString c) {
  resetCorrelations();
  corSource = "No correlations for one observable";
}
