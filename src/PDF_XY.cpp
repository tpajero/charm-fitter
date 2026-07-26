/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include <PDF_XY.h>

#include <CharmUtils.h>
#include <ParametersCharmCombo.h>

#include <Utils.h>

#include <RooFormulaVar.h>
#include <RooMultiVarGaussian.h>
#include <RooRealVar.h>

#include <TString.h>

#include <algorithm>
#include <iostream>

PDF_XY::PDF_XY(const TString measurement_id, const theory_config th_cfg) : PDF_Abs{2}, th_cfg{th_cfg} {
  name = "XY_" + measurement_id;

  TString label = measurement_id;
  if (measurement_id.EqualTo("BaBar_Kshh"))
    label = "BaBar #it{K}_{S}^{0}#it{h}^{+}#it{h}^{#minus}";
  else if (measurement_id.EqualTo("BaBar_pipipi0"))
    label = "BaBar #it{#pi}^{+}#it{#pi}^{#minus}#it{#pi}^{0}";
  else if (measurement_id.EqualTo("LHCb_KSpipi"))
    label = "LHCb #it{K}^{0}_{s}#it{#pi}^{+}#pi^{#minus}";
  else if (measurement_id.EqualTo("Belle_Belle2"))
    label = "Belle 1+2 #it{K}^{0}_{s}#it{#pi}^{+}#pi^{#minus}";
  initParameters();
  initRelations();
  initObservables(label);
  setObservables(measurement_id);
  setUncertainties(measurement_id);
  setCorrelations(measurement_id);
  build();
}

void PDF_XY::initParameters() {
  ParametersCharmCombo p;
  parameters = new RooArgList("parameters");
  switch (th_cfg) {
  case theory_config::phenomenological:
    parameters->add(*(p.get("x")));
    parameters->add(*(p.get("y")));
    break;
  case theory_config::theoretical:
    parameters->add(*(p.get("x12")));
    parameters->add(*(p.get("y12")));
    parameters->add(*(p.get("phiM")));
    parameters->add(*(p.get("phiG")));
    break;
  default:
    std::cout << "PDF_XY::initRelations : ERROR : "
                 "theory_config not supported."
              << std::endl;
    exit(1);
  }
}

void PDF_XY::initRelations() {
  theory = new RooArgList("theory");  ///< the order of this list must match that of the COR matrix!

  switch (th_cfg) {
  case theory_config::phenomenological:
    theory->add(*(Utils::makeTheoryVar("x_th", "x_th", "x", parameters)));
    theory->add(*(Utils::makeTheoryVar("y_th", "y_th", "y", parameters)));
    break;
  case theory_config::theoretical:
    theory->add(*(Utils::makeTheoryVar("x_th", "x_th", CharmUtils::x_to_theoretical, parameters)));
    theory->add(*(Utils::makeTheoryVar("y_th", "y_th", CharmUtils::y_to_theoretical, parameters)));
    break;
  default:
    std::cout << "PDF_XY::initRelations : ERROR : "
                 "theory_config not supported."
              << std::endl;
    exit(1);
  }
}

void PDF_XY::initObservables(const TString setName) {
  observables = new RooArgList("observables");  ///< the order of this list must match that of the COR matrix!
  observables->add(*(new RooRealVar("x_obs", setName + "   #it{x}", 0, -1e4, 1e4)));
  observables->add(*(new RooRealVar("y_obs", setName + "   #it{y}", 0, -1e4, 1e4)));
}

void PDF_XY::setObservables(const TString c) {
  if (c.EqualTo("truth"))
    setObservablesTruth();
  else if (c.EqualTo("toy"))
    setObservablesToy();
  else if (c.EqualTo("BaBar_Kshh")) {
    obsValSource = "https://inspirehep.net/literature/853279";
    setObservable("x_obs", 1.6e-3);
    setObservable("y_obs", 5.7e-3);
  } else if (c.EqualTo("BaBar_pipipi0")) {
    obsValSource = "https://inspirehep.net/literature/1441203";
    setObservable("x_obs", 15e-3);
    setObservable("y_obs", 2e-3);
  } else if (c.EqualTo("LHCb_KSpipi")) {
    obsValSource = "https://inspirehep.net/literature/1396327";
    setObservable("x_obs", -8.6e-3);
    setObservable("y_obs", 0.3e-3);
  } else if (c.EqualTo("Belle_Belle2")) {
    obsValSource = "https://arxiv.org/abs/2410.22961";
    setObservable("x_obs", 4.0e-3);
    setObservable("y_obs", 2.9e-3);
  } else {
    std::cout << "PDF_XY::setObservables() : ERROR : config " + c + " not found." << std::endl;
    exit(1);
  }
}

void PDF_XY::setUncertainties(const TString c) {
  if (c.EqualTo("BaBar_Kshh")) {
    obsErrSource = "https://inspirehep.net/literature/853279";
    StatErr[0] = sqrt(pow(2.3e-3, 2) + pow(1.2e-3, 2) + pow(0.8e-3, 2));  // x
    StatErr[1] = sqrt(pow(2.0e-3, 2) + pow(1.3e-3, 2) + pow(0.7e-3, 2));  // y
    std::ranges::fill(SystErr, 0.);
  } else if (c.EqualTo("BaBar_pipipi0")) {
    obsErrSource = "https://inspirehep.net/literature/1441203";
    StatErr[0] = sqrt(pow(12e-3, 2) + pow(6e-3, 2));  // x
    StatErr[1] = sqrt(pow(9e-3, 2) + pow(5e-3, 2));   // y
    std::ranges::fill(SystErr, 0.);
  } else if (c.EqualTo("LHCb_KSpipi")) {
    obsErrSource = "https://inspirehep.net/literature/1396327";
    StatErr[0] = sqrt(pow(5.3e-3, 2) + pow(1.7e-3, 2));  // x
    StatErr[1] = sqrt(pow(4.6e-3, 2) + pow(1.3e-3, 2));  // y
    std::ranges::fill(SystErr, 0.);
  } else if (c.EqualTo("Belle_Belle2")) {
    obsErrSource = "https://arxiv.org/abs/2410.22961";
    StatErr[0] = 1.7e-3;  // x
    StatErr[1] = 1.4e-3;  // y
    SystErr[0] = 0.4e-3;  // x
    SystErr[1] = 0.3e-3;  // y
  } else {
    std::cout << "PDF_XY::setUncertainties() : ERROR : config " + c + " not found." << std::endl;
    exit(1);
  }
}

void PDF_XY::setCorrelations(const TString c) {
  resetCorrelations();
  if (c.EqualTo("BaBar_Kshh")) {
    corSource = "https://inspirehep.net/literature/853279";
    corStatMatrix[1][0] = 0.0586;
  } else if (c.EqualTo("BaBar_pipipi0")) {
    corSource = "https://inspirehep.net/literature/1441203";
    corStatMatrix[1][0] = -0.006;
  } else if (c.EqualTo("LHCb_KSpipi")) {
    corSource = "https://inspirehep.net/literature/1396327";
    corStatMatrix[1][0] = 0.37;
  } else if (c.EqualTo("Belle_Belle2")) {
    // Correlations are negligible
    corSource = "https://arxiv.org/abs/2410.22961";
  } else {
    std::cout << "PDF_XY::setCorrelations() : ERROR : config " + c + " not found." << std::endl;
    exit(1);
  }
}

void PDF_XY::buildPdf() {
  pdf = new RooMultiVarGaussian("pdf_" + name, "pdf_" + name, *(RooArgSet*)observables, *(RooArgSet*)theory, covMatrix);
}
