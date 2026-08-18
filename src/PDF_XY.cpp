/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include <PDF_XY.h>

#include <CharmUtils.h>

#include <Utils.h>

#include <RooFormulaVar.h>
#include <RooRealVar.h>

#include <TString.h>

#include <algorithm>
#include <cmath>
#include <format>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

PDF_XY::PDF_XY(const TString measurement_id, const parametrisations::mix mix_param)
    : PDF_Charm{2}, mix_param{mix_param}, measurement_id{measurement_id} {
  name = "XY_" + measurement_id;
  initialise(measurement_id, measurement_id, measurement_id);
}

std::set<std::string> PDF_XY::getParameterNames() const {
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    return {"x", "y"};
  case mix::theo:
    return {"x12", "y12", "phiM", "phiG"};
  default:
    throw std::runtime_error(
        std::format("PDF_XY::getParameterNames ERROR Parametrisation {} not supported", utils::to_string(mix_param)));
  }
}

void PDF_XY::initRelations() {
  theory = new RooArgList("theory");  ///< the order of this list must match that of the COR matrix!

  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    theory->add(*(Utils::makeTheoryVar("x_th", "x", parameters)));
    theory->add(*(Utils::makeTheoryVar("y_th", "y", parameters)));
    break;
  case mix::theo:
    theory->add(*(Utils::makeTheoryVar("x_th", utils::x_expression(mix_param), parameters)));
    theory->add(*(Utils::makeTheoryVar("y_th", utils::y_expression(mix_param), parameters)));
    break;
  default:
    throw std::runtime_error(
        std::format("PDF_XY::initRelations ERROR Parametrisation {} not supported", utils::to_string(mix_param)));
  }
}

void PDF_XY::initObservables() {
  static const std::map<std::string, std::string> labels = {
      {"BaBar-KShh", "BaBar #it{K}_{S}^{0}#it{h}^{+}#it{h}^{#minus}"},
      {"BaBar-pipipi0", "BaBar #it{#pi}^{+}#it{#pi}^{#minus}#it{#pi}^{0}"},
      {"LHCb-KSpipi-2011-prompt", "LHCb #it{K}^{0}_{s}#it{#pi}^{+}#pi^{#minus}"},
      {"Belle-Belle2", "Belle 1+2 #it{K}^{0}_{s}#it{#pi}^{+}#pi^{#minus}"},
  };
  const auto it = labels.find(measurement_id.Data());
  const TString label = it != labels.end() ? it->second : measurement_id;

  observables = new RooArgList("observables");  ///< the order of this list must match that of the COR matrix!
  observables->add(*(new RooRealVar("x_obs", label + "   #it{x}", 0, -1e4, 1e4)));
  observables->add(*(new RooRealVar("y_obs", label + "   #it{y}", 0, -1e4, 1e4)));
}

void PDF_XY::setObservables(const TString c) {
  if (c.EqualTo("truth"))
    setObservablesTruth();
  else if (c.EqualTo("toy"))
    setObservablesToy();
  else if (c.EqualTo("BaBar-KShh")) {
    obsValSource = "https://inspirehep.net/literature/853279";
    setObservable("x_obs", 1.6e-3);
    setObservable("y_obs", 5.7e-3);
  } else if (c.EqualTo("BaBar-pipipi0")) {
    obsValSource = "https://inspirehep.net/literature/1441203";
    setObservable("x_obs", 15.0e-3);
    setObservable("y_obs", 1.9e-3);
  } else if (c.EqualTo("LHCb-KSpipi-2011-prompt")) {
    obsValSource = "https://inspirehep.net/literature/1396327";
    setObservable("x_obs", -8.6e-3);
    setObservable("y_obs", 0.3e-3);
  } else if (c.EqualTo("Belle-Belle2")) {
    obsValSource = "https://inspirehep.net/literature/2843831";
    setObservable("x_obs", 4.0e-3);
    setObservable("y_obs", 2.9e-3);
  } else {
    throw std::runtime_error(std::format("PDF_XY::setObservables ERROR config {} not found", c.Data()));
  }
}

void PDF_XY::setUncertainties(const TString c) {
  if (c.EqualTo("BaBar-KShh")) {
    obsErrSource = "https://inspirehep.net/literature/853279";
    StatErr = {std::hypot(2.3e-3, 1.2e-3, 0.8e-3),   // x
               std::hypot(2.0e-3, 1.3e-3, 0.7e-3)};  // y
    std::ranges::fill(SystErr, 0.0);
  } else if (c.EqualTo("BaBar-pipipi0")) {
    obsErrSource = "https://inspirehep.net/literature/1441203";
    StatErr = {std::hypot(11.7e-3, 5.6e-3),  // x
               std::hypot(8.9e-3, 4.6e-3)};  // y
    std::ranges::fill(SystErr, 0.0);
  } else if (c.EqualTo("LHCb-KSpipi-2011-prompt")) {
    obsErrSource = "https://inspirehep.net/literature/1396327";
    StatErr = {std::hypot(5.3e-3, 1.7e-3),   // x
               std::hypot(4.6e-3, 1.3e-3)};  // y
    std::ranges::fill(SystErr, 0.0);
  } else if (c.EqualTo("Belle-Belle2")) {
    obsErrSource = "https://inspirehep.net/literature/2843831";
    //         x       y
    StatErr = {1.7e-3, 1.4e-3};
    SystErr = {0.4e-3, 0.3e-3};
  } else {
    throw std::runtime_error(std::format("PDF_XY::setUncertainties ERROR config {} not found", c.Data()));
  }
}

void PDF_XY::setCorrelations(const TString c) {
  resetCorrelations();
  if (c.EqualTo("BaBar-KShh")) {
    corSource = "https://inspirehep.net/literature/853279";
    corStatMatrix[1][0] = 0.0586;
    /* --- Python ------------------------------------------------------------------------------------------------------
    xerr = math.hypot(2.3, 1.2, 0.8)
    yerr = math.hypot(2.0, 1.3, 0.7)
    cor  = (2.3 * 2.0 * 3.5e-2 + 1.2 * 1.3 * 16.0e-2 - 0.8 * 0.7 * 2.7e-2) / xerr / yerr
    ----------------------------------------------------------------------------------------------------------------- */
  } else if (c.EqualTo("BaBar-pipipi0")) {
    corSource = "https://inspirehep.net/literature/1441203";
    corStatMatrix[1][0] = -0.006;
  } else if (c.EqualTo("LHCb-KSpipi-2011-prompt")) {
    corSource = "https://inspirehep.net/literature/1396327";
    corStatMatrix[1][0] = 0.37;
  } else if (c.EqualTo("Belle-Belle2")) {
    corSource = "https://inspirehep.net/literature/2843831";
    // Correlations are negligible
  } else {
    throw std::runtime_error(std::format("PDF_XY::setCorrelations ERROR config {} not found", c.Data()));
  }
}
