/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include <PDF_Kpipi0.h>

#include <CharmUtils.h>

#include <Utils.h>

#include <RooFormulaVar.h>
#include <RooRealVar.h>

#include <algorithm>
#include <cmath>
#include <format>
#include <iostream>
#include <stdexcept>
#include <vector>

PDF_Kpipi0::PDF_Kpipi0(const TString measurement_id, const parametrisations::mix mix_param)
    : PDF_Charm{2}, mix_param{mix_param}, measurement_id{measurement_id} {
  name = measurement_id + "_Kpipi0";
  initialise(measurement_id, measurement_id, measurement_id);
}

std::set<std::string> PDF_Kpipi0::getParameterNames() const {
  std::set<std::string> names = {"Delta_Kpipi0"};
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    names.insert({"x", "y"});
    break;
  case mix::theo:
    names.insert({"phiG", "x12", "y12", "phiM"});
    break;
  default:
    throw std::runtime_error(std::format("PDF_Kpipi0::getParameterNames ERROR Parametrisation {} not supported",
                                         utils::to_string(mix_param)));
  }
  return names;
}

void PDF_Kpipi0::initRelations() {
  theory = new RooArgList("theory");
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    theory->add(*(Utils::makeTheoryVar("xpp_th", "x*cos(Delta_Kpipi0) - y*sin(Delta_Kpipi0)", parameters)));
    theory->add(*(Utils::makeTheoryVar("ypp_th", "y*cos(Delta_Kpipi0) + x*sin(Delta_Kpipi0)", parameters)));
    break;
  case mix::theo:
    theory->add(*(Utils::makeTheoryVar("xpp_th",
                                       "  x12*cos(Delta_Kpipi0) * cos(phiM)"
                                       "- y12*sin(Delta_Kpipi0) * cos(phiG)",
                                       parameters)));
    theory->add(*(Utils::makeTheoryVar("ypp_th",
                                       "  y12 * cos(Delta_Kpipi0) * cos(phiG)"
                                       "+ x12 * sin(Delta_Kpipi0) * cos(phiM)",
                                       parameters)));
    break;
  default:
    throw std::runtime_error(
        std::format("PDF_Kpipi0::initRelations ERROR Parametrisation {} not supported", utils::to_string(mix_param)));
  }
}

void PDF_Kpipi0::initObservables() {
  TString label;
  if (measurement_id.EqualTo("BaBar"))
    label = "BaBar #it{K}^{+}#pi^{#minus}#pi^{0}";
  else
    throw std::runtime_error(
        std::format("PDF_Kpipi0::initObservables ERROR Measurement ID {} not supported", measurement_id.Data()));

  observables = new RooArgList("observables");  ///< the order of this list must match that of the COR matrix!
  observables->add(*(new RooRealVar("xpp_obs", label + "   #it{x''}", 0., -1e4, 1e4)));
  observables->add(*(new RooRealVar("ypp_obs", label + "   #it{y''}", 0., -1e4, 1e4)));
}

void PDF_Kpipi0::setObservables(const TString c) {
  if (c.EqualTo("truth"))
    setObservablesTruth();
  else if (c.EqualTo("toy"))
    setObservablesToy();
  else if (c.EqualTo("BaBar")) {
    obsValSource = "https://inspirehep.net/literature/791715";
    setObservable("xpp_obs", 2.61e-2);
    setObservable("ypp_obs", -0.06e-2);
  } else {
    throw std::runtime_error(std::format("PDF_Kpipi0::setObservables ERROR config {} not found", c.Data()));
  }
}

void PDF_Kpipi0::setUncertainties(const TString c) {
  if (c.EqualTo("BaBar")) {
    obsErrSource = "https://inspirehep.net/literature/791715";
    StatErr = {std::hypot(0.625e-2, 0.39e-2), std::hypot(0.595e-2, 0.34e-2)};
    std::ranges::fill(SystErr, 0.0);
  } else {
    throw std::runtime_error(std::format("PDF_Kpipi0::setUncertainties ERROR config {} not found", c.Data()));
  }
}

void PDF_Kpipi0::setCorrelations(const TString c) {
  resetCorrelations();
  if (c.EqualTo("BaBar")) {
    obsErrSource = "https://inspirehep.net/literature/791715";
    std::vector<double> dataStat = {1., -0.75, 1.};
    corStatMatrix = Utils::buildCorMatrix(nObs, dataStat);
  } else {
    throw std::runtime_error(std::format("PDF_Kpipi0::setCorrelations ERROR config {} not found", c.Data()));
  }
}
