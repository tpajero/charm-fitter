/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include <PDF_CLEO_Kpi.h>

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

PDF_CLEO_Kpi::PDF_CLEO_Kpi(const TString measurement_id, const parametrisations::mix mix_param)
    : PDF_Charm{5}, mix_param{mix_param} {
  name = "CLEO";
  initialise(measurement_id, measurement_id, measurement_id);
}

std::set<std::string> PDF_CLEO_Kpi::getParameterNames() const {
  std::set<std::string> names = {"r_Kpi", "Delta_Kpi"};
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    names.insert({"x", "y"});
    break;
  case mix::theo:
    names.insert({"phiG", "phiM", "x12", "y12"});
    break;
  default:
    throw std::runtime_error(std::format("PDF_CLEO_Kpi::getParameterNames ERROR Parametrisation {} not supported",
                                         utils::to_string(mix_param)));
  }
  return names;
}

void PDF_CLEO_Kpi::initRelations() {
  theory = new RooArgList("theory");
  theory->add(*(Utils::makeTheoryVar("RD_th", "r_Kpi * r_Kpi", parameters)));
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    theory->add(*(Utils::makeTheoryVar("x2_th", "x*x", parameters)));
    break;
  case mix::theo:
    theory->add(*(Utils::makeTheoryVar("x2_th",
                                       "0.5 * (x12*x12 - y12*y12 + sqrt("
                                       "      TMath::Sq(x12*x12 + y12*y12) "
                                       "    - TMath::Sq(2 * x12 * y12 * sin(phiM - phiG))))",
                                       parameters)));
    break;
  default:
    throw std::runtime_error(
        std::format("PDF_CLEO_Kpi::initRelations ERROR Parametrisation {} not supported", utils::to_string(mix_param)));
  }
  theory->add(*(Utils::makeTheoryVar("y_th", utils::y_expression(mix_param), parameters)));
  theory->add(*(Utils::makeTheoryVar("cos_th", "cos(Delta_Kpi)", parameters)));
  theory->add(*(Utils::makeTheoryVar("sin_th", "-sin(Delta_Kpi)", parameters)));
}

void PDF_CLEO_Kpi::initObservables() {
  observables = new RooArgList("observables");  ///< the order of this list must match that of the COR matrix!
  observables->add(*(new RooRealVar("RD_obs", name + "   #it{R_{K#pi}}", 0., 0., 1e4)));
  observables->add(*(new RooRealVar("x2_obs", name + "   #it{x}^{2}", 0., -1e4, 1e4)));
  observables->add(*(new RooRealVar("y_obs", name + "   #it{y}", 0., -1e4, 1e4)));
  observables->add(*(new RooRealVar("cos_obs", name + "   cos#Delta_{#it{K#pi}}", 0., -1., 1.)));
  observables->add(*(new RooRealVar("sin_obs", name + "   #minussin#Delta_{#it{K#pi}}", 0., -1., 1.)));
}

void PDF_CLEO_Kpi::setObservables(const TString c) {
  if (c.EqualTo("truth"))
    setObservablesTruth();
  else if (c.EqualTo("toy"))
    setObservablesToy();
  else if (c.EqualTo("CLEO-c")) {
    obsValSource = "https://inspirehep.net/literature/1189182";
    setObservable("RD_obs", 5.33e-3);
    setObservable("x2_obs", 0.6e-3);
    setObservable("y_obs", 4.2e-2);
    setObservable("cos_obs", 0.81);
    setObservable("sin_obs", -0.01);
  } else {
    throw std::runtime_error(std::format("PDF_CLEO_Kpi::setObservables ERROR config {} not found", c.Data()));
  }
}

void PDF_CLEO_Kpi::setUncertainties(const TString c) {
  if (c.EqualTo("CLEO-c")) {
    obsErrSource = "https://inspirehep.net/literature/1189182";
    StatErr[0] = std::hypot(1.07e-3, 0.45e-3);  // RD
    StatErr[1] = std::hypot(2.3e-3, 1.1e-3);    // x2
    StatErr[2] = std::hypot(2e-2, 1e-2);        // y
    StatErr[3] = std::hypot(0.20, 0.06);        // cos
    StatErr[4] = std::hypot(0.41, 0.04);        // sin
    std::ranges::fill(SystErr, 0.0);
  } else {
    throw std::runtime_error(std::format("PDF_CLEO_Kpi::setUncertainties ERROR config {} not found", c.Data()));
  }
}

void PDF_CLEO_Kpi::setCorrelations(const TString c) {
  resetCorrelations();
  if (c.EqualTo("CLEO-c")) {
    corSource = "https://inspirehep.net/literature/1189182";
    std::vector<double> dataStat = {
        // clang-format off
        1., 0.,  0.,  -0.42,  0.01,  // RD
            1., -0.73, 0.39,  0.02,  // x2
                 1.,  -0.53, -0.03,  // y
                       1.,    0.04,  // cos
                              1.     // sin
        // clang-format on
    };
    corStatMatrix = Utils::buildCorMatrix(nObs, dataStat);
  } else {
    throw std::runtime_error(std::format("PDF_CLEO_Kpi::setCorrelations ERROR config {} not found", c.Data()));
  }
}
