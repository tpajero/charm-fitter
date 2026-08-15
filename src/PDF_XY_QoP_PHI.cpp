/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include <PDF_XY_QoP_PHI.h>

#include <CharmUtils.h>

#include <Utils.h>

#include <RooFormulaVar.h>
#include <RooRealVar.h>

#include <TString.h>

#include <algorithm>
#include <cmath>
#include <format>
#include <iostream>
#include <stdexcept>
#include <vector>

PDF_XY_QoP_PHI::PDF_XY_QoP_PHI(const TString measurement_id, const parametrisations::mix mix_param)
    : PDF_Charm{4}, mix_param{mix_param}, measurement_id{measurement_id} {
  name = "Kshh_" + measurement_id;
  initialise(measurement_id, measurement_id, measurement_id);
}

std::set<std::string> PDF_XY_QoP_PHI::getParameterNames() const {
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    return {"x", "y", "qop", "phi"};
  case mix::theo:
    return {"phiG", "x12", "y12", "phiM"};
  default:
    throw std::runtime_error(std::format("PDF_XY_QoP_PHI::getParameterNames ERROR Parametrisation {} not supported",
                                         utils::to_string(mix_param)));
  }
}

void PDF_XY_QoP_PHI::initRelations() {
  theory = new RooArgList("theory");
  theory->add(*(Utils::makeTheoryVar("x_th", utils::x_expression(mix_param), parameters)));
  theory->add(*(Utils::makeTheoryVar("y_th", utils::y_expression(mix_param), parameters)));
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    theory->add(*(Utils::makeTheoryVar("qop_th", "qop", parameters)));
    theory->add(*(Utils::makeTheoryVar("phi_th", "phi", parameters)));
    break;
  case mix::theo:
    theory->add(*(Utils::makeTheoryVar("qop_th",
                                       "sqrt(  (x12*x12 + y12*y12 + 2 * x12 * y12 * sin(phiM - phiG))"
                                       "     / sqrt(  TMath::Sq(x12*x12 + y12*y12)                       "
                                       "            - TMath::Sq(2 * x12 * y12 * sin(phiM - phiG))))  ",
                                       parameters)));
    theory->add(*(Utils::makeTheoryVar("phi_th",
                                       "-0.5 * TMath::ATan("
                                       "      (x12*x12 * sin(2*phiM) + y12*y12 * sin(2*phiG))"
                                       "    / (x12*x12 * cos(2*phiM) + y12*y12 * cos(2*phiG)))",
                                       parameters)));
    break;
  default:
    throw std::runtime_error(std::format("PDF_XY_QoP_PHI::initRelations ERROR Parametrisation {} not supported",
                                         utils::to_string(mix_param)));
  }
}

void PDF_XY_QoP_PHI::initObservables() {
  TString label = measurement_id;
  if (measurement_id == "Belle") label = "Belle #it{K}_{S}^{0}#it{#pi^{+}#pi^{#minus}}";

  observables = new RooArgList("observables");  ///< the order of this list must match that of the COR matrix!
  observables->add(*(new RooRealVar("x_obs", label + "   #it{x}", 0., -1e4, 1e4)));
  observables->add(*(new RooRealVar("y_obs", label + "   #it{y}", 0., -1e4, 1e4)));
  observables->add(*(new RooRealVar("qop_obs", label + "   |#it{q}/#it{p}|", 0., -1e4, 1e4)));
  observables->add(*(new RooRealVar("phi_obs", label + "   #it{#phi}_{2}", 0., -1e4, 1e4)));
}

void PDF_XY_QoP_PHI::setObservables(const TString c) {
  if (c.EqualTo("truth"))
    setObservablesTruth();
  else if (c.EqualTo("toy"))
    setObservablesToy();
  else if (c.EqualTo("Belle")) {
    obsValSource = "https://inspirehep.net/literature/1289224";
    setObservable("x_obs", 0.56e-2);
    setObservable("y_obs", 0.30e-2);
    setObservable("qop_obs", 0.90);
    setObservable("phi_obs", Utils::DegToRad(-6.));
  } else {
    throw std::runtime_error(std::format("PDF_XY_QoP_PHI::setObservables ERROR config {} not found", c.Data()));
  }
}

void PDF_XY_QoP_PHI::setUncertainties(const TString c) {
  if (c.EqualTo("Belle")) {
    obsErrSource = "https://inspirehep.net/literature/1289224";
    StatErr[0] = std::hypot(0.19e-2, (std::hypot(0.04e-2, 0.06e-2) + std::hypot(0.08e-2, 0.08e-2)) / 2.0);  // x
    StatErr[1] = std::hypot(0.15e-2, (std::hypot(0.04e-2, 0.03e-2) + std::hypot(0.05e-2, 0.07e-2)) / 2.0);  // y
    StatErr[2] = (std::hypot(0.16, 0.05, 0.06) + std::hypot(0.15, 0.04, 0.05)) / 2.0;                       // qop
    StatErr[3] = Utils::DegToRad(std::hypot(11.0, 3.0, 3.5));                                               // phi
    std::ranges::fill(SystErr, 0.0);
  } else {
    throw std::runtime_error(std::format("PDF_XY_QoP_PHI::setUncertainties ERROR config {} not found", c.Data()));
  }
}

void PDF_XY_QoP_PHI::setCorrelations(const TString c) {
  resetCorrelations();
  if (c.EqualTo("Belle")) {
    corSource = "hflav";
    std::vector<double> dataStat = {
        // clang-format off
        1., 0.054, -0.074, -0.031,  // x
            1.,     0.034, -0.019,  // y
                    1.,     0.044,  // qop
                            1.      // phi
        // clang-format on
    };
    corStatMatrix = Utils::buildCorMatrix(nObs, dataStat);
  } else {
    throw std::runtime_error(std::format("PDF_XY_QoP_PHI::setCorrelations ERROR config {} not found", c.Data()));
  }
}
