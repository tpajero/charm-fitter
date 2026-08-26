/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: June 2024
 **/

#include <PDF_K3pi.h>

#include <CharmUtils.h>

#include <Utils.h>

#include <RooFormulaVar.h>
#include <RooRealVar.h>

#include <algorithm>
#include <format>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
  using parametrisations::mix;
  std::map<std::string, std::map<mix, std::string>> theory_expressions = {
      {"c1",
       {
           {mix::pheno, "- k_K3pi * 0.5 * (      qop * (y*cos(Delta_K3pi - phi) + x*sin(Delta_K3pi - phi)) "
                        "                  + 1 / qop * (y*cos(Delta_K3pi + phi) + x*sin(Delta_K3pi + phi)))"},
           {mix::theo, "-k_K3pi * (y12 * cos(Delta_K3pi) * cos(phiG) + x12 * sin(Delta_K3pi) * cos(phiM))"},
       }},
      {"c2",
       {
           {mix::pheno, "(x * x + y * y) / 4"},
           {mix::theo, "0.25 * sqrt( "
                       "    + TMath::Sq(x12*x12 + y12*y12)"
                       "    - TMath::Sq(2 * x12 * y12 * sin(phiM - phiG)))"},
       }},
  };
}  // namespace

PDF_K3pi::PDF_K3pi(const TString measurement_id, const parametrisations::mix mix_param)
    : PDF_Charm{3}, mix_param{mix_param} {
  name = "K3pi_" + measurement_id;
  initialise(measurement_id, measurement_id, measurement_id);
}

std::set<std::string> PDF_K3pi::getParameterNames() const {
  std::set<std::string> names = {"r_K3pi", "k_K3pi", "Delta_K3pi"};
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    names.insert({"x", "y", "qop", "phi"});
    break;
  case mix::theo:
    names.insert({"x12", "y12", "phiM", "phiG"});
    break;
  default:
    throw std::runtime_error(
        std::format("PDF_K3pi::getParameterNames ERROR Parametrisation {} not supported", utils::to_string(mix_param)));
  }
  return names;
}

void PDF_K3pi::initRelations() {
  theory = new RooArgList("theory");  ///< the order of this list must match that of the COR matrix!
  theory->add(*(Utils::makeTheoryVar("r_K3pi_th", "r_K3pi", parameters)));
  try {
    theory->add(*(Utils::makeTheoryVar("c1_th", theory_expressions.at("c1").at(mix_param), parameters)));
    theory->add(*(Utils::makeTheoryVar("c2_th", theory_expressions.at("c2").at(mix_param), parameters)));
  } catch (const std::out_of_range& e) {
    throw std::runtime_error(
        std::format("PDF_K3pi::initRelations ERROR Parametrisation {} not handled for observable c1 or c2: {}",
                    utils::to_string(mix_param), e.what()));
  }
}

void PDF_K3pi::initObservables() {
  const TString label = "LHCb R1";
  observables = new RooArgList("observables");  ///< the order of this list must match that of the COR matrix!
  observables->add(*(new RooRealVar("r_K3pi_obs", label + "   #it{r_{K3#pi}}", 0, -1e4, 1e4)));
  observables->add(*(new RooRealVar("c1_obs", label + "   #it{#kappa_{K3#pi}y'}", 0, -1e4, 1e4)));
  observables->add(*(new RooRealVar("c2_obs", label + "   (#it{x}^{2}+#it{y}^{2})/4", 0, -1e4, 1e4)));
}

void PDF_K3pi::setObservables(const TString c) {
  if (c.EqualTo("truth"))
    setObservablesTruth();
  else if (c.EqualTo("toy"))
    setObservablesToy();
  else if (c.EqualTo("LHCb-R1")) {
    obsValSource = "https://arxiv.org/abs/1602.07224v2";
    setObservable("r_K3pi_obs", 5.67e-2);
    setObservable("c1_obs", 0.3e-3);
    setObservable("c2_obs", 4.8e-5);
  } else {
    throw std::runtime_error(std::format("PDF_K3pi::setObservables ERROR config {} not found", c.Data()));
  }
}

void PDF_K3pi::setUncertainties(const TString c) {
  if (c.EqualTo("LHCb-R1")) {
    obsErrSource = "https://arxiv.org/abs/1602.07224v2";
    StatErr = {0.12e-2, 1.8e-3, 1.8e-5};
    std::ranges::fill(SystErr, 0.0);
  } else {
    throw std::runtime_error(std::format("PDF_K3pi::setUncertainties ERROR config {} not found", c.Data()));
  }
}

void PDF_K3pi::setCorrelations(const TString c) {
  resetCorrelations();
  if (c.EqualTo("LHCb-R1")) {
    corSource = "https://arxiv.org/abs/1602.07224v2";
    std::vector<double> corrs = {
        // clang-format off
        1., 0.91, 0.80,  // r
            1.,   0.94,  // c1
                  1.,    // c2
        // clang-format on
    };
    corStatMatrix = Utils::buildCorMatrix(nObs, corrs);
  } else {
    throw std::runtime_error(std::format("PDF_K3pi::setCorrelations ERROR config {} not found", c.Data()));
  }
}
