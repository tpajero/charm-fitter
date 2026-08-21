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

namespace {
  using parametrisations::mix;
  const std::map<mix, std::map<std::string, std::string>> theory_expressions = {
      {mix::pheno,
       {
           {"xpp_th", "            x * cos(Delta_Kpipi0      ) - y * sin(Delta_Kpipi0       )"},
           {"ypp_th", "            y * cos(Delta_Kpipi0      ) + x * sin(Delta_Kpipi0       )"},
           {"xpp+_th", "    qop * (x * cos(Delta_Kpipi0 - phi) - y * sin(Delta_Kpipi0 - phi))"},
           {"ypp+_th", "    qop * (y * cos(Delta_Kpipi0 - phi) + x * sin(Delta_Kpipi0 - phi))"},
           {"xpp-_th", "1 / qop * (x * cos(Delta_Kpipi0 + phi) - y * sin(Delta_Kpipi0 + phi))"},
           {"ypp-_th", "1 / qop * (y * cos(Delta_Kpipi0 + phi) + x * sin(Delta_Kpipi0 + phi))"},
       }},
      {mix::theo,
       {
           {"xpp_th", "x12 * cos(Delta_Kpipi0) * cos(phiM) - y12 * sin(Delta_Kpipi0) * cos(phiG)"},
           {"ypp_th", "y12 * cos(Delta_Kpipi0) * cos(phiG) + x12 * sin(Delta_Kpipi0) * cos(phiM)"},
           {"xpp+_th", "x12 * cos(Delta_Kpipi0 + phiM) - y12 * sin(Delta_Kpipi0 + phiG)"},
           {"ypp+_th", "y12 * cos(Delta_Kpipi0 + phiG) + x12 * sin(Delta_Kpipi0 + phiM)"},
           {"xpp-_th", "x12 * cos(Delta_Kpipi0 - phiM) - y12 * sin(Delta_Kpipi0 - phiG)"},
           {"ypp-_th", "y12 * cos(Delta_Kpipi0 - phiG) + x12 * sin(Delta_Kpipi0 - phiM)"},
       }},
  };

  std::string get_formula(const std::string observable, const mix mix_param) {
    try {
      return theory_expressions.at(mix_param).at(observable);
    } catch (const std::out_of_range& e) {
      std::cerr << std::format("Out of range error, parametrisation {} not handled for observable {}: {}",
                               utils::to_string(mix_param), observable, e.what())
                << std::endl;
      throw;
    }
  }
}  // namespace

PDF_Kpipi0::PDF_Kpipi0(const TString measurement_id, const parametrisations::mix mix_param, const bool allow_cpv)
    : PDF_Charm{allow_cpv ? 4 : 2}, mix_param{mix_param}, measurement_id{measurement_id}, allow_cpv{allow_cpv} {
  name = measurement_id + "_Kpipi0";
  initialise(measurement_id, measurement_id, measurement_id);
}

std::set<std::string> PDF_Kpipi0::getParameterNames() const {
  std::set<std::string> names = {"Delta_Kpipi0"};
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    names.insert({"x", "y"});
    if (allow_cpv) names.insert({"qop", "phi"});
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
  if (allow_cpv) {
    theory->add(*(Utils::makeTheoryVar("xpp+_th", get_formula("xpp+_th", mix_param), parameters)));
    theory->add(*(Utils::makeTheoryVar("ypp+_th", get_formula("ypp+_th", mix_param), parameters)));
    theory->add(*(Utils::makeTheoryVar("xpp-_th", get_formula("xpp-_th", mix_param), parameters)));
    theory->add(*(Utils::makeTheoryVar("ypp-_th", get_formula("ypp-_th", mix_param), parameters)));
  } else {
    theory->add(*(Utils::makeTheoryVar("xpp_th", get_formula("xpp_th", mix_param), parameters)));
    theory->add(*(Utils::makeTheoryVar("ypp_th", get_formula("ypp_th", mix_param), parameters)));
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
  if (allow_cpv) {
    observables->add(*(new RooRealVar("xpp+_obs", label + "   #it{x''}^{+}", 0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("ypp+_obs", label + "   #it{y''}^{+}", 0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("xpp-_obs", label + "   #it{x''}^{#minus}", 0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("ypp-_obs", label + "   #it{y''}^{#minus}", 0., -1e4, 1e4)));
  } else {
    observables->add(*(new RooRealVar("xpp_obs", label + "   #it{x''}", 0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("ypp_obs", label + "   #it{y''}", 0., -1e4, 1e4)));
  }
}

void PDF_Kpipi0::setObservables(const TString c) {
  if (c.EqualTo("truth"))
    setObservablesTruth();
  else if (c.EqualTo("toy"))
    setObservablesToy();
  else if (c.EqualTo("BaBar")) {
    obsValSource = "https://inspirehep.net/literature/791715";
    if (allow_cpv) {
      setObservable("xpp+_obs", 2.53e-2);
      setObservable("ypp+_obs", -0.05e-2);
      setObservable("xpp-_obs", 3.55e-2);
      setObservable("ypp-_obs", -0.54e-2);
    } else {
      setObservable("xpp_obs", 2.61e-2);
      setObservable("ypp_obs", -0.06e-2);
    }
  } else {
    throw std::runtime_error(std::format("PDF_Kpipi0::setObservables ERROR config {} not found", c.Data()));
  }
}

void PDF_Kpipi0::setUncertainties(const TString c) {
  if (c.EqualTo("BaBar")) {
    obsErrSource = "https://inspirehep.net/literature/791715";
    if (allow_cpv) {
      StatErr = {
          std::hypot(0.585e-2, 0.39e-2),
          std::hypot(0.65e-2, 0.50e-2),
          std::hypot(0.78e-2, 0.65e-2),
          std::hypot(0.78e-2, 0.41e-2),
      };
    } else {
      StatErr = {std::hypot(0.625e-2, 0.39e-2), std::hypot(0.595e-2, 0.34e-2)};
    }
    std::ranges::fill(SystErr, 0.0);
  } else {
    throw std::runtime_error(std::format("PDF_Kpipi0::setUncertainties ERROR config {} not found", c.Data()));
  }
}

void PDF_Kpipi0::setCorrelations(const TString c) {
  resetCorrelations();
  if (c.EqualTo("BaBar")) {
    obsErrSource = "https://inspirehep.net/literature/791715";
    std::vector<double> dataStat;
    if (allow_cpv)
      dataStat = {
          // clang-format off
          1., -0.69, 0.,  0.,
               1.  , 0.,  0.,
                     1., -0.66,
                          1.
          // clang-format on
      };
    else
      dataStat = {1., -0.75, 1.};
    corStatMatrix = Utils::buildCorMatrix(nObs, dataStat);
  } else {
    throw std::runtime_error(std::format("PDF_Kpipi0::setCorrelations ERROR config {} not found", c.Data()));
  }
}
