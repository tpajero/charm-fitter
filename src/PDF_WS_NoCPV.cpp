/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include <PDF_WS_NoCPV.h>

#include <CharmUtils.h>

#include <Utils.h>

#include <RooFormulaVar.h>
#include <RooMultiVarGaussian.h>
#include <RooRealVar.h>

#include <algorithm>
#include <cmath>
#include <format>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

namespace {
  using parametrisations::mix;
  // Map containing the expressions for the observables in the various parametrisations
  std::map<std::string, std::map<mix, std::string>> theory_expressions = {
      {"y'",
       {
           {mix::pheno, "y*cos(Delta_Kpi) + x*sin(Delta_Kpi)"},
           {mix::theo, "  y12 * cos(Delta_Kpi) * TMath::Sign(1.,cos(phiG)) "
                       "+ x12 * sin(Delta_Kpi) * TMath::Sign(1.,cos(phiM))"},
           {mix::d0_to_kpi, "yp"},
       }},
      {"x'2",
       {
           {mix::pheno, "TMath::Sq(x*cos(Delta_Kpi) - y*sin(Delta_Kpi))"},
           {mix::theo, "TMath::Sq(- y12*sin(Delta_Kpi) * TMath::Sign(1.,cos(phiG))"
                       "    + x12*cos(Delta_Kpi) * TMath::Sign(1.,cos(phiM)))"},
           {mix::d0_to_kpi, "xp2"},
       }},
  };
}  // namespace

PDF_WS_NoCPV::PDF_WS_NoCPV(const TString measurement_id, const parametrisations::mix mix_param)
    : PDF_Charm{3}, mix_param{mix_param}, measurement_id{measurement_id} {
  name = measurement_id + "_WS_NoCPV";
  initialise(measurement_id, measurement_id, measurement_id);
}

std::set<std::string> PDF_WS_NoCPV::getParameterNames() const {
  using parametrisations::mix;
  std::set<std::string> names = {"R_Kpi"};
  if (mix_param != mix::d0_to_kpi) names.insert("Delta_Kpi");
  switch (mix_param) {
  case mix::pheno:
    names.insert({"x", "y", "qop", "phi"});
    break;
  case mix::theo:
    names.insert({"phiG", "x12", "y12", "phiM"});
    break;
  case mix::d0_to_kpi:
    names.insert({"yp", "xp2"});
    break;
  default:
    throw std::runtime_error(std::format("PDF_WS_NoCPV::getParameterNames ERROR Parametrisation {} not supported",
                                         utils::to_string(mix_param)));
  }
  return names;
}

void PDF_WS_NoCPV::initRelations() {
  theory = new RooArgList("theory");
  theory->add(*(Utils::makeTheoryVar("RD_th", "RD_th", "R_Kpi", parameters)));
  theory->add(*(Utils::makeTheoryVar("yp_th", "yp_th", theory_expressions["y'"][mix_param], parameters)));
  theory->add(*(Utils::makeTheoryVar("xp2_th", "xp2_th", theory_expressions["x'2"][mix_param], parameters)));
}

void PDF_WS_NoCPV::initObservables() {
  observables = new RooArgList("observables");  ///< the order of this list must match that of the COR matrix!
  observables->add(*(new RooRealVar("RD_obs", measurement_id + "   #it{R_{K#pi}}", 0., -1e4, 1e4)));
  observables->add(*(new RooRealVar("yp_obs", measurement_id + "   #it{y'}", 0., -1e4, 1e4)));
  observables->add(*(new RooRealVar("xp2_obs", measurement_id + "   #it{x'}^{2}", 0., -1e4, 1e4)));
}

void PDF_WS_NoCPV::setObservables(const TString c) {
  if (c.EqualTo("truth"))
    setObservablesTruth();
  else if (c.EqualTo("toy"))
    setObservablesToy();
  else if (c.EqualTo("CDF")) {
    obsValSource = "https://inspirehep.net/literature/1254229";
    setObservable("RD_obs", 0.00351);
    setObservable("yp_obs", 4.3e-3);
    setObservable("xp2_obs", 0.8e-4);
  } else if (c.EqualTo("BaBar")) {
    obsValSource = "https://inspirehep.net/literature/746245";
    setObservable("RD_obs", 0.00303);
    setObservable("yp_obs", 9.7e-3);
    setObservable("xp2_obs", -2.2e-4);
  } else if (c.EqualTo("Belle")) {
    obsValSource = "https://inspirehep.net/literature/1277238";
    setObservable("RD_obs", 0.00353);
    setObservable("yp_obs", 4.6e-3);
    setObservable("xp2_obs", 0.9e-4);
  } else {
    throw std::runtime_error(std::format("PDF_WS_NoCPV::setObservables ERROR config {} not found", c.Data()));
  }
}

void PDF_WS_NoCPV::setUncertainties(const TString c) {
  if (c.EqualTo("CDF")) {
    obsErrSource = "https://inspirehep.net/literature/1254229";
    StatErr[0] = 0.00035;  // RD
    StatErr[1] = 4.3e-3;   // y'
    StatErr[2] = 1.8e-4;   // x'2
    std::ranges::fill(SystErr, 0.0);
  } else if (c.EqualTo("BaBar")) {
    obsErrSource = "https://inspirehep.net/literature/746245";
    StatErr[0] = std::hypot(0.00016, 0.00010);  // RD
    StatErr[1] = std::hypot(4.4e-3, 3.1e-3);    // y'+
    StatErr[2] = std::hypot(3.0e-4, 2.1e-4);    // x'2+
    std::ranges::fill(SystErr, 0.0);
  } else if (c.EqualTo("Belle")) {
    obsErrSource = "https://inspirehep.net/literature/1277238";
    StatErr[0] = 0.00013;  // RD
    StatErr[1] = 3.4e-3;   // y'
    StatErr[2] = 2.2e-4;   // x'2
    std::ranges::fill(SystErr, 0.0);
  } else {
    throw std::runtime_error(std::format("PDF_WS_NoCPV::setUncertainties ERROR config {} not found", c.Data()));
  }
}

void PDF_WS_NoCPV::setCorrelations(const TString c) {
  resetCorrelations();
  if (c.EqualTo("CDF")) {
    corSource = "https://inspirehep.net/literature/1254229";
    std::cout << "INFO [PDF_WS_NoCPV]: The correlation matrix of https://inspirehep.net/literature/1254229\n"
                 "                     (CDF) is not positive definite. It has been modified to avoid non convergence.\n"
              << std::endl;
    std::vector<double> dataStat = {
        // clang-format off
        1., -0.967, 0.900,  // RD
             1.,   -0.975,  // y'
                    1.      // x'2
        // clang-format on
    };
    corStatMatrix = Utils::buildCorMatrix(nObs, dataStat);
  } else if (c.EqualTo("BaBar")) {
    corSource = "https://inspirehep.net/literature/746245";
    std::vector<double> dataStat = {
        // clang-format off
        1., -0.87, 0.77,  // RD
             1.,  -0.94,  // y'
                   1.     // x'2
        // clang-format on
    };
    corStatMatrix = Utils::buildCorMatrix(nObs, dataStat);
  } else if (c.EqualTo("Belle")) {
    corSource = "https://inspirehep.net/literature/1277238";
    std::vector<double> dataStat = {
        // clang-format off
        1., -0.865, 0.737,  // RD
             1.,   -0.948,  // y'
                    1.      // x'2
        // clang-format on
    };
    corStatMatrix = Utils::buildCorMatrix(nObs, dataStat);
  } else {
    throw std::runtime_error(std::format("PDF_WS_NoCPV::setCorrelations ERROR config {} not found", c.Data()));
  }
}

void PDF_WS_NoCPV::buildPdf() {
  pdf = new RooMultiVarGaussian("pdf_" + name, "pdf_" + name, *(RooArgSet*)observables, *(RooArgSet*)theory, covMatrix);
}
