/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: August 2026
 **/

#include <PDF_BES_Kpi.h>

#include <CharmUtils.h>

#include <Utils.h>

#include <RooFormulaVar.h>
#include <RooRealVar.h>

#include <format>
#include <iostream>
#include <stdexcept>
#include <string>

PDF_BES_Kpi::PDF_BES_Kpi(const parametrisations::mix mix_param) : PDF_Charm{1}, mix_param{mix_param} {
  name = "BES";
  initialise("BES", "BES", "BES");
}

std::set<std::string> PDF_BES_Kpi::getParameterNames() const {
  std::set<std::string> names = {"Delta_Kpi", "r_Kpi"};
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    names.insert({"x", "y"});
    break;
  case mix::theo:
    names.insert({"phiG", "x12", "y12", "phiM"});
    break;
  default:
    throw std::runtime_error(std::format("PDF_BES_Kpi::getParameterNames ERROR Parametrisation {} not supported",
                                         utils::to_string(mix_param)));
  }
  return names;
}

void PDF_BES_Kpi::initRelations() {
  theory = new RooArgList("theory");
  std::string a_kpi_formula =
      std::format("(2 * r_Kpi * cos(Delta_Kpi) + {0}) / (1 + r_Kpi * r_Kpi)", utils::y_expression(mix_param));
  theory->add(*(Utils::makeTheoryVar("A_kpi_th", a_kpi_formula, parameters)));
}

void PDF_BES_Kpi::initObservables() {
  observables = new RooArgList("observables");
  observables->add(*(new RooRealVar("A_kpi_obs", name + "   #it{A_{K#pi}^{CP}}", 0., -1e4, 1e4)));
}

void PDF_BES_Kpi::setObservables(const TString c) {
  if (c.EqualTo("truth"))
    setObservablesTruth();
  else if (c.EqualTo("toy"))
    setObservablesToy();
  else if (c.EqualTo("BES")) {
    obsValSource = "http://inspirehep.net/record/1291279";
    setObservable("A_kpi_obs", 12.7e-2);
  } else {
    throw std::runtime_error(std::format("PDF_BES_Kpi::setObservables ERROR config {} not found", c.Data()));
  }
}

void PDF_BES_Kpi::setUncertainties(const TString c) {
  if (c.EqualTo("BES")) {
    obsErrSource = "http://inspirehep.net/record/1291279";
    StatErr = {1.3e-2};
    SystErr = {0.7e-2};
  } else {
    throw std::runtime_error(std::format("PDF_BES_Kpi::setUncertainties ERROR config {} not found", c.Data()));
  }
}
