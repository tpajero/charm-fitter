/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include <PDF_BES_Kpi_1d.h>

#include <CharmUtils.h>
#include <ParametersCharmCombo.h>

#include <Utils.h>

#include <RooFormulaVar.h>
#include <RooMultiVarGaussian.h>
#include <RooRealVar.h>

#include <boost/algorithm/string.hpp>

#include <iostream>
#include <string>

PDF_BES_Kpi_1d::PDF_BES_Kpi_1d(const parametrisations::mix mix_param) : PDF_Abs{1}, mix_param{mix_param} {
  name = "BES";
  initParameters();
  initRelations();
  initObservables(name);
  setObservables("BES");
  setUncertainties("BES");
  setCorrelations("BES");
  build();
}

void PDF_BES_Kpi_1d::initParameters() {
  ParametersCharmCombo p;
  parameters = new RooArgList("parameters");
  parameters->add(*(p.get("Delta_Kpi")));
  parameters->add(*(p.get("R_Kpi")));

  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    parameters->add(*(p.get("x")));
    parameters->add(*(p.get("y")));
    break;
  case mix::theo:
    parameters->add(*(p.get("phiG")));
    parameters->add(*(p.get("x12")));
    parameters->add(*(p.get("y12")));
    parameters->add(*(p.get("phiM")));
    break;
  default:
    std::cout << "PDF_BES_Kpi_1d::initParameters : ERROR : "
                 "parametrisations::mix not supported."
              << std::endl;
    exit(1);
  }
}

void PDF_BES_Kpi_1d::initRelations() {
  theory = new RooArgList("theory");
  std::string a_kpi_formula = "(2 * sqrt(R_Kpi) * cos(Delta_Kpi) + y) / (1 + R_Kpi)";
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    break;
  case mix::theo:
    boost::replace_all(a_kpi_formula, "y", utils::y_expression(mix_param));
    break;
  default:
    std::cout << "PDF_BES_Kpi_1d::initRelations : ERROR : "
                 "parametrisations::mix not supported."
              << std::endl;
    exit(1);
  }
  theory->add(*(Utils::makeTheoryVar("A_kpi_th", "A_kpi_th", a_kpi_formula, parameters)));
}

void PDF_BES_Kpi_1d::initObservables(const TString setName) {
  observables = new RooArgList("observables");
  observables->add(*(new RooRealVar("A_kpi_obs", setName + "   #it{A_{K#pi}^{CP}}", 0., -1e4, 1e4)));
}

void PDF_BES_Kpi_1d::setObservables(const TString c) {
  if (c.EqualTo("truth"))
    setObservablesTruth();
  else if (c.EqualTo("toy"))
    setObservablesToy();
  else if (c.EqualTo("BES")) {
    obsValSource = "http://inspirehep.net/record/1291279";
    setObservable("A_kpi_obs", 12.7e-2);
  } else {
    std::cout << "PDF_BES_Kpi_1d::setObservables() : ERROR : config " + c + " not found." << std::endl;
    exit(1);
  }
}

void PDF_BES_Kpi_1d::setUncertainties(const TString c) {
  if (c.EqualTo("BES")) {
    obsErrSource = "http://inspirehep.net/record/1291279";
    StatErr[0] = 1.3e-2;
    SystErr[0] = 0.7e-2;
  } else {
    std::cout << "PDF_BES_Kpi_1d::setUncertainties() : ERROR : config " + c + " not found." << std::endl;
    exit(1);
  }
}

void PDF_BES_Kpi_1d::setCorrelations(const TString c) {
  resetCorrelations();
  corSource = "No correlations for one observable";
}

void PDF_BES_Kpi_1d::buildPdf() {
  pdf = new RooMultiVarGaussian("pdf_" + name, "pdf_" + name, *(RooArgSet*)observables, *(RooArgSet*)theory, covMatrix);
}
