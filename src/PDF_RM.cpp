/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include <PDF_RM.h>

#include <CharmUtils.h>
#include <ParametersCharmCombo.h>

#include <Utils.h>

#include <RooFormulaVar.h>
#include <RooMultiVarGaussian.h>
#include <RooRealVar.h>

#include <TString.h>

#include <algorithm>
#include <format>
#include <iostream>
#include <stdexcept>

PDF_RM::PDF_RM(const TString measurement_id, const parametrisations::mix mix_param) : PDF_Abs{1}, mix_param{mix_param} {
  name = "RM_" + measurement_id;
  initParameters();
  initRelations();
  initObservables(measurement_id);
  setObservables(measurement_id);
  setUncertainties(measurement_id);
  setCorrelations(measurement_id);
  build();
}

void PDF_RM::initParameters() {
  ParametersCharmCombo p;
  parameters = new RooArgList("parameters");

  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    parameters->add(*(p.get("x")));
    parameters->add(*(p.get("y")));
    break;
  case mix::theo:
    parameters->add(*(p.get("x12")));
    parameters->add(*(p.get("y12")));
    parameters->add(*(p.get("phiM")));
    parameters->add(*(p.get("phiG")));
    break;
  default:
    throw std::runtime_error(
        std::format("PDF_RM::initParameters ERROR Parametrisation {} not supported", utils::to_string(mix_param)));
  }
}

void PDF_RM::initRelations() {
  theory = new RooArgList("theory");
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    theory->add(*(Utils::makeTheoryVar("RM_th", "RM_th", "(pow(x,2) + pow(y,2))/2", parameters)));
    break;
  case mix::theo:
    theory->add(*(Utils::makeTheoryVar("RM_th", "RM_th",
                                       "0.5 * pow( "
                                       "    + pow(pow(x12,2) + pow(y12,2),2)"
                                       "    - pow(2 * x12 * y12 * sin(phiM - phiG),2), 0.5)",
                                       parameters)));
    break;
  default:
    throw std::runtime_error(
        std::format("PDF_RM::initRelations ERROR Parametrisation {} not supported", utils::to_string(mix_param)));
  }
}

void PDF_RM::initObservables(const TString setName) {
  observables = new RooArgList("observables");
  observables->add(*(new RooRealVar("RM_obs", setName + "   #it{R_{M}}", 0, 0., 1e4)));
}

void PDF_RM::setObservables(const TString c) {
  if (c.EqualTo("truth"))
    setObservablesTruth();
  else if (c.EqualTo("toy"))
    setObservablesToy();
  else if (c.EqualTo("HFLAV2016")) {
    obsValSource = "https://hflav-eos.web.cern.ch/hflav-eos/charm/CHARM21/results_mixing.html";
    setObservable("RM_obs", 1.30e-4);
  } else if (c.EqualTo("LHCb_K3pi_Run1")) {
    obsValSource = "https://inspirehep.net/literature/1423070";
    setObservable("RM_obs", 2 * 0.48e-4);
  } else {
    throw std::runtime_error(std::format("PDF_RM::setObservables ERROR config {} not found", c.Data()));
  }
}

void PDF_RM::setUncertainties(const TString c) {
  if (c.EqualTo("HFLAV2016")) {
    obsErrSource = "https://hflav-eos.web.cern.ch/hflav-eos/charm/CHARM21/results_mixing.html";
    StatErr[0] = 2.69e-4;
    std::ranges::fill(SystErr, 0.0);
  } else if (c.EqualTo("LHCb_K3pi_Run1")) {
    obsErrSource = "https://inspirehep.net/literature/1423070";
    StatErr[0] = 2 * 0.18e-4;
    std::ranges::fill(SystErr, 0.0);
  } else {
    throw std::runtime_error(std::format("PDF_RM::setUncertainties ERROR config {} not found", c.Data()));
  }
}

void PDF_RM::setCorrelations(const TString c) {
  resetCorrelations();
  corSource = "No correlations for one observable";
}

void PDF_RM::buildPdf() {
  pdf = new RooMultiVarGaussian("pdf_" + name, "pdf_" + name, *(RooArgSet*)observables, *(RooArgSet*)theory, covMatrix);
}
