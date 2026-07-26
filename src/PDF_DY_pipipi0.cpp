/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: June 2024
 **/

#include <PDF_DY_pipipi0.h>

#include <CharmUtils.h>
#include <ParametersCharmCombo.h>

#include <Utils.h>

#include <RooFormulaVar.h>
#include <RooMultiVarGaussian.h>
#include <RooRealVar.h>

#include <format>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

PDF_DY_pipipi0::PDF_DY_pipipi0(const TString measurement_id, const parametrisations::mix mix_param)
    : PDF_Abs{1}, mix_param{mix_param} {
  name = "DY_pipipi0_" + measurement_id;
  initParameters();
  initRelations();
  initObservables(measurement_id);
  setObservables(measurement_id);
  setUncertainties(measurement_id);
  setCorrelations(measurement_id);
  build();
}

void PDF_DY_pipipi0::initParameters() {
  std::vector<std::string> param_names = {"F_pipipi0"};
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    param_names.insert(param_names.end(), {"x", "y", "qop", "phi"});
    break;
  case mix::theo:
    param_names.insert(param_names.end(), {"x12", "y12", "phiM"});
    break;
  default:
    throw std::runtime_error(std::format("PDF_DY_pipipi0::initParameters ERROR Parametrisation {} not supported",
                                         utils::to_string(mix_param)));
  }
  ParametersCharmCombo p;
  parameters = new RooArgList("parameters");
  for (const auto& par : param_names) parameters->add(*(p.get(par)));
}

void PDF_DY_pipipi0::initRelations() {
  theory = new RooArgList("theory");
  theory->add(
      *(Utils::makeTheoryVar("DY_pipipi0_th", "DY_pipipi0_th",
                             std::format("-(2 * F_pipipi0 - 1) * ({})", utils::dy_expression(mix_param)), parameters)));
}

void PDF_DY_pipipi0::initObservables(const TString setName) {
  observables = new RooArgList("observables");
  observables->add(
      *(new RooRealVar("DY_pipipi0_obs", setName + "   #Delta#it{Y}(#pi^{+}#pi^{#minus}#pi^{0})", 0, -1e4, 1e4)));
}

void PDF_DY_pipipi0::setObservables(const TString c) {
  obsValSource = "https://arxiv.org/abs/2405.06556";
  if (c.EqualTo("truth"))
    setObservablesTruth();
  else if (c.EqualTo("toy"))
    setObservablesToy();
  else if (c.EqualTo("LHCb-R2"))
    setObservable("DY_pipipi0_obs", -1.21e-4);
  else {
    throw std::runtime_error(std::format(
        "PDF_DY_pipipi0::setObservables ERROR config {} not found for {} DY_pipipi0 observables", c.Data(), nObs));
  }
}

void PDF_DY_pipipi0::setUncertainties(const TString c) {
  obsErrSource = "https://arxiv.org/abs/2405.06556";
  if (c.EqualTo("LHCb-R2")) {
    StatErr[0] = 5.97e-4;
    SystErr[0] = 2.01e-4;  // Removed the sys. unc. for the time binning
  } else {
    throw std::runtime_error(std::format(
        "PDF_DY_pipipi0::setUncertainties ERROR config {} not found for {} DY_pipipi0 observables", c.Data(), nObs));
  }
}

void PDF_DY_pipipi0::setCorrelations(const TString c) {
  resetCorrelations();
  corSource = "No correlations for one observable";
}

void PDF_DY_pipipi0::buildPdf() {
  pdf = new RooMultiVarGaussian("pdf_" + name, "pdf_" + name, *(RooArgSet*)observables, *(RooArgSet*)theory, covMatrix);
}
