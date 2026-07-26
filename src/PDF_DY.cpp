/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include <PDF_DY.h>

#include <CharmUtils.h>
#include <ParametersCharmCombo.h>

#include <Utils.h>

#include <RooFormulaVar.h>
#include <RooMultiVarGaussian.h>
#include <RooRealVar.h>

#include <format>
#include <iostream>
#include <stdexcept>

PDF_DY::PDF_DY(const TString measurement_id, const parametrisations::mix mix_param,
               const parametrisations::dy_fsc dy_fsc_param)
    : PDF_Abs{dy_fsc_param == parametrisations::dy_fsc::none ? 1 : 2}, mix_param{mix_param},
      dy_fsc_param{dy_fsc_param} {
  name = "DY_" + measurement_id;
  initParameters();
  initRelations();
  initObservables(measurement_id);
  setObservables(measurement_id);
  setUncertainties(measurement_id);
  setCorrelations(measurement_id);
  build();
}

void PDF_DY::initParameters() {
  ParametersCharmCombo p;
  parameters = new RooArgList("parameters");
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    parameters->add(*(p.get("x")));
    parameters->add(*(p.get("y")));
    parameters->add(*(p.get("qop")));
    parameters->add(*(p.get("phi")));
    break;
  case mix::theo:
    parameters->add(*(p.get("x12")));
    parameters->add(*(p.get("y12")));
    parameters->add(*(p.get("phiM")));
    break;
  default:
    throw std::runtime_error(
        std::format("PDF_DY::initParameters ERROR Parametrisation {} not supported", utils::to_string(mix_param)));
  }
  using parametrisations::dy_fsc;
  switch (dy_fsc_param) {
  case dy_fsc::none:
    break;
  case dy_fsc::partial:
    parameters->add(*(p.get("Acp_KK")));
    parameters->add(*(p.get("Acp_PP")));
    break;
  case dy_fsc::full:
    parameters->add(*(p.get("Acp_KK")));
    parameters->add(*(p.get("Acp_PP")));
    parameters->add(*(p.get("cot_delta_KK")));
    parameters->add(*(p.get("cot_delta_PP")));
    break;
  }
}

void PDF_DY::initRelations() {
  theory = new RooArgList("theory");
  if (nObs == 1) {
    theory->add(*(Utils::makeTheoryVar("DY_th", "DY_th", utils::dy_hh_expression(mix_param), parameters)));
  } else if (nObs == 2) {
    theory->add(*(Utils::makeTheoryVar("DY_KK_th", "DY_KK_th", utils::dy_hh_expression(mix_param, dy_fsc_param, "KK"),
                                       parameters)));
    theory->add(*(Utils::makeTheoryVar("DY_PP_th", "DY_PP_th", utils::dy_hh_expression(mix_param, dy_fsc_param, "PP"),
                                       parameters)));
  }
}

void PDF_DY::initObservables(const TString setName) {
  observables = new RooArgList("observables");
  if (nObs == 1) {
    observables->add(*(new RooRealVar("DY_obs", setName + "   #Delta#it{Y}", 0, -1e4, 1e4)));
  } else if (nObs == 2) {
    observables->add(
        *(new RooRealVar("DY_KK_obs", setName + "   #Delta#it{Y}_{#it{K}^{+}#it{K}^{#minus}}", 0, -1e4, 1e4)));
    observables->add(
        *(new RooRealVar("DY_PP_obs", setName + "   #Delta#it{Y}_{#it{#pi}^{+}#it{#pi}^{#minus}}", 0, -1e4, 1e4)));
  }
}

void PDF_DY::setObservables(const TString c) {
  obsValSource = "https://github.com/tpajero/charm-fitter/tree/master/charmcombo/blue/DY.cpp";
  if (c.EqualTo("truth"))
    setObservablesTruth();
  else if (c.EqualTo("toy"))
    setObservablesToy();
  else if (nObs == 1 && c.EqualTo("WA2019"))
    setObservable("DY_obs", 3.2e-4);
  else if (nObs == 1 && c.EqualTo("WA2020"))
    setObservable("DY_obs", 3.1e-4);
  else if (nObs == 1 && c.EqualTo("WA2021"))
    setObservable("DY_obs", -0.92e-4);
  else if (nObs == 1 && c.EqualTo("Belle&BaBar"))
    setObservable("DY_obs", -1.68e-4);
  else if (nObs == 2 && c.EqualTo("WA2020")) {
    setObservable("DY_KK_obs", 4.99e-4);
    setObservable("DY_PP_obs", -2.40e-4);
  } else if (nObs == 2 && c.EqualTo("WA2021")) {
    setObservable("DY_KK_obs", -0.20e-4);
    setObservable("DY_PP_obs", -3.53e-4);
  } else {
    throw std::runtime_error(
        std::format("PDF_DY::setObservables ERROR config {} not found for {} DY observables", c.Data(), nObs));
  }
}

void PDF_DY::setUncertainties(const TString c) {
  obsErrSource = "https://github.com/tpajero/charm-fitter/tree/master/charmcombo/blue/DY.cpp";
  if (nObs == 1 && c.EqualTo("WA2019")) {
    StatErr[0] = 2.6e-4;
    SystErr[0] = 0.;
  } else if (nObs == 1 && c.EqualTo("WA2020")) {
    StatErr[0] = 2.0e-4;
    SystErr[0] = 0.;
  } else if (nObs == 1 && c.EqualTo("WA2021")) {
    StatErr[0] = 1.11e-4;
    SystErr[0] = 0.33e-4;
  } else if (nObs == 1 && c.EqualTo("Belle&BaBar")) {
    StatErr[0] = 15.75e-4;
    SystErr[0] = 4.81e-4;
  } else if (nObs == 2 && c.EqualTo("WA2020")) {
    StatErr[0] = 2.35e-4;
    SystErr[0] = 0.57e-4;
    StatErr[1] = 4.30e-4;
    SystErr[1] = 0.70e-4;
  } else if (nObs == 2 && c.EqualTo("WA2021")) {
    StatErr[0] = 1.28e-4;
    SystErr[0] = 0.32e-4;
    StatErr[1] = 2.36e-4;
    SystErr[1] = 0.39e-4;
  } else {
    throw std::runtime_error(
        std::format("PDF_DY::setUncertainties ERROR config {} not found for {} DY observables", c.Data(), nObs));
  }
}

void PDF_DY::setCorrelations(const TString c) {
  resetCorrelations();
  corSource = "https://github.com/tpajero/charm-fitter/tree/master/charmcombo/blue/DY.cpp";
  if (nObs == 1)
    corSource = "No correlations for one observable";
  else if (nObs == 2 && c.EqualTo("WA2020"))
    corSystMatrix[0][1] = 0.63;  // np.sum(np.square([0.05, 0.42, 0.10, 0.04, 0.23, 0.09])) / 0.57 / 0.70
  else if (nObs == 2 && c.EqualTo("WA2021"))
    corSystMatrix[0][1] = 0.68;  // np.sum(np.square([0.18, 0.21, 0.06, 0.01, 0.07])) / 0.32 / 0.39
  else {
    throw std::runtime_error(
        std::format("PDF_DY::setCorrelations ERROR config {} not found for {} DY observables", c.Data(), nObs));
  }
}

void PDF_DY::buildPdf() {
  pdf = new RooMultiVarGaussian("pdf_" + name, "pdf_" + name, *(RooArgSet*)observables, *(RooArgSet*)theory, covMatrix);
}
