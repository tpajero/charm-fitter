/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include <PDF_DY.h>

#include <CharmUtils.h>

#include <Utils.h>

#include <RooFormulaVar.h>
#include <RooRealVar.h>

#include <format>
#include <iostream>
#include <stdexcept>

PDF_DY::PDF_DY(const TString measurement_id, const hypotheses::dy_fsc dy_fsc_hypo,
               const parametrisations::acp acp_param, const parametrisations::mix mix_param)
    : PDF_Charm{dy_fsc_hypo == hypotheses::dy_fsc::none ? 1 : 2}, dy_fsc_hypo{dy_fsc_hypo}, acp_param{acp_param},
      mix_param{mix_param}, measurement_id{measurement_id} {
  name = "DY_" + measurement_id;
  initialise(measurement_id, measurement_id, measurement_id);
}

std::set<std::string> PDF_DY::getParameterNames() const {
  return utils::dy_hh_parameters_names(dy_fsc_hypo, acp_param, mix_param, {"KK", "PP"});
}

void PDF_DY::initRelations() {
  theory = new RooArgList("theory");
  if (nObs == 1) {
    theory->add(
        *(Utils::makeTheoryVar("DY_th", utils::dy_hh_expression(dy_fsc_hypo, acp_param, mix_param), parameters)));
  } else if (nObs == 2) {
    for (const auto& hh : {"KK", "PP"})
      theory->add(*(Utils::makeTheoryVar(std::format("DY_{}_th", hh),
                                         utils::dy_hh_expression(dy_fsc_hypo, acp_param, mix_param, hh), parameters)));
  }
}

void PDF_DY::initObservables() {
  observables = new RooArgList("observables");
  if (nObs == 1) {
    observables->add(*(new RooRealVar("DY_obs", measurement_id + "   #Delta#it{Y}", 0, -1e4, 1e4)));
  } else if (nObs == 2) {
    observables->add(
        *(new RooRealVar("DY_KK_obs", measurement_id + "   #Delta#it{Y}_{#it{K}^{+}#it{K}^{#minus}}", 0, -1e4, 1e4)));
    observables->add(*(
        new RooRealVar("DY_PP_obs", measurement_id + "   #Delta#it{Y}_{#it{#pi}^{+}#it{#pi}^{#minus}}", 0, -1e4, 1e4)));
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
    StatErr = {2.6e-4};
    SystErr = {0.0};
  } else if (nObs == 1 && c.EqualTo("WA2020")) {
    StatErr = {2.0e-4};
    SystErr = {0.0};
  } else if (nObs == 1 && c.EqualTo("WA2021")) {
    StatErr = {1.11e-4};
    SystErr = {0.33e-4};
  } else if (nObs == 1 && c.EqualTo("Belle&BaBar")) {
    StatErr = {15.75e-4};
    SystErr = {4.81e-4};
  } else if (nObs == 2 && c.EqualTo("WA2020")) {
    StatErr = {2.35e-4, 4.30e-4};
    SystErr = {0.57e-4, 0.70e-4};
  } else if (nObs == 2 && c.EqualTo("WA2021")) {
    StatErr = {1.28e-4, 2.36e-4};
    SystErr = {0.32e-4, 0.39e-4};
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
