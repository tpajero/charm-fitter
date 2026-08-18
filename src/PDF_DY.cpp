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
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

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
    observables->add(*(new RooRealVar("DY_obs", measurement_id + "   #Delta#it{Y}", 0.0, -1e4, 1e4)));
  } else if (nObs == 2) {
    observables->add(
        *(new RooRealVar("DY_KK_obs", measurement_id + "   #Delta#it{Y}_{#it{K}^{+}#it{K}^{#minus}}", 0.0, -1e4, 1e4)));
    observables->add(*(new RooRealVar("DY_PP_obs", measurement_id + "   #Delta#it{Y}_{#it{#pi}^{+}#it{#pi}^{#minus}}",
                                      0.0, -1e4, 1e4)));
  }
}

void PDF_DY::setObservables(const TString c) {
  obsValSource = "https://github.com/tpajero/charm-fitter/blob/main/BLUE/main/dy.cpp";
  if (c.EqualTo("truth")) {
    setObservablesTruth();
    return;
  }
  if (c.EqualTo("toy")) {
    setObservablesToy();
    return;
  }

  static const std::map<std::string, double> obs_hh = {
      {"Belle+BaBar", -1.68e-4}, {"WA-2019", 3.21e-4}, {"WA-2020", 3.23e-4},
      {"WA-2021", -0.96e-4},     {"LHCb-R1", 2.86e-4}, {"LHCb-R12", -1.08e-4},
  };
  static const std::map<std::string, std::pair<double, double>> obs_kk_pp = {
      {"WA-2019", {5.10e-4, -2.36e-4}}, {"WA-2020", {4.99e-4, -2.40e-4}},   {"WA-2021", {-0.20e-4, -3.53e-4}},
      {"LHCb-R1", {4.51e-4, -2.67e-4}}, {"LHCb-R12", {-0.35e-4, -3.61e-4}},
  };
  try {
    if (nObs == 1)
      setObservable("DY_obs", obs_hh.at(c.Data()));
    else {
      const auto [kk, pp] = obs_kk_pp.at(c.Data());
      setObservable("DY_KK_obs", kk);
      setObservable("DY_PP_obs", pp);
    }
  } catch (const std::out_of_range&) {
    throw std::runtime_error(
        std::format("PDF_DY::setObservables ERROR config \"{}\" not found for {} DY observables", c.Data(), nObs));
  }
}

void PDF_DY::setUncertainties(const TString c) {
  obsErrSource = "https://github.com/tpajero/charm-fitter/blob/main/BLUE/main/dy.cpp";
  using errs_t = std::pair<double, double>;
  static const std::map<std::string, errs_t> err_hh = {
      {"Belle+BaBar", {15.75e-4, 4.81e-4}}, {"WA-2019", {2.50e-4, 0.87e-4}}, {"WA-2020", {2.05e-4, 0.57e-4}},
      {"WA-2021", {1.11e-4, 0.33e-4}},      {"LHCb-R1", {2.59e-4, 0.93e-4}}, {"LHCb-R12", {1.12e-4, 0.32e-4}},
  };
  static const std::map<std::string, std::pair<errs_t, errs_t>> err_kk_pp = {
      {"WA-2019", {{2.90e-4, 0.87e-4}, {5.16e-4, 1.01e-4}}},  {"WA-2020", {{2.35e-4, 0.57e-4}, {4.30e-4, 0.70e-4}}},
      {"WA-2021", {{1.28e-4, 0.32e-4}, {2.36e-4, 0.39e-4}}},  {"LHCb-R1", {{2.96e-4, 0.91e-4}, {5.39e-4, 1.11e-4}}},
      {"LHCb-R12", {{1.28e-4, 0.32e-4}, {2.38e-4, 0.40e-4}}},
  };
  try {
    if (nObs == 1) {
      const auto& [stat, syst] = err_hh.at(c.Data());
      StatErr = {stat};
      SystErr = {syst};
    } else {
      const auto& [err_kk, err_pp] = err_kk_pp.at(c.Data());
      StatErr = {err_kk.first, err_pp.first};
      SystErr = {err_kk.second, err_pp.second};
    }
  } catch (const std::out_of_range&) {
    throw std::runtime_error(
        std::format("PDF_DY::setUncertainties ERROR config \"{}\" not found for {} DY observables", c.Data(), nObs));
  }
}

void PDF_DY::setCorrelations(const TString c) {
  resetCorrelations();
  if (nObs == 1) {
    corSource = "No correlations for one observable";
    return;
  }
  corSource = "https://github.com/tpajero/charm-fitter/blob/main/BLUE/main/dy.cpp";
  static const std::map<std::string, double> cor = {
      {"WA-2019", 0.65},   // np.sum(np.square([0.08, 0.66, 0.16, 0.05, 0.32])) / 0.87 / 1.01
      {"WA-2020", 0.60},   // np.sum(np.square([0.05, 0.42, 0.10, 0.04, 0.22])) / 0.57 / 0.70
      {"WA-2021", 0.71},   // np.sum(np.square([0.19, 0.21, 0.06, 0.01, 0.07])) / 0.32 / 0.39
      {"LHCb-R1", 0.62},   // np.sum(np.square([0.09, 0.68, 0.17, 0.06, 0.35])) / 0.91 / 1.11
      {"LHCb-R12", 0.70},  // np.sum(np.square([0.19, 0.21, 0.07, 0.01, 0.07])) / 0.32 / 0.40
  };
  try {
    corSystMatrix[0][1] = cor.at(c.Data());
  } catch (const std::out_of_range&) {
    throw std::runtime_error(
        std::format("PDF_DY::setCorrelations ERROR config \"{}\" not found for {} DY observables", c.Data(), nObs));
  }
}
