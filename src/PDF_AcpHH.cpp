/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: August 2026
 **/

#include <PDF_AcpHH.h>

#include <CharmUtils.h>

#include <Utils.h>

#include <RooArgList.h>
#include <RooFormulaVar.h>
#include <RooRealVar.h>

#include <format>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

PDF_AcpHH::PDF_AcpHH(const TString obs_id, const TString unc_id, const hypotheses::dy_fsc dy_fsc_hypo,
                     const parametrisations::acp acp_param, const parametrisations::mix mix_param)
    : PDF_Charm{2}, dy_fsc_hypo{dy_fsc_hypo}, acp_param{acp_param}, mix_param{mix_param} {
  name = "Charm_AcpHH";
  initialise(obs_id, unc_id, unc_id);
}

std::set<std::string> PDF_AcpHH::getParameterNames() const {
  std::set<std::string> names = utils::acp_hh_parameters_names(acp_param, {"KK", "PP"});
  names.merge(utils::dy_hh_parameters_names(dy_fsc_hypo, acp_param, mix_param, {"KK", "PP"}));
  return names;
}

void PDF_AcpHH::initRelations() {
  // Average decay times taken from https://cds.cern.ch/record/2799916/
  theory = new RooArgList("theory");
  add_acpkk(theory, "acp_d0_to_kk_th", 7.0e-13);
  add_dacp(theory, "dacp_th", 6.946e-13, 6.407e-13);
}

void PDF_AcpHH::initObservables() {
  observables = new RooArgList("observables");
  observables->add(*(new RooRealVar("acp_d0_to_kk_obs", "#it{A_{CP}(K^{+}K^{#minus})}", 0, -1, 1)));
  observables->add(*(new RooRealVar("dacp_obs", "#it{#DeltaA_{CP}}", 0, -1, 1)));
}

void PDF_AcpHH::setObservables(const TString obs_id) {
  if (obs_id.EqualTo("truth"))
    setObservablesTruth();
  else if (obs_id.EqualTo("toy"))
    setObservablesToy();
  else
    throw std::runtime_error(std::format("PDF_AcpHH::setObservables ERROR ID \"{}\" not found", obs_id.Data()));
}

void PDF_AcpHH::setUncertainties(const TString unc_id) {
  if (unc_id.EqualTo("LHCb-UI") || unc_id.EqualTo("LHCb-UII")) {
    obsErrSource = "charm-fitter";
    // Run 2 values times scale factor
    using constants::lhcb_extrapolations;
    const auto scale = lhcb_extrapolations.at(unc_id.Data()) / lhcb_extrapolations.at("LHCb-R2");
    StatErr[0] = 5.3e-4 * scale;  // AcpKK
    StatErr[1] = 2.9e-4 * scale;  // DeltaACP
    std::ranges::fill(SystErr, 0.0);
  } else {
    throw std::runtime_error(std::format("PDF_AcpHH::setUncertainties ERROR ID \"{}\" not found", unc_id.Data()));
  }
}

void PDF_AcpHH::setCorrelations(const TString unc_id) {
  if (unc_id.EqualTo("LHCb-UI") || unc_id.EqualTo("LHCb-UII")) {
    corSource = "";
    resetCorrelations();
    corStatMatrix[0][1] = 0.07;  // Run 2 value (averaged between CDP and CDS)
  } else {
    throw std::runtime_error(std::format("PDF_AcpHH::setCorrelations ERROR ID \"{}\" not found", unc_id.Data()));
  }
}

void PDF_AcpHH::add_acpkk(RooArgList* theory, const TString name, const double avg_time) {
  theory->add(*(Utils::makeTheoryVar(name,
                                     std::format("{} + {:.5e} * ({})", utils::acp_expression(acp_param, "KK"),
                                                 avg_time / constants::d0_lifetime,
                                                 utils::dy_hh_expression(dy_fsc_hypo, acp_param, mix_param, "KK")),
                                     parameters)));
}

void PDF_AcpHH::add_dacp(RooArgList* theory, const TString name, const double avg_time_kk, const double avg_time_pp) {
  theory->add(
      *(Utils::makeTheoryVar(name,
                             std::format("{} + {:.5e} * ({}) - ({}) - {:.5e} * ({})",
                                         utils::acp_expression(acp_param, "KK"), avg_time_kk / constants::d0_lifetime,
                                         utils::dy_hh_expression(dy_fsc_hypo, acp_param, mix_param, "KK"),
                                         utils::acp_expression(acp_param, "PP"), avg_time_pp / constants::d0_lifetime,
                                         utils::dy_hh_expression(dy_fsc_hypo, acp_param, mix_param, "PP")),
                             parameters)));
}
