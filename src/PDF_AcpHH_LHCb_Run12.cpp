/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: September 2022
 **/

#include <PDF_AcpHH_LHCb_Run12.h>

#include <CharmUtils.h>

#include <Utils.h>

#include <RooArgList.h>
#include <RooFormulaVar.h>
#include <RooMultiVarGaussian.h>
#include <RooRealVar.h>

#include <format>
#include <iostream>
#include <stdexcept>
#include <vector>

PDF_AcpHH_LHCb_Run12::PDF_AcpHH_LHCb_Run12(const parametrisations::mix mix_param,
                                           const parametrisations::dy_fsc dy_fsc_param)
    : PDF_Charm{8}, mix_param{mix_param}, dy_fsc_param{dy_fsc_param} {
  name = "Charm_AcpHH_LHCb_Run12_Run1-2";
  initialise("lhcb-run12", "lhcb-run12", "lhcb-run12");
}

std::set<std::string> PDF_AcpHH_LHCb_Run12::getParameterNames() const {
  std::set<std::string> names = {"Acp_KK", "Acp_PP"};
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    names.insert({"x", "y", "qop", "phi"});
    break;
  case mix::theo:
    names.insert({"x12", "y12", "phiM"});
    break;
  default:
    throw std::runtime_error(std::format(
        "PDF_AcpHH_LHCb_Run12::getParameterNames ERROR Parametrisation {} not supported", utils::to_string(mix_param)));
  }
  using parametrisations::dy_fsc;
  switch (dy_fsc_param) {
  case dy_fsc::none:
  case dy_fsc::partial:
    break;
  case dy_fsc::full:
    names.insert({"cot_delta_KK", "cot_delta_PP"});
    break;
  }
  return names;
}

void PDF_AcpHH_LHCb_Run12::initRelations() {
  // Average decay times taken from https://cds.cern.ch/record/2799916/
  theory = new RooArgList("theory");
  add_acpkk(theory, "acp_d0_to_kk_run1_mu_th", 4.310e-13);
  add_acpkk(theory, "acp_d0_to_kk_run1_prompt_th", 9.180e-13);
  add_acpkk(theory, "acp_d0_to_kk_run2_cdp_th", 7.315e-13);
  add_acpkk(theory, "acp_d0_to_kk_run2_cds_th", 6.868e-13);
  add_dacp(theory, "dacp_run1_mu_th", 4.437e-13, 4.379e-13);
  add_dacp(theory, "dacp_run1_prompt_th", 8.827e-13, 8.354e-13);
  add_dacp(theory, "dacp_run2_mu_th", 4.918e-13, 4.931e-13);
  add_dacp(theory, "dacp_run2_prompt_th", 6.946e-13, 6.407e-13);
}

void PDF_AcpHH_LHCb_Run12::initObservables() {
  observables = new RooArgList("observables");
  observables->add(
      *(new RooRealVar("acp_d0_to_kk_run1_mu_obs", "#it{A_{CP}(K^{+}K^{#minus})} Run 1 #it{#mu}", 0, -1, 1)));
  observables->add(
      *(new RooRealVar("acp_d0_to_kk_run1_prompt_obs", "#it{A_{CP}(K^{+}K^{#minus})} Run 1 #it{#pi}", 0, -1, 1)));
  observables->add(
      *(new RooRealVar("acp_d0_to_kk_run2_cdp_obs", "#it{A_{CP}(K^{+}K^{#minus})} Run 2 #it{C_{D^{+}}}", 0, -1, 1)));
  observables->add(*(
      new RooRealVar("acp_d0_to_kk_run2_cds_obs", "#it{A_{CP}(K^{+}K^{#minus})} Run 2 #it{C_{D^{+}_{s}}}", 0, -1, 1)));
  observables->add(*(new RooRealVar("dacp_run1_mu_obs", "#it{#DeltaA_{CP}} Run 1 #it{#mu}", 0, -1, 1)));
  observables->add(*(new RooRealVar("dacp_run1_prompt_obs", "#it{#DeltaA_{CP}} Run 1 #it{#pi}", 0, -1, 1)));
  observables->add(*(new RooRealVar("dacp_run2_mu_obs", "#it{#DeltaA_{CP}} Run 2 #it{#mu}", 0, -1, 1)));
  observables->add(*(new RooRealVar("dacp_run2_prompt_obs", "#it{#DeltaA_{CP}} Run 2 #it{#pi}", 0, -1, 1)));
}

void PDF_AcpHH_LHCb_Run12::setObservables(const TString c) {
  if (c.EqualTo("truth"))
    setObservablesTruth();
  else if (c.EqualTo("toy"))
    setObservablesToy();
  else {
    obsValSource = "https://cds.cern.ch/record/2799916/";
    setObservable("acp_d0_to_kk_run1_mu_obs", -6.0e-4);
    setObservable("acp_d0_to_kk_run1_prompt_obs", 14.0e-4);
    setObservable("acp_d0_to_kk_run2_cdp_obs", 13.6e-4);
    setObservable("acp_d0_to_kk_run2_cds_obs", 2.8e-4);
    setObservable("dacp_run1_mu_obs", 14.0e-4);
    setObservable("dacp_run1_prompt_obs", -10.0e-4);
    setObservable("dacp_run2_mu_obs", -9.0e-4);
    setObservable("dacp_run2_prompt_obs", -18.2e-4);
  }
}

void PDF_AcpHH_LHCb_Run12::setUncertainties(const TString c) {
  obsErrSource = "https://cds.cern.ch/record/2799916/";
  StatErr = {15.0e-4, 15.0e-4, 8.8e-4, 6.7e-4, 16.0e-4, 8.0e-4, 8.0e-4, 3.2e-4};
  SystErr = {10.0e-4, 10.0e-4, 1.6e-4, 2.0e-4, 8.0e-4, 3.0e-4, 5.0e-4, 0.9e-4};
}

void PDF_AcpHH_LHCb_Run12::setCorrelations(const TString c) {
  corSource = "https://cds.cern.ch/record/2799916/";

  std::vector<double> dataStat = {
      // clang-format off
      1., 0.36, 0., 0.,   0.23, 0.,   0., 0.,    // ACP(KK) Run 1 mu
          1.,   0., 0.,   0.,   0.24, 0., 0.,    // ACP(KK) Run 1 prompt
                1., 0.05, 0.,   0.,   0., 0.06,  // ACP(KK) Run 2 CDP
                    1.,   0.,   0.,   0., 0.08,  // ACP(KK) Run 2 CDS
                          1.,   0.,   0., 0.,    // DeltaACP Run 1 mu
                                1.,   0., 0.,    // DeltaACP Run 1 prompt
                                      1., 0.,    // DeltaACP Run 2 mu
                                          1.,    // DeltaACP Run 2 prompt
      // clang-format on
  };
  corStatMatrix = Utils::buildCorMatrix(nObs, dataStat);

  std::vector<double> dataSyst = {
      // clang-format off
      1., 1., 0., 0.,   0.40, 0., 0., 0.,  // ACP(KK) Run 1 mu
          1., 0., 0.,   0.,   0., 0., 0.,  // ACP(KK) Run 1 prompt
              1., 0.28, 0.,   0., 0., 0.,  // ACP(KK) Run 2 CDP
                  1.,   0.,   0., 0., 0.,  // ACP(KK) Run 2 CDS
                        1.,   0., 0., 0.,  // DeltaACP Run 1 mu
                              1., 0., 0.,  // DeltaACP Run 1 prompt
                                  1., 0.,  // DeltaACP Run 2 mu
                                      1.,  // DeltaACP Run 2 prompt
      // clang-format on
  };
  corSystMatrix = Utils::buildCorMatrix(nObs, dataSyst);
}

void PDF_AcpHH_LHCb_Run12::buildPdf() {
  pdf = new RooMultiVarGaussian("pdf_" + name, "pdf_" + name, *(RooArgSet*)observables, *(RooArgSet*)theory, covMatrix);
}

void PDF_AcpHH_LHCb_Run12::add_acpkk(RooArgList* theory, TString name, double avg_time) {
  theory->add(*(Utils::makeTheoryVar(name, name,
                                     std::format("Acp_KK + {:.5e} * ({})", avg_time / constants::d0_lifetime,
                                                 utils::dy_hh_expression(mix_param, dy_fsc_param, "KK")),
                                     parameters)));
}

void PDF_AcpHH_LHCb_Run12::add_dacp(RooArgList* theory, TString name, double avg_time_kk, double avg_time_pp) {
  theory->add(*(Utils::makeTheoryVar(
      name, name,
      std::format("Acp_KK + {:.5e} * ({}) - Acp_PP - {:.5e} * ({})", avg_time_kk / constants::d0_lifetime,
                  utils::dy_hh_expression(mix_param, dy_fsc_param, "KK"), avg_time_pp / constants::d0_lifetime,
                  utils::dy_hh_expression(mix_param, dy_fsc_param, "PP")),
      parameters)));
}
