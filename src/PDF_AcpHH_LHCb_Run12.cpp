/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: September 2022
 **/

#include <PDF_AcpHH_LHCb_Run12.h>

#include <CharmUtils.h>
#include <ParametersCharmCombo.h>

#include <Utils.h>

#include <RooArgList.h>
#include <RooFormulaVar.h>
#include <RooMultiVarGaussian.h>
#include <RooRealVar.h>

#include <format>
#include <iostream>
#include <vector>

PDF_AcpHH_LHCb_Run12::PDF_AcpHH_LHCb_Run12(const parametrisations::mix mix_param,
                                           const parametrisations::dy_fsc dy_fsc_param)
    : PDF_Abs{8}, mix_param{mix_param}, dy_fsc_param{dy_fsc_param} {
  name = "Charm_AcpHH_LHCb_Run12_Run1-2";
  initParameters();
  initRelations();
  initObservables();
  setObservables("lhcb-run12");
  setUncertainties("lhcb-run12");
  setCorrelations("lhcb-run12");
  build();
}

void PDF_AcpHH_LHCb_Run12::initParameters() {
  ParametersCharmCombo p;
  parameters = new RooArgList("parameters");
  parameters->add(*(p.get("Acp_KK")));
  parameters->add(*(p.get("Acp_PP")));
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
    std::cout << "PDF_AcpHH_LHCb_Run12::initParameters : ERROR : "
                 "parametrisations::mix not supported."
              << std::endl;
    exit(1);
  }
  using parametrisations::dy_fsc;
  switch (dy_fsc_param) {
  case dy_fsc::none:
  case dy_fsc::partial:
    break;
  case dy_fsc::full:
    parameters->add(*(p.get("cot_delta_KK")));
    parameters->add(*(p.get("cot_delta_PP")));
    break;
  }
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

void PDF_AcpHH_LHCb_Run12::setUncertainties(const TString c) {
  obsErrSource = "https://cds.cern.ch/record/2799916/";

  StatErr[0] = 15.0e-4;
  StatErr[1] = 15.0e-4;
  StatErr[2] = 8.8e-4;
  StatErr[3] = 6.7e-4;
  StatErr[4] = 16.0e-4;
  StatErr[5] = 8.0e-4;
  StatErr[6] = 8.0e-4;
  StatErr[7] = 3.2e-4;

  SystErr[0] = 10.0e-4;
  SystErr[1] = 10.0e-4;
  SystErr[2] = 1.6e-4;
  SystErr[3] = 2.0e-4;
  SystErr[4] = 8.0e-4;
  SystErr[5] = 3.0e-4;
  SystErr[6] = 5.0e-4;
  SystErr[7] = 0.9e-4;
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
