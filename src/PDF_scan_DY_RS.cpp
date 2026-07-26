/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include <PDF_scan_DY_RS.h>

#include <CharmUtils.h>
#include <ParametersCharmCombo.h>

#include <Utils.h>

#include <RooFormulaVar.h>
#include <RooMultiVarGaussian.h>
#include <RooRealVar.h>

#include <format>
#include <iostream>
#include <stdexcept>

PDF_scan_DY_RS::PDF_scan_DY_RS(const parametrisations::mix mix_param) : PDF_Abs{1}, mix_param{mix_param} {
  name = "scan_DY_RS";
  initParameters();
  initRelations();
  initObservables();
  setObservables();
  setUncertainties();
  setCorrelations();
  build();
}

void PDF_scan_DY_RS::initParameters() {
  ParametersCharmCombo p;
  parameters = new RooArgList("parameters");
  parameters->add(*(p.get("R_Kpi")));
  parameters->add(*(p.get("Acp_KP")));
  parameters->add(*(p.get("Delta_Kpi")));
  parameters->add(*(p.get("DY_RS")));

  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    parameters->add(*(p.get("x")));
    parameters->add(*(p.get("y")));
    parameters->add(*(p.get("qop")));
    parameters->add(*(p.get("phi")));
    break;
  case mix::theo:
    parameters->add(*(p.get("phiG")));
    parameters->add(*(p.get("x12")));
    parameters->add(*(p.get("y12")));
    parameters->add(*(p.get("phiM")));
    break;
  default:
    throw std::runtime_error(std::format("PDF_scan_DY_RS::initParameters ERROR Parametrisation {} not supported",
                                         utils::to_string(mix_param)));
  }
}

void PDF_scan_DY_RS::initRelations() {
  theory = new RooArgList("theory");
  theory->add(*(Utils::makeTheoryVar("DY_RS_scan_th", "DY_RS_scan_th",
                                     std::format("DY_RS - abs({})", utils::dy_kp_expression(mix_param)), parameters)));
}

void PDF_scan_DY_RS::initObservables() {
  observables = new RooArgList("observables");
  observables->add(*(new RooRealVar("DY_RS_scan_obs", "scan   #Delta#it{Y}_{#it{K^{#minus}#pi^{+}}}", 0, -1e4, 1e4)));
}

void PDF_scan_DY_RS::setObservables() { setObservable("DY_RS_scan_obs", 0.); }

void PDF_scan_DY_RS::setUncertainties() {
  StatErr[0] = 5e-7;
  SystErr[0] = 0.;
}

void PDF_scan_DY_RS::setCorrelations() {
  resetCorrelations();
  corSource = "No correlations for one observable";
}

void PDF_scan_DY_RS::buildPdf() {
  pdf = new RooMultiVarGaussian("pdf_" + name, "pdf_" + name, *(RooArgSet*)observables, *(RooArgSet*)theory, covMatrix);
}
