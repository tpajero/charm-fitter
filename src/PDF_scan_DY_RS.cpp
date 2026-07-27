/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include <PDF_scan_DY_RS.h>

#include <CharmUtils.h>

#include <Utils.h>

#include <RooFormulaVar.h>
#include <RooMultiVarGaussian.h>
#include <RooRealVar.h>

#include <algorithm>
#include <format>
#include <iostream>
#include <stdexcept>

PDF_scan_DY_RS::PDF_scan_DY_RS(const parametrisations::mix mix_param) : PDF_Charm{1}, mix_param{mix_param} {
  name = "scan_DY_RS";
  initialise("", "", "");
}

std::set<std::string> PDF_scan_DY_RS::getParameterNames() const {
  std::set<std::string> names = {"r_Kpi", "Acp_KP", "Delta_Kpi", "DY_RS"};
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    names.insert({"x", "y", "qop", "phi"});
    break;
  case mix::theo:
    names.insert({"phiG", "x12", "y12", "phiM"});
    break;
  default:
    throw std::runtime_error(std::format("PDF_scan_DY_RS::getParameterNames ERROR Parametrisation {} not supported",
                                         utils::to_string(mix_param)));
  }
  return names;
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

void PDF_scan_DY_RS::setObservables(const TString) { setObservable("DY_RS_scan_obs", 0.); }

void PDF_scan_DY_RS::setUncertainties(const TString) {
  StatErr[0] = 5e-7;
  std::ranges::fill(SystErr, 0.0);
}

void PDF_scan_DY_RS::setCorrelations(const TString) {
  resetCorrelations();
  corSource = "No correlations for one observable";
}

void PDF_scan_DY_RS::buildPdf() {
  pdf = new RooMultiVarGaussian("pdf_" + name, "pdf_" + name, *(RooArgSet*)observables, *(RooArgSet*)theory, covMatrix);
}
