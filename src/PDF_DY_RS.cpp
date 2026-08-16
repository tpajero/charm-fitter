/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include <PDF_DY_RS.h>

#include <CharmUtils.h>

#include <Utils.h>

#include <RooFormulaVar.h>
#include <RooMultiVarGaussian.h>
#include <RooRealVar.h>

#include <format>
#include <set>
#include <stdexcept>
#include <string>

PDF_DY_RS::PDF_DY_RS(const TString measurement_id, const parametrisations::mix mix_param)
    : PDF_Charm{1}, mix_param{mix_param}, measurement_id{measurement_id} {
  name = "DY_RS_" + measurement_id;
  initialise(measurement_id, measurement_id, measurement_id);
}

std::set<std::string> PDF_DY_RS::getParameterNames() const {
  std::set<std::string> names = {"r_Kpi", "Acp_KP", "Delta_Kpi"};
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    names.insert({"x", "y", "qop", "phi"});
    break;
  case mix::theo:
    names.insert({"phiG", "x12", "y12", "phiM"});
    break;
  default:
    throw std::runtime_error(std::format("PDF_DY_RS::getParameterNames ERROR Parametrisation {} not supported",
                                         utils::to_string(mix_param)));
  }
  return names;
}

void PDF_DY_RS::initRelations() {
  theory = new RooArgList("theory");
  theory->add(*(Utils::makeTheoryVar("DY_RS_th", utils::dy_kp_expression(mix_param), parameters)));
}

void PDF_DY_RS::initObservables() {
  observables = new RooArgList("observables");
  observables->add(*(new RooRealVar("DY_RS_obs", measurement_id + "   #it{#Delta Y}^{#it{K#pi}}", 0, -1e4, 1e4)));
}

void PDF_DY_RS::setObservables(const TString c) {
  if (c.EqualTo("truth"))
    setObservablesTruth();
  else if (c.EqualTo("toy"))
    setObservablesToy();
  else if (c.EqualTo("LHCb2021")) {
    obsValSource = "https://inspirehep.net/literature/1864385";
    setObservable("DY_RS_obs", -0.36e-4);
  } else {
    throw std::runtime_error(std::format("PDF_DY_RS::setObservables ERROR config {} not found", c.Data()));
  }
}

void PDF_DY_RS::setUncertainties(const TString c) {
  if (c.EqualTo("LHCb2021")) {
    obsErrSource = "https://inspirehep.net/literature/1864385";
    StatErr = {0.50e-4};
    SystErr = {0.23e-4};
  } else {
    throw std::runtime_error(std::format("PDF_DY_RS::setUncertainties ERROR config {} not found", c.Data()));
  }
}

void PDF_DY_RS::setCorrelations(const TString c) {
  resetCorrelations();
  corSource = "No correlations for one observable";
}

void PDF_DY_RS::buildPdf() {
  pdf = new RooMultiVarGaussian("pdf_" + name, "pdf_" + name, *(RooArgSet*)observables, *(RooArgSet*)theory, covMatrix);
}
