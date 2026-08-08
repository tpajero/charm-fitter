/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include <PDF_RM.h>

#include <CharmUtils.h>

#include <Utils.h>

#include <RooFormulaVar.h>
#include <RooMultiVarGaussian.h>
#include <RooRealVar.h>

#include <TString.h>

#include <format>
#include <iostream>
#include <stdexcept>

PDF_RM::PDF_RM(const TString measurement_id, const parametrisations::mix mix_param)
    : PDF_Charm{1}, mix_param{mix_param}, measurement_id{measurement_id} {
  name = "RM_" + measurement_id;
  initialise(measurement_id, measurement_id, measurement_id);
}

std::set<std::string> PDF_RM::getParameterNames() const {
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    return {"x", "y"};
  case mix::theo:
    return {"x12", "y12", "phiM", "phiG"};
  default:
    throw std::runtime_error(
        std::format("PDF_RM::getParameterNames ERROR Parametrisation {} not supported", utils::to_string(mix_param)));
  }
}

void PDF_RM::initRelations() {
  theory = new RooArgList("theory");
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    theory->add(*(Utils::makeTheoryVar("RM_th", "RM_th", "(x*x + y*y)/2", parameters)));
    break;
  case mix::theo:
    theory->add(*(Utils::makeTheoryVar("RM_th", "RM_th",
                                       "0.5 * sqrt( "
                                       "    + TMath::Sq(x12*x12 + y12*y12)"
                                       "    - TMath::Sq(2 * x12 * y12 * sin(phiM - phiG)))",
                                       parameters)));
    break;
  default:
    throw std::runtime_error(
        std::format("PDF_RM::initRelations ERROR Parametrisation {} not supported", utils::to_string(mix_param)));
  }
}

void PDF_RM::initObservables() {
  observables = new RooArgList("observables");
  observables->add(*(new RooRealVar("RM_obs", measurement_id + "   #it{R_{M}}", 0, 0., 1e4)));
}

void PDF_RM::setObservables(const TString c) {
  if (c.EqualTo("truth"))
    setObservablesTruth();
  else if (c.EqualTo("toy"))
    setObservablesToy();
  else if (c.EqualTo("HFLAV2016")) {
    obsValSource = "https://hflav-eos.web.cern.ch/hflav-eos/charm/CHARM21/results_mixing.html";
    setObservable("RM_obs", 1.30e-4);
  } else if (c.EqualTo("LHCb_K3pi_Run1")) {
    obsValSource = "https://inspirehep.net/literature/1423070";
    setObservable("RM_obs", 2 * 0.48e-4);
  } else {
    throw std::runtime_error(std::format("PDF_RM::setObservables ERROR config {} not found", c.Data()));
  }
}

void PDF_RM::setUncertainties(const TString c) {
  if (c.EqualTo("HFLAV2016")) {
    obsErrSource = "https://hflav-eos.web.cern.ch/hflav-eos/charm/CHARM21/results_mixing.html";
    StatErr = {2.69e-4};
    SystErr = {0.0};
  } else if (c.EqualTo("LHCb_K3pi_Run1")) {
    obsErrSource = "https://inspirehep.net/literature/1423070";
    StatErr = {2 * 0.18e-4};
    SystErr = {0.0};
  } else {
    throw std::runtime_error(std::format("PDF_RM::setUncertainties ERROR config {} not found", c.Data()));
  }
}

void PDF_RM::setCorrelations(const TString c) {
  resetCorrelations();
  corSource = "No correlations for one observable";
}

void PDF_RM::buildPdf() {
  pdf = new RooMultiVarGaussian("pdf_" + name, "pdf_" + name, *(RooArgSet*)observables, *(RooArgSet*)theory, covMatrix);
}
