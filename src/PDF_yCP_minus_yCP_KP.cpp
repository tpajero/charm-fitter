/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include <PDF_yCP_minus_yCP_KP.h>

#include <CharmUtils.h>

#include <Utils.h>

#include <RooFormulaVar.h>
#include <RooMultiVarGaussian.h>
#include <RooRealVar.h>

#include <TString.h>

#include <format>
#include <iostream>
#include <stdexcept>

PDF_yCP_minus_yCP_KP::PDF_yCP_minus_yCP_KP(const TString measurement_id, const parametrisations::mix mix_param)
    : PDF_Charm{1}, mix_param{mix_param}, measurement_id{measurement_id} {
  name = "yCP_minus_yCP_KP_" + measurement_id;
  initialise(measurement_id, measurement_id, measurement_id);
}

std::set<std::string> PDF_yCP_minus_yCP_KP::getParameterNames() const {
  std::set<std::string> names = {"r_Kpi", "Delta_Kpi"};
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    names.insert({"x", "y", "qop", "phi"});
    break;
  case mix::theo:
    names.insert({"phiG", "x12", "y12"});
    break;
  default:
    throw std::runtime_error(std::format(
        "PDF_yCP_minus_yCP_KP::getParameterNames ERROR Parametrisation {} not supported", utils::to_string(mix_param)));
  }
  return names;
}

void PDF_yCP_minus_yCP_KP::initRelations() {
  theory = new RooArgList("theory");
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    theory->add(*(Utils::makeTheoryVar("yCP_minus_yCP_KP_th", "yCP_minus_yCP_KP_th",
                                       " 0.5*( "
                                       "       y*(qop + 1/qop)*cos(phi)"
                                       "     - x*(qop - 1/qop)*sin(phi))"
                                       " + r_Kpi * cos(Delta_Kpi) * ("
                                       "      y * (qop + 1/qop) * cos(phi)"
                                       "    - x * (qop - 1/qop) * sin(phi))",
                                       parameters)));
    break;
  case mix::theo:
    theory->add(*(Utils::makeTheoryVar("yCP_minus_yCP_KP_th", "yCP_minus_yCP_KP_th",
                                       "y12*cos(phiG)"
                                       "+ 2 * r_Kpi * y12 * cos(Delta_Kpi) * cos(phiG)",
                                       parameters)));
    break;
  default:
    throw std::runtime_error(std::format("PDF_yCP_minus_yCP_KP::initRelations ERROR Parametrisation {} not supported",
                                         utils::to_string(mix_param)));
  }
}

void PDF_yCP_minus_yCP_KP::initObservables() {
  observables = new RooArgList("observables");
  observables->add(*(
      new RooRealVar("yCP_minus_yCP_KP_obs", measurement_id + "   #it{y_{CP}}#minus#it{y_{CP}^{K#pi}}", 0, -1e4, 1e4)));
}

void PDF_yCP_minus_yCP_KP::setObservables(const TString c) {
  if (c.EqualTo("truth"))
    setObservablesTruth();
  else if (c.EqualTo("toy"))
    setObservablesToy();
  else if (c.EqualTo("WA2020")) {
    obsValSource = "https://cds.cern.ch/record/2747731";
    setObservable("yCP_minus_yCP_KP_obs", 0.732e-2);
  } else {
    throw std::runtime_error(std::format("PDF_yCP_minus_yCP_KP::setObservables ERROR config {} not found", c.Data()));
  }
}

void PDF_yCP_minus_yCP_KP::setUncertainties(const TString c) {
  if (c.EqualTo("WA2020")) {
    obsErrSource = "https://cds.cern.ch/record/2747731";
    StatErr = {3.068e-2};
    SystErr = {0.0};
  } else {
    throw std::runtime_error(std::format("PDF_yCP_minus_yCP_KP::setUncertainties ERROR config {} not found", c.Data()));
  }
}

void PDF_yCP_minus_yCP_KP::setCorrelations(const TString c) {
  resetCorrelations();
  corSource = "No correlations for one observable";
}

void PDF_yCP_minus_yCP_KP::buildPdf() {
  pdf = new RooMultiVarGaussian("pdf_" + name, "pdf_" + name, *(RooArgSet*)observables, *(RooArgSet*)theory, covMatrix);
}
