/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include <PDF_yCP.h>

#include <CharmUtils.h>

#include <Utils.h>

#include <RooFormulaVar.h>
#include <RooRealVar.h>

#include <TString.h>

#include <format>
#include <iostream>
#include <stdexcept>

PDF_yCP::PDF_yCP(const TString measurement_id, const parametrisations::mix mix_param)
    : PDF_Charm{1}, mix_param{mix_param}, measurement_id{measurement_id} {
  name = "yCP_" + measurement_id;
  initialise(measurement_id, measurement_id, measurement_id);
}

std::set<std::string> PDF_yCP::getParameterNames() const {
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    return {"x", "y", "qop", "phi"};
  case mix::theo:
    return {"phiG", "y12"};
  default:
    throw std::runtime_error(
        std::format("PDF_yCP::getParameterNames ERROR Parametrisation {} not supported", utils::to_string(mix_param)));
  }
}

void PDF_yCP::initRelations() {
  theory = new RooArgList("theory");
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    theory->add(*(Utils::makeTheoryVar("yCP_th",
                                       "0.5*(  y * (qop + 1/qop) * cos(phi)"
                                       "     - x * (qop - 1/qop) * sin(phi))",
                                       parameters)));
    break;
  case mix::theo:
    theory->add(*(Utils::makeTheoryVar("yCP_th", "y12*cos(phiG)", parameters)));
    break;
  default:
    throw std::runtime_error(
        std::format("PDF_yCP::initRelations ERROR Parametrisation {} not supported", utils::to_string(mix_param)));
  }
}

void PDF_yCP::initObservables() {
  observables = new RooArgList("observables");
  observables->add(*(new RooRealVar("yCP_obs", measurement_id + "   #it{y_{CP}}", 0, -1e4, 1e4)));
}

void PDF_yCP::setObservables(const TString c) {
  if (c.EqualTo("truth"))
    setObservablesTruth();
  else if (c.EqualTo("toy"))
    setObservablesToy();
  else if (c.EqualTo("WA2020")) {
    obsValSource = "https://cds.cern.ch/record/2747731";
    setObservable("yCP_obs", -3.70e-3);
  } else if (c.EqualTo("WA2020_biased")) {
    obsValSource = "HFLAV";
    setObservable("yCP_obs", 7.19e-3);
  } else if (c.EqualTo("LHCb2022_biased")) {
    obsValSource = "https://inspirehep.net/literature/2035063";
    setObservable("yCP_obs", 6.96e-3);
  } else {
    throw std::runtime_error(std::format("PDF_yCP::setObservables ERROR config {} not found", c.Data()));
  }
}

void PDF_yCP::setUncertainties(const TString c) {
  if (c.EqualTo("WA2020")) {
    obsErrSource = "https://cds.cern.ch/record/2747731";
    StatErr = {7.04e-3};
    SystErr = {0.0};
  } else if (c.EqualTo("WA2020_biased")) {
    obsErrSource = "HFLAV";
    StatErr = {1.13e-3};
    SystErr = {0.0};
  } else if (c.EqualTo("LHCb2022_biased")) {
    obsErrSource = "https://inspirehep.net/literature/2035063";
    StatErr = {0.26e-3};
    SystErr = {0.13e-3};
  } else {
    throw std::runtime_error(std::format("PDF_yCP::setUncertainties ERROR config {} not found", c.Data()));
  }
}

void PDF_yCP::setCorrelations(const TString c) {
  resetCorrelations();
  corSource = "No correlations for one observable";
}
