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
  obsValSource = "https://github.com/tpajero/charm-fitter/blob/main/BLUE/main/ycp.cpp";
  if (c.EqualTo("truth"))
    setObservablesTruth();
  else if (c.EqualTo("toy"))
    setObservablesToy();
  else if (c.EqualTo("Belle")) {
    obsValSource = "https://inspirehep.net/literature/821323";
    setObservable("yCP_obs", 1.1e-3);
  } else if (c.EqualTo("WA-2015")) {
    setObservable("yCP_obs", -3.70e-3);
    // Biased world averages including also yCP +/- yCP(RS/KP) measurements --------------------------------------------
  } else if (c.EqualTo("WA-biased-2015-01")) {
    setObservable("yCP_obs", 6.56e-3);
  } else if (c.EqualTo("WA-biased-2015-09")) {
    setObservable("yCP_obs", 8.41e-3);
  } else if (c.EqualTo("WA-biased-2018")) {
    setObservable("yCP_obs", 7.11e-3);
  } else if (c.EqualTo("WA-biased-2019")) {
    setObservable("yCP_obs", 7.14e-3);
  } else if (c.EqualTo("WA-biased-2022")) {
    setObservable("yCP_obs", 6.97e-3);
    // Biased world averages including also yCP +/- yCP(RS/KP) measurements, but without LHCb --------------------------
  } else if (c.EqualTo("WA-biased-no-LHCb-2019")) {
    setObservable("yCP_obs", 8.55e-3);
  } else {
    throw std::runtime_error(std::format("PDF_yCP::setObservables ERROR config {} not found", c.Data()));
  }
}

void PDF_yCP::setUncertainties(const TString c) {
  obsErrSource = "https://github.com/tpajero/charm-fitter/blob/main/BLUE/main/ycp.cpp";
  if (c.EqualTo("Belle")) {
    obsErrSource = "https://inspirehep.net/literature/821323";
    StatErr = {6.1e-3};
    SystErr = {5.2e-3};
  } else if (c.EqualTo("WA-2015")) {
    StatErr = {5.56e-3};
    SystErr = {4.32e-3};
    // Biased world averages including also yCP +/- yCP(RS/KP) measurements --------------------------------------------
  } else if (c.EqualTo("WA-biased-2015-01")) {
    StatErr = {1.64e-3};
    SystErr = {1.10e-3};
  } else if (c.EqualTo("WA-biased-2015-09")) {
    StatErr = {1.32e-3};
    SystErr = {0.75e-3};
  } else if (c.EqualTo("WA-biased-2018")) {
    StatErr = {0.93e-3};
    SystErr = {0.58e-3};
  } else if (c.EqualTo("WA-biased-2019")) {
    StatErr = {0.92e-3};
    SystErr = {0.58e-3};
  } else if (c.EqualTo("WA-biased-2022")) {
    StatErr = {0.25e-3};
    SystErr = {0.13e-3};
    // Biased world averages including also yCP +/- yCP(RS/KP) measurements, but without LHCb --------------------------
  } else if (c.EqualTo("WA-biased-no-LHCb-2019")) {
    StatErr = {1.34e-3};
    SystErr = {0.75e-3};
  } else {
    throw std::runtime_error(std::format("PDF_yCP::setUncertainties ERROR config {} not found", c.Data()));
  }
}
