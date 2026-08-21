/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include <PDF_yCP_plus_yCP_RS.h>

#include <CharmUtils.h>

#include <Utils.h>

#include <RooFormulaVar.h>
#include <RooRealVar.h>

#include <cmath>
#include <format>
#include <iostream>
#include <stdexcept>

PDF_yCP_plus_yCP_RS::PDF_yCP_plus_yCP_RS(const TString measurement_id, const parametrisations::mix mix_param)
    : PDF_Charm{1}, mix_param{mix_param}, measurement_id{measurement_id} {
  name = "yCP_plus_yCP_RS_" + measurement_id;
  initialise(measurement_id, measurement_id, measurement_id);
}

std::set<std::string> PDF_yCP_plus_yCP_RS::getParameterNames() const {
  std::set<std::string> names = {"r_Kpi", "Delta_Kpi"};
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    names.insert({"x", "y", "qop", "phi"});
    break;
  case mix::theo:
    names.insert({"phiG", "phiM", "x12", "y12"});
    break;
  default:
    throw std::runtime_error(std::format(
        "PDF_yCP_plus_yCP_RS::getParameterNames ERROR Parametrisation {} not supported", utils::to_string(mix_param)));
  }
  return names;
}

void PDF_yCP_plus_yCP_RS::initRelations() {
  theory = new RooArgList("theory");
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    theory->add(*(Utils::makeTheoryVar("yCP_plus_yCP_RS_th",
                                       "0.5*( "
                                       "      y*(qop + 1/qop)*cos(phi)"
                                       "    - x*(qop - 1/qop)*sin(phi)"
                                       " + r_Kpi * ("
                                       "    - (y * cos(Delta_Kpi) - x * sin(Delta_Kpi)) * (qop + 1/qop) * cos(phi)"
                                       "    + (x * cos(Delta_Kpi) + y * sin(Delta_Kpi)) * (qop - 1/qop) * sin(phi)))",
                                       parameters)));
    break;
  case mix::theo:
    theory->add(*(Utils::makeTheoryVar("yCP_plus_yCP_RS_th",
                                       " y12 * cos(phiG)"
                                       " + r_Kpi * ("
                                       "     - y12 * cos(Delta_Kpi) * cos(phiG)"
                                       "     + x12 * sin(Delta_Kpi) * cos(phiM))",
                                       parameters)));
    break;
  default:
    throw std::runtime_error(std::format("PDF_yCP_plus_yCP_RS::initRelations ERROR Parametrisation {} not supported",
                                         utils::to_string(mix_param)));
  }
}

void PDF_yCP_plus_yCP_RS::initObservables() {
  observables = new RooArgList("observables");
  observables->add(*(new RooRealVar("yCP_plus_yCP_RS_obs",
                                    measurement_id + "   #it{y_{CP}}+#it{y_{CP}^{K^{#minus}#pi^{+}}}", 0, -1e4, 1e4)));
}

void PDF_yCP_plus_yCP_RS::setObservables(const TString c) {
  if (c.EqualTo("truth"))
    setObservablesTruth();
  else if (c.EqualTo("toy"))
    setObservablesToy();
  else if (c.EqualTo("Belle")) {
    obsValSource = "https://inspirehep.net/literature/1772245";
    setObservable("yCP_plus_yCP_RS_obs", 9.6e-3);
  } else {
    throw std::runtime_error(std::format("PDF_yCP_plus_yCP_RS::setObservables ERROR config {} not found", c.Data()));
  }
}

void PDF_yCP_plus_yCP_RS::setUncertainties(const TString c) {
  if (c.EqualTo("Belle")) {
    obsErrSource = "https://inspirehep.net/literature/1772245";
    StatErr = {9.1e-3};
    SystErr = {std::hypot(6.2e-3, 1.7e-3)};  // The second uncertainty is not halved for symmetrisation.
  } else {
    throw std::runtime_error(std::format("PDF_yCP_plus_yCP_RS::setUncertainties ERROR config {} not found", c.Data()));
  }
}
