/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include <PDF_yCP.h>

#include <CharmUtils.h>
#include <ParametersCharmCombo.h>

#include <Utils.h>

#include <RooFormulaVar.h>
#include <RooMultiVarGaussian.h>
#include <RooRealVar.h>

#include <TString.h>

#include <format>
#include <iostream>
#include <stdexcept>

PDF_yCP::PDF_yCP(const TString measurement_id, const parametrisations::mix mix_param)
    : PDF_Abs{1}, mix_param{mix_param} {
  name = "yCP_" + measurement_id;
  initParameters();
  initRelations();
  initObservables(measurement_id);
  setObservables(measurement_id);
  setUncertainties(measurement_id);
  setCorrelations(measurement_id);
  build();
}

void PDF_yCP::initParameters() {
  ParametersCharmCombo p;
  parameters = new RooArgList("parameters");

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
    parameters->add(*(p.get("y12")));
    break;
  default:
    throw std::runtime_error(
        std::format("PDF_yCP::initParameters ERROR Parametrisation {} not supported", utils::to_string(mix_param)));
  }
}

void PDF_yCP::initRelations() {
  theory = new RooArgList("theory");
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    theory->add(*(Utils::makeTheoryVar("yCP_th", "yCP_th",
                                       "0.5*(  y * (qop + 1/qop) * cos(phi)"
                                       "     - x * (qop - 1/qop) * sin(phi))",
                                       parameters)));
    break;
  case mix::theo:
    theory->add(*(Utils::makeTheoryVar("yCP_th", "yCP_th", "y12*cos(phiG)", parameters)));
    break;
  default:
    throw std::runtime_error(
        std::format("PDF_yCP::initRelations ERROR Parametrisation {} not supported", utils::to_string(mix_param)));
  }
}

void PDF_yCP::initObservables(const TString setName) {
  observables = new RooArgList("observables");
  observables->add(*(new RooRealVar("yCP_obs", setName + "   #it{y_{CP}}", 0, -1e4, 1e4)));
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
    StatErr[0] = 7.04e-3;
    SystErr[0] = 0;
  } else if (c.EqualTo("WA2020_biased")) {
    obsValSource = "HFLAV";
    StatErr[0] = 1.13e-3;
    SystErr[0] = 0;
  } else if (c.EqualTo("LHCb2022_biased")) {
    obsValSource = "https://inspirehep.net/literature/2035063";
    StatErr[0] = 0.26e-3;
    SystErr[0] = 0.13e-3;
  } else {
    throw std::runtime_error(std::format("PDF_yCP::setUncertainties ERROR config {} not found", c.Data()));
  }
}

void PDF_yCP::setCorrelations(const TString c) {
  resetCorrelations();
  corSource = "No correlations for one observable";
}

void PDF_yCP::buildPdf() {
  pdf = new RooMultiVarGaussian("pdf_" + name, "pdf_" + name, *(RooArgSet*)observables, *(RooArgSet*)theory, covMatrix);
}
