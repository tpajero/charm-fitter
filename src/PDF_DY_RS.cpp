/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include <PDF_DY_RS.h>

#include <CharmUtils.h>
#include <ParametersCharmCombo.h>

#include <Utils.h>

#include <RooFormulaVar.h>
#include <RooMultiVarGaussian.h>
#include <RooRealVar.h>

#include <format>
#include <stdexcept>

PDF_DY_RS::PDF_DY_RS(const TString measurement_id, const parametrisations::mix mix_param)
    : PDF_Abs{1}, mix_param{mix_param} {
  name = "DY_RS_" + measurement_id;
  initParameters();
  initRelations();
  initObservables(measurement_id);
  setObservables(measurement_id);
  setUncertainties(measurement_id);
  setCorrelations(measurement_id);
  build();
}

void PDF_DY_RS::initParameters() {
  ParametersCharmCombo p;
  parameters = new RooArgList("parameters");
  parameters->add(*(p.get("R_Kpi")));
  parameters->add(*(p.get("Acp_KP")));
  parameters->add(*(p.get("Delta_Kpi")));

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
    throw std::runtime_error(
        std::format("PDF_DY_RS::initParameters ERROR Parametrisation {} not supported", utils::to_string(mix_param)));
  }
}

void PDF_DY_RS::initRelations() {
  theory = new RooArgList("theory");  ///< the order of this list must match that of the COR matrix!
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    theory->add(*(Utils::makeTheoryVar("DY_RS_th", "DY_RS_th",
                                       "0.5 * pow(R_Kpi, 0.5) * "
                                       "(  (y*cos(Delta_Kpi) - x*sin(Delta_Kpi))*(qop - 1/qop - Acp_KP)*cos(phi)"
                                       " - (x*cos(Delta_Kpi) + y*sin(Delta_Kpi))*(qop + 1/qop         )*sin(phi))",
                                       parameters)));
    break;
  case mix::theo:
    theory->add(*(Utils::makeTheoryVar("DY_RS_th", "DY_RS_th",
                                       "pow(R_Kpi, 0.5) * "
                                       "(  (-y12*cos(Delta_Kpi)*cos(phiG) + x12*sin(Delta_Kpi)*cos(phiM))*Acp_KP*0.5"
                                       " + ( y12*sin(Delta_Kpi)*sin(phiG) + x12*cos(Delta_Kpi)*sin(phiM))           )",
                                       parameters)));
    break;
  default:
    throw std::runtime_error(
        std::format("PDF_DY_RS::initRelations ERROR Parametrisation {} not supported", utils::to_string(mix_param)));
  }
}

void PDF_DY_RS::initObservables(const TString setName) {
  observables = new RooArgList("observables");
  observables->add(*(new RooRealVar("DY_RS_obs", setName + "   #it{A}_{#Gamma}^{#it{K#pi}}", 0, -1e4, 1e4)));
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
    StatErr[0] = 0.50e-4;
    SystErr[0] = 0.23e-4;
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
