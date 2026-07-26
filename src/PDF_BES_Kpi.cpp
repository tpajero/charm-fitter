/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: September 2022
 **/

#include <PDF_BES_Kpi.h>

#include <CharmUtils.h>
#include <ParametersCharmCombo.h>

#include <Utils.h>

#include <RooFormulaVar.h>
#include <RooMultiVarGaussian.h>
#include <RooRealVar.h>

#include <boost/algorithm/string.hpp>

#include <cmath>
#include <format>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

PDF_BES_Kpi::PDF_BES_Kpi(const parametrisations::mix mix_param) : PDF_Abs{4}, mix_param{mix_param} {
  name = "charm-bes-kpi";
  initParameters();
  initRelations();
  initObservables();
  setObservables("3fb");
  setUncertainties("3fb");
  setCorrelations("3fb");
  build();
}

void PDF_BES_Kpi::initParameters() {
  ParametersCharmCombo p;
  parameters = new RooArgList("parameters");
  parameters->add(*(p.get("R_Kpi")));
  parameters->add(*(p.get("Delta_Kpi")));
  parameters->add(*(p.get("F_pipipi0")));
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    parameters->add(*(p.get("y")));
    break;
  case mix::theo:
    parameters->add(*(p.get("phiG")));
    parameters->add(*(p.get("x12")));
    parameters->add(*(p.get("y12")));
    parameters->add(*(p.get("phiM")));
    break;
  default:
    throw std::runtime_error(
        std::format("PDF_BES_Kpi::initParameters ERROR Parametrisation {} not supported", utils::to_string(mix_param)));
  }
}

void PDF_BES_Kpi::initRelations() {
  const std::string y = utils::y_expression(mix_param);
  std::string a_kpi_formula = std::format("(2 * sqrt(R_Kpi) * cos(Delta_Kpi) + {0}) / (1 + R_Kpi)", y);
  std::string a_kpi_pipipi0_formula =
      std::format("F_pipipi0 * (2 * sqrt(R_Kpi) * cos(Delta_Kpi) + {0}) "
                  " / (1 + R_Kpi + (1 - F_pipipi0) * (-2 * sqrt(R_Kpi) * cos(Delta_Kpi) + {0}))",
                  y);
  using Utils::makeTheoryVar;
  theory = new RooArgList("theory");
  theory->add(*(makeTheoryVar("A_kpi_th", "A_kpi_th", a_kpi_formula, parameters)));
  theory->add(*(makeTheoryVar("A_kpi_pipipi0_th", "A_kpi_pipipi0_th", a_kpi_pipipi0_formula, parameters)));
  theory->add(*(makeTheoryVar("rcos_th", "rcos_th", "-sqrt(R_Kpi)*cos(Delta_Kpi)", parameters)));
  theory->add(*(makeTheoryVar("rsin_th", "rsin_th", " sqrt(R_Kpi)*sin(Delta_Kpi)", parameters)));
}

void PDF_BES_Kpi::initObservables() {
  observables = new RooArgList("observables");
  observables->add(*(new RooRealVar("A_kpi_obs", "#it{A_{K#pi}}", 0., -100., 100.)));
  observables->add(*(new RooRealVar("A_kpi_pipipi0_obs", "#it{A_{K#pi}^{#pi#pi#pi^{0}}}", 0., -150., 150)));
  observables->add(*(new RooRealVar("rcos_obs", "#minus#it{r_{D}^{K#pi}}cos#it{#Delta_{D}^{K#pi}}", 0., -20., 20.)));
  observables->add(*(new RooRealVar("rsin_obs", "#it{r_{D}^{K#pi}}sin#it{#Delta_{D}^{K#pi}}", 0., -20., 20.)));
}

void PDF_BES_Kpi::setObservables(const TString c) {
  if (c.EqualTo("3fb")) {
    obsValSource = "https://arxiv.org/pdf/2208.09402v2.pdf";
    setObservable("A_kpi_obs", 13.2e-2);
    setObservable("A_kpi_pipipi0_obs", 13.0e-2);
    setObservable("rcos_obs", -5.62e-2);
    setObservable("rsin_obs", -1.1e-2);
  } else {
    throw std::runtime_error(std::format("PDF_BES_Kpi::setObservables ERROR obs config {} not found", c.Data()));
  }
}

void PDF_BES_Kpi::setUncertainties(const TString c) {
  if (c.EqualTo("3fb")) {
    obsErrSource = "https://arxiv.org/pdf/2208.09402v2.pdf";
    StatErr[0] = 1.1e-2;
    SystErr[0] = 0.7e-2;
    StatErr[1] = 1.2e-2;
    SystErr[1] = 0.8e-2;
    // The stat. errs. for rcos and rsin include both the statistical and systematic components
    StatErr[2] = std::sqrt(std::pow(0.81e-2, 2) + std::pow(0.50e-2, 2) + std::pow(0.10e-2, 2));
    SystErr[2] = 0.;
    StatErr[3] = std::sqrt(std::pow(1.2e-2, 2) + std::pow(0.7e-2, 2) + std::pow(0.3e-2, 2));
    SystErr[3] = 0.;
  } else {
    throw std::runtime_error(std::format("PDF_BES_Kpi::setUncertainties ERROR err config {} not found", c.Data()));
  }
}

void PDF_BES_Kpi::setCorrelations(const TString c) {
  resetCorrelations();
  if (c.EqualTo("3fb")) {
    corSource = "https://arxiv.org/pdf/2208.09402v2.pdf";
    std::vector<double> dataStat = {
        // clang-format off
        1., 0.38, 0., 0.,    // A_Kpi
            1.,   0., 0.,    // A_Kpi^pipipi0
                  1., 0.02,  // rcos
                      1.,    // rsin
        // clang-format on
    };
    corStatMatrix = Utils::buildCorMatrix(nObs, dataStat);
    std::vector<double> dataSyst = {
        // clang-format off
        1., 0.16, 0., 0.,  // A_Kpi
            1.,   0., 0.,  // A_Kpi^pipipi0
                  1., 0.,  // rcos
                      1.,  // rsin
        // clang-format on
    };
    corSystMatrix = Utils::buildCorMatrix(nObs, dataSyst);
  } else {
    throw std::runtime_error(std::format("PDF_BES_Kpi::setCorrelations ERROR cor config {} not found", c.Data()));
  }
}

void PDF_BES_Kpi::buildPdf() {
  pdf = new RooMultiVarGaussian("pdf_" + name, "pdf_" + name, *(RooArgSet*)observables, *(RooArgSet*)theory, covMatrix);
}
