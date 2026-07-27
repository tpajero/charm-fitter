/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: October 2025
 **/

#include <PDF_BES_Kpi_7d.h>

#include <CharmUtils.h>

#include <Utils.h>

#include <RooArgList.h>
#include <RooFormulaVar.h>
#include <RooMultiVarGaussian.h>
#include <RooRealVar.h>

#include <boost/algorithm/string.hpp>
#include <cmath>
#include <format>
#include <stdexcept>
#include <vector>

PDF_BES_Kpi_7d::PDF_BES_Kpi_7d(const parametrisations::mix mix_param) : PDF_Charm{7}, mix_param{mix_param} {
  name = "charm-bes-kpi-2025";
  initialise("3+7fb", "3+7fb", "3+7fb");
}

std::set<std::string> PDF_BES_Kpi_7d::getParameterNames() const {
  std::set<std::string> names = {"R_Kpi", "Delta_Kpi", "F_pipipi0"};
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    names.insert("y");
    break;
  case mix::theo:
    names.insert({"phiG", "x12", "y12", "phiM"});
    break;
  default:
    throw std::runtime_error(std::format("PDF_BES_Kpi_7d::getParameterNames ERROR Parametrisation {} not supported",
                                         utils::to_string(mix_param)));
  }
  return names;
}

void PDF_BES_Kpi_7d::initRelations() {
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
  theory->add(*(makeTheoryVar("rcos_3fb_th", "rcos_3fb_th", "-sqrt(R_Kpi)*cos(Delta_Kpi)", parameters)));
  theory->add(*(makeTheoryVar("rsin_3fb_th", "rsin_3fb_th", " sqrt(R_Kpi)*sin(Delta_Kpi)", parameters)));
  theory->add(*(makeTheoryVar("rcos_7fbCP_th", "rcos_7fbCP_th", "-sqrt(R_Kpi)*cos(Delta_Kpi)", parameters)));
  theory->add(*(makeTheoryVar("rcos_7fb_th", "rcos_7fb_th", "-sqrt(R_Kpi)*cos(Delta_Kpi)", parameters)));
  theory->add(*(makeTheoryVar("rsin_7fb_th", "rsin_7fb_th", " sqrt(R_Kpi)*sin(Delta_Kpi)", parameters)));
}

void PDF_BES_Kpi_7d::initObservables() {
  observables = new RooArgList("observables");
  observables->add(*(new RooRealVar("A_kpi_obs", "#it{A_{K#pi}}", 0., -100., 100.)));
  observables->add(*(new RooRealVar("A_kpi_pipipi0_obs", "#it{A_{K#pi}^{#pi#pi#pi^{0}}}", 0., -150., 150)));
  observables->add(
      *(new RooRealVar("rcos_3fb_obs", "#minus#it{r_{D}^{K#pi}}cos#it{#Delta_{D}^{K#pi}}", 0., -20., 20.)));
  observables->add(*(new RooRealVar("rsin_3fb_obs", "#it{r_{D}^{K#pi}}sin#it{#Delta_{D}^{K#pi}}", 0., -20., 20.)));
  observables->add(
      *(new RooRealVar("rcos_7fbCP_obs", "#minus#it{r_{D}^{K#pi}}cos#it{#Delta_{D}^{K#pi}}", 0., -20., 20.)));
  observables->add(
      *(new RooRealVar("rcos_7fb_obs", "#minus#it{r_{D}^{K#pi}}cos#it{#Delta_{D}^{K#pi}}", 0., -20., 20.)));
  observables->add(*(new RooRealVar("rsin_7fb_obs", "#it{r_{D}^{K#pi}}sin#it{#Delta_{D}^{K#pi}}", 0., -20., 20.)));
}

void PDF_BES_Kpi_7d::setObservables(const TString c) {
  if (c.EqualTo("3+7fb")) {
    obsValSource = "https://arxiv.org/pdf/2208.09402v2.pdf, https://arxiv.org/pdf/2506.07907";
    setObservable("A_kpi_obs", 13.2e-2);
    setObservable("A_kpi_pipipi0_obs", 13.0e-2);
    setObservable("rcos_3fb_obs", -5.62e-2);
    setObservable("rsin_3fb_obs", -1.1e-2);
    setObservable("rcos_7fbCP_obs", -7.0e-2);
    setObservable("rcos_7fb_obs", -4.4e-2);
    setObservable("rsin_7fb_obs", -2.2e-2);
  } else {
    throw std::runtime_error(std::format("PDF_BES_Kpi_7d::setObservables ERROR obs config {} not found", c.Data()));
  }
}

void PDF_BES_Kpi_7d::setUncertainties(const TString c) {
  if (c.EqualTo("3+7fb")) {
    obsErrSource = "https://arxiv.org/pdf/2208.09402v2.pdf, https://arxiv.org/pdf/2506.07907";
    StatErr[0] = 1.1e-2;
    SystErr[0] = 0.7e-2;
    StatErr[1] = 1.2e-2;
    SystErr[1] = 0.8e-2;
    StatErr[2] = 0.81e-2;
    SystErr[2] = std::hypot(0.50e-2, 0.10e-2);
    StatErr[3] = 1.2e-2;
    SystErr[3] = std::hypot(0.7e-2, 0.3e-2);
    StatErr[4] = 0.8e-2;
    SystErr[4] = 0.15e-2;
    StatErr[5] = 1.4e-2;
    SystErr[5] = 0.18e-2;
    StatErr[6] = 1.7e-2;
    SystErr[6] = 0.31e-2;
  } else {
    throw std::runtime_error(std::format("PDF_BES_Kpi_7d::setUncertainties ERROR err config {} not found", c.Data()));
  }
}

void PDF_BES_Kpi_7d::setCorrelations(const TString c) {
  resetCorrelations();
  if (c.EqualTo("3+7fb")) {
    corSource = "https://arxiv.org/pdf/2208.09402v2.pdf, https://arxiv.org/pdf/2506.07907";
    std::vector<double> dataStat = {
        // clang-format off
            1.  , 0.38, 0.  , 0.  , 0.  , 0.  , 0.  , // A_Kpi
                  1.  , 0.  , 0.  , 0.  , 0.  , 0.  , // A_Kpi^pipipi0
                        1.  , 0.02, 0.  , 0.  , 0.  , // rcos
                              1.  , 0.  , 0.  , 0.  , // rsin
                                    1.  , 0.  , 0.  , // rcosCP7
                                          1.  , 0.03, // rcos7
                                                1.    // rsin7
        // clang-format on
    };
    corStatMatrix = Utils::buildCorMatrix(nObs, dataStat);
    std::vector<double> dataSyst = {
        // clang-format off
            1.  , 0.16,  0.  ,  0.  , 0.  , 0.  , 0.  , // A_Kpi
                  1.  ,  0.  ,  0.  , 0.  , 0.  , 0.  , // A_Kpi^pipipi0
                         1.  ,  0.  , 0.  , 0.09, 0.34, // rcos
                                1.  , 0.  , 0.  , 0.  , // rsin
                                      1.  , 0.18, 0.  , // rcosCP7
                                            1.  , 0.19, // rcos7
                                                  1.    // rsin7
        // clang-format on
    };
    corSystMatrix = Utils::buildCorMatrix(nObs, dataSyst);
  } else {
    throw std::runtime_error(std::format("PDF_BES_Kpi_7d::setCorrelations ERROR cor config {} not found", c.Data()));
  }
}

void PDF_BES_Kpi_7d::buildPdf() {
  pdf = new RooMultiVarGaussian("pdf_" + name, "pdf_" + name, *(RooArgSet*)observables, *(RooArgSet*)theory, covMatrix);
}
