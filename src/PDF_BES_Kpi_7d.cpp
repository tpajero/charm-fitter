/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: October 2025
 **/

#include <PDF_BES_Kpi_7d.h>

#include <CharmUtils.h>
#include <ParametersCharmCombo.h>

#include <Utils.h>

#include <RooArgList.h>
#include <RooFormulaVar.h>
#include <RooMultiVarGaussian.h>
#include <RooRealVar.h>

#include <boost/algorithm/string.hpp>
#include <vector>

PDF_BES_Kpi_7d::PDF_BES_Kpi_7d(const theory_config& th_cfg) : PDF_Abs{7}, th_cfg{th_cfg} {
  name = "charm-bes-kpi-2025";
  initParameters();
  initRelations();
  initObservables();
  setObservables("3+7fb");
  setUncertainties("3+7fb");
  setCorrelations("3+7fb");
  build();
}

void PDF_BES_Kpi_7d::initParameters() {
  ParametersCharmCombo p;
  parameters = new RooArgList("parameters");
  parameters->add(*(p.get("R_Kpi")));
  parameters->add(*(p.get("Delta_Kpi")));
  parameters->add(*(p.get("F_pipipi0")));
  switch (th_cfg) {
  case theory_config::phenomenological:
    parameters->add(*(p.get("y")));
    break;
  case theory_config::theoretical:
    parameters->add(*(p.get("phiG")));
    [[fallthrough]];
  case theory_config::superweak:
    parameters->add(*(p.get("x12")));
    parameters->add(*(p.get("y12")));
    parameters->add(*(p.get("phiM")));
    break;
  default:
    std::cout << "PDF_BES_Kpi_7d::initParameters : ERROR : theory_config not supported." << std::endl;
    std::exit(1);
  }
}

void PDF_BES_Kpi_7d::initRelations() {
  theory = new RooArgList("theory");
  std::string a_kpi_formula = "(2 * 10 * sqrt(R_Kpi) * cos(Delta_Kpi) + y) / (1 + R_Kpi/100)";
  std::string a_kpi_pipipi0_formula =
      "F_pipipi0 * (2 * 10 * sqrt(R_Kpi) * cos(Delta_Kpi) + y)"
      "/ (1 + R_Kpi/100. + (1 - F_pipipi0) * (-2 * sqrt(R_Kpi/100.) * cos(Delta_Kpi) + y/100.))";
  switch (th_cfg) {
  case theory_config::phenomenological:
    break;
  case theory_config::theoretical:
    boost::replace_all(a_kpi_formula, "y", CharmUtils::y_to_theoretical);
    boost::replace_all(a_kpi_pipipi0_formula, "y", CharmUtils::y_to_theoretical);
    break;
  case theory_config::superweak:
    boost::replace_all(a_kpi_formula, "y", CharmUtils::y_to_superweak);
    boost::replace_all(a_kpi_pipipi0_formula, "y", CharmUtils::y_to_superweak);
    break;
  default:
    std::cout << "PDF_BES_Kpi_7d::initRelations : ERROR : theory_config not supported." << std::endl;
    std::exit(1);
  }
  using Utils::makeTheoryVar;
  theory->add(*(makeTheoryVar("A_kpi_th", "A_kpi_th", a_kpi_formula, parameters)));
  theory->add(*(makeTheoryVar("A_kpi_pipipi0_th", "A_kpi_pipipi0_th", a_kpi_pipipi0_formula, parameters)));
  theory->add(*(makeTheoryVar("rcos_3fb_th", "rcos_3fb_th", "-10 * sqrt(R_Kpi)*cos(Delta_Kpi)", parameters)));
  theory->add(*(makeTheoryVar("rsin_3fb_th", "rsin_3fb_th", " 10 * sqrt(R_Kpi)*sin(Delta_Kpi)", parameters)));
  theory->add(*(makeTheoryVar("rcos_7fbCP_th", "rcos_7fbCP_th", "-10 * sqrt(R_Kpi)*cos(Delta_Kpi)", parameters)));
  theory->add(*(makeTheoryVar("rcos_7fb_th", "rcos_7fb_th", "-10 * sqrt(R_Kpi)*cos(Delta_Kpi)", parameters)));
  theory->add(*(makeTheoryVar("rsin_7fb_th", "rsin_7fb_th", " 10 * sqrt(R_Kpi)*sin(Delta_Kpi)", parameters)));
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

void PDF_BES_Kpi_7d::setObservables(TString c) {
  if (c.EqualTo("3+7fb")) {
    obsValSource = "https://arxiv.org/pdf/2208.09402v2.pdf, https://arxiv.org/pdf/2506.07907";
    setObservable("A_kpi_obs", 13.2);
    setObservable("A_kpi_pipipi0_obs", 13.0);
    setObservable("rcos_3fb_obs", -5.62);
    setObservable("rsin_3fb_obs", -1.1);
    setObservable("rcos_7fbCP_obs", -7.0);
    setObservable("rcos_7fb_obs", -4.4);
    setObservable("rsin_7fb_obs", -2.2);
  } else {
    std::cout << "PDF_BES_Kpi_7d::setObservables() : ERROR : obs config " << c << " not found." << std::endl;
    std::exit(1);
  }
}

void PDF_BES_Kpi_7d::setUncertainties(TString c) {
  if (c.EqualTo("3+7fb")) {
    obsErrSource = "https://arxiv.org/pdf/2208.09402v2.pdf, https://arxiv.org/pdf/2506.07907";
    StatErr[0] = 1.1;
    SystErr[0] = 0.7;
    StatErr[1] = 1.2;
    SystErr[1] = 0.8;
    StatErr[2] = 0.81;
    SystErr[2] = std::sqrt(std::pow(0.50, 2) + std::pow(0.10, 2));
    StatErr[3] = 1.2;
    StatErr[3] = std::sqrt(std::pow(0.7, 2) + std::pow(0.3, 2));
    StatErr[4] = 0.8;
    SystErr[4] = 0.15;
    StatErr[5] = 1.4;
    SystErr[5] = 0.18;
    StatErr[6] = 1.7;
    SystErr[6] = 0.31;
  } else {
    std::cout << "PDF_BES_Kpi_7d::setObservables() : ERROR : err config " << c << " not found." << std::endl;
    std::exit(1);
  }
}

void PDF_BES_Kpi_7d::setCorrelations(TString c) {
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
    std::cout << "PDF_BES_Kpi_7d::setCorrelations() : ERROR : cor config " << c << " not found." << std::endl;
    std::exit(1);
  }
}

void PDF_BES_Kpi_7d::buildPdf() {
  pdf = new RooMultiVarGaussian("pdf_" + name, "pdf_" + name, *(RooArgSet*)observables, *(RooArgSet*)theory, covMatrix);
}
