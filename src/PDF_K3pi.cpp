/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: June 2024
 **/

#include <PDF_K3pi.h>

#include <CharmUtils.h>
#include <ParametersCharmCombo.h>

#include <Utils.h>

#include <RooFormulaVar.h>
#include <RooMultiVarGaussian.h>
#include <RooRealVar.h>

#include <algorithm>
#include <map>
#include <string>
#include <vector>

namespace {
  using parametrisations::mix;
  std::map<std::string, std::map<mix, std::string>> theory_expressions = {
      {"c1",
       {
           {mix::pheno, "- k_K3pi * 0.5 * (      qop * (y*cos(Delta_K3pi - phi) + x*sin(Delta_K3pi - phi)) "
                        "                  + 1 / qop * (y*cos(Delta_K3pi + phi) + x*sin(Delta_K3pi + phi)))"},
           {mix::theo, "-k_K3pi * (y12 * cos(Delta_K3pi) * cos(phiG) + x12 * sin(Delta_K3pi) * cos(phiM))"},
       }},
      {"c2",
       {
           {mix::pheno, "(x * x + y * y) / 4"},
           {mix::theo, "(x12 * x12 + y12 * y12) / 4"},
       }},
  };
}  // namespace

PDF_K3pi::PDF_K3pi(const TString measurement_id, const parametrisations::mix mix_param)
    : PDF_Abs{3}, mix_param{mix_param} {
  name = "K3pi_" + measurement_id;
  initParameters();
  initRelations();
  initObservables("LHCb R1");
  setObservables(measurement_id);
  setUncertainties(measurement_id);
  setCorrelations(measurement_id);
  build();
}

void PDF_K3pi::initParameters() {
  std::vector<std::string> param_names = {"r_K3pi", "k_K3pi", "Delta_K3pi"};
  using parametrisations::mix;
  if (mix_param == mix::pheno)
    param_names.insert(param_names.end(), {"x", "y", "qop", "phi"});
  else
    param_names.insert(param_names.end(), {"x12", "y12", "phiM", "phiG"});
  ParametersCharmCombo p;
  parameters = new RooArgList("parameters");
  for (const auto& par : param_names) parameters->add(*(p.get(par)));
}

void PDF_K3pi::initRelations() {
  theory = new RooArgList("theory");  ///< the order of this list must match that of the COR matrix!
  theory->add(*(Utils::makeTheoryVar("r_K3pi_th", "r_K3pi_th", "r_K3pi", parameters)));
  theory->add(*(Utils::makeTheoryVar("c1_th", "c1_th", theory_expressions["c1"][mix_param], parameters)));
  theory->add(*(Utils::makeTheoryVar("c2_th", "c2_th", theory_expressions["c2"][mix_param], parameters)));
}

void PDF_K3pi::initObservables(const TString label) {
  observables = new RooArgList("observables");  ///< the order of this list must match that of the COR matrix!
  observables->add(*(new RooRealVar("r_K3pi_obs", label + "   #it{r_{K3#pi}}", 0, -1e4, 1e4)));
  observables->add(*(new RooRealVar("c1_obs", label + "   #it{#kappa_{K3#pi}y'}", 0, -1e4, 1e4)));
  observables->add(*(new RooRealVar("c2_obs", label + "   (#it{x}^{2}+#it{y}^{2})/4", 0, -1e4, 1e4)));
}

void PDF_K3pi::setObservables(const TString c) {
  if (c.EqualTo("truth"))
    setObservablesTruth();
  else if (c.EqualTo("toy"))
    setObservablesToy();
  else if (c.EqualTo("LHCb-run1")) {
    obsValSource = "https://arxiv.org/abs/1602.07224v2";
    setObservable("r_K3pi_obs", 5.67e-2);
    setObservable("c1_obs", 3e-4);
    setObservable("c2_obs", 4.8e-5);
  } else {
    std::cout << "PDF_K3pi::setObservables() : ERROR : config " + c + " not found." << std::endl;
    exit(1);
  }
}

void PDF_K3pi::setUncertainties(const TString c) {
  if (c.EqualTo("LHCb-run1")) {
    obsErrSource = "https://arxiv.org/abs/1602.07224v2";
    StatErr[0] = 0.12e-2;
    StatErr[1] = 1.8e-3;
    StatErr[2] = 1.8e-5;
    std::ranges::fill(SystErr, 0);
  } else {
    std::cout << "PDF_K3pi::setUncertainties() : ERROR : config " + c + " not found." << std::endl;
    exit(1);
  }
}

void PDF_K3pi::setCorrelations(const TString c) {
  resetCorrelations();
  if (c.EqualTo("LHCb-run1")) {
    corSource = "https://arxiv.org/abs/1602.07224v2";
    std::vector<double> corrs = {
        // clang-format off
        1., 0.91, 0.80,  // r
            1.,   0.94,  // c1
                  1.,    // c2
        // clang-format on
    };
    corStatMatrix = Utils::buildCorMatrix(nObs, corrs);
  } else {
    std::cout << "PDF_K3pi::setCorrelations() : ERROR : config " + c + " not found." << std::endl;
    exit(1);
  }
}

void PDF_K3pi::buildPdf() {
  pdf = new RooMultiVarGaussian("pdf_" + name, "pdf_" + name, *(RooArgSet*)observables, *(RooArgSet*)theory, covMatrix);
}
