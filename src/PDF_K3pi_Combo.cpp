/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: November 2026
 **/

#include <PDF_K3pi_Combo.h>

#include <CharmUtils.h>

#include <Utils.h>

#include <RooFormulaVar.h>
#include <RooRealVar.h>

#include <TMath.h>

#include <algorithm>
#include <format>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
  const TString lhcb_paper_2026_src =
      "https://arxiv.org/abs/2611.26449; "
      "https://gitlab.cern.ch/lhcb-b2oc/CommonAnalysisTools/k3pi-combiner/-/tags/lhcb-gamma-charm-2026";
}  // namespace

PDF_K3pi_Combo::PDF_K3pi_Combo(const TString measurement_id, const parametrisations::mix mix_param)
    : PDF_Charm{4}, mix_param{mix_param} {
  name = "K3pi_Combo_" + measurement_id;
  initialise(measurement_id, measurement_id, measurement_id);
}

std::set<std::string> PDF_K3pi_Combo::getParameterNames() const {
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    return {"x", "y", "qop", "phi"};
  case mix::theo:
    return {"x12", "y12", "phiM", "phiG"};
  default:
    throw std::runtime_error(std::format("PDF_K3pi_Combo::getParameterNames ERROR Parametrisation {} not supported",
                                         utils::to_string(mix_param)));
  }
}

void PDF_K3pi_Combo::initRelations() {
  const auto x = utils::x_expression(mix_param);
  const auto y = utils::y_expression(mix_param);
  theory = new RooArgList("theory");  ///< the order of this list must match that of the COR matrix!
  theory->add(*(Utils::makeTheoryVar("rMix_th", std::format("sqrt(TMath::Sq({}) + TMath::Sq({}))", x, y), parameters)));
  theory->add(*(Utils::makeTheoryVar("pMix_th", std::format("TMath::ATan2({}, {})", x, y), parameters)));
  theory->add(*(Utils::makeTheoryVar("dx_th", utils::dx_expression(mix_param), parameters)));
  theory->add(*(Utils::makeTheoryVar("dy_th", utils::dy_expression(mix_param), parameters)));
}

void PDF_K3pi_Combo::initObservables() {
  const TString label = "LHCbR2+CF K3#pi 2026";
  observables = new RooArgList("observables");  ///< the order of this list must match that of the COR matrix!
  observables->add(*(new RooRealVar("rMix_obs", label + "   #it{r}_{mix}", 0., -1e4, 1e4)));
  observables->add(*(new RooRealVar("pMix_obs", label + "   #it{p}_{mix}", 0., -1e4, 1e4)));
  observables->add(*(new RooRealVar("dx_obs", label + "   #it{#Deltax}", 0., -1e4, 1e4)));
  observables->add(*(new RooRealVar("dy_obs", label + "   #it{#Deltay}", 0., -1e4, 1e4)));
}

void PDF_K3pi_Combo::setObservables(const TString c) {
  if (c.EqualTo("truth"))
    setObservablesTruth();
  else if (c.EqualTo("toy"))
    setObservablesToy();
  else if (c.EqualTo("LHCb-PAPER-2026-058")) {
    obsValSource = lhcb_paper_2026_src;
    setObservable("rMix_obs", 7.48771e-3);
    setObservable("pMix_obs", 0.833236);
    setObservable("dx_obs", 0.141428e-4);
    setObservable("dy_obs", -0.166178e-4);
  } else {
    throw std::runtime_error(std::format("PDF_K3pi_Combo::setObservables ERROR config {} not found", c.Data()));
  }
}

void PDF_K3pi_Combo::setUncertainties(const TString c) {
  if (c.EqualTo("LHCb-PAPER-2026-058")) {
    obsErrSource = lhcb_paper_2026_src;
    StatErr = {0.721005e-3, 0.272758, 2.86472e-4, 1.94352e-4};  // total uncertainties
    std::ranges::fill(SystErr, 0.0);
  } else {
    throw std::runtime_error(std::format("PDF_K3pi_Combo::setUncertainties ERROR config {} not found", c.Data()));
  }
}

void PDF_K3pi_Combo::setCorrelations(const TString c) {
  resetCorrelations();
  if (c.EqualTo("LHCb-PAPER-2026-058")) {
    corSource = lhcb_paper_2026_src;  // total correlations
    std::vector<double> dataStat = {
        // clang-format off
        //   rMix       pMix         dx         dy
         1.      , -0.515566, -0.088650,  0.107954,  // rMix
                    1.      ,  0.143738, -0.099128,  // pMix
                               1.      , -0.587988,  // dx
                                          1.      ,  // dy
        // clang-format on
    };
    corStatMatrix = Utils::buildCorMatrix(nObs, dataStat);
  } else {
    throw std::runtime_error(std::format("PDF_K3pi_Combo::setCorrelations ERROR config {} not found", c.Data()));
  }
}
