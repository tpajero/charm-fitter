/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include <PDF_BinFlip.h>

#include <CharmUtils.h>

#include <Utils.h>

#include <RooFormulaVar.h>
#include <RooRealVar.h>

#include <algorithm>
#include <cmath>
#include <format>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

PDF_BinFlip::PDF_BinFlip(const TString measurement_id, const parametrisations::mix mix_param)
    : PDF_BinFlip{measurement_id, measurement_id, mix_param} {}

PDF_BinFlip::PDF_BinFlip(const TString obs_id, const TString unc_id, const parametrisations::mix mix_param)
    : PDF_Charm{4}, mix_param{mix_param}, measurement_id{obs_id}, unc_id{unc_id} {
  name = "BinFlip_" + unc_id;
  initialise(obs_id, unc_id, unc_id);
}

std::set<std::string> PDF_BinFlip::getParameterNames() const {
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    return {"x", "y", "qop", "phi"};
  case mix::theo:
    return {"phiG", "x12", "y12", "phiM"};
  default:
    throw std::runtime_error(std::format("PDF_BinFlip::getParameterNames ERROR Parametrisation {} not supported",
                                         utils::to_string(mix_param)));
  }
}

void PDF_BinFlip::initRelations() {
  theory = new RooArgList("theory");  ///< the order of this list must match that of the COR matrix!
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    theory->add(*(Utils::makeTheoryVar("xCP_th",
                                       "0.5*(  x*cos(phi)*(qop + 1/qop)"
                                       "     + y*sin(phi)*(qop - 1/qop))",
                                       parameters)));
    theory->add(*(Utils::makeTheoryVar("yCP_th",
                                       "0.5*(  y*cos(phi)*(qop + 1./qop)"
                                       "     - x*sin(phi)*(qop - 1./qop))",
                                       parameters)));
    break;
  case mix::theo:
    theory->add(*(Utils::makeTheoryVar("xCP_th", " x12*cos(phiM)", parameters)));
    theory->add(*(Utils::makeTheoryVar("yCP_th", " y12*cos(phiG)", parameters)));
    break;
  default:
    throw std::runtime_error(
        std::format("PDF_BinFlip::initRelations ERROR Parametrisation {} not supported", utils::to_string(mix_param)));
  }
  theory->add(*(Utils::makeTheoryVar("dx_th", utils::dx_expression(mix_param), parameters)));
  theory->add(*(Utils::makeTheoryVar("dy_th", utils::dy_expression(mix_param), parameters)));
}

void PDF_BinFlip::initObservables() {
  static const std::map<std::string, std::string> labels = {
      {"LHCb-R1", "LHCb Binflip Run 1"},
      {"LHCb-R2-prompt", "LHCb Binflip Run 2 (prompt)"},
      {"LHCb-R2-SL", "LHCb Binflip Run 2 (muon-tagged)"},
      {"LHCb-R2", "LHCb Binflip Run 2"},
  };
  const auto it = labels.find(measurement_id.Data());
  const TString label = it != labels.end() ? it->second : name;

  observables = new RooArgList("observables");  ///< the order of this list must match that of the COR matrix!
  observables->add(*(new RooRealVar("xCP_obs", label + "   #it{x_{CP}}", 0., -1e4, 1e4)));
  observables->add(*(new RooRealVar("yCP_obs", label + "   #it{y_{CP}}", 0., -1e4, 1e4)));
  observables->add(*(new RooRealVar("dx_obs", label + "   #it{#Deltax}", 0., -1e4, 1e4)));
  observables->add(*(new RooRealVar("dy_obs", label + "   #it{#Deltay}", 0., -1e4, 1e4)));
}

void PDF_BinFlip::setObservables(const TString c) {
  if (c.EqualTo("truth"))
    setObservablesTruth();
  else if (c.EqualTo("toy"))
    setObservablesToy();
  else if (c.EqualTo("LHCb-R1")) {
    obsValSource = "https://inspirehep.net/literature/1724179";
    setObservable("xCP_obs", 2.7e-3);
    setObservable("yCP_obs", 7.4e-3);
    setObservable("dx_obs", -0.53e-3);
    setObservable("dy_obs", 0.6e-3);
  } else if (c.EqualTo("LHCb-R2-prompt")) {
    obsValSource = "https://inspirehep.net/literature/1867376; https://cds.cern.ch/record/2713310 Table 28";
    setObservable("xCP_obs", 3.973e-3);
    setObservable("yCP_obs", 4.589e-3);
    setObservable("dx_obs", -0.271e-3);
    setObservable("dy_obs", 0.203e-3);
  } else if (c.EqualTo("LHCb-R2-SL")) {
    obsValSource = "https://inspirehep.net/literature/2135966";
    setObservable("xCP_obs", 4.29e-3);
    setObservable("yCP_obs", 12.61e-3);
    setObservable("dx_obs", -0.77e-3);
    setObservable("dy_obs", 3.01e-3);
  } else if (c.EqualTo("LHCb-R2")) {
    obsValSource = "https://inspirehep.net/literature/2135966";
    setObservable("xCP_obs", 4.01e-3);
    setObservable("yCP_obs", 5.51e-3);
    setObservable("dx_obs", -0.29e-3);
    setObservable("dy_obs", 0.31e-3);
  } else {
    throw std::runtime_error(std::format("PDF_BinFlip::setObservables ERROR config {} not found", c.Data()));
  }
}

void PDF_BinFlip::setUncertainties(const TString c) {
  const std::vector<double> stat_lhcb_run2 = {0.45e-3, 1.16e-3, 0.18e-3, 0.35e-3};
  const std::vector<double> syst_lhcb_run2 = {0.195e-3, 0.594e-3, 0.013e-3, 0.128e-3};
  // x, y, dx, dy
  if (c.EqualTo("LHCb-R1")) {
    obsErrSource = "https://inspirehep.net/literature/1724179";
    StatErr = {1.6e-3, 3.6e-3, 0.70e-3, 1.6e-3};
    SystErr = {0.4e-3, 1.1e-3, 0.22e-3, 0.3e-3};
  } else if (c.EqualTo("LHCb-R2-prompt")) {
    obsErrSource = "https://inspirehep.net/literature/1867376; https://cds.cern.ch/record/2713310 Table 28";
    StatErr = {0.459e-3, 1.198e-3, 0.182e-3, 0.365e-3};
    SystErr = {0.29e-3, 0.85e-3, 0.01e-3, 0.11e-3};
  } else if (c.EqualTo("LHCb-R2-SL")) {
    obsErrSource = "https://inspirehep.net/literature/2135966";
    StatErr = {1.48e-3, 3.12e-3, 0.93e-3, 1.92e-3};
    SystErr = {0.26e-3, 0.83e-3, 0.28e-3, 0.26e-3};
  } else if (c.EqualTo("LHCb-R2")) {
    obsErrSource = "https://inspirehep.net/literature/2135966";
    StatErr = stat_lhcb_run2;
    SystErr = syst_lhcb_run2;
  } else if (c.EqualTo("LHCb-UI") || c.EqualTo("LHCb-UII")) {
    obsErrSource = "charm-fitter";
    // Run 2 values times scale factor
    using constants::lhcb_extrapolations;
    const auto scale = lhcb_extrapolations.at(c.Data()) / lhcb_extrapolations.at("LHCb-R2");
    std::transform(stat_lhcb_run2.begin(), stat_lhcb_run2.end(), StatErr.begin(),
                   [scale](double x) { return x * scale; });
    std::ranges::fill(SystErr, 0.0);
  } else {
    throw std::runtime_error(std::format("PDF_BinFlip::setUncertainties ERROR config {} not found", c.Data()));
  }
}

void PDF_BinFlip::setCorrelations(const TString c) {
  resetCorrelations();
  if (c.EqualTo("LHCb-R1")) {
    corSource = "https://inspirehep.net/literature/1724179";
    std::vector<double> dataStat = {
        // clang-format off
        1., -0.17, 0.04, -0.02,  // x
             1.,  -0.03,  0.01,  // y
                   1.,   -0.13,  // dx
                          1.     // dy
        // clang-format on
    };
    corStatMatrix = Utils::buildCorMatrix(nObs, dataStat);
    std::vector<double> dataSyst = {
        // clang-format off
        1., 0.15,  0.01, -0.02,  // x
            1.,   -0.05, -0.03,  // y
                   1.,    0.14,  // dx
                          1.     // dy
        // clang-format on
    };
    corSystMatrix = Utils::buildCorMatrix(nObs, dataSyst);
  } else if (c.EqualTo("LHCb-R2-prompt")) {
    corSource = "https://inspirehep.net/literature/1867376";
    std::vector<double> dataStat = {
        // clang-format off
        1., 0.11, -0.02, -0.01,  // x
            1.,   -0.01, -0.05,  // y
                   1.,    0.08,  // dx
                          1.     // dy
        // clang-format on
    };
    corStatMatrix = Utils::buildCorMatrix(nObs, dataStat);
    std::vector<double> dataSyst = {
        // clang-format off
        1., 0.13,  0.01, 0.01,  // x
            1.,   -0.02, 0.01,  // y
                   1.,   0.31,  // dx
                         1.     // dy
        // clang-format on
    };
    corSystMatrix = Utils::buildCorMatrix(nObs, dataSyst);
  } else if (c.EqualTo("LHCb-R2-SL")) {
    corSource = "https://inspirehep.net/literature/2135966";
    std::vector<double> dataStat = {
        // clang-format off
        1., 0.085, -0.011, -0.009,  // x
            1.,    -0.001, -0.050,  // y
                    1.,     0.070,  // dx
                            1.      // dy
        // clang-format on
    };
    corStatMatrix = Utils::buildCorMatrix(nObs, dataStat);
    std::vector<double> dataSyst = {
        // clang-format off
        1., 0.11,  -0.25, -0.02,  // x
            1.,    -0.05, -0.20,  // y
                    1.,    0.11,  // dx
                           1.     // dy
        // clang-format on
    };
    corSystMatrix = Utils::buildCorMatrix(nObs, dataSyst);
  } else if (c.EqualTo("LHCb-R2")) {
    corSource = "https://inspirehep.net/literature/2135966";
    std::vector<double> dataStat = {
        // clang-format off
        1., 0.121,  -0.018, -0.016,  // x
            1.,     -0.012, -0.058,  // y
                     1.,     0.069,  // dx
                             1.      // dy
        // clang-format on
    };
    corStatMatrix = Utils::buildCorMatrix(nObs, dataStat);
    std::vector<double> dataSyst = {
        // clang-format off
        1., 0.08,  0.,   -0.01,  // x
            1.,   -0.02, -0.04,  // y
                   1.,    0.33,  // dx
                          1.     // dy
        // clang-format on
    };
    corSystMatrix = Utils::buildCorMatrix(nObs, dataSyst);
  } else if (c.EqualTo("LHCb-UII")) {
    corSource = "";
    resetCorrelations();
    // Run 2 values
    std::vector<double> dataStat = {
        // clang-format off
        1., 0.121,  -0.018, -0.016,  // x
            1.,     -0.012, -0.058,  // y
                     1.,     0.069,  // dx
                             1.      // dy
        // clang-format on
    };
    corStatMatrix = Utils::buildCorMatrix(nObs, dataStat);
  } else {
    throw std::runtime_error(std::format("PDF_BinFlip::setCorrelations ERROR config {} not found", c.Data()));
  }
}
