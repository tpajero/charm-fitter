/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include <PDF_BinFlip.h>

#include <CharmUtils.h>

#include <Utils.h>

#include <RooFormulaVar.h>
#include <RooMultiVarGaussian.h>
#include <RooRealVar.h>

#include <algorithm>
#include <cmath>
#include <format>
#include <iostream>
#include <stdexcept>
#include <vector>

PDF_BinFlip::PDF_BinFlip(const TString measurement_id, const parametrisations::mix mix_param)
    : PDF_Charm{4}, mix_param{mix_param}, measurement_id{measurement_id} {
  name = "BinFlip_" + measurement_id;
  initialise(measurement_id, measurement_id, measurement_id);
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
    theory->add(*(Utils::makeTheoryVar("x_th", "x_th",
                                       "0.5*(  x*cos(phi)*(qop + 1/qop)"
                                       "     + y*sin(phi)*(qop - 1/qop))",
                                       parameters)));
    theory->add(*(Utils::makeTheoryVar("y_th", "y_th",
                                       "0.5*(  y*cos(phi)*(qop + 1./qop)"
                                       "     - x*sin(phi)*(qop - 1./qop))",
                                       parameters)));
    break;
  case mix::theo:
    theory->add(*(Utils::makeTheoryVar("x_th", "x_th", " x12*cos(phiM)", parameters)));
    theory->add(*(Utils::makeTheoryVar("y_th", "y_th", " y12*cos(phiG)", parameters)));
    break;
  default:
    throw std::runtime_error(
        std::format("PDF_BinFlip::initRelations ERROR Parametrisation {} not supported", utils::to_string(mix_param)));
  }
  theory->add(*(Utils::makeTheoryVar("dx_th", "dx_th", utils::dx_expression(mix_param), parameters)));
  theory->add(*(Utils::makeTheoryVar("dy_th", "dy_th", utils::dy_expression(mix_param), parameters)));
}

void PDF_BinFlip::initObservables() {
  TString label;
  if (measurement_id.EqualTo("LHCb_Run1"))
    label = "LHCb Binflip Run 1";
  else if (measurement_id.EqualTo("LHCb_Run2_prompt"))
    label = "LHCb Binflip Run 2 (prompt)";
  else if (measurement_id.EqualTo("LHCb_Run2_sl"))
    label = "LHCb Binflip Run 2 (muon-tagged)";
  else if (measurement_id.EqualTo("LHCb_Run2"))
    label = "LHCb Binflip Run 2";

  observables = new RooArgList("observables");  ///< the order of this list must match that of the COR matrix!
  observables->add(*(new RooRealVar("x_obs", label + "   #it{x_{CP}}", 0., -1e4, 1e4)));
  observables->add(*(new RooRealVar("y_obs", label + "   #it{y_{CP}}", 0., -1e4, 1e4)));
  observables->add(*(new RooRealVar("dx_obs", label + "   #it{#Deltax}", 0., -1e4, 1e4)));
  observables->add(*(new RooRealVar("dy_obs", label + "   #it{#Deltay}", 0., -1e4, 1e4)));
}

void PDF_BinFlip::setObservables(const TString c) {
  if (c.EqualTo("truth"))
    setObservablesTruth();
  else if (c.EqualTo("toy"))
    setObservablesToy();
  else if (c.EqualTo("LHCb_Run1")) {
    obsValSource = "https://inspirehep.net/literature/1724179";
    setObservable("x_obs", 2.7e-3);
    setObservable("y_obs", 7.4e-3);
    setObservable("dx_obs", -0.53e-3);
    setObservable("dy_obs", 0.6e-3);
  } else if (c.EqualTo("LHCb_Run2_prompt")) {
    obsValSource = "https://inspirehep.net/literature/1867376";
    setObservable("x_obs", 3.973e-3);
    setObservable("y_obs", 4.589e-3);
    setObservable("dx_obs", -0.271e-3);
    setObservable("dy_obs", 0.203e-3);
  } else if (c.EqualTo("LHCb_Run2_sl")) {
    obsValSource = "https://inspirehep.net/literature/2135966";
    setObservable("x_obs", 4.29e-3);
    setObservable("y_obs", 12.61e-3);
    setObservable("dx_obs", -0.77e-3);
    setObservable("dy_obs", 3.01e-3);
  } else if (c.EqualTo("LHCb_Run2")) {
    obsValSource = "https://inspirehep.net/literature/2135966";
    setObservable("x_obs", 4.00e-3);
    setObservable("y_obs", 5.51e-3);
    setObservable("dx_obs", -0.29e-3);
    setObservable("dy_obs", 0.31e-3);
  } else {
    throw std::runtime_error(std::format("PDF_BinFlip::setObservables ERROR config {} not found", c.Data()));
  }
}

void PDF_BinFlip::setUncertainties(const TString c) {
  if (c.EqualTo("LHCb_Run1")) {
    obsErrSource = "https://inspirehep.net/literature/1724179";
    StatErr[0] = 1.6e-3;   // x
    StatErr[1] = 3.6e-3;   // y
    StatErr[2] = 0.7e-3;   // dx
    StatErr[3] = 1.6e-3;   // dy
    SystErr[0] = 0.4e-3;   // x
    SystErr[1] = 1.1e-3;   // y
    SystErr[2] = 0.22e-3;  // dx
    SystErr[3] = 0.3e-3;   // dy
  } else if (c.EqualTo("LHCb_Run2_prompt")) {
    obsErrSource = "https://inspirehep.net/literature/1867376";
    StatErr[0] = std::hypot(0.459e-3, 0.29e-3);  // x
    StatErr[1] = std::hypot(1.198e-3, 0.85e-3);  // y
    StatErr[2] = std::hypot(0.182e-3, 0.01e-3);  // dx
    StatErr[3] = std::hypot(0.365e-3, 0.11e-3);  // dy
    std::ranges::fill(SystErr, 0.0);
  } else if (c.EqualTo("LHCb_Run2_sl")) {
    obsErrSource = "https://inspirehep.net/literature/2135966";
    StatErr[0] = 1.48e-3;  // x
    StatErr[1] = 3.12e-3;  // y
    StatErr[2] = 0.93e-3;  // dx
    StatErr[3] = 1.92e-3;  // dy
    SystErr[0] = 0.26e-3;  // x
    SystErr[1] = 0.83e-3;  // y
    SystErr[2] = 0.28e-3;  // dx
    SystErr[3] = 0.26e-3;  // dy
  } else if (c.EqualTo("LHCb_Run2")) {
    obsErrSource = "https://inspirehep.net/literature/2135966";
    StatErr[0] = 0.45e-3;   // x
    StatErr[1] = 1.16e-3;   // y
    StatErr[2] = 0.18e-3;   // dx
    StatErr[3] = 0.35e-3;   // dy
    SystErr[0] = 0.195e-3;  // x
    SystErr[1] = 0.594e-3;  // y
    SystErr[2] = 0.013e-3;  // dx
    SystErr[3] = 0.128e-3;  // dy
  } else {
    throw std::runtime_error(std::format("PDF_BinFlip::setUncertainties ERROR config {} not found", c.Data()));
  }
}

void PDF_BinFlip::setCorrelations(const TString c) {
  resetCorrelations();
  if (c.EqualTo("LHCb_Run1")) {
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
  } else if (c.EqualTo("LHCb_Run2_prompt")) {
    corSource = "https://inspirehep.net/literature/1867376";
    std::vector<double> dataStat = {
        // clang-format off
        1., 0.111,  -0.017, -0.010,  // x
            1.,     -0.011, -0.051,  // y
                     1.,     0.077,  // dx
                             1.      // dy
        // clang-format on
    };
    corStatMatrix = Utils::buildCorMatrix(nObs, dataStat);
  } else if (c.EqualTo("LHCb_Run2_sl")) {
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
  } else if (c.EqualTo("LHCb_Run2")) {
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
  } else {
    throw std::runtime_error(std::format("PDF_BinFlip::setCorrelations ERROR config {} not found", c.Data()));
  }
}

void PDF_BinFlip::buildPdf() {
  pdf = new RooMultiVarGaussian("pdf_" + name, "pdf_" + name, *(RooArgSet*)observables, *(RooArgSet*)theory, covMatrix);
}
