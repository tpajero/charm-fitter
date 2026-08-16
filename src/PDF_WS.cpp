/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include <PDF_WS.h>

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

namespace {
  using parametrisations::mix;
  // Expressions for the observables in the various parametrisations
  const std::map<std::string, std::map<mix, std::string>> theory_expressions = {
      {"y'+",
       {
           {mix::pheno, "qop*(  y * cos(Delta_Kpi - phi)"
                        "     + x * sin(Delta_Kpi - phi))"},
           {mix::theo, "  y12 * cos(Delta_Kpi + phiG)"
                       "+ x12 * sin(Delta_Kpi + phiM)"},
           {mix::d0_to_kpi, "yp + dyp"},
       }},
      {"y'-",
       {
           {mix::pheno, "1/qop*(  y * cos(Delta_Kpi + phi)"
                        "       + x * sin(Delta_Kpi + phi))"},
           {mix::theo, "  y12 * cos(Delta_Kpi - phiG)"
                       "+ x12 * sin(Delta_Kpi - phiM)"},
           {mix::d0_to_kpi, "yp - dyp"},
       }},
      {"x'2+",
       {
           {mix::pheno, "TMath::Sq(qop*(  x * cos(Delta_Kpi - phi)"
                        "               - y * sin(Delta_Kpi - phi)))"},
           {mix::theo, "TMath::Sq(- y12 * sin(Delta_Kpi + phiG)"
                       "          + x12 * cos(Delta_Kpi + phiM))"},
           {mix::d0_to_kpi, "xp2 + dxp2"},
       }},
      {"x'2-",
       {
           {mix::pheno, "TMath::Sq(1/qop*(  x * cos(Delta_Kpi + phi)"
                        "                 - y * sin(Delta_Kpi + phi)))"},
           {mix::theo, "TMath::Sq(- y12 * sin(Delta_Kpi - phiG)"
                       "          + x12 * cos(Delta_Kpi - phiM))"},
           {mix::d0_to_kpi, "xp2 - dxp2"},
       }},
      {"c",
       {
           {mix::pheno, "0.5 * (  qop     * (y*cos(Delta_Kpi - phi) + x*sin(Delta_Kpi - phi)) "
                        "       + 1 / qop * (y*cos(Delta_Kpi + phi) + x*sin(Delta_Kpi + phi)))"},
           {mix::theo, "y12 * cos(Delta_Kpi) * cos(phiG) + x12 * sin(Delta_Kpi) * cos(phiM)"},
           {mix::d0_to_kpi, "yp"},
       }},
      {"c'",
       {
           {mix::pheno, "0.125 * (x*x + y*y) * (qop*qop + 1 / (qop*qop))"},
           {mix::theo, "0.25 * (x12*x12 + y12*y12)"
                       "+ 0.25 * r_Kpi * r_Kpi * (y12*y12 - x12*x12)"},  // 2nd order corrections
           {mix::d0_to_kpi, "(yp*yp + xp2) / 4"},
       }},
      {"dc",
       {
           {mix::pheno, "0.5 * (      qop*(  y*cos(Delta_Kpi - phi) + x*sin(Delta_Kpi - phi)) "
                        "       - 1 / qop*(  y*cos(Delta_Kpi + phi) + x*sin(Delta_Kpi + phi)))"},
           {mix::theo, "  x12 * cos(Delta_Kpi) * sin(phiM)"
                       "- y12 * sin(Delta_Kpi) * sin(phiG)"},
           {mix::d0_to_kpi, "dyp"},
       }},
      {"dc'",
       {
           {mix::pheno, "1 / 8 * (x*x + y*y) * (qop*qop - 1 / (qop*qop))"},
           {mix::theo, "0.5 * x12 * y12 * sin(phiM - phiG)"},
           {mix::d0_to_kpi, "(2 * yp * dyp + dyp*dyp + dxp2) / 4"},
       }},
  };

  std::string get_formula(const std::string observable, const mix mix_param) {
    try {
      return theory_expressions.at(observable).at(mix_param);
    } catch (const std::out_of_range& e) {
      std::cerr << std::format("Out of range error, parametrisation {} not handled for observable {}: {}",
                               utils::to_string(mix_param), observable, e.what())
                << std::endl;
      throw;
    }
  }

  const std::map<std::string, std::string> labels = {
      {"BaBar", "WS/RS BaBar CPV"},
      {"Belle", "WS/RS Belle CPV"},
      {"LHCb_DT_Run1", "WS/RS LHCb Run 1 DT"},
      {"LHCb_Run1", "WS/RS LHCb Run 1"},
      {"LHCb_Prompt_2011_2016", "WS/RS LHCb (11/16) prompt"},
      {"LHCb_Prompt_Run12_sec9", "WS/RS LHCb Run 1+2 prompt Sec. 9"},
      {"LHCb_Prompt_Run12_appB", "WS/RS LHCb Run 1+2 prompt App. B"},
      {"LHCb_DT_Run2", "WS/RS LHCb Run 2 DT"},
      {"LHCb_DT_Run12", "WS/RS LHCb Run 1+2 DT"},
  };
}  // namespace

PDF_WS::PDF_WS(const TString measurement_id, const parametrisations::mix mix_param, parametrisations::kpi p)
    : PDF_Charm{measurement_id.EqualTo("LHCb_Prompt_Run12_appB") ? 9 : 6}, mix_param{mix_param}, ws_param{p} {
  try {
    label = labels.at(measurement_id.Data());
  } catch (const std::out_of_range&) {
    throw std::runtime_error(
        std::format("PDF_WS::PDF_WS ERROR Measurement ID {} not supported", measurement_id.Data()));
  }

  if (ws_param == parametrisations::kpi::ccprime && !measurement_id.BeginsWith("LHCb_Prompt_Run12")) {
    throw std::runtime_error(
        "PDF_WS::PDF_WS ERROR The c/c' parametrisation was introduced only with the LHCb Run 2 measurement");
  }

  name = "WS_" + measurement_id;
  initialise(measurement_id, measurement_id, measurement_id);
}

PDF_WS::PDF_WS(const TString val, TString err, const parametrisations::mix mix_param)
    : PDF_Charm{6}, mix_param{mix_param} {
  if (err.EqualTo("LHCb_Run12"))
    label = "WS/RS LHCb prompt (Run 1+2)";
  else
    throw std::runtime_error(std::format("PDF_WS::PDF_WS ERROR Measurement ID {} not supported", err.Data()));

  name = "WS_" + err;
  initialise(val, err, err);
}

std::set<std::string> PDF_WS::getParameterNames() const {
  using parametrisations::mix;
  std::set<std::string> names = {"r_Kpi", "Acp_KP"};
  if (mix_param != mix::d0_to_kpi) names.insert("Delta_Kpi");
  if (nObs == 9) names.insert("Acp_KK");
  switch (mix_param) {
  case mix::pheno:
    names.insert({"x", "y", "qop", "phi"});
    break;
  case mix::theo:
    names.insert({"phiG", "x12", "y12", "phiM"});
    break;
  case mix::d0_to_kpi:
    names.insert({"yp", "dyp", "xp2", "dxp2"});
    break;
  default:
    throw std::runtime_error(
        std::format("PDF_WS::getParameterNames ERROR Parametrisation {} not supported", utils::to_string(mix_param)));
  }
  return names;
}

void PDF_WS::initRelations() {
  using parametrisations::kpi;
  switch (ws_param) {
  case kpi::raxy:
    initRelationsRAXY();
    break;
  case kpi::rrxy:
    initRelationsRRXY();
    break;
  case kpi::ccprime:
    initRelationsCCPrime();
    break;
  default:
    throw std::runtime_error(
        std::format("PDF_WS::initRelations ERROR WS parametrisation {} not supported", static_cast<int>(ws_param)));
  }
}

void PDF_WS::initRelationsCCPrime() {
  using parametrisations::dy_fsc;
  theory = new RooArgList("theory");
  theory->add(*(Utils::makeTheoryVar("RD_th", "r_Kpi * r_Kpi", parameters)));
  theory->add(*(Utils::makeTheoryVar("c_th", get_formula("c", mix_param), parameters)));
  theory->add(*(Utils::makeTheoryVar("c'_th", get_formula("c'", mix_param), parameters)));
  theory->add(*(Utils::makeTheoryVar("AD_th", "Acp_KP", parameters)));
  theory->add(*(Utils::makeTheoryVar("dc_th", get_formula("dc", mix_param), parameters)));
  theory->add(*(Utils::makeTheoryVar("dc'_th", get_formula("dc'", mix_param), parameters)));
  if (nObs == 9) {
    theory->add(*(Utils::makeTheoryVar("ADt_th", "Acp_KP - 2 * Acp_KK", parameters)));
    theory->add(*(Utils::makeTheoryVar(
        "dc~_th",
        std::format("{} - 2 * r_Kpi * ({}) - Acp_KK * ({})", get_formula("dc", mix_param),
                    utils::dy_hh_expression(mix_param, dy_fsc::none, "KK"), get_formula("c", mix_param)),
        parameters)));
    theory->add(*(Utils::makeTheoryVar("dc'~_th",
                                       std::format("{} - 2 * r_Kpi * ({}) * ({}) - 2 * Acp_KK * ({})",
                                                   get_formula("dc'", mix_param), get_formula("c", mix_param),
                                                   utils::dy_hh_expression(mix_param, dy_fsc::none, "KK"),
                                                   get_formula("c'", mix_param)),
                                       parameters)));
  }
}

void PDF_WS::initRelationsRAXY() {
  theory = new RooArgList("theory");
  theory->add(*(Utils::makeTheoryVar("RD_th", "r_Kpi * r_Kpi", parameters)));
  theory->add(*(Utils::makeTheoryVar("y'+_th", get_formula("y'+", mix_param), parameters)));
  theory->add(*(Utils::makeTheoryVar("x'2+_th", get_formula("x'2+", mix_param), parameters)));
  theory->add(*(Utils::makeTheoryVar("AD_th", "Acp_KP", parameters)));
  theory->add(*(Utils::makeTheoryVar("y'-_th", get_formula("y'-", mix_param), parameters)));
  theory->add(*(Utils::makeTheoryVar("x'2-_th", get_formula("x'2-", mix_param), parameters)));
}

void PDF_WS::initRelationsRRXY() {
  theory = new RooArgList("theory");
  theory->add(*(Utils::makeTheoryVar("RD_p_th", "r_Kpi * r_Kpi * (1 + Acp_KP)", parameters)));
  theory->add(*(Utils::makeTheoryVar("y'+_th", get_formula("y'+", mix_param), parameters)));
  theory->add(*(Utils::makeTheoryVar("x'2+_th", get_formula("x'2+", mix_param), parameters)));
  theory->add(*(Utils::makeTheoryVar("RD_m_th", "r_Kpi * r_Kpi * (1 - Acp_KP)", parameters)));
  theory->add(*(Utils::makeTheoryVar("y'-_th", get_formula("y'-", mix_param), parameters)));
  theory->add(*(Utils::makeTheoryVar("x'2-_th", get_formula("x'2-", mix_param), parameters)));
}

void PDF_WS::initObservables() {
  observables = new RooArgList("observables");  // the order of this list must match that of the COR matrix!
  using parametrisations::kpi;
  switch (ws_param) {
  case kpi::raxy:
    observables->add(*(new RooRealVar("RD_obs", label + "   #it{R_{K#pi}}", 0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("y'+_obs", label + "   #it{y'^{+}}", 0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("x'2+_obs", label + "   #it{x'^{+2}}", 0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("AD_obs", label + "   #it{A_{K#pi}}", 0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("y'-_obs", label + "   #it{y'}^{#minus}", 0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("x'2-_obs", label + "   #it{x'}^{#minus2}", 0., -1e4, 1e4)));
    break;
  case kpi::rrxy:
    observables->add(*(new RooRealVar("RD_p_obs", label + "   #it{R_{K#pi}^{+}}", 0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("y'+_obs", label + "   #it{y'^{+}}", 0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("x'2+_obs", label + "   #it{x'^{+2}}", 0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("RD_m_obs", label + "   #it{R_{K#pi}^{#minus}}", 0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("y'-_obs", label + "   #it{y'}^{#minus}", 0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("x'2-_obs", label + "   #it{x'}^{#minus2}", 0., -1e4, 1e4)));
    break;
  case kpi::ccprime:
    observables->add(*(new RooRealVar("RD_obs", label + "   #it{R_{K#pi}}", 0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("c_obs", label + "   #it{c_{K#pi}}", 0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("c'_obs", label + "   #it{c'_{K#pi}}", 0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("AD_obs", label + "   #it{A_{K#pi}}", 0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("dc_obs", label + "   #it{#Deltac_{K#pi}}", 0., -1e4, 1e4)));
    observables->add(*(new RooRealVar("dc'_obs", label + "   #it{#Deltac'_{K#pi}}", 0., -1e4, 1e4)));
    if (nObs == 9) {
      observables->add(*(new RooRealVar("ADt_obs", label + "   #it{#tilde{A}_{K#pi}}", 0., -1e4, 1e4)));
      observables->add(*(new RooRealVar("dc~_obs", label + "   #it{#Delta#tilde{c}_{K#pi}}", 0., -1e4, 1e4)));
      observables->add(*(new RooRealVar("dc'~_obs", label + "   #it{#Delta#tilde{c}'_{K#pi}}", 0., -1e4, 1e4)));
    }
    break;
  default:
    throw std::runtime_error(
        std::format("PDF_WS::initObservables ERROR WS parametrisation {} not supported", static_cast<int>(ws_param)));
  }
}

void PDF_WS::setObservables(const TString c) {
  if (c.EqualTo("truth"))
    setObservablesTruth();
  else if (c.EqualTo("toy"))
    setObservablesToy();
  else if (c.EqualTo("BaBar")) {
    obsValSource = "https://inspirehep.net/literature/746245";
    // setObservable("RD_obs", 0.00303); TODO
    setObservable("RD_p_obs", 2.97e-3);
    setObservable("y'+_obs", 9.8e-3);
    setObservable("x'2+_obs", -2.4e-4);
    // setObservable("AD_obs", -2.1e-2); TODO
    setObservable("RD_m_obs", 3.09e-3);
    setObservable("y'-_obs", 9.6e-3);
    setObservable("x'2-_obs", -2.0e-4);
  } else if (c.EqualTo("Belle")) {
    obsValSource = "http://belle.kek.jp/belle/theses/doctor/lmzhang06/phd-mix-400.ps.gz";
    setObservable("RD_p_obs", 3.73e-3);
    setObservable("y'+_obs", -1.2e-3);
    setObservable("x'2+_obs", 3.2e-4);
    setObservable("RD_m_obs", 3.56e-3);
    setObservable("y'-_obs", 2e-3);
    setObservable("x'2-_obs", 0.6e-4);
  } else if (c.EqualTo("LHCb_DT_Run1")) {
    obsValSource = "https://inspirehep.net/literature/1499047";
    setObservable("RD_p_obs", 3.38e-3);
    setObservable("y'+_obs", 5.81e-3);
    setObservable("x'2+_obs", -1.9e-5);
    setObservable("RD_m_obs", 3.60e-3);
    setObservable("y'-_obs", 3.32e-3);
    setObservable("x'2-_obs", 7.9e-5);
  } else if (c.EqualTo("LHCb_Run1")) {
    obsValSource = "https://inspirehep.net/literature/1499047";
    setObservable("RD_p_obs", 3.474e-3);
    setObservable("y'+_obs", 5.97e-3);
    setObservable("x'2+_obs", 1.1e-5);
    setObservable("RD_m_obs", 3.591e-3);
    setObservable("y'-_obs", 4.50e-3);
    setObservable("x'2-_obs", 6.1e-5);
  } else if (c.EqualTo("LHCb_Prompt_2011_2016")) {
    obsValSource = "https://inspirehep.net/literature/1642234";
    setObservable("RD_p_obs", 3.454e-3);
    setObservable("y'+_obs", 5.01e-3);
    setObservable("x'2+_obs", 6.1e-5);
    setObservable("RD_m_obs", 3.454e-3);
    setObservable("y'-_obs", 5.54e-3);
    setObservable("x'2-_obs", 1.6e-5);
  } else if (c.EqualTo("LHCb_Prompt_Run12_sec9")) {
    obsValSource = "https://inspirehep.net/literature/2871248 Tab III";
    setObservable("RD_obs", 3.427e-3);
    setObservable("c_obs", 5.28e-3);
    setObservable("c'_obs", 1.20e-5);
    setObservable("AD_obs", -0.66e-2);
    setObservable("dc_obs", 2.0e-4);
    setObservable("dc'_obs", -0.7e-6);
  } else if (c.EqualTo("LHCb_Prompt_Run12_appB")) {
    obsValSource = "https://inspirehep.net/literature/2871248 Tab IV";
    setObservable("RD_obs", 3.427e-3);
    setObservable("c_obs", 5.28e-3);
    setObservable("c'_obs", 1.20e-5);
    setObservable("AD_obs", -0.9e-2);
    setObservable("dc_obs", -0.1e-3);
    setObservable("dc'_obs", 4.6e-6);
    setObservable("ADt_obs", -0.82e-2);
    setObservable("dc~_obs", 3.2e-4);
    setObservable("dc'~_obs", -2.0e-6);
  } else if (c.EqualTo("LHCb_DT_Run2")) {
    obsValSource = "https://inspirehep.net/literature/2871248";
    setObservable("RD_p_obs", 3.55e-3);
    setObservable("y'+_obs", 3.56e-3);
    setObservable("x'2+_obs", 1.086e-4);
    setObservable("RD_m_obs", 3.39e-3);
    setObservable("y'-_obs", 8.11e-3);
    setObservable("x'2-_obs", -1.129e-4);
  } else if (c.EqualTo("LHCb_DT_Run12")) {
    obsValSource = "https://inspirehep.net/literature/2871248";
    setObservable("RD_p_obs", 3.50e-3);
    setObservable("y'+_obs", 4.14e-3);
    setObservable("x'2+_obs", 7.84e-5);
    setObservable("RD_m_obs", 0.00344);
    setObservable("y'-_obs", 6.81e-3);
    setObservable("x'2-_obs", -4.86e-5);
  } else {
    throw std::runtime_error(std::format("PDF_WS::setObservables ERROR config {} not found", c.Data()));
  }
}

void PDF_WS::setUncertainties(const TString c) {
  if (c.EqualTo("BaBar")) {
    obsErrSource = "https://inspirehep.net/literature/746245";
    // StatErr[0] = 0.000189; // RD TODO
    // StatErr[3] = 5.4e-2;    // AD TODO
    //         RD+       y'+     x'2+    RD-       y'-     x'2-
    StatErr = {0.267e-3, 7.8e-3, 5.2e-4, 0.267e-3, 7.5e-3, 5.0e-4};
    std::ranges::fill(SystErr, 0.0);
  } else if (c.EqualTo("Belle")) {
    obsErrSource = "http://belle.kek.jp/belle/theses/doctor/lmzhang06/phd-mix-400.ps.gz";
    //         RD+      y'+     x'2+    RD-      y'-     x'2-
    StatErr = {0.24e-3, 5.7e-3, 3.1e-4, 0.24e-3, 5.4e-3, 2.9e-4};
    std::ranges::fill(SystErr, 0.0);
  } else if (c.EqualTo("LHCb_DT_Run1")) {
    obsErrSource = "https://inspirehep.net/literature/1499047";
    StatErr = {std::hypot(1.5e-3, 0.6e-3),     // RD+
               std::hypot(5.25e-3, 0.32e-3),   // y'+
               std::hypot(4.46e-4, 0.31e-4),   // x'2+
               std::hypot(1.5e-3, 0.7e-3),     // RD-
               std::hypot(5.21e-3, 0.40e-3),   // y'-
               std::hypot(4.31e-4, 0.38e-4)};  // x'2-
    std::ranges::fill(SystErr, 0.0);
  } else if (c.EqualTo("LHCb_Run1")) {
    obsErrSource = "https://inspirehep.net/literature/1499047";
    //         RD+      y'+      x'2+    RD-      y'-      x'2-
    StatErr = {0.081e-3, 1.25e-3, 6.5e-5, 0.081e-3, 1.21e-3, 6.1e-5};
    std::ranges::fill(SystErr, 0.0);
  } else if (c.EqualTo("LHCb_Prompt_2011_2016")) {
    obsErrSource = "https://inspirehep.net/literature/1642234";
    //         RD+       y'+     x'2+    RD-       y'-     x'2-
    StatErr = {0.045e-3, 7.4e-4, 3.7e-5, 0.045e-3, 7.4e-4, 3.9e-5};
    std::ranges::fill(SystErr, 0.0);
  } else if (c.EqualTo("LHCb_Prompt_Run12_sec9")) {
    obsErrSource = "https://inspirehep.net/literature/2871248 Tab III";
    //         RD        c       c'      AD       dc      dc'
    StatErr = {0.019e-3, 3.3e-4, 3.5e-6, 0.57e-2, 3.4e-4, 3.6e-6};
    std::ranges::fill(SystErr, 0.0);
  } else if (c.EqualTo("LHCb_Prompt_Run12_appB")) {
    obsErrSource = "https://inspirehep.net/literature/2871248 Tab IV";
    //         RD        c       c'      AD      dc      dc'     ADt      dc~     dc'~
    StatErr = {0.019e-3, 3.3e-4, 3.5e-6, 2.0e-2, 1.0e-3, 9.8e-6, 0.59e-2, 3.6e-4, 3.8e-6};
    std::ranges::fill(SystErr, 0.0);
  } else if (c.EqualTo("LHCb_DT_Run2")) {
    obsErrSource = "https://inspirehep.net/literature/2871248";
    //         RD+       y'+     x'2+      RD-      y'-      x'2-
    StatErr = {0.08e-3, 2.25e-3, 1.623e-4, 0.08e-3, 2.36e-3, 1.859e-4};
    std::ranges::fill(SystErr, 0.0);
  } else if (c.EqualTo("LHCb_DT_Run12")) {
    obsErrSource = "https://inspirehep.net/literature/2871248";
    //        RD+       y'+     x'2+      RD-      y'-      x'2-
    StatErr = {0.07e-3, 2.04e-3, 1.522e-4, 0.07e-3, 2.11e-3, 1.665e-4};
    std::ranges::fill(SystErr, 0.0);
  } else {
    throw std::runtime_error(std::format("PDF_WS::setUncertainties ERROR config {} not found", c.Data()));
  }
}

void PDF_WS::setCorrelations(const TString c) {
  resetCorrelations();
  if (c.EqualTo("BaBar")) {  // TODO
    corSource = "https://hflav-eos.web.cern.ch/hflav-eos/charm/CKM23/results_mix_cpv.html";
    std::vector<double> data = {
        // clang-format off
        // RD+   y'+    x'2+   RD-  y'-   x'2-
           1.,  -0.87,  0.77,  0.,  0.,   0.,    // RD+
                 1.,   -0.94,  0.,  0.,   0.,    // y'+
                        1.,    0.,  0.,   0.,    // x'2+
                               1., -0.87, 0.77,  // RD-
                                    1.,  -0.94,  // y'-
                                          1.     // x'2-
        // clang-format on
    };
    corStatMatrix = Utils::buildCorMatrix(nObs, data);
  } else if (c.EqualTo("Belle")) {
    corSource = "http://belle.kek.jp/belle/theses/doctor/lmzhang06/phd-mix-400.ps.gz";
    std::vector<double> data = {
        // clang-format off
        // RD+  y'+    x'2+   RD-  y'-    x'2-
        1.,    -0.834, 0.655, 0.,  0.,    0.,     // RD+
                1.,   -0.909, 0.,  0.,    0.,     // y'+
                       1.,    0.,  0.,    0.,     // x'2+
                              1., -0.834, 0.655,  // RD-
                                   1.,   -0.909,  // y'-
                                          1.      // x'2-
        // clang-format on
    };
    corStatMatrix = Utils::buildCorMatrix(nObs, data);
  } else if (c.EqualTo("LHCb_DT_Run1")) {
    corSource = "https://inspirehep.net/literature/1499047";
    std::vector<double> data = {
        // clang-format off
        // RD+  y'+    x'2+    RD-    y'-    x'2-
        1.,    -0.732, 0.625, -0.008, 0.,    0.,     // RD+
                1.,   -0.963,  0.,    0.,    0.,     // y'+
                       1.,     0.,    0.,    0.,     // x'2+
                               1.,   -0.707, 0.602,  // RD-
                                      1.,   -0.958,  // y'-
                                             1.      // x'2-
        // clang-format on
    };
    corStatMatrix = Utils::buildCorMatrix(nObs, data);
  } else if (c.EqualTo("LHCb_Run1")) {
    corSource = "https://inspirehep.net/literature/1499047";
    std::vector<double> data = {
        // clang-format off
        //  RD+     y'+    x'2+     RD-     y'-    x'2-
        1., -0.920, 0.823, -0.007, -0.010,  0.008,  // RD+
             1.,   -0.962, -0.011,  0.000, -0.002,  // y'+
                    1.,     0.009, -0.002,  0.004,  // x'2+
                            1.,    -0.918,  0.812,  // RD-
                                    1.,    -0.956,  // y'-
                                            1.      // x'2-
        // clang-format on
    };
    corStatMatrix = Utils::buildCorMatrix(nObs, data);
  } else if (c.EqualTo("LHCb_Prompt_2011_2016")) {
    corSource = "https://inspirehep.net/literature/1642234";
    std::vector<double> data = {
        // clang-format off
        //  RD+     y'+    x'2+     RD-     y'-    x'2-
        1., -0.935, 0.843, -0.012, -0.003, -0.002,  // RD+
             1.,   -0.963, -0.003,  0.004, -0.003,  // y'+
                    1.,     0.002, -0.003,  0.003,  // x'2+
                            1.,    -0.935,  0.846,  // RD-
                                    1.,    -0.964,  // y'-
                                            1.      // x'2-
        // clang-format on
    };
    corStatMatrix = Utils::buildCorMatrix(nObs, data);
  } else if (c.EqualTo("LHCb_Prompt_Run12_sec9")) {
    corSource = "https://inspirehep.net/literature/2871248 Tab III";
    std::vector<double> data = {
        // clang-format off
        // RD  c      c'      AD      c       c'
        1.,   -0.927, 0.803,  0.009, -0.007,  0.002,  // RD
               1.,   -0.942, -0.013,  0.012, -0.007,  // c
                      1.,     0.007, -0.007,  0.002,  // c'
                              1.,    -0.919,  0.797,  // AD
                                      1.,    -0.941,  // dc
                                              1.      // dc'
        // clang-format on
    };
    corStatMatrix = Utils::buildCorMatrix(nObs, data);
  } else if (c.EqualTo("LHCb_Prompt_Run12_appB")) {
    corSource = "https://inspirehep.net/literature/2871248 Tab IV";
    std::vector<double> data = {
        // clang-format off
        // RD  c      c'      AD      c       c'      ADt     ct      c't
        1.,   -0.927, 0.803,  0.003, -0.002,  0.002,  0.008, -0.007,  0.000,  // RD
               1.,   -0.943, -0.005,  0.004, -0.004, -0.014,  0.013, -0.006,  // c
                      1.,     0.003, -0.003,  0.003,  0.007, -0.006,  0.000,  // c'
                              1.,    -0.938,  0.811,  0.,     0.,     0.,     // AD
                                      1.,    -0.943,  0.,     0.,     0.,     // dc
                                              1.,     0.,     0.,     0.,     // dc'
                                                      1.,    -0.934,  0.810,  // ADt
                                                              1.,    -0.943,  // dc~
                                                                      1.      // dc'~
        // clang-format on
    };
    corStatMatrix = Utils::buildCorMatrix(nObs, data);
  } else if (c.EqualTo("LHCb_DT_Run2")) {
    corSource = "https://inspirehep.net/literature/2871248";
    std::vector<double> data = {
        // clang-format off
        // RD+  y'+    x'2+    RD-     y'-    x'2-
        1.,    -0.768, 0.639, -0.015, -0.001, 0.,     // RD+
                1.,   -0.940,  0.,     0.,    0.,     // y'+
                       1.,     0.,     0.,    0.,     // x'2+
                               1.,    -0.769, 0.653,  // RD-
                                       1.,   -0.950,  // y'-
                                              1.      // x'2-
        // clang-format on
    };
    corStatMatrix = Utils::buildCorMatrix(nObs, data);
  } else if (c.EqualTo("LHCb_DT_Run12")) {
    corSource = "https://inspirehep.net/literature/2871248";
    std::vector<double> data = {
        // clang-format off
        //  RD+  y'+    x'2+    RD-    y'-    x'2-
        1.,     -0.759, 0.631, -0.013, 0.,    0.,     // RD+
                 1.,   -0.945,  0.,    0.,    0.,     // y'+
                        1.,     0.,    0.,    0.,     // x'2+
                                1.,   -0.756, 0.637,  // RD-
                                       1.,   -0.949,  // y'-
                                              1.      // x'2-
        // clang-format on
    };
    corStatMatrix = Utils::buildCorMatrix(nObs, data);
  } else {
    throw std::runtime_error(std::format("PDF_WS::setCorrelations ERROR config {} not found", c.Data()));
  }
}
