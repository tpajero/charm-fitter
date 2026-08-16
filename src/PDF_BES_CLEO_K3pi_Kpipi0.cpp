/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: June 2024
 *
 * Note that this part of the combination is slightly sloppy, as the CLEO + BESIII combination:
 *   1. employs external inputs for the mixing parameters (the effect on this combination is probably negligible,
 *      as the relative precision coming from K3pi and Kpipi0 is much poorer than than of the other decay channels.
 *   2. the PDF is not Gaussian, and this effect cannot be taken into account as the full likelihood is not public.
 **/

#include <PDF_BES_CLEO_K3pi_Kpipi0.h>

#include <CharmUtils.h>

#include <Utils.h>

#include <RooArgList.h>
#include <RooFormulaVar.h>
#include <RooRealVar.h>

#include <algorithm>
#include <map>
#include <string>
#include <vector>

using Utils::DegToRad;

PDF_BES_CLEO_K3pi_Kpipi0::PDF_BES_CLEO_K3pi_Kpipi0(const TString measurement_id) : PDF_Charm{6} {
  name = "K3pi_" + measurement_id;
  initialise(measurement_id, measurement_id, measurement_id);
}

std::set<std::string> PDF_BES_CLEO_K3pi_Kpipi0::getParameterNames() const {
  return {"r_K3pi", "k_K3pi", "Delta_K3pi", "r_Kpipi0", "k_Kpipi0", "Delta_Kpipi0"};
}

void PDF_BES_CLEO_K3pi_Kpipi0::initRelations() {
  theory = new RooArgList("theory");  ///< the order of this list must match that of the COR matrix!
  using Utils::makeTheoryVar;
  theory->add(*(makeTheoryVar("k_K3pi_th", "k_K3pi", parameters)));
  theory->add(*(makeTheoryVar("Delta_K3pi_th", "Delta_K3pi", parameters)));
  theory->add(*(makeTheoryVar("k_Kpipi0_th", "k_Kpipi0", parameters)));
  theory->add(*(makeTheoryVar("Delta_Kpipi0_th", "Delta_Kpipi0", parameters)));
  theory->add(*(makeTheoryVar("r_K3pi_th", "r_K3pi", parameters)));
  theory->add(*(makeTheoryVar("r_Kpipi0_th", "r_Kpipi0", parameters)));
}

void PDF_BES_CLEO_K3pi_Kpipi0::initObservables() {
  const TString label = "BES3 + CLEO";
  observables = new RooArgList("observables");  ///< the order of this list must match that of the COR matrix!
  observables->add(*(new RooRealVar("k_K3pi_obs", label + "   #it{#kappa_{K3#pi}}", 1, -2, 2)));
  observables->add(
      *(new RooRealVar("Delta_K3pi_obs", label + "   #it{#delta_{K3#pi}}", 0, DegToRad(-180), DegToRad(180))));
  observables->add(*(new RooRealVar("k_Kpipi0_obs", label + "   #it{#kappa_{K#pi#pi^{0}}}", 1, -2, 2)));
  observables->add(
      *(new RooRealVar("Delta_Kpipi0_obs", label + "   #it{#delta_{K#pi#pi^{0}}}", 0, DegToRad(-180), DegToRad(180))));
  observables->add(*(new RooRealVar("r_K3pi_obs", label + "   #it{r_{K3#pi}}", 1, -1e4, 1e4)));
  observables->add(*(new RooRealVar("r_Kpipi0_obs", label + "   #it{r_{K#pi#pi^{0}}}", 1, -1e4, 1e4)));
}

void PDF_BES_CLEO_K3pi_Kpipi0::setObservables(const TString c) {
  if (c.EqualTo("truth"))
    setObservablesTruth();
  else if (c.EqualTo("toy"))
    setObservablesToy();
  else if (c.EqualTo("BES3-CLEO")) {
    obsValSource = "BES+CLEO, arXiv:2103.05988";
    setObservable("k_K3pi_obs", 0.49);
    setObservable("Delta_K3pi_obs", DegToRad(26));
    setObservable("k_Kpipi0_obs", 0.79);
    setObservable("Delta_Kpipi0_obs", DegToRad(-16));
    setObservable("r_K3pi_obs", 5.46e-2);
    setObservable("r_Kpipi0_obs", 4.41e-2);
  } else {
    std::cout << "PDF_BES_CLEO_K3pi_Kpipi0::setObservables() : ERROR : config " + c + " not found." << std::endl;
    exit(1);
  }
}

void PDF_BES_CLEO_K3pi_Kpipi0::setUncertainties(const TString c) {
  if (c.EqualTo("BES3-CLEO")) {
    obsErrSource = "BES+CLEO, arXiv:2103.05988";
    // Values are the average of the upper and lower asymmetric uncertainties
    StatErr[0] = 0.105;         // k_K3pi
    StatErr[1] = DegToRad(18);  // Delta_K3pi
    StatErr[2] = 0.04;          // k_Kpipi0
    StatErr[3] = DegToRad(11);  // Delta_Kpipi0
    StatErr[4] = 0.08e-2;       // r_K3pi
    StatErr[5] = 0.11e-2;       // r_Kpipi0

    std::ranges::fill(SystErr, 0.0);
  } else {
    std::cout << "PDF_BES_CLEO_K3pi_Kpipi0::setUncertainties() : ERROR : config " + c + " not found." << std::endl;
    exit(1);
  }
}

void PDF_BES_CLEO_K3pi_Kpipi0::setCorrelations(const TString c) {
  resetCorrelations();
  if (c.EqualTo("BES3-CLEO")) {
    corSource = "BES+CLEO, arXiv:2103.05988";
    std::vector<double> dataStat = {
        // clang-format off
        // k_K3pi Delta_K3pi k_Kpipi0 Delta_Kpipi0 r_K3pi r_Kpipi0
        1.,       0.78,       0.04,    0.07,       0.50,  -0.06,  // k_K3pi
                  1.,        -0.15,   -0.15,       0.34,   0.04,  // Delta_K3pi
                              1.,     -0.23,       0.02,   0.05,  // k_Kpipi0
                                       1.,         0.11,   0.11,  // Delta_Kpipi0
                                                   1.,    -0.02,  // r_K3pi
                                                           1.,    // r_Kpipi0
        // clang-format on
    };
    corStatMatrix = Utils::buildCorMatrix(nObs, dataStat);
  } else {
    std::cout << "PDF_BES_CLEO_K3pi_Kpipi0::setCorrelations() : ERROR : config " + c + " not found." << std::endl;
    exit(1);
  }
}
