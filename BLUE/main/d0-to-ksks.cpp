/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: September 2024
 *
 * Perform the average of ACP(D0 -> KS KS) measurements.
 *
 * Print usage and available options with
 *
 *    ./bin/BLUE/ksks -h
 *
 **/

#include <BLUE/Blue.h>
#include <BLUE/Utils.h>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

const BLUE::Combinations combinations = {
    {0, {"LHCb average 2025", {1, 3, 7}}},
    {1, {"World average 2025", {0, 1, 3, 4, 6, 7}}},
};

const std::vector<TString> names = {
    "CLEO      ",  //  0  https://inspirehep.net/literature/539090
    "LHCb 2015 ",  //  1  https://inspirehep.net/literature/1389705
    "Belle 2017",  //  2  https://inspirehep.net/literature/1599959
    "LHCb 2021 ",  //  3  https://inspirehep.net/literature/1861934
    "CMS       ",  //  4  https://inspirehep.net/literature/2788405
    "Belle 2024",  //  5  https://inspirehep.net/literature/2914970
    "Belle 2025",  //  6  https://inspirehep.net/literature/2844984
    "LHCb 2025 ",  //  7  https://inspirehep.net/literature/3070434
};

/**
 * Perform the combination of a given set of measurements.
 *
 * @param flag Steers which set of measurements should be employed (according to the values of "combinations" above).
 */
void run(const int flag) {

  const std::vector<TString> names_obs = {" ACP(D0 -> KS KS)"};
  const auto num_est = names.size();

  // Array of estimates and uncertainties
  const std::vector<TString> names_unc = {"     Stat", "    Sys1", "    Sys2"};
  const std::vector<double> x_est = {
      // clang-format off
      //       0       1      2
      // Val   Stat    Sys1   Sys2
      -0.23  , 0.19  , 0.    , 0.     ,  // CLEO
      -0.029 , 0.052 , 0.022 , 0.     ,  // LHCb 2015
      -0.0002, 0.0153, 0.0002, 0.0017 ,  // Belle 2017
      -0.031 , 0.012 , 0.004 , 0.00056,  // LHCb 2021 updated with https://arxiv.org/abs/2209.03179
       0.063 , 0.030 , 0.002 , 0.     ,  // CMS
      -0.014 , 0.013 , 0.001 , 0.     ,  // Belle 2024
      -0.006 , 0.011 , 0.001 , 0.     ,  // Belle 2025
       0.0186, 0.0104, 0.0041, 0.     ,  // LHCb 2025
      // clang-format on
  };
  const auto num_unc = names_unc.size();
  if (x_est.size() != num_est * (num_unc + 1)) {
    std::cout << "The size of the uncertainty matrix is inconsistent with the number of estimates and the number of "
                 "uncertainties per estimate"
              << std::endl;
    exit(EXIT_FAILURE);
  }
  auto inp_est = std::make_unique<const TMatrixD>(num_est, num_unc + 1, &x_est[0]);

  // Statistical uncertainties are not correlated
  // All systematic uncertainties are uncorrelated
  const std::vector<double> rho_val = {0., 0., 0.};
  if (rho_val.size() == num_unc + 1) {
    std::cout << "The size of correlation vector among uncertainties is inconsistent with the number of uncertainties"
              << std::endl;
    exit(EXIT_FAILURE);
  }

  // Statistical precision in systematic uncertainties
  const std::vector<double> s_unc(num_est * num_unc, 0.);
  auto inp_sta = std::make_unique<const TMatrixD>(num_est, num_unc, &s_unc[0]);

  // Format for the output
  const TString for_val = "%+1.4f";
  const TString for_unc = "%2.4f";
  const TString for_wei = "%2.4f";
  const TString for_rho = "%2.4f";
  const TString for_pul = "%2.4f";
  const TString for_chi = "%2.4f";
  const TString for_uni = "";

  // Initialise the combiner
  auto my_blue = std::make_unique<Blue>(num_est, num_unc);
  my_blue->SetFormat(for_val, for_unc, for_wei, for_rho, for_pul, for_chi, for_uni);
  my_blue->FillNamEst(&names[0]);
  my_blue->FillNamUnc(&names_unc[0]);
  my_blue->FillNamObs(&names_obs[0]);
  my_blue->FillEst(inp_est.get());
  my_blue->FillSta(inp_sta.get());
  for (int i = 0; i < num_unc; ++i) my_blue->FillCor(i, rho_val[i]);

  // Perform the combination
  const auto [title, inputs] = combinations.at(flag);
  BLUE::print_banner("ACP(D0 -> KS KS)", title);
  for (auto i = 0; i < num_est; ++i) my_blue->SetInActiveEst(i);
  for (auto i : inputs) my_blue->SetActiveEst(i);
  my_blue->FixInp();
  my_blue->PrintEst();
  my_blue->Solve();
  my_blue->PrintResult();
}

int main(int argc, char** argv) {
  const int flag = BLUE::parse_args(argc, argv, combinations, names);
  run(flag);
  return EXIT_SUCCESS;
}
