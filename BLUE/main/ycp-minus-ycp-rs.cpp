/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: August 2026
 *
 * Perform the average of yCP - yCP^(K-pi+) measurements. Values are expressed in units of 10^-3.
 *
 * Print usage and available options with
 *
 *    ./bin/BLUE/ycp-minus-ycp-rs -h
 *
 **/

#include <BLUE/Blue.h>
#include <BLUE/Utils.h>

#include <TMatrixD.h>
#include <TString.h>

#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <vector>

const BLUE::Combinations combinations = {
    {0, {"World average 2016", {0, 1, 2, 3, 4}}},
    {1, {"World average 2018", {0, 1, 2, 3, 4, 5}}},
    {2, {"World average 2022", {0, 1, 2, 3, 4, 5, 6}}},
};

// number of estimates, uncertainties, observables
const std::vector<TString> names = {
    "E791  1999 ",  //  0
    "FOCUS 2000 ",  //  1
    "CLEO  2002 ",  //  2
    "BaBar 2012 ",  //  3
    "Belle 2016 ",  //  4
    "LHCb  2018 ",  //  5
    "LHCb  2022 ",  //  6
};

/**
 * Perform the combination of a given set of measurements.
 *
 * @param flag Steers which set of measurements should be employed (according to the values of "combinations" above).
 */
void run(int flag) {

  const std::vector<TString> names_obs = {"yCP - yCP(Kpi)"};
  const auto num_est = names.size();

  // Array of estimates and uncertainties
  const std::vector<TString> names_unc = {"  Stat", "  Syst"};

  const std::vector<double> x_est = {
      // clang-format off
      //      0       1
      // Val  Stat     Syst
        7.32 , 28.91 , 10.32,  // E791  1999  https://inspirehep.net/literature/496510
       34.2  , 13.9  ,  7.4 ,  // FOCUS 2000  https://inspirehep.net/literature/526693
      -12.0  , 25.0  , 14.0 ,  // CLEO  2002  https://inspirehep.net/literature/565839
        7.2  ,  1.8  ,  1.2 ,  // BaBar 2012  https://inspirehep.net/literature/1186384
       11.1  ,  2.2  ,  0.9 ,  // Belle 2016  https://inspirehep.net/literature/1395100
        5.7  ,  1.3  ,  0.9 ,  // LHCb  2018  https://inspirehep.net/literature/1698962
        6.96 ,  0.26 ,  0.13,  // LHCb  2022  https://inspirehep.net/literature/2035063
      // clang-format on
  };
  const auto num_unc = names_unc.size();
  if (x_est.size() != num_est * (num_unc + 1))
    throw std::runtime_error("The size of the uncertainty matrix is inconsistent with the number of estimates and the "
                             "number of uncertainties per estimate");
  auto inp_est = std::make_unique<const TMatrixD>(num_est, num_unc + 1, &x_est[0]);

  // Statistical and systematic uncertainties are not correlated among the different experiments
  const std::vector<double> rho_val = {0., 0.};
  if (rho_val.size() != num_unc)
    throw std::runtime_error("The size of the correlation vector is inconsistent with the number of uncertainties");

  // Statistical precision in systematic uncertainties
  const std::vector<double> s_unc(num_est * num_unc, 0.);
  auto inp_sta = std::make_unique<const TMatrixD>(num_est, num_unc, &s_unc[0]);

  // Format for the output
  const TString for_val = "%+6.2f";
  const TString for_unc = "%5.2f";
  const TString for_wei = "%2.2f";
  const TString for_rho = "%2.2f";
  const TString for_pul = "%2.2f";
  const TString for_chi = "%2.2f";
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
  BLUE::print_banner("yCP - yCP(K-pi+)", title);
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
