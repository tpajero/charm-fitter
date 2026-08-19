/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: August 2026
 *
 * Perform the average of yCP measurements, based on CP-eveen final states like K+ K-, pi+ pi-, KS K+ K-.
 *
 * Print usage and available options with
 *
 *    ./bin/BLUE/ycp -h
 *
 */

#include <BLUE/Blue.h>
#include <BLUE/Utils.h>

#include <TMatrixD.h>
#include <TString.h>

#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <vector>

const BLUE::Combinations combinations = {
    {0, {"World average 2015", {0, 1}}},
};

const std::vector<TString> names = {
    "Belle 2009 ",  //  0  https://inspirehep.net/literature/821323
    "BESIII 2015",  //  1  https://inspirehep.net/literature/1337473
};

/**
 * Perform the combination of a given set of measurements.
 *
 * @param flag Steers which set of measurements should be employed (according to the values of "combinations" above).
 */
void run(const int flag) {

  const std::vector<TString> names_obs = {"    yCP"};
  const auto num_est = names.size();

  // Array of estimates and uncertainties
  const std::vector<TString> names_unc = {"  Stat", "  Syst"};

  const std::vector<double> x_est = {
      // clang-format off
      // Val  Stat    Syst
        1.1 ,  6.1 ,  5.2,  // Belle 2009, D0 -> KS0 K+ K-
      -20.0 , 13.0 ,  7.0,  // BESIII 2015
      // clang-format on
  };
  const auto num_unc = names_unc.size();
  if (x_est.size() != num_est * (num_unc + 1))
    throw std::runtime_error("The size of the uncertainty matrix is inconsistent with the number of estimates and the "
                             "number of uncertainties per estimate");
  auto inp_est = std::make_unique<const TMatrixD>(num_est, num_unc + 1, &x_est[0]);

  // Statistical and systematic uncertainties are not correlated among the two measurements
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
  BLUE::print_banner("yCP (CP-even)", title);
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
