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

#include <BLUE/Utils.h>

#include <TString.h>

#include <vector>

const BLUE::Estimates estimates = {
    // clang-format off
    //                Val    Stat    Syst
    {"Belle 2009",  {  1.1 ,  6.1 ,  5.2}},  // https://inspirehep.net/literature/821323   D0 -> KS0 K+ K-
    {"BESIII 2015", {-20.0 , 13.0 ,  7.0}},  // https://inspirehep.net/literature/1337473
    // clang-format on
};
const std::vector<double> rho = {0., 0.};
const std::vector<TString> names_obs = {"    yCP"};
const std::vector<TString> names_unc = {"  Stat", "  Syst"};

const BLUE::Combinations combinations = {
    {0, {"World average 2015", {"Belle 2009", "BESIII 2015"}}},
};

/**
 * Perform the combination of a set of yCP (CP-even) measurements.
 */
int main(int argc, char** argv) {
  const int flag = BLUE::parse_args(argc, argv, combinations);

  BLUE::run_combination(flag, "yCP (CP-even)", combinations, estimates, names_obs, names_unc, rho,
                        BLUE::OutputFormat{});
  return EXIT_SUCCESS;
}
