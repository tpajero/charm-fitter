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

#include <BLUE/Utils.h>

#include <TString.h>

#include <vector>

const BLUE::Estimates estimates = {
    // clang-format off
    //               Val    Stat     Syst
    {"E791 1999",   { 7.32 , 28.91 , 10.32}},  // https://inspirehep.net/literature/496510
    {"FOCUS 2000",  {34.2  , 13.9  ,  7.4 }},  // https://inspirehep.net/literature/526693
    {"CLEO 2002",   {-12.0 , 25.0  , 14.0 }},  // https://inspirehep.net/literature/565839
    {"BaBar 2012",  { 7.2  ,  1.8  ,  1.2 }},  // https://inspirehep.net/literature/1186384
    {"Belle 2016",  {11.1  ,  2.2  ,  0.9 }},  // https://inspirehep.net/literature/1395100
    {"LHCb 2018",   { 5.7  ,  1.3  ,  0.9 }},  // https://inspirehep.net/literature/1698962
    {"LHCb 2022",   { 6.96 ,  0.26 ,  0.13}},  // https://inspirehep.net/literature/2035063
    // clang-format on
};
const std::vector<double> rho = {0., 0.};
const std::vector<TString> names_obs = {"yCP - yCP(Kpi)"};
const std::vector<TString> names_unc = {"Stat", "Syst"};

const BLUE::Combinations combinations = {
    {0, {"World average 2016", {"E791 1999", "FOCUS 2000", "CLEO 2002", "BaBar 2012", "Belle 2016"}}},
    {1, {"World average 2018", {"E791 1999", "FOCUS 2000", "CLEO 2002", "BaBar 2012", "Belle 2016", "LHCb 2018"}}},
    {2,
     {"World average 2022",
      {"E791 1999", "FOCUS 2000", "CLEO 2002", "BaBar 2012", "Belle 2016", "LHCb 2018", "LHCb 2022"}}},
};

/**
 * Perform the combination of a set of yCP - yCP^(K-pi+) measurements. Values are expressed in units of 10^-3.
 */
int main(int argc, char** argv) {
  const int flag = BLUE::parse_args(argc, argv, combinations);

  BLUE::run_combination(flag, "yCP - yCP(K-pi+)", combinations, estimates, names_obs, names_unc, rho,
                        BLUE::OutputFormat{});
  return EXIT_SUCCESS;
}
