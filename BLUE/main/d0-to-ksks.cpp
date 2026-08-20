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

#include <BLUE/Utils.h>

#include <TString.h>

#include <vector>

const BLUE::Estimates estimates = {
    // clang-format off
    //                 Val       Stat     Sys1    Sys2
    {"CLEO",          {-0.23  , 0.19  , 0.    , 0.     }},  // https://inspirehep.net/literature/539090
    {"LHCb 2015",     {-0.029 , 0.052 , 0.022 , 0.     }},  // https://inspirehep.net/literature/1389705
    {"Belle 2017",    {-0.0002, 0.0153, 0.0002, 0.0017 }},  // https://inspirehep.net/literature/1599959
    {"LHCb 2021",     {-0.031 , 0.012 , 0.004 , 0.00056}},  // https://inspirehep.net/literature/1861934  updated with https://arxiv.org/abs/2209.03179
    {"CMS",           { 0.063 , 0.030 , 0.002 , 0.     }},  // https://inspirehep.net/literature/2788405
    {"Belle 2024",    {-0.014 , 0.013 , 0.001 , 0.     }},  // https://inspirehep.net/literature/2914970
    {"Belle 2025",    {-0.006 , 0.011 , 0.001 , 0.     }},  // https://inspirehep.net/literature/2844984
    {"LHCb 2025",     { 0.0186, 0.0104, 0.0041, 0.     }},  // https://inspirehep.net/literature/3070434
    // clang-format on
};
const std::vector<double> rho = {0., 0., 0.};
const std::vector<TString> names_obs = {" ACP(D0 -> KS KS)"};
const std::vector<TString> names_unc = {"Stat", "Sys1", "Sys2"};

const BLUE::Combinations combinations = {
    {0, {"LHCb average 2025", {"LHCb 2015", "LHCb 2021", "LHCb 2025"}}},
    {1, {"World average 2025", {"CLEO", "LHCb 2015", "LHCb 2021", "CMS", "Belle 2025", "LHCb 2025"}}},
};

/**
 * Perform the combination of a set of ACP(D0 -> KS KS) measurements.
 */
int main(int argc, char** argv) {
  const int flag = BLUE::parse_args(argc, argv, combinations);

  BLUE::run_combination(flag, "ACP(D0 -> KS KS)", combinations, estimates, names_obs, names_unc, rho,
                        BLUE::OutputFormat{"%+1.4f", "%2.4f", "%2.4f", "%2.4f", "%2.4f", "%2.4f", ""});
  return EXIT_SUCCESS;
}
