/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: October 2022
 *
 * Perform the average of DeltaY (a.k.a. -A_Gamma) measurements. Values are expressed in units of 10^-4.
 *
 * Print usage and available options with
 *
 *    ./bin/BLUE/dy -h
 *
 */

#include <BLUE/Utils.h>

#include <TString.h>

#include <vector>

// Inverse of Fp_pipipi0 = 0.942554 +/- 0.00404575
constexpr auto corr_3pi = 1.061;

/**
 * Vector of all estimates, listed in chronological order.
 *
 * The scale factor to account for dilution in the mu-tagged measurement by LHCb (Run2_mu) is taken from
 * LHCb-ANA-2019-021, and gets contributions from:
 *  - Fig. 11, where the average from the deviation from unity of the slope is -0.0554
 *  - Fig. 12, where the average from the deviation from unity of the slope is -0.0521
 * The scale factor is thus equal to 1 / (1 - 0.0554) / (1 - 0.0521) = 1.117.
 */
const BLUE::Estimates estimates = {
    // clang-format off
    //                           Val               Stat              m(hhh) Sec   m(KK) m(PP) Weight TimeRes Mistag Run1Mu Other
    {"BaBar",                   {-8.8            , 25.5            , 0.  ,  0.  , 0.  , 0.  , 0.  ,  0.  ,   0.  ,  0.  ,  5.8            }},  // https://inspirehep.net/literature/1186384  2012
    {"CDF KK",                  {19.             , 15.             , 0.  ,  0.  , 0.  , 0.  , 0.  ,  0.  ,   0.  ,  0.  ,  0.             }},  // https://inspirehep.net/literature/1323066  2014
    {"CDF PP",                  { 1.             , 18.             , 0.  ,  0.  , 0.  , 0.  , 0.  ,  0.  ,   0.  ,  0.  ,  0.             }},  // https://inspirehep.net/literature/1323066  2014
    {"LHCb Run 1 mu KK",        {13.40           ,  7.70           , 0.  ,  0.  , 0.  , 0.  , 0.  ,  0.49,   0.  ,  2.55,  0.             }},  // https://inspirehep.net/literature/1341286  2015-01
    {"LHCb Run 1 mu PP",        { 9.2            , 14.5            , 0.  ,  0.  , 0.  , 0.  , 0.  ,  0.42,   0.  ,  2.48,  0.             }},  // https://inspirehep.net/literature/1341286  2015-01
    {"Belle",                   { 3.             , 20.             , 0.  ,  0.  , 0.  , 0.  , 0.  ,  0.  ,   0.  ,  0.  ,  7.             }},  // https://inspirehep.net/literature/1395100  2015-09
    {"LHCb Run 1 prompt KK",    { 3.0            ,  3.2            , 0.1 ,  0.8 , 0.5 , 0.  , 0.2 ,  0.  ,   0.  ,  0.  ,  0.             }},  // https://inspirehep.net/literature/1514549  2017
    {"LHCb Run 1 prompt PP",    {-4.6            ,  5.8            , 0.1 ,  1.2 , 0.  , 0.  , 0.2 ,  0.  ,   0.  ,  0.  ,  0.             }},  // https://inspirehep.net/literature/1514549  2017
    {"LHCb Run 2 mu KK",        { 4.8            ,  4.0            , 0.  ,  0.  , 0.  , 0.  , 0.  ,  0.  ,   0.  ,  0.  ,  0.3            }},  // https://inspirehep.net/literature/1762838  2020 (scale factor of 1.12 to account for dilution)
    {"LHCb Run 2 mu PP",        {-2.5            ,  7.8            , 0.  ,  0.  , 0.  , 0.  , 0.  ,  0.  ,   0.  ,  0.  ,  0.3            }},  // https://inspirehep.net/literature/1762838  2020 (scale factor of 1.12 to account for dilution)
    {"LHCb Run 2 prompt KK",    {-2.321          ,  1.524          , 0.24,  0.13, 0.06, 0.  , 0.05,  0.  ,   0.  ,  0.  ,  0.             }},  // https://inspirehep.net/literature/1864385  2021
    {"LHCb Run 2 prompt PP",    {-4.014          ,  2.814          , 0.34,  0.13, 0.  , 0.03, 0.05,  0.  ,   0.  ,  0.  ,  0.             }},  // https://inspirehep.net/literature/1864385  2021
    {"LHCb Run 2 pi+ pi- pi0",  {-1.21 * corr_3pi,  5.97 * corr_3pi, 0.  ,  0.  , 0.  , 0.  , 0.  ,  0.  ,   0.  ,  0.  ,  2.00 * corr_3pi}},  // https://inspirehep.net/literature/2785424  2024 (remove sys. unc. from time binning)
    // clang-format on
};
const std::vector<double> rho = {0., 1., 1., 1., 1., 1., 1., 1., 1., 0.};
const std::vector<TString> names_obs = {"DY"};
const std::vector<TString> names_unc = {"Stat",   "m(hhh)", "Sec",    "m(KK)",  "m(PP)",
                                        "Weight", "TimRes", "Mistag", "Run1Mu", " Other"};

const BLUE::Combinations combinations = {
    // LHCb only
    {0, {"LHCb Run 1", {"LHCb Run 1 mu KK", "LHCb Run 1 mu PP", "LHCb Run 1 prompt KK", "LHCb Run 1 prompt PP"}}},
    {1, {"LHCb Run 1 (K+ K-)", {"LHCb Run 1 mu KK", "LHCb Run 1 prompt KK"}}},
    {2, {"LHCb Run 1 (pi+ pi-)", {"LHCb Run 1 mu PP", "LHCb Run 1 prompt PP"}}},
    {3,
     {"LHCb Run 1+2 (h+ h-)",
      {"LHCb Run 1 mu KK", "LHCb Run 1 mu PP", "LHCb Run 1 prompt KK", "LHCb Run 1 prompt PP", "LHCb Run 2 mu KK",
       "LHCb Run 2 mu PP", "LHCb Run 2 prompt KK", "LHCb Run 2 prompt PP"}}},
    {4,
     {"LHCb Run 1+2 (K+ K-)",
      {"LHCb Run 1 mu KK", "LHCb Run 1 prompt KK", "LHCb Run 2 mu KK", "LHCb Run 2 prompt KK"}}},
    {5,
     {"LHCb Run 1+2 (pi+ pi-)",
      {"LHCb Run 1 mu PP", "LHCb Run 1 prompt PP", "LHCb Run 2 mu PP", "LHCb Run 2 prompt PP"}}},
    {6,
     {"LHCb Run 1+2",
      {"LHCb Run 1 mu KK", "LHCb Run 1 mu PP", "LHCb Run 1 prompt KK", "LHCb Run 1 prompt PP", "LHCb Run 2 mu KK",
       "LHCb Run 2 mu PP", "LHCb Run 2 prompt KK", "LHCb Run 2 prompt PP", "LHCb Run 2 pi+ pi- pi0"}}},
    // World averages
    {100,
     {"World average 2019",
      {"BaBar", "CDF KK", "CDF PP", "LHCb Run 1 mu KK", "LHCb Run 1 mu PP", "Belle", "LHCb Run 1 prompt KK",
       "LHCb Run 1 prompt PP"}}},
    {101, {"World average 2019 (K+ K-)", {"CDF KK", "LHCb Run 1 mu KK", "LHCb Run 1 prompt KK"}}},
    {102, {"World average 2019 (pi+ pi-)", {"CDF PP", "LHCb Run 1 mu PP", "LHCb Run 1 prompt PP"}}},
    {103,
     {"World average 2020",
      {"BaBar", "CDF KK", "CDF PP", "LHCb Run 1 mu KK", "LHCb Run 1 mu PP", "Belle", "LHCb Run 1 prompt KK",
       "LHCb Run 1 prompt PP", "LHCb Run 2 mu KK", "LHCb Run 2 mu PP"}}},
    {104, {"World average 2020 (K+ K-)", {"CDF KK", "LHCb Run 1 mu KK", "LHCb Run 1 prompt KK", "LHCb Run 2 mu KK"}}},
    {105, {"World average 2020 (pi+ pi-)", {"CDF PP", "LHCb Run 1 mu PP", "LHCb Run 1 prompt PP", "LHCb Run 2 mu PP"}}},
    {106,
     {"World average 2021",
      {"BaBar", "CDF KK", "CDF PP", "LHCb Run 1 mu KK", "LHCb Run 1 mu PP", "Belle", "LHCb Run 1 prompt KK",
       "LHCb Run 1 prompt PP", "LHCb Run 2 mu KK", "LHCb Run 2 mu PP", "LHCb Run 2 prompt KK",
       "LHCb Run 2 prompt PP"}}},
    {107,
     {"World average 2021 (K+ K-)",
      {"CDF KK", "LHCb Run 1 mu KK", "LHCb Run 1 prompt KK", "LHCb Run 2 mu KK", "LHCb Run 2 prompt KK"}}},
    {108,
     {"World average 2021 (pi+ pi-)",
      {"CDF PP", "LHCb Run 1 mu PP", "LHCb Run 1 prompt PP", "LHCb Run 2 mu PP", "LHCb Run 2 prompt PP"}}},
    {109,
     {"World average 2024",
      {"BaBar", "CDF KK", "CDF PP", "LHCb Run 1 mu KK", "LHCb Run 1 mu PP", "Belle", "LHCb Run 1 prompt KK",
       "LHCb Run 1 prompt PP", "LHCb Run 2 mu KK", "LHCb Run 2 mu PP", "LHCb Run 2 prompt KK", "LHCb Run 2 prompt PP",
       "LHCb Run 2 pi+ pi- pi0"}}},
    // B-factories (the only D0 -> h+ h- measurements not separated by K+ K- / pi+ pi-)
    {200, {"BaBar + Belle", {"BaBar", "Belle"}}},
    // Non-LHCb averages
    {300, {"BaBar + CDF + Belle", {"BaBar", "CDF KK", "CDF PP", "Belle"}}},
};

/**
 * Perform the combination of a set of DeltaY (a.k.a. -A_Gamma) measurements. Values are expressed in units of 10^-4.
 */
int main(int argc, char** argv) {
  const int flag = BLUE::parse_args(argc, argv, combinations);

  BLUE::run_combination(flag, "DeltaY(h- h+ (pi0))", combinations, estimates, names_obs, names_unc, rho,
                        BLUE::OutputFormat{"%+6.2f", "%5.2f", "%2.2f", "%2.2f", "%2.2f", "%2.2f", ""});
  return EXIT_SUCCESS;
}
