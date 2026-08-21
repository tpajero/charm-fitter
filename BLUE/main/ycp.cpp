/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: August 2026
 *
 * Perform the average of yCP measurements, adopting the sign convention of CP-even final states (like K+ K-, pi+ pi-,
 * KS K+ K-). Different combinations allow for grouping the measurements based on the control-channel correction.
 *
 * Print usage and available options with
 *
 *    ./bin/BLUE/ycp -h
 *
 */

#include <BLUE/Utils.h>

#include <TString.h>

#include <vector>

// Chronologically ordered, according to the arXiv submission date.
// All values are expressed in units of 10^-3.
const BLUE::Estimates estimates = {
    // clang-format off
    //                Val    Stat    Syst
    {"E791 1999",   { 7.32 , 28.91 , 10.32}},  // https://inspirehep.net/literature/496510   yCP - yCP(KP)    |
    {"FOCUS 2000",  {34.2  , 13.9  ,  7.4 }},  // https://inspirehep.net/literature/526693   yCP - yCP(KP/RS) | KP likely more abundant
    {"CLEO 2001",   {-12.0 , 25.0  , 14.0 }},  // https://inspirehep.net/literature/565839   yCP - yCP(RS)    |
    {"Belle 2009",  {  1.1 ,  6.1  ,  5.2 }},  // https://inspirehep.net/literature/821323   yCP              | D0 -> KS0 K+ K-
    {"LHCb 2011",   {  5.5 ,  6.3  ,  4.1 }},  // https://inspirehep.net/literature/1082327  yCP - yCP(RS)    |
    {"BaBar 2012",  {  7.2 ,  1.8  ,  1.2 }},  // https://inspirehep.net/literature/1186384  yCP - yCP(KP)    |
    {"BESIII 2015", {-20.0 , 13.0  ,  7.0 }},  // https://inspirehep.net/literature/1337473  yCP              |
    {"Belle 2015",  { 11.1 ,  2.2  ,  0.9 }},  // https://inspirehep.net/literature/1395100  yCP - yCP(RS)    |
    {"LHCb 2018",   {  5.7 ,  1.3  ,  0.9 }},  // https://inspirehep.net/literature/1698962  yCP - yCP(RS)    |
    {"Belle 2019",  {  9.6 ,  9.1  ,  6.43}},  // https://inspirehep.net/literature/1772245  yCP + yCP(RS)    | D0 -> KS omega
    {"LHCb 2022",   {  6.96,  0.26 ,  0.13}},  // https://inspirehep.net/literature/2035063  yCP - yCP(RS)    |
    // clang-format on
};
const std::vector<double> rho = {0., 0.};
const std::vector<TString> names_obs = {"yCP"};
const std::vector<TString> names_unc = {"Stat", "Syst"};

const BLUE::Combinations combinations = {
    // yCP - yCP(RS) ---------------------------------------------------------------------------------------------------
    {0, {"yCP - yCP(RS) - World average 2015", {"CLEO 2001", "LHCb 2011", "Belle 2015"}}},
    {1, {"yCP - yCP(RS) - World average 2018", {"CLEO 2001", "LHCb 2011", "Belle 2015", "LHCb 2018"}}},
    {2, {"yCP - yCP(RS) - World average 2022", {"CLEO 2001", "LHCb 2011", "Belle 2015", "LHCb 2018", "LHCb 2022"}}},
    // yCP + yCP(RS) ---------------------------------------------------------------------------------------------------
    // yCP - yCP(KP) ---------------------------------------------------------------------------------------------------
    {200, {"yCP - yCP(KP) - World average 2015", {"E791 1999", "FOCUS 2000", "BaBar 2012"}}},
    // yCP -------------------------------------------------------------------------------------------------------------
    {300, {"yCP - World average 2015", {"Belle 2009", "BESIII 2015"}}},
    // All yCP-like measurements ---------------------------------------------------------------------------------------
    {400,
     {"All - World average Jan 2015",
      {"E791 1999", "FOCUS 2000", "CLEO 2001", "Belle 2009", "LHCb 2011", "BaBar 2012", "BESIII 2015"}}},
    {401,
     {"All - World average Sep 2015",
      {"E791 1999", "FOCUS 2000", "CLEO 2001", "Belle 2009", "LHCb 2011", "BaBar 2012", "BESIII 2015", "Belle 2015"}}},
    {402,
     {"All - World average 2018",
      {"E791 1999", "FOCUS 2000", "CLEO 2001", "Belle 2009", "LHCb 2011", "BaBar 2012", "BESIII 2015", "Belle 2015",
       "LHCb 2018"}}},
    {403,
     {"All - World average 2019",
      {"E791 1999", "FOCUS 2000", "CLEO 2001", "Belle 2009", "LHCb 2011", "BaBar 2012", "BESIII 2015", "Belle 2015",
       "LHCb 2018", "Belle 2019"}}},
    {404,
     {"All - World average 2022",
      {"E791 1999", "FOCUS 2000", "CLEO 2001", "Belle 2009", "LHCb 2011", "BaBar 2012", "BESIII 2015", "Belle 2015",
       "LHCb 2018", "Belle 2019", "LHCb 2022"}}},
    // yCP - yCP(RS) without LHCb --------------------------------------------------------------------------------------
    {500, {"yCP - yCP(RS) - World average no LHCb 2015", {"CLEO 2001", "Belle 2015"}}},
    // All yCP-like measurements without LHCb --------------------------------------------------------------------------
    {600,
     {"All - World average no LHCb 2019",
      {"E791 1999", "FOCUS 2000", "CLEO 2001", "Belle 2009", "BaBar 2012", "BESIII 2015", "Belle 2015", "Belle 2019"}}},
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
