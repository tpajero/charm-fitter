#pragma once

#include <TString.h>

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace BLUE {

  /**
   * Representation of a measurement, identified by a name and containing a vector containing:
   * - the estimate of the observable;
   * - its statistical uncertainty;
   * - the systematic uncertainties (uncertainties that are correlated between different estimates must be listed
   *   separately).
   */
  struct Estimate {
    std::string name;
    std::vector<double> vals;
  };

  using Estimates = std::vector<Estimate>;

  /// Named combination: human-readable title plus the names of active estimates.
  using Combination = std::pair<std::string, std::vector<std::string>>;

  /// Map combining integer flags to named combinations.
  using Combinations = std::map<int, Combination>;

  struct OutputFormat {
    TString for_val = "%+4.2f";
    TString for_unc = "%2.2f";
    TString for_wei = "%2.2f";
    TString for_rho = "%2.2f";
    TString for_pul = "%2.2f";
    TString for_chi = "%2.2f";
    TString for_uni = "";
  };

  /**
   * Parse the arguments of the main function of BLUE executables.
   *
   * Looks for a single integer argument, which is interpreted as the combination flag, and for the "-h" or "--help",
   * which instead prints all available combinations, including flag, title, and the names of the measurements that
   * it includes.
   */
  int parse_args(int argc, char** argv, const Combinations& combinations);

  /**
   * Run the combination of a set of 1D estimates.
   *
   * @param flag The integer flag identifying the combination to be run.
   * @param combo_category The human-readable title of the category of combinations run by the executable, to be
   * printed in the banner.
   * @param combinations The map of all available combinations.
   * @param estimates The vector of all estimates that the combinations are based on.
   * @param names_obs The vector containing the name of the observable.
   * @param names_unc The vector of names of the uncertainties, in the order they appear in the `vals` vector of each
   * estimate.
   * @param rho The vector of correlation coefficients for the uncertainties, in the order they appear in `names_unc`.
   * @param format The output format for printing the results.
   */
  void run_combination(const int flag, const std::string& combo_category, const Combinations& combinations,
                       const Estimates& estimates, const std::vector<TString>& names_obs,
                       const std::vector<TString>& names_unc, const std::vector<double>& rho,
                       const OutputFormat& format);

}  // namespace BLUE
