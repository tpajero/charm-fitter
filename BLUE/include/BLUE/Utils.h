#pragma once

#include <TString.h>

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace BLUE {

  /// Named combination: human-readable title plus the indices of active estimates.
  using Combination = std::pair<std::string, std::vector<int>>;

  /// Map combining integer flags to named combinations.
  using Combinations = std::map<int, Combination>;

  /// Print a banner with the ID of the current executable, plus the title for the combination being run.
  void print_banner(const std::string& combiner_title, const std::string& avg_title);

  /**
   * Parse the arguments of the main function of BLUE executables.
   *
   * Looks for a single integer argument, which is interpreted as the combination flag, and for the "-h" or "--help",
   * which instead prints all available combinations, including flag, title, and the names of the measurements that
   * it includes.
   */
  int parse_args(int argc, char** argv, const Combinations& combinations, const std::vector<TString>& names);

}  // namespace BLUE
