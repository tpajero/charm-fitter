#include <BLUE/Blue.h>
#include <BLUE/Utils.h>

#include <TMatrixD.h>
#include <TString.h>

#include <algorithm>
#include <cstdlib>
#include <format>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
  /// Get the names of a vector of estimates, padded to the same width for aligned printing.
  std::vector<TString> get_names(const BLUE::Estimates& estimates) {
    size_t name_width = 0;
    for (const auto& e : estimates) name_width = std::max(name_width, e.name.size());

    std::vector<TString> names;
    for (const auto& e : estimates) { names.push_back(e.name + std::string(name_width - e.name.size(), ' ')); }
    return names;
  }

  void print_help(const BLUE::Combinations& combinations) {
    std::cout << "Available combinations:\n\n";
    for (const auto& [flag, combination] : combinations) {
      const auto& [title, inputs] = combination;
      std::cout << std::format("  {:4d}: {}\n", flag, title);
      for (const auto input : inputs) std::cout << std::format("        - {}\n", input);
      std::cout << "\n";
    }
  }

  /// Build a name -> index map for the estimates, throwing on a duplicate name.
  std::map<std::string, int> index_estimates(const BLUE::Estimates& estimates) {
    std::map<std::string, int> index;
    for (size_t i = 0; i < estimates.size(); ++i)
      if (!index.emplace(estimates[i].name, static_cast<int>(i)).second)
        throw std::runtime_error("Duplicate estimate name '" + estimates[i].name + "'");
    return index;
  }

  /// Check that every estimate name referenced by `combinations` exists in `index`.
  void validate(const BLUE::Combinations& combinations, const std::map<std::string, int>& index) {
    for (const auto& [flag, combination] : combinations) {
      const auto& [title, inputs] = combination;
      for (const auto& input : inputs)
        if (!index.count(input))
          throw std::runtime_error(
              std::format("Combination {} ('{}') references unknown estimate '{}'", flag, title, input));
    }
  }

  /// Print a banner with a title identifying the BLUE combination executable, plus the title of the run combination.
  void print_banner(const std::string& combiner_title, const std::string& avg_title) {
    const int line_length = 120;
    std::cout << std::string(line_length, '-') << '\n'
              << combiner_title << " - " << avg_title << "\n"
              << std::string(line_length, '-') << std::endl;
  }

}  // namespace

namespace BLUE {

  int parse_args(const int argc, char** argv, const Combinations& combinations) {
    if (argc != 2) {
      std::cerr << std::format("Usage: {0} <n-combination>\n"
                               "       {0} --help\n",
                               argv[0]);
      std::exit(EXIT_FAILURE);
    }
    if (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help") {
      std::cerr << std::format("\nUsage: {} <n-combination>\n\n", argv[0]);
      print_help(combinations);
      std::exit(EXIT_SUCCESS);
    }
    return std::atoi(argv[1]);
  }

  void run_combination(const int flag, const std::string& combo_category, const Combinations& combinations,
                       const Estimates& estimates, const std::vector<TString>& names_obs,
                       const std::vector<TString>& names_unc, const std::vector<double>& rho,
                       const OutputFormat& format) {

    const auto num_est = estimates.size();
    const auto num_unc = names_unc.size();

    // Catch inconsistent configuration of estimates and combinations at startup, even if the combination/estimate
    // is not selected for the current flag.
    const auto index = index_estimates(estimates);
    validate(combinations, index);

    const auto names = get_names(estimates);

    // Transform the input estimate values into the TMatrixD format expected by the BLUE class
    std::vector<double> x_est;
    for (const auto& e : estimates) {
      if (e.vals.size() != num_unc + 1)
        throw std::runtime_error("Estimate '" + std::string(e.name) + "' has " + std::to_string(e.vals.size()) +
                                 " values, expected " + std::to_string(num_unc + 1));
      x_est.insert(x_est.end(), e.vals.begin(), e.vals.end());
    }

    // Consistency checks
    if (names_obs.size() != 1) throw std::runtime_error("Only 1D combinations are currently supported");
    if (names_unc.size() != num_unc)
      throw std::runtime_error("The number of uncertainty names is inconsistent with the number of uncertainties");
    if (rho.size() != num_unc)
      throw std::runtime_error("The size of the correlation vector is inconsistent with the number of uncertainties");

    auto inp_est = std::make_unique<const TMatrixD>(num_est, num_unc + 1, &x_est[0]);

    // Statistical precision in systematic uncertainties
    const std::vector<double> s_unc(num_est * num_unc, 0.0);
    auto inp_sta = std::make_unique<const TMatrixD>(num_est, num_unc, &s_unc[0]);

    // Initialise the combiner
    auto my_blue = std::make_unique<Blue>(num_est, num_unc);
    my_blue->SetFormat(format.for_val, format.for_unc, format.for_wei, format.for_rho, format.for_pul, format.for_chi,
                       format.for_uni);
    my_blue->FillNamEst(&names[0]);
    my_blue->FillNamUnc(&names_unc[0]);
    my_blue->FillNamObs(&names_obs[0]);
    my_blue->FillEst(inp_est.get());
    my_blue->FillSta(inp_sta.get());
    for (int i = 0; i < num_unc; ++i) my_blue->FillCor(i, rho[i]);

    // Perform the combination
    const auto [title, inputs] = combinations.at(flag);
    print_banner(combo_category, title);
    for (auto i = 0; i < num_est; ++i) my_blue->SetInActiveEst(i);
    for (const auto& input : inputs) my_blue->SetActiveEst(index.at(input));
    my_blue->FixInp();
    my_blue->PrintEst();
    my_blue->Solve();
    my_blue->PrintResult();
  }

}  // namespace BLUE
