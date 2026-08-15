#include <BLUE/Utils.h>

#include <TString.h>

#include <cstdlib>
#include <format>
#include <iostream>
#include <vector>

namespace BLUE {

  void print_banner(const std::string& combiner_title, const std::string& avg_title) {
    const int line_length = 120;
    std::cout << std::string(line_length, '-') << '\n'
              << combiner_title << " - " << avg_title << "\n"
              << std::string(line_length, '-') << std::endl;
  }

  namespace {

    void print_help(const Combinations& combinations, const std::vector<TString>& names) {
      std::cout << "Available combinations:\n\n";
      for (const auto& [flag, combination] : combinations) {
        const auto& [title, inputs] = combination;
        std::cout << std::format("  {:4d}: {}\n", flag, title);
        for (const auto i : inputs) std::cout << std::format("        - {}\n", names.at(i).Data());
        std::cout << "\n";
      }
    }

  }  // namespace

  int parse_args(const int argc, char** argv, const Combinations& combinations, const std::vector<TString>& names) {
    if (argc != 2) {
      std::cerr << std::format(
          "Usage: {0} <n-combination>\n"
          "       {0} --help\n",
          argv[0]);
      std::exit(EXIT_FAILURE);
    }
    if (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help") {
      std::cerr << std::format("\nUsage: {} <n-combination>\n\n", argv[0]);
      print_help(combinations, names);
      std::exit(EXIT_SUCCESS);
    }
    return std::atoi(argv[1]);
  }

}  // namespace BLUE
