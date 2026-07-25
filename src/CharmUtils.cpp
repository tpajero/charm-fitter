#include <CharmUtils.h>

#include <boost/algorithm/string.hpp>

#include <format>
#include <iostream>
#include <map>
#include <ostream>
#include <stdexcept>
#include <string>

namespace {
  std::string str_repr(const parametrisations::mix par, const bool id) {
    using parametrisations::mix;
    switch (par) {
    case mix::pheno:
      return id ? "pheno" : "(x, y, q/p, phi)";
    case mix::theo:
      return id ? "theo" : "(x12, y12, phiM, phiG)";
    default:
      throw std::runtime_error(
          std::format("ERROR parametrisations::mix {} not supported by \"str_repr\"", static_cast<int>(par)));
    }
  }

  std::string str_repr(const parametrisations::dy_fsc par, const bool id) {
    using parametrisations::dy_fsc;
    switch (par) {
    case dy_fsc::none:
      return id ? "no_dyfsc" : "(no final-state correction)";
    case dy_fsc::partial:
      return id ? "partial_dyfsc" : "(y12 x a_CP term only)";
    case dy_fsc::full:
      return id ? "full_dyfsc" : "(full correction including delta_HH dependence)";
    default:
      throw std::runtime_error(
          std::format("ERROR parametrisations::dy_fsc {} not supported by \"str_repr\"", static_cast<int>(par)));
    }
  }
}  // namespace

std::string utils::x_expression(const parametrisations::mix mix_param) {
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    return "x";
  case mix::theo:
    return "pow(2,-0.5) * sqrt("
           "    pow(x12,2) - pow(y12,2) + sqrt("
           "       + pow(pow(x12,2) + pow(y12,2),2)"
           "       - pow(2 * x12 * y12 * sin(phiM - phiG),2)"
           "    )"
           ") * TMath::Sign(1., cos(phiM - phiG))";
  default:
    throw std::runtime_error(
        std::format("utils::x_expression ERROR Parametrisation {} not supported", utils::to_string(mix_param)));
  }
}

std::string utils::y_expression(const parametrisations::mix mix_param) {
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    return "y";
  case mix::theo:
    return "pow(2,-0.5) * sqrt("
           "    pow(y12,2) - pow(x12,2) + sqrt("
           "       + pow(pow(x12,2) + pow(y12,2),2)"
           "       - pow(2 * x12 * y12 * sin(phiM - phiG),2)"
           "    )"
           ")";
  default:
    throw std::runtime_error(
        std::format("utils::y_expression ERROR Parametrisation {} not supported", utils::to_string(mix_param)));
  }
}

std::string utils::dx_expression(const parametrisations::mix mix_param) {
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    return "0.5 * (x * cos(phi) * (qop - 1/qop) + y * sin(phi) * (qop + 1/qop))";
  case mix::theo:
    return "-y12 * sin(phiG)";
  default:
    throw std::runtime_error(
        std::format("utils::dx_expression ERROR Parametrisation {} not supported", utils::to_string(mix_param)));
  }
}

std::string utils::dy_expression(const parametrisations::mix mix_param) {
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    return "0.5 * (y * cos(phi) * (qop - 1/qop) - x * sin(phi) * (qop + 1/qop))";
  case mix::theo:
    return "x12 * sin(phiM)";
  default:
    throw std::runtime_error(
        std::format("utils::dy_expression ERROR Parametrisation {} not supported", utils::to_string(mix_param)));
  }
}

std::string utils::dy_hh_expression(const parametrisations::mix mix_param, const parametrisations::dy_fsc dy_fsc_param,
                                    const std::string fs) {
  const std::string fs_independent = std::format("- {}", utils::dy_expression(mix_param));

  std::string fs_dependent;
  using parametrisations::dy_fsc;
  switch (dy_fsc_param) {
  case dy_fsc::none:
    break;
  case dy_fsc::partial:
    fs_dependent = "y * Acp_HH";
    break;
  case dy_fsc::full:
    fs_dependent = "y * Acp_HH * (1 + x / y * cot_delta_HH)";
    break;
  default:
    throw std::runtime_error(
        std::format("utils::dy_expression ERROR {} not supported", utils::to_string(dy_fsc_param)));
  }
  boost::replace_all(fs_dependent, "HH", fs);
  if (mix_param == parametrisations::mix::theo) {
    boost::replace_all(fs_dependent, "x", "x12");
    boost::replace_all(fs_dependent, "y", "y12");
  }

  return fs_dependent.empty() ? fs_independent : std::format("{} + {}", fs_independent, fs_dependent);
}

std::string utils::get_id(const parametrisations::mix par) { return str_repr(par, true); }

std::string utils::get_id(const parametrisations::dy_fsc par) { return str_repr(par, true); }

std::string utils::to_string(const parametrisations::mix par) { return str_repr(par, false); }

std::string utils::to_string(const parametrisations::dy_fsc par) { return str_repr(par, false); }

std::ostream& operator<<(std::ostream& os, const parametrisations::mix param) {
  os << utils::to_string(param);
  return os;
}

std::ostream& operator<<(std::ostream& os, const parametrisations::dy_fsc param) {
  os << utils::to_string(param);
  return os;
}
