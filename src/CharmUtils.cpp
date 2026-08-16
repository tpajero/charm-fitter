#include <CharmUtils.h>

#include <boost/algorithm/string.hpp>

#include <format>
#include <iostream>
#include <map>
#include <ostream>
#include <set>
#include <stdexcept>
#include <string>

namespace {
  std::string str_repr(const hypotheses::dy_fsc hypo, const bool id) {
    using hypotheses::dy_fsc;
    switch (hypo) {
    case dy_fsc::none:
      return id ? "no" : "(no final-state correction)";
    case dy_fsc::partial:
      return id ? "partial" : "(y12 x a_CP term only)";
    case dy_fsc::full:
      return id ? "full" : "(full correction including delta_HH dependence)";
    default:
      throw std::runtime_error(
          std::format("ERROR hypotheses::dy_fsc {} not supported by \"str_repr\"", static_cast<int>(hypo)));
    }
  }

  std::string str_repr(const parametrisations::acp par, const bool id) {
    using parametrisations::acp;
    switch (par) {
    case acp::acp_dy:
      return id ? "acp-dy" : "(aCP, DeltaY_HH)";
    case acp::acp_cot:
      return id ? "acp-cot" : "(aCP, cot(delta_HH))";
    case acp::r_delta:
      return id ? "r-delta" : "(r_HH, delta_HH)";
    default:
      throw std::runtime_error(
          std::format("ERROR parametrisations::acp {} not supported by \"str_repr\"", static_cast<int>(par)));
    }
  }

  std::string str_repr(const parametrisations::mix par, const bool id) {
    using parametrisations::mix;
    switch (par) {
    case mix::pheno:
      return id ? "pheno" : "(x, y, q/p, phi)";
    case mix::theo:
      return id ? "theo" : "(x12, y12, phiM, phiG)";
    case mix::d0_to_kpi:
      return id ? "d0_to_kpi" : "(y', x'^2, dy', dx'2)";
    default:
      throw std::runtime_error(
          std::format("ERROR parametrisations::mix {} not supported by \"str_repr\"", static_cast<int>(par)));
    }
  }

  std::set<std::string> r_delta_parameter_names(const std::set<std::string>& fs) {
    std::set<std::string> names;
    for (const auto& hh : fs) names.insert({std::format("r_{}", hh), std::format("delta_{}", hh)});
    return names;
  }

}  // namespace

std::string utils::x_expression(const parametrisations::mix mix_param) {
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    return "x";
  case mix::theo:
    return "1/sqrt(2) * sqrt("
           "    x12*x12 - y12*y12 + sqrt("
           "       + TMath::Sq(x12*x12 + y12*y12)"
           "       - TMath::Sq(2 * x12 * y12 * sin(phiM - phiG))"
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
    return "1/sqrt(2) * sqrt("
           "    y12*y12 - x12*x12 + sqrt("
           "       + TMath::Sq(x12*x12 + y12*y12)"
           "       - TMath::Sq(2 * x12 * y12 * sin(phiM - phiG))"
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

std::string utils::dy_hh_expression(const hypotheses::dy_fsc dy_fsc_hypo, const parametrisations::acp acp_param,
                                    const parametrisations::mix mix_param, const std::string fs) {
  check_compatibility(dy_fsc_hypo, acp_param);

  const std::string fs_independent = std::format("- {}", utils::dy_expression(mix_param));
  std::string fs_dependent;

  using hypotheses::dy_fsc;
  switch (dy_fsc_hypo) {
  case dy_fsc::none:
    return fs_independent;
  case dy_fsc::partial:
    fs_dependent = "y * Acp_HH";
    break;
  case dy_fsc::full:
    switch (acp_param) {
    case parametrisations::acp::acp_dy:
      return std::format("DY_{}", fs);
    case parametrisations::acp::acp_cot:
      fs_dependent = "y * Acp_HH * (1 + x / y * cot_delta_HH)";
      break;
    case parametrisations::acp::r_delta:
      fs_dependent = std::format("{} ({}) * r_HH * (x * cos(delta_HH) + y * sin(delta_HH))", fs == "KK" ? "+" : "-",
                                 constants::acp_prefix);
      break;
    default:
      throw std::runtime_error(
          std::format("utils::dy_hh_expression ERROR {} not supported for final-state correction {}",
                      utils::to_string(acp_param), utils::to_string(dy_fsc_hypo)));
    }
    break;
  default:
    throw std::runtime_error(
        std::format("utils::dy_hh_expression ERROR {} not supported", utils::to_string(dy_fsc_hypo)));
  }

  boost::replace_all(fs_dependent, "HH", fs);
  if (mix_param == parametrisations::mix::theo) {
    boost::replace_all(fs_dependent, "x", "x12");
    boost::replace_all(fs_dependent, "y", "y12");
  }

  return std::format("{} + {}", fs_independent, fs_dependent);
}

std::string utils::dy_kp_expression(const parametrisations::mix mix_param) {
  using parametrisations::mix;
  switch (mix_param) {
  case mix::pheno:
    return "0.5 * r_Kpi * (  (y*cos(Delta_Kpi) - x*sin(Delta_Kpi))*(qop - 1/qop - Acp_KP)*cos(phi) "
           "               - (x*cos(Delta_Kpi) + y*sin(Delta_Kpi))*(qop + 1/qop         )*sin(phi))";
    break;
  case mix::theo:
    return "r_Kpi * (  (-y12*cos(Delta_Kpi)*cos(phiG) + x12*sin(Delta_Kpi)*cos(phiM))*Acp_KP*0.5 "
           "         + ( y12*sin(Delta_Kpi)*sin(phiG) + x12*cos(Delta_Kpi)*sin(phiM))           )";
    break;
  default:
    throw std::runtime_error(
        std::format("utils::dy_kp_expression ERROR {} not supported", utils::to_string(mix_param)));
  }
}

std::set<std::string> utils::acp_hh_parameters_names(const parametrisations::acp acp_param,
                                                     const std::set<std::string>& fs) {
  using parametrisations::acp;

  std::set<std::string> names;
  switch (acp_param) {
  case acp::acp_dy:
    [[fallthrough]];
  case acp::acp_cot:
    for (const auto& hh : fs) names.insert(std::format("Acp_{}", hh));
    break;
  case acp::r_delta:
    names = r_delta_parameter_names(fs);
    break;
  default:
    throw std::runtime_error(std::format("utils::acp_hh_parameters_names ERROR ACP parametrisation {} not supported",
                                         utils::to_string(acp_param)));
  }
  return names;
}

std::set<std::string> utils::dy_hh_parameters_names(const hypotheses::dy_fsc dy_fsc_hypo,
                                                    const parametrisations::acp acp_param,
                                                    const parametrisations::mix mix_param,
                                                    const std::set<std::string>& fs) {
  using hypotheses::dy_fsc;
  using parametrisations::acp;
  using parametrisations::mix;

  check_compatibility(dy_fsc_hypo, acp_param);

  std::set<std::string> names;
  switch (dy_fsc_hypo) {
  case dy_fsc::none:
    break;
  case dy_fsc::partial:
    for (const auto& hh : fs) names.insert(std::format("Acp_{}", hh));
    break;
  case dy_fsc::full:
    switch (acp_param) {
    case acp::acp_dy:
      for (const auto& hh : fs) names.insert(std::format("DY_{}", hh));
      return names;
    case acp::acp_cot:
      for (const auto& hh : fs) names.insert({std::format("Acp_{}", hh), std::format("cot_delta_{}", hh)});
      break;
    case acp::r_delta:
      names = r_delta_parameter_names(fs);
      break;
    default:
      throw std::runtime_error(std::format("utils::dy_hh_parameters_names ERROR ACP parametrisation {} not supported",
                                           utils::to_string(acp_param)));
    }
    break;
  }
  switch (mix_param) {
  case mix::pheno:
    names.insert({"x", "y", "qop", "phi"});
    break;
  case mix::theo:
    names.insert({"x12", "phiM"});
    if (dy_fsc_hypo != dy_fsc::none) names.insert("y12");
    break;
  default:
    throw std::runtime_error(std::format("utils::dy_hh_parameters_names ERROR Parametrisation {} not supported",
                                         utils::to_string(mix_param)));
  }
  return names;
}

void utils::check_compatibility(const hypotheses::dy_fsc dy_fsc_hypo, const parametrisations::acp acp_param) {
  using hypotheses::dy_fsc;
  using parametrisations::acp;
  switch (dy_fsc_hypo) {
  case dy_fsc::none:
    [[fallthrough]];
  case dy_fsc::partial:
    if (acp_param == acp::r_delta)
      throw std::runtime_error(std::format("utils::check_compatibility ERROR {} parametrisation not compatible with {} "
                                           "final-state correction to DeltaY(h- h+)",
                                           utils::to_string(acp_param), utils::to_string(dy_fsc_hypo)));
    [[fallthrough]];
  case dy_fsc::full:
    return;
  }
}

std::string utils::get_id(const hypotheses::dy_fsc par) { return str_repr(par, true); }
std::string utils::get_id(const parametrisations::acp par) { return str_repr(par, true); }
std::string utils::get_id(const parametrisations::mix par) { return str_repr(par, true); }

std::string utils::to_string(const hypotheses::dy_fsc par) { return str_repr(par, false); }
std::string utils::to_string(const parametrisations::acp par) { return str_repr(par, false); }
std::string utils::to_string(const parametrisations::mix par) { return str_repr(par, false); }

std::ostream& operator<<(std::ostream& os, const hypotheses::dy_fsc param) {
  os << utils::to_string(param);
  return os;
}

std::ostream& operator<<(std::ostream& os, const parametrisations::acp param) {
  os << utils::to_string(param);
  return os;
}

std::ostream& operator<<(std::ostream& os, const parametrisations::mix param) {
  os << utils::to_string(param);
  return os;
}

std::string utils::acp_expression(const parametrisations::acp acp_param, const std::string final_state) {
  if (final_state != "KK" && final_state != "PP")
    throw std::runtime_error(std::format("utils::acp_expression ERROR final state {} not supported", final_state));

  using parametrisations::acp;
  switch (acp_param) {
  case acp::acp_dy:
    [[fallthrough]];
  case acp::acp_cot:
    return std::format("Acp_{}", final_state);
  case acp::r_delta:
    return std::format("{} {:.3e} * sin(delta_{})", final_state == "KK" ? "+" : "-", constants::acp_prefix,
                       final_state);
  default:
    throw std::runtime_error("utils::acp_expression ERROR acp_param not supported");
  }
}
