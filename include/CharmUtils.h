#ifndef CHARM_UTILS_H
#define CHARM_UTILS_H

#include <ostream>
#include <set>
#include <string>

namespace constants {
  constexpr double d0_lifetime = 4.103e-13;

  /**
   * CKM constant responsible for CP violation in charm decays:
   *
   *    acp_prefix = - 2 * Imag((V*cb Vub) / (V*cs Vus)).
   *
   * The CP asymmetries in the decay can then be expressed as:
   *
   *    ACP_KK ~ + acp_prefix * sin(delta_KK),
   *    ACP_PP ~ - acp_prefix * sin(delta_PP).
   */
  constexpr double acp_prefix = 1.3e-3;
}  // namespace constants

namespace hypotheses {
  /**
   * Sets what final-state dependent correction to DeltaY in D0 -> h- h+ decays should be applied.
   *
   * The slope of the time-dependent asymmetry of the decay widths of D0 -> h- h+ decays can be written as
   *
   *     DeltaY = - x12 * sin(phiM) + y12 * aCP_HH * (1 + (x12 / y12) * cot(delta_HH)),
   *
   * see e.g. eq. (A.10) of https://doi.org/10.1103/PhysRevD.104.072010.
   *
   * The options are:
   *   - none:    Neglect all final-state dependent corrections, i.e. DeltaY = - x12 * sin(phiM);
   *   - partial: Correct for the final-state dependence that is precisely known experimentally,
   *              i.e. use the known value of aCP_HH but set cot(delta_HH) = 0. This is the default choice of HFLAV,
   *              but is not fully consistent, since the term proportional to cot(delta_HH) is a priori of the same
   *              size of (and potentially much larger than) y12 * aCP_HH.
   *   - full:    Full final-state correction. Unless one is interested in determining the value of delta_HH,
   *              if one neglects the small dependence of the recent measurements of the WS/RS ratio in D0 -> K+ pi-
   *              decays on delta_KK and DeltaY(KK) (cf. Appendix B of https://doi.org/10.1103/PhysRevD.111.012001),
   *              this is equivalent to decoupling the (aCP, DeltaY) fit from that of all other observables, i.e.
   *              DeltaY(h- h+) cannot be used to improve the knowledge of phiM any longer.
   *
   * The default choice is "none", to improve the knowledge of aCP(KK) under the assumption that r_HH
   * cannot be much larger than unity (one would start to see significant corrections at LHCb Upgrade II only if
   * r_HH = O(10)). On the other hand, "full" should be used for unbiased estimates.
   */
  enum class dy_fsc { none, partial, full };
}  // namespace hypotheses

namespace parametrisations {
  /**
   * Sets how to parametrise the CP asymmetry in the decay.
   *
   * The options are:
   *   - acp_dy: (aCP, DeltaY_HH), ie the direct CP asymmetry and the slope of the time-dependent asymmetry;
   *   - acp_cot: (aCP, cot(delta_HH)), ie the direct CP asymmetry and the cotangent of the relative strong phase
   *       between the penguin and tree amplitudes;
   *   - r_delta: (r_HH, delta_HH), ie the magnitude and the strong phase of the penguin-to-tree ratio.
   */
  enum class acp { acp_dy, acp_cot, r_delta };

  /**
   * Sets how to parametrise mixing observables in the fit.
   *
   * The options are:
   *   - pheno:     (x, y, phi, |q/p| - 1);
   *   - theo:      (x12, y12, phiM, phiG);
   *   - d0_to_kpi: (y', x'^2, dy', dx'2) (to be used for WS/RS measurements of D0 -> K+ pi- only).
   */
  enum class mix { pheno, theo, d0_to_kpi };
}  // namespace parametrisations

namespace utils {
  // Helper functions to get the expressions of observables in different parametrisations

  /// Mixing parameter x of the phenomenological parametrisation.
  std::string x_expression(parametrisations::mix);
  /// Mixing parameter y of the phenomenological parametrisation.
  std::string y_expression(parametrisations::mix);
  /// Observable DeltaX of the bin-flip parametrisation.
  std::string dx_expression(parametrisations::mix);
  /// Observable DeltaY of the bin-flip parametrisation.
  std::string dy_expression(parametrisations::mix);
  /// DeltaY(D0 -> h- h+).
  std::string dy_hh_expression(hypotheses::dy_fsc, parametrisations::acp, parametrisations::mix, std::string fs = "");
  /// DeltaY(D0 -> K- pi+).
  std::string dy_kp_expression(parametrisations::mix);
  /// A_CP(h+ h-), as a function of either Acp_HH or (r_HH, delta_HH), depending on the acp_param choice.
  std::string acp_expression(parametrisations::acp, std::string fs);

  /// Get the theory parameter names for describing aCP(D0 -> h- h+), for a set of final states h- h+.
  std::set<std::string> acp_hh_parameters_names(parametrisations::acp, const std::set<std::string>& fs = {});

  /// Get the theory parameter names for describing DeltaY(D0 -> h- h+), for a set of final states h- h+.
  std::set<std::string> dy_hh_parameters_names(hypotheses::dy_fsc, parametrisations::acp, parametrisations::mix,
                                               const std::set<std::string>& fs = {});

  /// Check whether a final-state correction to DeltaY is compatible with the a parametrisation of aCP.
  void check_compatibility(hypotheses::dy_fsc, parametrisations::acp);

  std::string get_id(hypotheses::dy_fsc);
  std::string get_id(parametrisations::mix);
  std::string get_id(parametrisations::acp);
  std::string to_string(hypotheses::dy_fsc);
  std::string to_string(parametrisations::mix);
  std::string to_string(parametrisations::acp);
}  // namespace utils

std::ostream& operator<<(std::ostream& os, hypotheses::dy_fsc par);
std::ostream& operator<<(std::ostream& os, parametrisations::acp par);
std::ostream& operator<<(std::ostream& os, parametrisations::mix par);

#endif
