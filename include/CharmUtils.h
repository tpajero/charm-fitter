#ifndef CHARM_UTILS_H
#define CHARM_UTILS_H

#include <ostream>
#include <string>

namespace constants {
  constexpr double d0_lifetime = 4.103e-13;
}  // namespace constants

namespace parametrisations {
  /**
   * Sets how to parametrise mixing observables in the fit.
   *
   * The options are:
   *   - pheno:     (x, y, phi, |q/p| - 1);
   *   - theo:      (x12, y12, phiM, phiG);
   *   - d0_to_kpi: (y', x'^2, dy', dx'2) (to be used for WS/RS measurements of D0 -> K+ pi- only).
   */
  enum class mix { pheno, theo, d0_to_kpi };

  /**
   * Sets how to parametrise the final-state dependent correction to DeltaY in D0 -> h- h+ decays.
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
   *              The final-state correction is parametrised as a function of cot(delta_HH).
   *
   * The default choice is "none", to improve the knowledge of aCP(KK) under the assumption that r_HH
   * cannot be much larger than unity (one would start to see significant corrections at LHCb Upgrade II only if
   * r_HH = O(10)). On the other hand, "full" should be used for completely unbiased estimates.
   */
  enum class dy_fsc { none, partial, full };

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
  std::string dy_hh_expression(parametrisations::mix, parametrisations::dy_fsc = parametrisations::dy_fsc::none,
                               std::string fs = "");
  /// DeltaY(D0 -> K- pi+).
  std::string dy_kp_expression(parametrisations::mix);

  std::string get_id(parametrisations::mix par);
  std::string get_id(parametrisations::dy_fsc par);
  std::string to_string(parametrisations::mix par);
  std::string to_string(parametrisations::dy_fsc par);
}  // namespace utils

std::ostream& operator<<(std::ostream& os, parametrisations::mix par);
std::ostream& operator<<(std::ostream& os, parametrisations::dy_fsc par);

#endif
