#ifndef CHARM_UTILS_H
#define CHARM_UTILS_H

#include <ostream>
#include <string>

/**
 * Options for the final-state dependent correction to DeltaY:
 *   - none;
 *   - partial (correct for the y aCP term, but ignore the term
 *              depending on the strong-phase difference);
 *   - full (full correction including the term depending on the
 *           strong-phase difference)
 */
enum class FSC { none, partial, full };

std::ostream& operator<<(std::ostream& os, const FSC& fsc);

/**
 * Enumeration of the possible parametrisations to be used in the fit.
 */
enum class theory_config {
    phenomenological,
    theoretical,
    superweak,
    d0_to_kpi
};

std::ostream& operator<<(std::ostream& os, const theory_config& th_cfg);

namespace CharmUtils {

    std::string get_dy_expression(const theory_config th_cfg, const FSC fsc=FSC::none, const std::string fs="");

    // Helper strings to go from phenomenological to phenomenological parametrisation
    const std::string x_to_superweak =
        "pow(2,-0.5) * sqrt("
        "    pow(x12,2) - pow(y12,2) + sqrt("
        "       + pow(pow(x12,2) + pow(y12,2),2)"
        "       - pow(2 * x12 * y12 * sin(phiM),2)"
        "    )"
        ") * TMath::Sign(1., cos(phiM))";
    const std::string y_to_superweak =
        "pow(2,-0.5) * sqrt("
        "    pow(y12,2) - pow(x12,2) + sqrt("
        "       + pow(pow(x12,2) + pow(y12,2),2)"
        "       - pow(2 * x12 * y12 * sin(phiM),2)"
        "    )"
        ")";

    // Helper strings to go from phenomenological to theoretical parametrisation
    const std::string x_to_theoretical =
        "pow(2,-0.5) * sqrt("
        "    pow(x12,2) - pow(y12,2) + sqrt("
        "       + pow(pow(x12,2) + pow(y12,2),2)"
        "       - pow(2 * x12 * y12 * sin(phiM - phiG),2)"
        "    )"
        ") * TMath::Sign(1., cos(phiM - phiG))";
    const std::string y_to_theoretical =
        "pow(2,-0.5) * sqrt("
        "    pow(y12,2) - pow(x12,2) + sqrt("
        "       + pow(pow(x12,2) + pow(y12,2),2)"
        "       - pow(2 * x12 * y12 * sin(phiM - phiG),2)"
        "    )"
        ")";
}

#endif
