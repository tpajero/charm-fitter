#include <string>

namespace CharmUtils {
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
