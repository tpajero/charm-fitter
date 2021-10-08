/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#ifndef ParametersCharmCombo_h
#define ParametersCharmCombo_h

#include "ParametersAbs.h"
#include "Utils.h"

using namespace std;
using namespace Utils;


enum theory_config {
    phenomenological,
    theoretical,
    superweak
};


class ParametersCharmCombo : public ParametersAbs
{
public:
    ParametersCharmCombo(theory_config parametrisation,
                         bool predict_dy_rs=kFALSE);

protected:
    void defineParameters(theory_config parametrisation,
                          bool predict_dy_rs);
};

#endif
