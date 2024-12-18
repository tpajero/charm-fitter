/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#ifndef ParametersCharmCombo_h
#define ParametersCharmCombo_h

#include "ParametersAbs.h"

class ParametersCharmCombo : public ParametersAbs {
 public:
  ParametersCharmCombo();

 private:
  void defineParameters();
};

#endif
