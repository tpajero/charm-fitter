/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#pragma once

#include <ParametersAbs.h>

class CharmParameters : public ParametersAbs {
 public:
  CharmParameters();

 private:
  void defineParameters();
};
