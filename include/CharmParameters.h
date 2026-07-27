/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#ifndef CharmParameters_h
#define CharmParameters_h

#include <ParametersAbs.h>

class CharmParameters : public ParametersAbs {
 public:
  CharmParameters();

 private:
  void defineParameters();
};

#endif
