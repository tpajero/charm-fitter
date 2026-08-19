/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#pragma once

#include "CharmUtils.h"

#include <PDF_Charm.h>

#include <TString.h>

#include <set>
#include <string>

/**
 * Implements a constraint with very small uncertainty (where "very small" is currently 5e-7 and should be tuned to the
 * precision of the combination) to set an upper bound on the value of |DeltaY(D0 -> K- pi+)|.
 */
class PDF_scan_DY_RS : public PDF_Charm {
 public:
  PDF_scan_DY_RS(parametrisations::mix mix_param);
  void initObservables() override;
  void initRelations() override;
  void setObservables(TString c) override;
  void setUncertainties(TString c) override;

 private:
  std::set<std::string> getParameterNames() const override;
  const parametrisations::mix mix_param;
};
