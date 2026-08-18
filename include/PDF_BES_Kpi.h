/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: September 2022
 **/

#pragma once

#include "CharmUtils.h"

#include <PDF_Charm.h>

#include <TString.h>

#include <set>
#include <string>

/**
 * Models the first measurement of strong parameters of D0 -> K+ pi- decays by the BESIII collaboration.
 *
 * References:
 *   - http://inspirehep.net/record/1291279
 */
class PDF_BES_Kpi : public PDF_Charm {
 public:
  PDF_BES_Kpi(parametrisations::mix mix_param);
  void initObservables() override;
  void initRelations() override;
  void setObservables(TString c) override;
  void setUncertainties(TString c) override;

 private:
  std::set<std::string> getParameterNames() const override;
  const parametrisations::mix mix_param;
};
