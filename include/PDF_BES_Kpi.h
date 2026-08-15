/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#pragma once

#include "CharmUtils.h"

#include <PDF_Charm.h>

#include <TString.h>

#include <set>
#include <string>

/**
 * Models the measurement of strong parameters of D0 -> K+ pi- decays by the BESIII collaboration.
 */
class PDF_BES_Kpi : public PDF_Charm {
 public:
  PDF_BES_Kpi(parametrisations::mix mix_param);
  void initObservables() override;
  void initRelations() override;
  void setCorrelations(TString measurement_id) override;
  void setObservables(TString measurement_id) override;
  void setUncertainties(TString measurement_id) override;

 private:
  std::set<std::string> getParameterNames() const override;
  const parametrisations::mix mix_param;
};
