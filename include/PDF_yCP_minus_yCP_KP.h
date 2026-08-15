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
 * Model for measurements of the yCP parameter in CP-even final states, using untagged D -> K- pi+ decays as reference
 * channel.
 *
 * Currently used only for the E791 measurement of D0 -> K- K+ decays.
 */
class PDF_yCP_minus_yCP_KP : public PDF_Charm {
 public:
  PDF_yCP_minus_yCP_KP(TString measurement_id, parametrisations::mix mix_param);
  void initObservables() override;
  void initRelations() override;
  void setCorrelations(TString measurement_id) override;
  void setObservables(TString measurement_id) override;
  void setUncertainties(TString measurement_id) override;

 private:
  std::set<std::string> getParameterNames() const override;
  const parametrisations::mix mix_param;
  const TString measurement_id;
};
