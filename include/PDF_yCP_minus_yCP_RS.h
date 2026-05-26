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
 * Model for measurements of the yCP parameter in CP-even final states, using D0 -> K- pi+ decays as reference channel.
 */
class PDF_yCP_minus_yCP_RS : public PDF_Charm {
 public:
  PDF_yCP_minus_yCP_RS(TString measurement_id, parametrisations::mix mix_param);
  PDF_yCP_minus_yCP_RS(TString obs_id, TString unc_id, parametrisations::mix mix_param);
  void initObservables() override;
  void initRelations() override;
  void setObservables(TString measurement_id) override;
  void setUncertainties(TString measurement_id) override;

 private:
  std::set<std::string> getParameterNames() const override;
  const parametrisations::mix mix_param;
  const TString measurement_id;
};
