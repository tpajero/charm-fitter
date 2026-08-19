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
 * Models measurements of the mixing parameter RM = (x^2 + y^2) / 2, as defined by the HFLAV.
 *
 * This is mostly used for measurements of semileptonic decays, but can also be employed to model part of the
 * information of the LHCb Run 1 measurement of the WS/RS ratio with D0 -> K- pi- pi+ pi+ decays (NOT recommended).
 */
class PDF_RM : public PDF_Charm {
 public:
  PDF_RM(TString measurement_id, parametrisations::mix mix_param);
  void initObservables() override;
  void initRelations() override;
  void setObservables(TString measurement_id) override;
  void setUncertainties(TString measurement_id) override;

 private:
  std::set<std::string> getParameterNames() const override;
  const parametrisations::mix mix_param;
  const TString measurement_id;
};
