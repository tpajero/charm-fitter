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
 * Implements any measurement using the bin-flip variables (xCP, yCP, DeltaX, DeltaY).
 *
 * Mostly useful for analyses of D0 -> KS h- h+ decays.
 */
class PDF_BinFlip : public PDF_Charm {
 public:
  PDF_BinFlip(TString measurement_id, parametrisations::mix mix_param);
  PDF_BinFlip(TString obs_id, TString unc_id, parametrisations::mix mix_param);
  void initObservables() override;
  void initRelations() override;
  void setCorrelations(TString measurement_id) override;
  void setObservables(TString measurement_id) override;
  void setUncertainties(TString measurement_id) override;

 private:
  std::set<std::string> getParameterNames() const override;
  const parametrisations::mix mix_param;
  const TString measurement_id;
  const TString unc_id;
};
