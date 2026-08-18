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
 * Models measurements of the WS/RS ratio of D0 -> K+ pi- decays, without allowing for CP violation.
 *
 * Should be mostly equivalent to using PDF_WS with CP violation fixed to zero, up to non-Gaussian effects in the
 * likelihood of the various measurements.
 */
class PDF_WS_NoCPV : public PDF_Charm {
 public:
  PDF_WS_NoCPV(TString measurement_id, parametrisations::mix mix_param);
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
