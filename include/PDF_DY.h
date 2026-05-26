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
 * Models measurements of DeltaY(D0 -> h- h+).
 *
 * Depending on the value of `dy_fsc_hypo`, separate measurements for the K- K+ and pi- pi+ final states are taken into
 * account, or a single average is employed.
 */
class PDF_DY : public PDF_Charm {
 public:
  PDF_DY(TString measurement_id, hypotheses::dy_fsc, parametrisations::acp, parametrisations::mix);
  PDF_DY(TString obs_id, TString unc_id, hypotheses::dy_fsc, parametrisations::acp, parametrisations::mix);
  void initObservables() override;
  void initRelations() override;
  void setCorrelations(TString measurement_id) override;
  void setObservables(TString measurement_id) override;
  void setUncertainties(TString measurement_id) override;

 private:
  std::set<std::string> getParameterNames() const override;
  const hypotheses::dy_fsc dy_fsc_hypo;
  const parametrisations::acp acp_param;
  const parametrisations::mix mix_param;
  const TString measurement_id;
};
