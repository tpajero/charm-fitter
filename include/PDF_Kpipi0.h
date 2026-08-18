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
 * Models the measurements of strong-phase rotated mixing parameters in D0 -> K+ pi- pi0 decays.
 */
class PDF_Kpipi0 : public PDF_Charm {
 public:
  /**
   * Constructor.
   *
   * @param allow_cpv Sets whether to use results allowing for CP violation or not. While one would like to allow for
   *     CP violation to obtain unbiased results, the likelihood of these inputs are less Gaussian.
   */
  PDF_Kpipi0(TString measurement_id, parametrisations::mix mix_param, bool allow_cpv = true);
  void initObservables() override;
  void initRelations() override;
  void setCorrelations(TString measurement_id) override;
  void setObservables(TString measurement_id) override;
  void setUncertainties(TString measurement_id) override;

 private:
  std::set<std::string> getParameterNames() const override;
  const parametrisations::mix mix_param;
  const TString measurement_id;
  const bool allow_cpv;
};
