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

namespace parametrisations {
  enum class kpi { raxy, rrxy, ccprime };
}  // namespace parametrisations

/**
 * Models measurements of the WS/RS ratio of D0 -> K+ pi- decays.
 */
class PDF_WS : public PDF_Charm {
 public:
  PDF_WS(TString measurement_id, parametrisations::mix, parametrisations::kpi p = parametrisations::kpi::rrxy,
         hypotheses::dy_fsc dy_fsc_hypo = hypotheses::dy_fsc::none,
         parametrisations::acp acp_param = parametrisations::acp::acp_dy);
  PDF_WS(TString val, TString err, parametrisations::mix, parametrisations::kpi p = parametrisations::kpi::rrxy,
         hypotheses::dy_fsc dy_fsc_hypo = hypotheses::dy_fsc::none,
         parametrisations::acp acp_param = parametrisations::acp::acp_dy);
  void initObservables() override;
  void initRelations() override;
  void setCorrelations(TString measurement_id) override;
  void setObservables(TString measurement_id) override;
  void setUncertainties(TString measurement_id) override;

 private:
  void initRelationsCCPrime();
  void initRelationsRAXY();
  void initRelationsRRXY();

  /// Add the y'+_th and x'2+_th theory relations.
  void addYpXp2Plus();
  /// Add the y'-_th and x'2-_th theory relations.
  void addYpXp2Minus();

  std::set<std::string> getParameterNames() const override;

  const parametrisations::mix mix_param;
  const parametrisations::kpi kpi_param;
  const hypotheses::dy_fsc dy_fsc_hypo;
  const parametrisations::acp acp_param;
  TString label;
};
