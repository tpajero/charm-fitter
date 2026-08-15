/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: October 2025
 **/

#pragma once

#include "CharmUtils.h"

#include <PDF_Charm.h>

#include <TString.h>

#include <set>
#include <string>

/**
 * Models the simultaneous measurements of strong parameters of D0 -> K+ pi- and D0 -> pi+ pi- pi0 decays by the BESIII
 * collaboration.
 */
class PDF_BES_Kpi_pipipi0 : public PDF_Charm {
 public:
  PDF_BES_Kpi_pipipi0(TString, parametrisations::mix mix_param);
  void initObservables() override;
  void initRelations() override;
  void setCorrelations(TString c) override;
  void setObservables(TString c) override;
  void setUncertainties(TString c) override;

 private:
  std::set<std::string> getParameterNames() const override;
  const parametrisations::mix mix_param;
  const TString measurement_id;
};
