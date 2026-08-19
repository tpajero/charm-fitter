/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: June 2024
 **/

#pragma once

#include <PDF_Charm.h>

#include <TString.h>

#include <set>
#include <string>

/**
 * Models the combination of BESIII and CLEO-c measurements of the coherence factors and strong-phase differences
 * in D0 -> K+ pi- pi- pi+ and D0 -> K+ pi- pi0 decays.
 */
class PDF_BES_CLEO_K3pi_Kpipi0 : public PDF_Charm {
 public:
  PDF_BES_CLEO_K3pi_Kpipi0(TString measurement_id);
  void initObservables() override;
  void initRelations() override;
  void setCorrelations(TString measurement_id) override;
  void setObservables(TString measurement_id) override;
  void setUncertainties(TString measurement_id) override;

 private:
  std::set<std::string> getParameterNames() const override;
};
