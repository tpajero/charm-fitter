/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: September 2022
 **/

#pragma once

#include <PDF_Charm.h>

#include <TString.h>

#include <set>
#include <string>

/**
 * Models measurements of the CP-even fraction of D0 -> pi- pi+ pi0 decays.
 */
class PDF_Fp_pipipi0 : public PDF_Charm {
 public:
  PDF_Fp_pipipi0(TString measurement_id);
  void initObservables() override;
  void initRelations() override;
  void setCorrelations(TString measurement_id) override;
  void setObservables(TString measurement_id) override;
  void setUncertainties(TString measurement_id) override;

 private:
  std::set<std::string> getParameterNames() const override;
  const TString measurement_id;
};
