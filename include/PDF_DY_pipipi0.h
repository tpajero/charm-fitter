/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: June 2024
 **/

#pragma once

#include "CharmUtils.h"

#include <PDF_Charm.h>

#include <TString.h>

#include <set>
#include <string>

/**
 * Models measurements of the DeltaY(D0 -> pi- pi+ pi0) observable.
 *
 * Kept separated from `PDF_DY.h` since it depends on the CP-even fraction of the D0 -> pi- pi+ pi0 decay.
 */
class PDF_DY_pipipi0 : public PDF_Charm {
 public:
  PDF_DY_pipipi0(TString measurement_id, parametrisations::mix mix_param);
  void initObservables() override;
  void initRelations() override;
  void setObservables(TString measurement_id) override;
  void setUncertainties(TString measurement_id) override;

 private:
  std::set<std::string> getParameterNames() const override;
  const parametrisations::mix mix_param;
  const TString measurement_id;
};
