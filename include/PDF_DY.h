/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#ifndef PDF_DY_h
#define PDF_DY_h

#include "CharmUtils.h"

#include <PDF_Charm.h>

#include <TString.h>

#include <set>
#include <string>

class PDF_DY : public PDF_Charm {
 public:
  PDF_DY(TString measurement_id, parametrisations::mix mix_param, parametrisations::dy_fsc dy_fsc_param);
  void initObservables() override;
  void initRelations() override;
  void setCorrelations(TString measurement_id) override;
  void setObservables(TString measurement_id) override;
  void setUncertainties(TString measurement_id) override;

 private:
  std::set<std::string> getParameterNames() const override;
  const parametrisations::dy_fsc dy_fsc_param;
  const parametrisations::mix mix_param;
  const TString measurement_id;
};

#endif
