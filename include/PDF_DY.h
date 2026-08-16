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
  PDF_DY(TString measurement_id, hypotheses::dy_fsc, parametrisations::acp, parametrisations::mix);
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

#endif
