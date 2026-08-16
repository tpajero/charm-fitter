/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: October 2025
 **/

#ifndef PDF_BES_Kpi_7d_h
#define PDF_BES_Kpi_7d_h

#include "CharmUtils.h"

#include <PDF_Charm.h>

#include <TString.h>

#include <set>
#include <string>

class PDF_BES_Kpi_7d : public PDF_Charm {
 public:
  PDF_BES_Kpi_7d(parametrisations::mix mix_param);
  void initObservables() override;
  void initRelations() override;
  void setCorrelations(TString c) override;
  void setObservables(TString c) override;
  void setUncertainties(TString c) override;

 private:
  std::set<std::string> getParameterNames() const override;
  const parametrisations::mix mix_param;
};

#endif
