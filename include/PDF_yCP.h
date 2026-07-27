/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#ifndef PDF_yCP_h
#define PDF_yCP_h

#include "CharmUtils.h"

#include <PDF_Charm.h>

#include <TString.h>

#include <set>
#include <string>

class PDF_yCP : public PDF_Charm {
 public:
  PDF_yCP(TString measurement_id, parametrisations::mix mix_param);
  void buildPdf() override;
  void initObservables() override;
  void initRelations() override;
  void setCorrelations(TString measurement_id) override;
  void setObservables(TString measurement_id) override;
  void setUncertainties(TString measurement_id) override;

 private:
  std::set<std::string> getParameterNames() const override;
  const parametrisations::mix mix_param;
  const TString measurement_id;
};

#endif
