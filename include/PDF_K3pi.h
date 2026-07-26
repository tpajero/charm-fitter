/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: June 2024
 **/

#ifndef PDF_K3pi_h
#define PDF_K3pi_h

#include "CharmUtils.h"

#include <PDF_Abs.h>

#include <TString.h>

class PDF_K3pi : public PDF_Abs {
 public:
  PDF_K3pi(TString measurement_id, parametrisations::mix mix_param);
  void buildPdf() override;
  void initObservables(TString setName);
  void initParameters() override;
  void initRelations() override;
  void setCorrelations(TString measurement_id) override;
  void setObservables(TString measurement_id) override;
  void setUncertainties(TString measurement_id) override;

 private:
  const parametrisations::mix mix_param;
};

#endif
