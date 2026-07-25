/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#ifndef PDF_WSMixingNoCPV_h
#define PDF_WSMixingNoCPV_h

#include "CharmUtils.h"

#include <PDF_Abs.h>

#include <TString.h>

class PDF_WS_NoCPV : public PDF_Abs {
 public:
  PDF_WS_NoCPV(TString measurement_id, parametrisations::mix mix_param);
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
