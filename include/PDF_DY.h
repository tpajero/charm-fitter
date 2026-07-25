/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#ifndef PDF_DY_h
#define PDF_DY_h

#include "CharmUtils.h"

#include <PDF_Abs.h>

#include <TString.h>

class PDF_DY : public PDF_Abs {
 public:
  PDF_DY(TString measurement_id, parametrisations::mix mix_param, parametrisations::dy_fsc dy_fsc_param);
  void buildPdf() override;
  void initObservables(TString setName);
  void initParameters() override;
  void initRelations() override;
  void setCorrelations(TString measurement_id) override;
  void setObservables(TString measurement_id) override;
  void setUncertainties(TString measurement_id) override;

 private:
  const parametrisations::dy_fsc dy_fsc_param;
  const parametrisations::mix mix_param;
};

#endif
