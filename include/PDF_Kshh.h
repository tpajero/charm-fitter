/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#ifndef PDF_Kshh_h
#define PDF_Kshh_h

#include "CharmUtils.h"

#include <PDF_Abs.h>

#include <TString.h>

class PDF_Kshh : public PDF_Abs {
 public:
  PDF_Kshh(TString measurement_id, theory_config th_cf);
  void buildPdf() override;
  void initObservables(TString setName);
  void initParameters() override;
  void initRelations() override;
  void setCorrelations(TString measurement_id) override;
  void setObservables(TString measurement_id) override;
  void setUncertainties(TString measurement_id) override;

 private:
  const theory_config th_cfg;
};

#endif
