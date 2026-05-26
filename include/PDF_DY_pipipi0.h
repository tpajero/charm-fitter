/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: June 2024
 **/

#ifndef PDF_DY_pipipi0_h
#define PDF_DY_pipipi0_h

#include "CharmUtils.h"

#include <PDF_Abs.h>

#include <TString.h>

class PDF_DY_pipipi0 : public PDF_Abs {
 public:
  PDF_DY_pipipi0(TString measurement_id, theory_config th_cf);
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
