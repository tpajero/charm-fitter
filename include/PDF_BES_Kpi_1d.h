/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#ifndef PDF_BES_Kpi_1d_h
#define PDF_BES_Kpi_1d_h

#include "CharmUtils.h"

#include <PDF_Abs.h>

#include <TString.h>

class PDF_BES_Kpi_1d : public PDF_Abs {
 public:
  PDF_BES_Kpi_1d(theory_config th_cfg);
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
