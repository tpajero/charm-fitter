/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: October 2025
 **/

#ifndef PDF_BES_Kpi_7d_h
#define PDF_BES_Kpi_7d_h

#include <TString.h>

#include <PDF_Abs.h>

#include "CharmUtils.h"

class PDF_BES_Kpi_7d : public PDF_Abs {
 public:
  PDF_BES_Kpi_7d(const theory_config& th_cfg);
  void buildPdf() override;
  void initObservables() override;
  void initParameters() override;
  void initRelations() override;
  void setCorrelations(TString c) override;
  void setObservables(TString c) override;
  void setUncertainties(TString c) override;

 private:
  const theory_config th_cfg;
};

#endif
