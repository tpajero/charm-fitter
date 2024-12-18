/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#ifndef PDF_BES_Kpi_1d_h
#define PDF_BES_Kpi_1d_h

#include <TString.h>

#include <PDF_Abs.h>

#include "CharmUtils.h"

class PDF_BES_Kpi_1d : public PDF_Abs {
 public:
  PDF_BES_Kpi_1d(const theory_config& th_cfg);
  ~PDF_BES_Kpi_1d();
  void buildPdf() override;
  void initObservables(const TString& setName);
  virtual void initParameters() override;
  virtual void initRelations() override;
  void setCorrelations(TString measurement_id) override;
  void setObservables(TString measurement_id) override;
  void setUncertainties(TString measurement_id) override;

 private:
  const theory_config th_cfg;
};

#endif
