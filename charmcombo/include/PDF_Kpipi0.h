/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#ifndef PDF_Kpipi0_h
#define PDF_Kpipi0_h

#include <TString.h>

#include <PDF_Abs.h>

#include "CharmUtils.h"

class PDF_Kpipi0 : public PDF_Abs {
 public:
  PDF_Kpipi0(TString measurement_id, const theory_config& th_cfg);
  ~PDF_Kpipi0();
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
