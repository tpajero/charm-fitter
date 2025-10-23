/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#ifndef PDF_scan_DY_RS_h
#define PDF_scan_DY_RS_h

#include <TString.h>

#include <PDF_Abs.h>

#include "CharmUtils.h"

class PDF_scan_DY_RS : public PDF_Abs {
 public:
  PDF_scan_DY_RS(const theory_config& th_cfg);
  void buildPdf() override;
  virtual void initObservables() override;
  virtual void initParameters() override;
  virtual void initRelations() override;
  void setCorrelations();
  void setObservables();
  void setUncertainties();

 private:
  const theory_config th_cfg;
};

#endif
