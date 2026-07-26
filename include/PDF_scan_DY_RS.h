/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#ifndef PDF_scan_DY_RS_h
#define PDF_scan_DY_RS_h

#include "CharmUtils.h"

#include <PDF_Abs.h>

#include <TString.h>

class PDF_scan_DY_RS : public PDF_Abs {
 public:
  PDF_scan_DY_RS(parametrisations::mix mix_param);
  void buildPdf() override;
  void initObservables() override;
  void initParameters() override;
  void initRelations() override;
  void setCorrelations();
  void setObservables();
  void setUncertainties();

 private:
  const parametrisations::mix mix_param;
};

#endif
