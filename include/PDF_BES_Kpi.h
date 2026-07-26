/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: September 2022
 **/

#ifndef PDF_BES_Kpi_h
#define PDF_BES_Kpi_h

#include "CharmUtils.h"

#include <PDF_Abs.h>

#include <TString.h>

class PDF_BES_Kpi : public PDF_Abs {
 public:
  PDF_BES_Kpi(parametrisations::mix mix_param);
  void buildPdf() override;
  void initObservables() override;
  void initParameters() override;
  void initRelations() override;
  void setCorrelations(TString c) override;
  void setObservables(TString c) override;
  void setUncertainties(TString c) override;

 private:
  const parametrisations::mix mix_param;
};

#endif
