/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: June 2024
 **/

#ifndef PDF_BES_CLEO_K3pi_Kpipi0_h
#define PDF_BES_CLEO_K3pi_Kpipi0_h

#include <PDF_Abs.h>

#include <TString.h>

class PDF_BES_CLEO_K3pi_Kpipi0 : public PDF_Abs {
 public:
  PDF_BES_CLEO_K3pi_Kpipi0(TString measurement_id);
  void buildPdf() override;
  void initObservables(TString measurement_label);
  void initParameters() override;
  void initRelations() override;
  void setCorrelations(TString measurement_id) override;
  void setObservables(TString measurement_id) override;
  void setUncertainties(TString measurement_id) override;
};

#endif
