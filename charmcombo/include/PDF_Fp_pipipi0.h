/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: September 2022
 **/

#ifndef PDF_Fp_pipipi0_h
#define PDF_Fp_pipipi0_h

#include <TString.h>

#include <PDF_Abs.h>

class PDF_Fp_pipipi0 : public PDF_Abs {
 public:
  PDF_Fp_pipipi0(TString measurement_id);
  void buildPdf() override;
  void initObservables(TString measurement_id);
  void initParameters() override;
  void initRelations() override;
  void setCorrelations(TString measurement_id) override;
  void setObservables(TString measurement_id) override;
  void setUncertainties(TString measurement_id) override;
};

#endif
