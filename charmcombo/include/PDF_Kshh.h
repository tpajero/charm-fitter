/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#ifndef PDF_Kshh_h
#define PDF_Kshh_h

#include <TString.h>

#include <PDF_Abs.h>

#include "CharmUtils.h"

class PDF_Kshh : public PDF_Abs {
 public:
  PDF_Kshh(TString measurement_id, const theory_config& th_cf);
  ~PDF_Kshh();
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
