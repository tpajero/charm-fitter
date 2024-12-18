/**
 * Charm Combination
 * Author: tommaso.pajero@cern.ch
 * Date: June 2024
 **/

#ifndef PDF_DY_pipipi0_h
#define PDF_DY_pipipi0_h

#include <TString.h>

#include <PDF_Abs.h>

#include "CharmUtils.h"

class PDF_DY_pipipi0 : public PDF_Abs {
 public:
  PDF_DY_pipipi0(TString measurement_id, const theory_config th_cf);
  ~PDF_DY_pipipi0();
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
