#ifndef PDF_AcpHH_Nico_h
#define PDF_AcpHH_Nico_h

#include <TString.h>

#include <PDF_Abs.h>

#include "CharmUtils.h"

class PDF_AcpHH_Nico : public PDF_Abs {
 public:
  PDF_AcpHH_Nico(const theory_config th_cf, const FSC fsc);
  void buildPdf() override;
  void initObservables() override;
  void initParameters() override;
  void initRelations() override;
  void setCorrelations(TString c) override;
  void setObservables(TString c) override;
  void setUncertainties(TString c) override;

 private:
  // Helper functions to avoid boilerplate code
  void add_acpkk(RooArgList* theory, TString name, double avg_time);
  void add_dacp(RooArgList* theory, TString name, double avg_time_kk, double avg_time_pipi);
  const double d0_lifetime = 4.103e-13;
  const FSC fsc;
  const theory_config th_cfg;
};

#endif
