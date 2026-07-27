/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: June 2024
 **/

#ifndef PDF_BES_CLEO_K3pi_Kpipi0_h
#define PDF_BES_CLEO_K3pi_Kpipi0_h

#include <PDF_Charm.h>

#include <TString.h>

#include <set>
#include <string>

class PDF_BES_CLEO_K3pi_Kpipi0 : public PDF_Charm {
 public:
  PDF_BES_CLEO_K3pi_Kpipi0(TString measurement_id);
  void buildPdf() override;
  void initObservables() override;
  void initRelations() override;
  void setCorrelations(TString measurement_id) override;
  void setObservables(TString measurement_id) override;
  void setUncertainties(TString measurement_id) override;

 private:
  std::set<std::string> getParameterNames() const override;
};

#endif
