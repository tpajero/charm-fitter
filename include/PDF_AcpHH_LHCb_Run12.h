/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: September 2022
 **/

#ifndef PDF_AcpHH_LHCb_Run12_h
#define PDF_AcpHH_LHCb_Run12_h

#include "CharmUtils.h"

#include <PDF_Charm.h>

#include <TString.h>

#include <set>
#include <string>

class RooArgList;

class PDF_AcpHH_LHCb_Run12 : public PDF_Charm {
 public:
  PDF_AcpHH_LHCb_Run12(parametrisations::mix mix_param, parametrisations::dy_fsc dy_fsc_param);
  void initObservables() override;
  void initRelations() override;
  void setCorrelations(TString c) override;
  void setObservables(TString c) override;
  void setUncertainties(TString c) override;

 private:
  // Helper functions to avoid boilerplate code
  void add_acpkk(RooArgList* theory, TString name, double avg_time);
  void add_dacp(RooArgList* theory, TString name, double avg_time_kk, double avg_time_pipi);
  std::set<std::string> getParameterNames() const override;
  const parametrisations::dy_fsc dy_fsc_param;
  const parametrisations::mix mix_param;
};

#endif
