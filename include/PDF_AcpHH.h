/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: August 2026
 **/

#ifndef PDF_AcpHH_h
#define PDF_AcpHH_h

#include "CharmUtils.h"

#include "PDF_Charm.h"

#include <TString.h>

#include <set>
#include <string>

class RooArgList;

class PDF_AcpHH : public PDF_Charm {
 public:
  PDF_AcpHH(TString obs_id, TString unc_id, hypotheses::dy_fsc, parametrisations::acp, parametrisations::mix);
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
  const hypotheses::dy_fsc dy_fsc_hypo;
  const parametrisations::acp acp_param;
  const parametrisations::mix mix_param;
};

#endif
