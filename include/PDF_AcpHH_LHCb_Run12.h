/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: September 2022
 **/

#pragma once

#include "CharmUtils.h"

#include <PDF_Charm.h>

#include <TString.h>

#include <set>
#include <string>

class RooArgList;

/**
 * Models the LHCb measurements of DeltaACP and of ACP(K- K+) performed during Run 1 and 2.
 *
 * This requires a specific PDF owing to their non-negligible correlations. On the other hand, measurements of
 * DeltaY(h- h+) are assumed to be uncorrelated (due, in particular, to the kinematic equalisation of these last
 * measurements).
 */
class PDF_AcpHH_LHCb_Run12 : public PDF_Charm {
 public:
  PDF_AcpHH_LHCb_Run12(hypotheses::dy_fsc, parametrisations::acp, parametrisations::mix);
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
