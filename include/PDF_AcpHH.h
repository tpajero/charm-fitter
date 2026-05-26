/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: September 2026
 **/

#pragma once

#include "CharmUtils.h"

#include "PDF_Charm.h"

#include <TString.h>

#include <set>
#include <string>

class RooArgList;

class PDF_AcpHH : public PDF_Charm {
 public:
  /**
   * @param deferInitialise Skip calling initialise() in the constructor. For use by derived classes that override
   *                        PDF_AcpHH's virtuals.
   *                        See https://isocpp.org/wiki/faq/strange-inheritance#calling-virtuals-from-ctors for details.
   */
  PDF_AcpHH(TString obs_id, TString unc_id, hypotheses::dy_fsc, parametrisations::acp, parametrisations::mix,
            int nObs = 2, bool deferInitialise = false);
  void initObservables() override;
  void initRelations() override;
  void setCorrelations(TString c) override;
  void setObservables(TString c) override;
  void setUncertainties(TString c) override;

 protected:
  // Helper functions to avoid boilerplate code
  void add_acpkk(RooArgList* theory, TString name, double avg_time);
  void add_dacp(RooArgList* theory, TString name, double avg_time_kk, double avg_time_pipi);

 private:
  std::set<std::string> getParameterNames() const override;
  const hypotheses::dy_fsc dy_fsc_hypo;
  const parametrisations::acp acp_param;
  const parametrisations::mix mix_param;
};
