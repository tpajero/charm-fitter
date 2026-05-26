/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: September 2022
 **/

#pragma once

#include "CharmUtils.h"

#include "PDF_AcpHH.h"

#include <TString.h>

#include <set>
#include <string>

/**
 * Models the LHCb measurements of DeltaACP and of ACP(K- K+) performed during Run 1 and 2.
 *
 * This requires a specific PDF owing to their non-negligible correlations. On the other hand, measurements of
 * DeltaY(h- h+) are assumed to be uncorrelated (due, in particular, to the kinematic equalisation performed in these
 * last measurements).
 */
class PDF_AcpHH_LHCb_Run12 : public PDF_AcpHH {
 public:
  PDF_AcpHH_LHCb_Run12(hypotheses::dy_fsc, parametrisations::acp, parametrisations::mix);
  void initObservables() override;
  void initRelations() override;
  void setCorrelations(TString c) override;
  void setObservables(TString c) override;
  void setUncertainties(TString c) override;
};
