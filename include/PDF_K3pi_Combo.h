/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: November 2026
 **/

#pragma once

#include "CharmUtils.h"

#include <PDF_Charm.h>

#include <TString.h>

#include <set>
#include <string>

/**
 * Charm-only output of the combined analysis of D -> K+ pi- pi- pi+ decays (LHCb + charm-factory inputs).
 *
 * The four observables are
 *
 *   rMix = sqrt(x^2 + y^2),
 *   pMix = atan2(x, y)  (in radians),
 *   dx,
 *   dy,
 *
 * with (dx, dy) following the same conventions of the bin-flip variables (see PDF_BinFlip). The hadronic parameters of
 * D -> K+ pi- pi- pi+ decays are profiled in the combination, so this PDF should not be combined with other
 * measurements of the same decay that are included in the combination.
 *
 * Tags:
 * - LHCb-PAPER-2026-058: built from
 *     - the binned LHCb Run 2 measurement https://arxiv.org/abs/2510.04963
 *     - CLEO inputs from https://arxiv.org/abs/1909.10196
 *     - BESIII inputs from https://arxiv.org/abs/2602.13002
 */
class PDF_K3pi_Combo : public PDF_Charm {
 public:
  PDF_K3pi_Combo(TString measurement_id, parametrisations::mix mix_param);
  void initObservables() override;
  void initRelations() override;
  void setCorrelations(TString measurement_id) override;
  void setObservables(TString measurement_id) override;
  void setUncertainties(TString measurement_id) override;

 private:
  std::set<std::string> getParameterNames() const override;
  const parametrisations::mix mix_param;
};
