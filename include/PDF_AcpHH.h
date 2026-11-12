/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: November 2026
 **/

#pragma once

#include "CharmUtils.h"

#include <PDF_Charm.h>

#include <TString.h>

#include <map>
#include <set>
#include <string>

/**
 * Models sets of correlated measurements of ACP(K- K+), ACP(pi- pi+) and DeltaACP.
 *
 * The measurements are stored in a dictionary of datasets, defined in PDF_AcpHH.cpp and identified by a string ID
 * (e.g. "LHCb-R12"). Each dataset lists its observables (with the average decay times entering the theory
 * expressions), their central values and uncertainties, and the statistical and systematic correlation matrices.
 *
 * Central values, uncertainties and correlations can be taken from different datasets, provided that these define
 * the same observables (e.g. to use the true values of the observables with the uncertainties expected from a given
 * dataset). The observables of the PDF are those of the dataset providing the uncertainties.
 */
class PDF_AcpHH final : public PDF_Charm {
 public:
  /// Take central values, uncertainties and correlations from the same dataset.
  PDF_AcpHH(const TString& id, hypotheses::dy_fsc, parametrisations::acp, parametrisations::mix);
  /**
   * @param obs_id ID of the dataset providing the central values, or "truth"/"toy".
   * @param unc_id ID of the dataset providing the observables and their uncertainties.
   * @param cor_id ID of the dataset providing the correlations.
   */
  PDF_AcpHH(const TString& obs_id, const TString& unc_id, const TString& cor_id, hypotheses::dy_fsc,
            parametrisations::acp, parametrisations::mix);

  void initObservables() override;
  void initRelations() override;
  void setCorrelations(TString cor_id) override;
  void setObservables(TString obs_id) override;
  void setUncertainties(TString unc_id) override;

 private:
  struct Dataset;  // defined in PDF_AcpHH.cpp

  /// Dictionary of all the available datasets, indexed by their ID.
  static const std::map<std::string, Dataset>& datasets();
  /// Get the dataset with the given ID.
  static const Dataset& dataset(const TString& id);
  /// Get the dataset with the given ID, checking that it defines the same observables as the PDF.
  const Dataset& compatible_dataset(const TString& id) const;

  std::set<std::string> getParameterNames() const override;

  const Dataset& structure;  ///< Dataset defining the observables of the PDF
  const hypotheses::dy_fsc dy_fsc_hypo;
  const parametrisations::acp acp_param;
  const parametrisations::mix mix_param;
};
