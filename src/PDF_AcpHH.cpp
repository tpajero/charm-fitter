/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: November 2026
 **/

#include <PDF_AcpHH.h>

#include <CharmUtils.h>

#include <Utils.h>

#include <RooArgList.h>
#include <RooFormulaVar.h>
#include <RooRealVar.h>

#include <TString.h>

#include <algorithm>
#include <format>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
  enum class ObsType {
    acp_kk,  ///< ACP(K- K+); one average decay time, {t(KK)}
    acp_pp,  ///< ACP(pi- pi+); one average decay time, {t(pipi)}
    dacp,    ///< DeltaACP = ACP(K- K+) - ACP(pi- pi+); two average decay times, {t(KK), t(pipi)}
  };

  struct Measurement {
    ObsType type;
    std::string name;  ///< Base name; theory and observable variables are called `<name>_th` and `<name>_obs`
    std::string title;
    std::vector<double> avg_times;  ///< Average decay times, e.g. `{431.0 * units::fs}`; see ObsType
    double value;
    double stat;
    double syst;
  };
}  // namespace

struct PDF_AcpHH::Dataset {
  std::string source;  ///< Reference for central values, uncertainties and correlations
  std::vector<Measurement> measurements;
  std::vector<double> cor_stat;  ///< Upper triangle (including the diagonal) of the stat. correlation matrix, by row
  std::vector<double> cor_syst;  ///< Upper triangle (including the diagonal) of the syst. correlation matrix, by row
};

const std::map<std::string, PDF_AcpHH::Dataset>& PDF_AcpHH::datasets() {
  // Built on first use (rather than as a global variable) to avoid depending on the initialisation order of globals
  static const std::map<std::string, Dataset> datasets = [] {
    std::map<std::string, Dataset> d;

    d["LHCb-R12"] = {
        // Values, uncertainties, correlations and average decay times from https://cds.cern.ch/record/2799916/
        .source = "https://cds.cern.ch/record/2799916/",
        .measurements =
            {
                // clang-format off
                // type            name                         title                                                     avg. times                                  value      stat      syst
                {ObsType::acp_kk,  "acp_d0_to_kk_run1_mu",      "#it{A_{CP}(K^{+}K^{#minus})} Run 1 #it{#mu}",            {431.0 * units::fs},                      -6.0e-4,  15.0e-4,   9.7e-4},
                {ObsType::acp_kk,  "acp_d0_to_kk_run1_prompt",  "#it{A_{CP}(K^{+}K^{#minus})} Run 1 #it{#pi}",            {918.0 * units::fs},                      14.0e-4,  15.0e-4,  10.0e-4},
                {ObsType::acp_kk,  "acp_d0_to_kk_run2_cdp",     "#it{A_{CP}(K^{+}K^{#minus})} Run 2 #it{C_{D^{+}}}",      {731.5 * units::fs},                      13.6e-4,   8.8e-4,   1.6e-4},
                {ObsType::acp_kk,  "acp_d0_to_kk_run2_cds",     "#it{A_{CP}(K^{+}K^{#minus})} Run 2 #it{C_{D^{+}_{s}}}",  {686.8 * units::fs},                       2.8e-4,   6.7e-4,   2.0e-4},
                {ObsType::dacp,    "dacp_run1_mu",              "#it{#DeltaA_{CP}} Run 1 #it{#mu}",                       {443.7 * units::fs, 437.9 * units::fs},   14.0e-4,  16.0e-4,   8.0e-4},
                {ObsType::dacp,    "dacp_run1_prompt",          "#it{#DeltaA_{CP}} Run 1 #it{#pi}",                       {882.7 * units::fs, 835.4 * units::fs},  -10.0e-4,   8.0e-4,   3.0e-4},
                {ObsType::dacp,    "dacp_run2_mu",              "#it{#DeltaA_{CP}} Run 2 #it{#mu}",                       {491.8 * units::fs, 493.1 * units::fs},   -9.0e-4,   8.0e-4,   5.0e-4},
                {ObsType::dacp,    "dacp_run2_prompt",          "#it{#DeltaA_{CP}} Run 2 #it{#pi}",                       {694.6 * units::fs, 640.7 * units::fs},  -18.2e-4,   3.2e-4,   0.9e-4},
                // clang-format on
            },
        .cor_stat =
            {
                // clang-format off
                1., 0.36, 0., 0.,   0.23, 0.,   0., 0.,    // ACP(KK) Run 1 mu
                    1.,   0., 0.,   0.,   0.24, 0., 0.,    // ACP(KK) Run 1 prompt
                          1., 0.05, 0.,   0.,   0., 0.06,  // ACP(KK) Run 2 CDP
                              1.,   0.,   0.,   0., 0.08,  // ACP(KK) Run 2 CDS
                                    1.,   0.,   0., 0.,    // DeltaACP Run 1 mu
                                          1.,   0., 0.,    // DeltaACP Run 1 prompt
                                                1., 0.,    // DeltaACP Run 2 mu
                                                    1.,    // DeltaACP Run 2 prompt
                // clang-format on
            },
        /**
         * The correlation of the two ACP(KK) Run 1 measurements is calculated from:
         * - Tab. 43 of https://cds.cern.ch/record/1556318;
         * - Tab. 27 of https://cds.cern.ch/record/2126060.
         *
         * The following uncertainties are taken to be 100% correlated,
         *                 ACP(KK) mu  ACP(KK) prompt
         * - Weighting     5.3e-4      6.2e-4
         * - Weighting toy 3.0e-4      5.4e-4
         * - K0 asymmetry  1.4e-4      1.4e-4
         * while the others are assumed to be uncorrelated.
         *
         * The correlation between the two ACP(KK) measurements is then
         * (1.4**2 + 5.3 * 6.2 + 3.0 * 5.4) / (9.7 * 10.0) = 0.53
         */
        .cor_syst =
            {
                // clang-format off
                1., 0.53, 0., 0.,   0.40, 0., 0., 0.,  // ACP(KK) Run 1 mu
                    1.  , 0., 0.,   0.,   0., 0., 0.,  // ACP(KK) Run 1 prompt
                          1., 0.28, 0.,   0., 0., 0.,  // ACP(KK) Run 2 CDP
                              1.,   0.,   0., 0., 0.,  // ACP(KK) Run 2 CDS
                                    1.,   0., 0., 0.,  // DeltaACP Run 1 mu
                                          1., 0., 0.,  // DeltaACP Run 1 prompt
                                              1., 0.,  // DeltaACP Run 2 mu
                                                  1.,  // DeltaACP Run 2 prompt
                // clang-format on
            },
    };

    // Check the internal consistency of all datasets
    for (const auto& [id, dataset] : d) {
      for (const auto& m : dataset.measurements) {
        const std::size_t n_times = m.type == ObsType::dacp ? 2 : 1;
        if (m.avg_times.size() != n_times) {
          throw std::logic_error(
              std::format("PDF_AcpHH ERROR in dataset \"{}\": measurement {} needs {} average decay time(s), {} given",
                          id, m.name, n_times, m.avg_times.size()));
        }
      }
      const std::size_t n_obs = dataset.measurements.size();
      const std::size_t n_cor = n_obs * (n_obs + 1) / 2;
      if (dataset.cor_stat.size() != n_cor || dataset.cor_syst.size() != n_cor) {
        throw std::logic_error(
            std::format("PDF_AcpHH ERROR in dataset \"{}\": correlation matrices need {} elements for {} "
                        "measurements, got {} (stat) and {} (syst)",
                        id, n_cor, n_obs, dataset.cor_stat.size(), dataset.cor_syst.size()));
      }
    }
    return d;
  }();
  return datasets;
}

const PDF_AcpHH::Dataset& PDF_AcpHH::dataset(const TString& id) {
  const auto it = datasets().find(id.Data());
  if (it == datasets().end())
    throw std::runtime_error(std::format("PDF_AcpHH ERROR dataset \"{}\" not found", id.Data()));
  return it->second;
}

const PDF_AcpHH::Dataset& PDF_AcpHH::compatible_dataset(const TString& id) const {
  const Dataset& d = dataset(id);
  if (!std::ranges::equal(d.measurements, structure.measurements, {}, &Measurement::name, &Measurement::name)) {
    throw std::runtime_error(std::format(
        "PDF_AcpHH ERROR dataset \"{}\" does not define the same observables as PDF {}", id.Data(), name.Data()));
  }
  return d;
}

PDF_AcpHH::PDF_AcpHH(const TString& id, hypotheses::dy_fsc dy_fsc_hypo, parametrisations::acp acp_param,
                     parametrisations::mix mix_param)
    : PDF_AcpHH{id, id, id, dy_fsc_hypo, acp_param, mix_param} {}

PDF_AcpHH::PDF_AcpHH(const TString& obs_id, const TString& unc_id, const TString& cor_id,
                     hypotheses::dy_fsc dy_fsc_hypo, parametrisations::acp acp_param, parametrisations::mix mix_param)
    : PDF_Charm{static_cast<int>(dataset(unc_id).measurements.size())}, structure{dataset(unc_id)},
      dy_fsc_hypo{dy_fsc_hypo}, acp_param{acp_param}, mix_param{mix_param} {
  name = "Charm_AcpHH_" + obs_id;
  if (unc_id != obs_id || cor_id != obs_id) name += "_" + unc_id + "_" + cor_id;
  initialise(obs_id, unc_id, cor_id);
}

std::set<std::string> PDF_AcpHH::getParameterNames() const {
  std::set<std::string> names = utils::acp_hh_parameters_names(acp_param, {"KK", "PP"});
  names.merge(utils::dy_hh_parameters_names(dy_fsc_hypo, acp_param, mix_param, {"KK", "PP"}));
  return names;
}

void PDF_AcpHH::initRelations() {
  const auto acp = [&](const std::string& fs) { return utils::acp_expression(acp_param, fs); };
  const auto dy = [&](const std::string& fs) { return utils::dy_hh_expression(dy_fsc_hypo, acp_param, mix_param, fs); };

  theory = new RooArgList("theory");
  for (const auto& m : structure.measurements) {
    std::string formula;
    switch (m.type) {
    case ObsType::acp_kk:
      formula = std::format("{} + {:.5e} * ({})", acp("KK"), m.avg_times[0] / constants::d0_lifetime, dy("KK"));
      break;
    case ObsType::acp_pp:
      formula = std::format("{} + {:.5e} * ({})", acp("PP"), m.avg_times[0] / constants::d0_lifetime, dy("PP"));
      break;
    case ObsType::dacp:
      formula =
          std::format("{} + {:.5e} * ({}) - ({}) - {:.5e} * ({})", acp("KK"), m.avg_times[0] / constants::d0_lifetime,
                      dy("KK"), acp("PP"), m.avg_times[1] / constants::d0_lifetime, dy("PP"));
      break;
    }
    theory->add(*(Utils::makeTheoryVar(m.name + "_th", formula, parameters)));
  }
}

void PDF_AcpHH::initObservables() {
  observables = new RooArgList("observables");
  for (const auto& m : structure.measurements)
    observables->add(*(new RooRealVar((m.name + "_obs").c_str(), m.title.c_str(), 0., -1., 1.)));
}

void PDF_AcpHH::setObservables(const TString obs_id) {
  if (obs_id.EqualTo("truth"))
    setObservablesTruth();
  else if (obs_id.EqualTo("toy"))
    setObservablesToy();
  else {
    const Dataset& d = compatible_dataset(obs_id);
    obsValSource = d.source;
    for (const auto& m : d.measurements) setObservable(m.name + "_obs", m.value);
  }
}

void PDF_AcpHH::setUncertainties(const TString unc_id) {
  const Dataset& d = compatible_dataset(unc_id);
  obsErrSource = d.source;
  StatErr.clear();
  SystErr.clear();
  for (const auto& m : d.measurements) {
    StatErr.push_back(m.stat);
    SystErr.push_back(m.syst);
  }
}

void PDF_AcpHH::setCorrelations(const TString cor_id) {
  const Dataset& d = compatible_dataset(cor_id);
  corSource = d.source;
  corStatMatrix = Utils::buildCorMatrix(nObs, d.cor_stat);
  corSystMatrix = Utils::buildCorMatrix(nObs, d.cor_syst);
}
