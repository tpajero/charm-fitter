// Core
#include <Combiner.h>
#include <GammaComboEngine.h>

// CharmFitter
#include <CharmUtils.h>
#include <PDF_AcpHH_LHCb_Run12.h>
#include <PDF_BES_CLEO_K3pi_Kpipi0.h>
#include <PDF_BES_Kpi.h>
#include <PDF_BES_Kpi_pipipi0.h>
#include <PDF_BinFlip.h>
#include <PDF_CLEO_Kpi.h>
#include <PDF_DY.h>
#include <PDF_DY_RS.h>
#include <PDF_DY_pipipi0.h>
#include <PDF_Fp_pipipi0.h>
#include <PDF_K3pi.h>
#include <PDF_Kpipi0.h>
#include <PDF_RM.h>
#include <PDF_WS.h>
#include <PDF_WS_NoCPV.h>
#include <PDF_XY.h>
#include <PDF_XY_QoP_PHI.h>
#include <PDF_scan_DY_RS.h>
#include <PDF_yCP.h>
#include <PDF_yCP_minus_yCP_KP.h>
#include <PDF_yCP_minus_yCP_RS.h>
#include <PDF_yCP_plus_yCP_RS.h>

#include <format>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
  std::vector<int> lhcb_pdfs(const std::string period, const hypotheses::dy_fsc dy_fsc_hypo) {
    std::vector<int> pdfs;
    if (period == "R1") {  // LHCb Run 1
      pdfs = {
          140,  // yCP - yCP(RS)
          190,  // DeltaY(h- h+)
          220,  // DeltaACP, ACP(K- K+)
          281,  // D0 -> K+ pi-          Prompt
          282,  // D0 -> K+ pi-          DT
          322,  // D0 -> KS pi- pi+      Prompt 2011 (AmAna)
          323,  // D0 -> KS pi- pi+      Binflip
          371,  // D0 -> K+ pi- pi- pi+
      };
    } else if (period == "R12") {  // LHCb Run 1 + Run 2
      pdfs = {
          140,  // yCP - yCP(RS)         Run 1
          141,  // yCP - yCP(RS)         Run 2
          191,  // DeltaY(h- h+)         Run 1+2
          222,  // DeltaACP, ACP(K- K+)  Run 1+2
          286,  // D0 -> K+ pi-          Run 1+2 prompt (App. B)
          288,  // D0 -> K+ pi-          Run 1+2 DT
          // 310,  // DeltaY(K- pi+)        Run 2
          322,  // D0 -> KS pi- pi+      Prompt 2011 (AmAna)
          323,  // D0 -> KS pi- pi+      Binflip Run 1
          326,  // D0 -> KS pi- pi+      Binflip Run 2
          371,  // D0 -> K+ pi- pi- pi+  Run 1
                // TODO add D0 -> K+ pi- pi- pi+ Run 2
      };
      if (dy_fsc_hypo == hypotheses::dy_fsc::none) {
        pdfs.push_back(351);  // DY_pipipi0 Run 2
      }
    } else {
      throw std::runtime_error(
          std::format("lhcb_pdfs ERROR The list of the LHCb results from the period '{}' is not supported", period));
    }
    return pdfs;
  }

  using hypotheses::dy_fsc;
  using parametrisations::acp;
  using parametrisations::mix;

  /// Auxiliary inputs from Charm Factories.
  const std::vector<int> cf_pdfs = {
      1,   // Delta_Kpi            Cleo-c
      21,  // Delta_Kpi F+_pipipi0 BESIII 3+7 fb
      40,  // F+_pipipi0           Cleo-c
      41,  // F+_pipipi0           BESIII
      60,  // K3pi Kpipi0          Cleo-c + BESIII
  };

  struct ParsedArgs {
    dy_fsc dy_fsc_hypo;
    acp acp_param;
    mix mix_param;
    bool dcs_cpv;
    bool help;
    std::vector<char*> combiner_argv;
  };

  const std::set<dy_fsc> supported_dyfsc{dy_fsc::none, dy_fsc::partial, dy_fsc::full};
  const std::set<acp> supported_acp{acp::acp_dy, acp::acp_cot, acp::r_delta};
  const std::set<mix> supported_mix{mix::pheno, mix::theo};

  /// Join the id strings of a set of enum values with "|", e.g. "no|partial|full" for dy_fsc.
  template <typename Enum>
  std::string join_ids(const std::set<Enum>& values) {
    std::string result;
    for (auto it = values.begin(); it != values.end(); ++it) {
      if (it != values.begin()) result += "|";
      result += utils::get_id(*it);
    }
    return result;
  }

  /// Print the command-line options that are specific to charm-combo.
  void print_charm_help() {
    std::cout << "charm-combo: Global fit of charm-quark mixing and CP-violation measurements\n\n"
              << "CharmFitter-specific options (in addition to all GammaCombo options listed below):\n\n"
              << "  --dy-fsc [" << join_ids(supported_dyfsc) << "]  (default: " << utils::get_id(dy_fsc::none) << ")\n"
              << "      Choose the final-state correction to be applied to DeltaY(h- h+).\n"
              << "      Changes the combiner name to <combiner-name>_<dy_fsc>-dyfsc.\n\n"
              << "  --acp [" << join_ids(supported_acp) << "]  (default: " << utils::get_id(acp::acp_dy) << ")\n"
              << "      Choose the parametrisation for aCP(h- h+).\n"
              << "      Changes the combiner name to <combiner-name>_<acp>.\n\n"
              << "  --mix [" << join_ids(supported_mix) << "]  (default: " << utils::get_id(mix::theo) << ")\n"
              << "      Choose the mixing parametrisation.\n"
              << "      Changes the combiner name to <combiner-name>_<mix>.\n\n"
              << "  --dcs-cpv\n"
              << "      Allow for CP violation in doubly Cabibbo-suppressed D0 -> K+ pi- decays.\n"
              << "      Changes the combiner name to <combiner-name>_dcs-cpv. If not set, `--fix Acp_KP=0` is\n"
              << "      automatically passed to GammaComboEngine.\n\n"
              << "-------------------------------------------------------------------------------------------"
              << std::endl;
  }

  /**
   * Parse a `--<flag> <value>` pair out of argv, matching `value` against the id string of one of `supported`, and
   * write the result to `out`. Marks both tokens for removal from the argv that gets forwarded to GammaComboEngine.
   */
  template <typename Enum>
  Enum parse_enum_option(const int argc, char* argv[], const int i, const char* flag, const std::set<Enum>& supported,
                         std::set<int>& to_remove) {
    if (i == argc - 1) throw std::runtime_error(std::format("main ERROR Option \"{}\" requires an argument", flag));
    to_remove.insert({i, i + 1});
    for (auto val : supported) {
      if (!strcmp(argv[i + 1], utils::get_id(val).c_str())) { return val; }
    }
    throw std::runtime_error(std::format("main ERROR Option \"{}\" is not supported by \"{}\"", argv[i + 1], flag));
  }

  /**
   * Parse the command-line arguments of the main executable that are specific to charm-fitter and will not be parsed
   * by GammaComboEngine.
   */
  ParsedArgs parse_args(int argc, char* argv[]) {
    dy_fsc dy_fsc_hypo = dy_fsc::none;
    acp acp_param = acp::acp_dy;
    mix mix_param = mix::theo;
    bool dcs_cpv = false;
    bool help = false;

    std::set<int> to_remove;
    for (int i = 1; i < argc; ++i) {
      if (!strcmp(argv[i], "--dy-fsc")) {
        dy_fsc_hypo = parse_enum_option(argc, argv, i, "--dy-fsc", supported_dyfsc, to_remove);
      } else if (!strcmp(argv[i], "--acp")) {
        acp_param = parse_enum_option(argc, argv, i, "--acp", supported_acp, to_remove);
      } else if (!strcmp(argv[i], "--mix")) {
        mix_param = parse_enum_option(argc, argv, i, "--mix", supported_mix, to_remove);
      } else if (!strcmp(argv[i], "--dcs-cpv")) {
        dcs_cpv = true;
        to_remove.insert(i);
      } else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
        help = true;
      }
    }
    if (!help) utils::check_compatibility(dy_fsc_hypo, acp_param);

    // Prepare the arguments to pass to GammaComboEngine
    std::vector<const char*> extra_args;
    if (!dcs_cpv) { extra_args.insert(extra_args.end(), {"--fix", "Acp_KP=0"}); }

    std::vector<char*> combiner_argv = {};
    for (int i = 0; i < argc; ++i) {
      if (!to_remove.contains(i)) combiner_argv.emplace_back(argv[i]);
    }
    for (auto arg : extra_args) combiner_argv.emplace_back(const_cast<char*>(arg));

    return {dy_fsc_hypo, acp_param, mix_param, dcs_cpv, help, std::move(combiner_argv)};
  }
}  // namespace

/**
 * Run the combination of measurements of charm mixing and CP violation.
 *
 * Accepts all command-line arguments of GammaComboEngine, and in addition:
 *   --dy-fsc [no|partial|full] Choose the final-state correction to be applied to DeltaY(h- h+).
 *       Changes the combiner name to <combiner_name>_<dy_fsc>-dyfsc, where <dy_fsc> is one of `no`, `partial`,
 *       or `full`.
 *   --acp [acp-dy|acp-cot|r-delta] Choose the parametrisation for aCP(h- h+) (default: `acp-dy`).
 *       Changes the combiner name to <combiner_name>_<acp>, where <acp> is one of `acp-dy`, `acp-cot` or `r-delta`.
 *   --mix [pheno|theo] Choose the mixing parametrisation (default: `theo`).
 *       Changes the combiner name to <combiner_name>_<mix>, where <mix> is one of `pheno` or `theo`.
 *   --dcs-cpv Do allow for CP violation in doubly Cabibbo-suppressed D0 -> K+ pi- decays.
 *       It changes the combiner name to `<combiner_name>_dcs-cpv`. If not set, the argument `--fix Acp_KP=0` is
 *       automatically passed to GammaComboEngine.
 *
 * Passing "-h" or "--help" prints the options above, followed by the full list of GammaCombo options.
 */
int main(int argc, char* argv[]) {
  auto parsed_args = parse_args(argc, argv);
  const auto dy_fsc_hypo = parsed_args.dy_fsc_hypo;
  const auto acp_param = parsed_args.acp_param;
  const auto mix_param = parsed_args.mix_param;
  const bool dcs_cpv = parsed_args.dcs_cpv;
  std::vector<char*> combiner_argv = std::move(parsed_args.combiner_argv);

  if (parsed_args.help) {
    print_charm_help();
  } else {
    std::cout << "INFO The combination will be run with the following configuration:\n"
              << "     - DeltaY(h- h+) final-state correction:          " << dy_fsc_hypo << "\n"
              << "     - aCP(h- h+) asymmetry parametrisation:          " << acp_param << "\n"
              << "     - Mixing parametrisation:                        " << mix_param << "\n"
              << "     - Allow for CP violation in D0 -> K+ pi- decays: " << std::boolalpha << dcs_cpv
              << std::noboolalpha << std::endl;
  }

  std::string combiner_name = std::format("charm-combo_{}-dyfsc_{}_{}", utils::get_id(dy_fsc_hypo),
                                          utils::get_id(acp_param), utils::get_id(mix_param));
  if (dcs_cpv) { combiner_name += "_dcs-cpv"; }
  GammaComboEngine gc(combiner_name, combiner_argv.size(), &combiner_argv[0]);

  ///////////////////////////////////////////////////
  //
  // define PDFs
  //
  ///////////////////////////////////////////////////

  using hypotheses::dy_fsc;

  // clang-format off

  // --- Inputs from Charm Factories, grouped by final state ---

  gc.addPdf(1, new PDF_CLEO_Kpi("CLEO-c", mix_param),                            "Delta_Kpi    CLEO-c                         ");
  gc.addPdf(2, new PDF_BES_Kpi(mix_param),                                       "Delta_Kpi    BESIII   3fb      [A_kpi only] ");  // Superseded by 20, 21

  gc.addSubsetPdf(20, new PDF_BES_Kpi_pipipi0("3fb", mix_param), 0, 1, 2, 3,     "Kpi+pipipi0  BESIII   3fb                   ");  // Supersedes 2
  gc.addPdf(21, new PDF_BES_Kpi_pipipi0("3+7fb", mix_param),                     "Kpi+pipipi0  BESIII   3+7fb                 ");  // Supersedes 2, 20

  gc.addPdf(40, new PDF_Fp_pipipi0("CLEO-c"),                                    "Fpipipi0     CLEO-c                         ");
  gc.addPdf(41, new PDF_Fp_pipipi0("BESIII-8fb"),                                "Fpipipi0     BESIII   8fb                   ");  // Assumed uncorrelated from 20, 21 due to much larger sample (8 vs. 3 fb)

  gc.addPdf(60, new PDF_BES_CLEO_K3pi_Kpipi0("BESIII--CLEO-c"),                  "K3pi-Kpipi0  BESIII + CLEO-c                ");

  // --- Charm measurements ---

  // SL Mixing
  gc.addPdf(100, new PDF_RM("HFLAV-2008", mix_param),                            "R_M          HFLAV    2008                  ");

  // yCP
  gc.addPdf(110, new PDF_yCP("WA-2015", mix_param),                              "yCP          WA       2015                  ");

  gc.addPdf(120, new PDF_yCP_plus_yCP_RS("Belle", mix_param),                    "yCP+yCP(RS)  Belle                          ");

  gc.addPdf(130, new PDF_yCP_minus_yCP_KP("WA-2015", mix_param),                 "yCP-yCP(KP)  WA       2015                  ");

  gc.addPdf(140, new PDF_yCP_minus_yCP_RS("LHCb-R1", mix_param),                 "yCP-yCP(RS)  LHCb     Run 1    [B -> D* mu] ");  // Included in WA >= 2018
  gc.addPdf(141, new PDF_yCP_minus_yCP_RS("LHCb-R2", mix_param),                 "yCP-yCP(RS)  LHCb     Run 2    [D* -> D0 pi]");  // Included in WA >= 2022

  gc.addPdf(160, new PDF_yCP_minus_yCP_RS("WA-2015", mix_param),                 "yCP-yCP(RS)  WA       2015                  ");
  gc.addPdf(161, new PDF_yCP_minus_yCP_RS("WA-2018", mix_param),                 "yCP-yCP(RS)  WA       2018                  ");
  gc.addPdf(162, new PDF_yCP_minus_yCP_RS("WA-2022", mix_param),                 "yCP-yCP(RS)  WA       2022                  ");

  // D0 -> h+ h-
  if (dy_fsc_hypo == dy_fsc::none) {
    gc.addPdf(170, new PDF_DY("Belle+BaBar", dy_fsc_hypo, acp_param, mix_param), "DY           B-Factories                    ");  // Included in WA >= 2019
  }
  gc.addPdf(171, new PDF_DY("NonLHCb-2015", dy_fsc_hypo, acp_param, mix_param),  "DY           NonLHCb  2015                  ");  // Superset of 170; included in WA >= 2019

  gc.addPdf(190, new PDF_DY("LHCb-R1",  dy_fsc_hypo, acp_param, mix_param),      "DY           LHCb     Run 1                 ");  // Included in WA >= 2019
  gc.addPdf(191, new PDF_DY("LHCb-R12", dy_fsc_hypo, acp_param, mix_param),      "DY           LHCb     Run 1+2               ");  // Included in WA >= 2021

  gc.addPdf(200, new PDF_DY("WA-2019", dy_fsc_hypo, acp_param, mix_param),       "DY           WA       2019                  ");
  gc.addPdf(201, new PDF_DY("WA-2020", dy_fsc_hypo, acp_param, mix_param),       "DY           WA       2020                  ");
  gc.addPdf(202, new PDF_DY("WA-2021", dy_fsc_hypo, acp_param, mix_param),       "DY           WA       2021                  ");

  gc.addSubsetPdf(220, new PDF_AcpHH_LHCb_Run12(dy_fsc_hypo, acp_param, mix_param), 0, 1, 4, 5,
                                                                                 "ACP(KK/PP)   LHCb     Run1                  ");  // Superseded by 221, 222
  gc.addSubsetPdf(221, new PDF_AcpHH_LHCb_Run12(dy_fsc_hypo, acp_param, mix_param), 0, 1, 4, 5, 6, 7,
                                                                                 "ACP(KK/PP)   LHCb     DeltaACP R12, ACPKK R1");  // Supersedes 220, superseded by 222
  gc.addPdf(222, new PDF_AcpHH_LHCb_Run12(dy_fsc_hypo, acp_param, mix_param),    "ACP(KK/PP)   LHCb     Run1+2                ");  // Supersedes 220, 221

  // D0 -> K+ pi-

  gc.addPdf(250, new PDF_WS_NoCPV("CDF", mix_param),                             "WS/RS        CDF      No CPV                ");
  gc.addPdf(251, new PDF_WS_NoCPV("BaBar", mix_param),                           "WS/RS        BaBar    No CPV                ");  // Alternative to 260
  gc.addPdf(252, new PDF_WS_NoCPV("Belle", mix_param),                           "WS/RS        Belle    No CPV                ");  // Incompatible with 261 and based on a larger sample

  gc.addPdf(260, new PDF_WS("BaBar", mix_param),                                 "WS/RS        BaBar                          ");  // Alternative to 251
  gc.addPdf(261, new PDF_WS("Belle", mix_param),                                 "WS/RS        Belle                          ");  // Incompatible with 252 and based on a smaller sample

  gc.addPdf(280, new PDF_WS("LHCb-R1-prompt-biased", mix_param),                 "WS/RS        LHCb     Run 1    [D* -> D0 pi, biased]");  // Neglects ghost pions
  gc.addPdf(281, new PDF_WS("LHCb-R1-prompt", mix_param),                        "WS/RS        LHCb     Run 1    [D* -> D0 pi]");  // Supersedes 280
  gc.addPdf(282, new PDF_WS("LHCb-R1-DT", mix_param),                            "WS/RS        LHCb     Run 1    [B -> D* mu] ");
  gc.addPdf(283, new PDF_WS("LHCb-R1-biased", mix_param),                        "WS/RS        LHCb     Run 1    [biased]     ");  // Prompt measurements neglects ghost pions
  gc.addPdf(284, new PDF_WS("LHCb-2011-2016-prompt", mix_param),                 "WS/RS        LHCb     2011-6   [D* -> D0 pi]");  // Supersedes 281 and, partially, 283
  gc.addPdf(285, new PDF_WS("LHCb-R12-prompt-sec9", mix_param),                  "WS/RS        LHCb     Run 1+2  [D* -> D0 pi, sec 9]");  // Supersedes 284
  gc.addPdf(286, new PDF_WS("LHCb-R12-prompt-appB", mix_param, dy_fsc_hypo, acp_param),
                                                                                 "WS/RS        LHCb     Run 1+2  [D* -> D0 pi, app B]");  // Supersedes 284
  gc.addPdf(287, new PDF_WS("LHCb-R2-DT", mix_param),                            "WS/RS        LHCb     Run 2    [B -> D* mu] ");
  gc.addPdf(288, new PDF_WS("LHCb-R12-DT", mix_param),                           "WS/RS        LHCb     Run 1+2  [B -> D* mu] ");  // Supersedes 282, 287

  // D0 -> K- pi+

  gc.addPdf(310, new PDF_DY_RS("LHCb-R2-prompt", mix_param),                     "DY(RS)       LHCb     Run 2    [D* -> D0 pi] ");

  // D0 -> KS pi+ pi-

  gc.addPdf(320, new PDF_XY("BaBar-KShh", mix_param),                            "KS h+ h-     BaBar    2010                  ");
  gc.addPdf(321, new PDF_XY_QoP_PHI("Belle", mix_param),                         "KS h+ h-     Belle    2014                  ");  // Superseded by 327
  gc.addPdf(322, new PDF_XY("LHCb-KSpipi-2011-prompt", mix_param),               "KS pi+ pi-   LHCb     2011     [D* -> D0 pi]");
  gc.addPdf(323, new PDF_BinFlip("LHCb-R1", mix_param),                          "Bin-flip     LHCb     Run 1                 ");
  gc.addPdf(324, new PDF_BinFlip("LHCb-R2-prompt", mix_param),                   "Bin-flip     LHCb     Run 2    [D* -> D0 pi]");
  gc.addPdf(325, new PDF_BinFlip("LHCb-R2-SL", mix_param),                       "Bin-flip     LHCb     Run 2    [B -> D0 mu] ");
  gc.addPdf(326, new PDF_BinFlip("LHCb-R2", mix_param),                          "Bin-flip     LHCb     Run 2                 ");
  gc.addPdf(327, new PDF_XY("Belle-Belle2", mix_param),                          "KS pi+ pi-   Belle+Belle2 (951+408 fb-1)    ");  // Supersedes 321

  // D0 -> pi+ pi- pi0

  gc.addPdf(350, new PDF_XY("BaBar-pipipi0", mix_param),                         "pi+ pi- pi0  BaBar                          ");
  if (dy_fsc_hypo == dy_fsc::none) {
    gc.addPdf(351, new PDF_DY_pipipi0("LHCb-R2", mix_param),                     "DY(pipipi0)  LHCb     Run 2                 ");
  }

  // D0 -> K+ pi- pi0

  gc.addPdf(360, new PDF_Kpipi0("BaBar", mix_param, false),                      "Kpipi0       BaBar                          ");

  // D0 -> K+ pi- pi- pi+

  gc.addPdf(370, new PDF_RM("LHCb-K3pi-R1", mix_param),                          "R_M  K3pi    LHCb     R1                    ");  // Subset of 371
  gc.addPdf(371, new PDF_K3pi("LHCb-R1", mix_param),                             "K3pi         LHCb     R1                    ");  // Superset of 370

  // --- Nuisance parameters ---

  gc.addPdf(390, new PDF_scan_DY_RS(mix_param),                                      "ScanDYRS     This is just a nuisance parameter");

  // clang-format on

  ///////////////////////////////////////////////////
  //
  // Define combinations
  //
  ///////////////////////////////////////////////////

  gc.newCombiner(0, "empty", "empty");

  // World averages ----------------------------------------------------------------------------------------------------

  /*
  // New LHCb measurement of the WS/RS ratio in D0 -> K+ pi- decays with LHCb 2011-2016 prompt data
  gc.cloneCombiner(, , "WA-2018-02", "World average (Feb 2018)");
  gc.getCombiner()->delPdf(gc[283]);  // LHCb Run 1 prompt + DT
  gc.getCombiner()->addPdf(gc[282]);  // LHCb Run 1 DT
  gc.getCombiner()->addPdf(gc[284]);  // LHCb 2011-2016 prompt

  // New LHCb measurement of yCP(h- h+) - yCP(RS) with Run 1 SL data (CKM)
  gc.cloneCombiner(, , "WA-2018-09", "World average (Sep 2018)");
  gc.getCombiner()->delPdf(gc[160]);
  gc.getCombiner()->addPdf(gc[161]);

  // New LHCb bin-flip measurement of D0 -> KS pi- pi+ with Run 1 data (La Thuile)
  gc.cloneCombiner(, , "WA-2019-03", "World average (Mar 2019)");
  gc.getCombiner()->addPdf(gc[323]);
  */

  gc.newCombiner(1, "base", "base",
                 {1, 2, 100, 110, 120, 130, 161, 200, 250, 252, 260, 282, 284, 320, 321, 323, 350, 360, 371});

  // New LHCb measurement of DeltaY(h- h+) with Run 2 SL data (Brookhaven Forum)
  gc.cloneCombiner(2, 1, "WA-2019-09", "World average (Sep 2019)");
  gc.getCombiner(2)->delPdf(gc[200]);
  gc.getCombiner(2)->addPdf(gc[201]);

  // New LHCb measurement of DeltaY(h- h+) with Run 2 prompt data (LHC seminar)
  gc.cloneCombiner(3, 2, "WA-2021-02", "World average (Feb 2021)");
  gc.getCombiner(3)->delPdf(gc[201]);
  gc.getCombiner(3)->addPdf(gc[202]);

  // New LHCb bin-flip measurement of D0 -> KS pi- pi+ with Run 2 prompt data (Moriond QCD)
  // New BESIII deterrmination of strong parameters for D0 -> K+ pi- pi- pi+ and D0 -> K+ pi- pi0
  gc.cloneCombiner(4, 3, "WA-2021-03", "World average (Mar 2021)");
  gc.getCombiner(4)->addPdf(gc[324]);
  gc.getCombiner(4)->addPdf(gc[60]);  // BESIII + CLEO K3pi, Kpipi0

  // New LHCb bin-flip measurement of D0 -> KS pi- pi+ with Run 2 SL data (PANIC 2021)
  gc.cloneCombiner(5, 4, "WA-2021-09", "World average (Sep 2021)");
  gc.getCombiner(5)->delPdf(gc[324]);  // Bin-flip LHCb Run 2 prompt
  gc.getCombiner(5)->addPdf(gc[326]);  // Bin-flip LHCb Run 2

  // New LHCb measurement of yCP - yCP(RS) with Run 2 data (La Thuile 2022)
  gc.cloneCombiner(6, 5, "WA-2022-03", "World average (Mar 2022)");
  gc.getCombiner(6)->delPdf(gc[161]);
  gc.getCombiner(6)->addPdf(gc[162]);

  // New LHCb measurement of ACP(K- K+) (ICHEP 2022)
  gc.cloneCombiner(7, 6, "WA-2022-07", "World average (Jul 2022)");
  gc.getCombiner(7)->delPdf(gc[221]);  // LHCb DeltaACP Run 1+2, ACP(K- K+) Run 1
  gc.getCombiner(7)->addPdf(gc[222]);  // LHCb DeltaACP Run 1+2, ACP(K- K+) Run 1+2

  // New BESIII measurement of Delta_Kpi and F+(pi+ pi- pi0) with 3 fb
  gc.cloneCombiner(8, 7, "WA-2022-08", "World average (Aug 2022)");
  gc.getCombiner(8)->delPdf(gc[2]);   // BESIII delta_Kpi
  gc.getCombiner(8)->addPdf(gc[20]);  // BESIII delta_Kpi + pipipi0 3fb

  // New LHCb measurements for Moriond EW:
  //   - WS/RS ratio in D0 -> K+ pi- decays with LHCb Run 2 prompt data
  //   - DeltaY(pi- pi+ pi0) with 2012 + Run 2 data
  gc.cloneCombiner(9, 8, "WA-2024-03", "World average (Mar 2024)");
  gc.getCombiner(9)->delPdf(gc[284]);  // LHCb WS/RS in D0 -> K+ pi+ 2011-2016 prompt
  gc.getCombiner(9)->addPdf(gc[286]);  // LHCb WS/RS in D0 -> K+ pi+ Run 1+2   prompt
  if (dy_fsc_hypo == dy_fsc::none) { gc.getCombiner(9)->addPdf(gc[351]); }

  // New measurement of F+(pi+pi-pi0) by BESIII
  gc.cloneCombiner(10, 9, "WA-2024-09", "World average (Sep 2024)");
  gc.getCombiner(10)->addPdf(gc[41]);  // F+_pipipi0 BESIII 8fb

  // - New measurement of the WS/RS ratio in D0 -> K+ pi- decays with LHCb Run 2 DT data (LHCb Implications Workshop)
  // - New Belle + Belle2 bin-flip measurement of D0 -> KS pi- pi+
  gc.cloneCombiner(11, 10, "WA-2024-10", "World average (Oct 2024)");
  gc.getCombiner(11)->delPdf(gc[282]);  // WS/RS in D0 -> Kpi, LHCb Run 1   DT
  gc.getCombiner(11)->addPdf(gc[288]);  // WS/RS in D0 -> Kpi, LHCb Run 1+2 DT
  gc.getCombiner(11)->delPdf(gc[321]);  // D0 -> KS h- h+   from Belle
  gc.getCombiner(11)->addPdf(gc[327]);  // D0 -> KS pi- pi+ from Belle + Belle 2

  // New BESIII measurement of Delta_Kpi with 8 fb
  gc.cloneCombiner(12, 11, "WA-2025-06", "World average (Jun 2025)");
  gc.getCombiner(12)->delPdf(gc[20]);  // BESIII delta_Kpi + pipipi0 3   fb
  gc.getCombiner(12)->addPdf(gc[21]);  // BESIII delta_Kpi + pipipi0 3+7 fb

  // New LHCb measurements of WS/RS in D0 -> K+ pi- pi- pi+ with Run 2 prompt data (CKM 2025) TODO
  // gc.cloneCombiner(13, 12, "WA-2025-09", "World average (Sep 2025)");
  // gc.getCombiner(13)->addPdf(gc[]);

  // LHCb-only averages ------------------------------------------------------------------------------------------------

  // All averages rely on the latest determinations of strong parameters at Charm factories.

  gc.newCombiner(300, "LHCb-R1", "LHCb Run 1");
  for (const auto pdf : lhcb_pdfs("R1", dy_fsc_hypo)) gc.getCombiner(300)->addPdf(gc[pdf]);
  for (const auto pdf : cf_pdfs) gc.getCombiner(300)->addPdf(gc[pdf]);

  gc.newCombiner(301, "LHCb-R12", "LHCb Run 1+2");
  for (const auto pdf : lhcb_pdfs("R12", dy_fsc_hypo)) gc.getCombiner(301)->addPdf(gc[pdf]);
  for (const auto pdf : cf_pdfs) gc.getCombiner(301)->addPdf(gc[pdf]);

  // Impact of LHCb upgrades -------------------------------------------------------------------------------------------

  ///////////////////////////////////////////////////
  //
  // Run
  //
  ///////////////////////////////////////////////////

  gc.run();

  return 0;
}
