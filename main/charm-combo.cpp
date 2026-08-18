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
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
  std::vector<int> get_lhcb_pdfs(const std::string run, const hypotheses::dy_fsc dy_fsc_hypo) {
    std::vector<int> pdfs;
    if (run == "run12") {
      pdfs = {
          3,   // XY               KSpipi Prompt Run 1
          11,  // RM               K3pi Run 1
          21,  // BinFlip          Run 1
          24,  // BinFlip          Run 2
          35,  // WS               DT Run 1
          39,  // WS               Prompt Run 1+2
          64,  // yCP_minus_yCP_RS WA 2022
          72,  // DY average
          90,  // AcpHH
      };
      if (dy_fsc_hypo == hypotheses::dy_fsc::none) {
        pdfs.push_back(85);  // DY_pipipi0 Run 2
      }
    } else {
      throw std::runtime_error(
          std::format("get_lhcb_pdfs ERROR The list of the LHCb results from the period `{}` is not supported", run));
    }
    return pdfs;
  }

  using hypotheses::dy_fsc;
  using parametrisations::acp;
  using parametrisations::mix;

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
              << "     DeltaY(h- h+) final-state correction: " << dy_fsc_hypo << "\n"
              << "     aCP(h- h+) asymmetry parametrisation: " << acp_param << "\n"
              << "     Mixing parametrisation: " << mix_param << "\n"
              << "     Allow for CP violation in DCS D0 -> K+ pi- decays: " << dcs_cpv << std::endl;
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
  using parametrisations::kpi;

  // clang-format off
  gc.addPdf(1, new PDF_XY("BaBar_Kshh", mix_param),                                  "XY KShh      BaBar                          ");
  gc.addPdf(2, new PDF_XY("BaBar_pipipi0", mix_param),                               "XY pipipi0   BaBar                          ");
  gc.addPdf(3, new PDF_XY("LHCb_KSpipi", mix_param),                                 "XY KSpipi    LHCb     2011     [D* -> D0 pi]");
  gc.addPdf(4, new PDF_Kpipi0("BaBar", mix_param, false),                            "Kpipi0       BaBar                          ");
  gc.addPdf(5, new PDF_K3pi("LHCb-run1", mix_param),                                 "K3pi         LHCb     Run 1                 ");
  gc.addPdf(6, new PDF_XY("Belle_Belle2", mix_param),                                "XY KSpipi    Belle+Belle2 (951+408 fb-1)    ");

  gc.addPdf(10, new PDF_RM("HFLAV-2008", mix_param),                                 "R_M          HFLAV    2008                  ");
  gc.addPdf(11, new PDF_RM("LHCb-K3pi-R1", mix_param),                               "K3pi         LHCb     Run 1    [RM only]    ");

  gc.addPdf(20, new PDF_XY_QoP_PHI("Belle", mix_param),                              "KShh         Belle                          ");
  gc.addPdf(21, new PDF_BinFlip("LHCb_Run1", mix_param),                             "Bin-flip     LHCb     Run 1                 ");
  gc.addPdf(22, new PDF_BinFlip("LHCb_Run2_prompt", mix_param),                      "Bin-flip     LHCb     Run 2    [D* -> D0 pi]");
  gc.addPdf(23, new PDF_BinFlip("LHCb_Run2_sl", mix_param),                          "Bin-flip     LHCb     Run 2    [B -> D0 mu] ");
  gc.addPdf(24, new PDF_BinFlip("LHCb_Run2", mix_param),                             "Bin-flip     LHCb     Run 2                 ");

  gc.addPdf(30, new PDF_WS_NoCPV("CDF", mix_param),                                  "WS/RS        CDF                            ");
  gc.addPdf(31, new PDF_WS_NoCPV("BaBar", mix_param),                                "WS/RS        BaBar    no CPV                ");
  gc.addPdf(32, new PDF_WS_NoCPV("Belle", mix_param),                                "WS/RS        Belle    no CPV                ");
  gc.addPdf(33, new PDF_WS("BaBar", mix_param),                                      "WS/RS        BaBar                          ");
  gc.addPdf(34, new PDF_WS("Belle", mix_param),                                      "WS/RS        Belle                          ");
  gc.addPdf(35, new PDF_WS("LHCb_DT_Run1", mix_param),                               "WS/RS        LHCb     Run 1    [B -> D* mu] ");
  gc.addPdf(36, new PDF_WS("LHCb_Run1", mix_param),                                  "WS/RS        LHCb     Run 1                 ");
  gc.addPdf(37, new PDF_WS("LHCb_Prompt_2011_2016", mix_param),                      "WS/RS        LHCb     2011-6   [D* -> D0 pi]");
  gc.addPdf(38, new PDF_WS("LHCb_Prompt_Run12_sec9", mix_param, kpi::ccprime),       "WS/RS        LHCb     Run 1+2  [D* -> D0 pi]");
  gc.addPdf(39, new PDF_WS("LHCb_Prompt_Run12_appB", mix_param, kpi::ccprime, dy_fsc_hypo, acp_param), "WS/RS        LHCb     Run 1+2  [D* -> D0 pi]");
  gc.addPdf(40, new PDF_WS("LHCb_DT_Run2", mix_param),                               "WS/RS        LHCb     Run 2    [B -> D* mu] ");
  gc.addPdf(41, new PDF_WS("LHCb_DT_Run12", mix_param),                              "WS/RS        LHCb     Run 1-2  [B -> D* mu] ");

  gc.addPdf(50, new PDF_CLEO_Kpi("Cleo-c", mix_param),                               "Delta_Kpi    Cleo-c                         ");
  gc.addPdf(51, new PDF_BES_Kpi(mix_param),                                          "Delta_Kpi    BES      3fb      [A_kpi only] ");
  gc.addSubsetPdf(52, new PDF_BES_Kpi_pipipi0("3fb", mix_param), 0, 1, 2, 3,         "Kpi+pipipi0  BES      3fb                   ");
  gc.addPdf(53, new PDF_Fp_pipipi0("Cleo-c"),                                        "Fpipipi0     Cleo-c                         ");
  gc.addPdf(54, new PDF_BES_CLEO_K3pi_Kpipi0("BES3-CLEO"),                           "K3pi-Kpipi0  BES3 + Cleo                    ");
  gc.addPdf(55, new PDF_Fp_pipipi0("BESIII"),                                        "Fpipipi0     BES3                           ");
  gc.addPdf(56, new PDF_BES_Kpi_pipipi0("3+7fb", mix_param),                         "Kpi+pipipi0  BES3     3+7fb                 ");

  gc.addPdf(60, new PDF_yCP("WA-2015", mix_param),                                   "yCP          WA       2015                  ");
  gc.addPdf(61, new PDF_yCP_minus_yCP_RS("WA-2018", mix_param),                      "yCP-yCP(RS)  WA       2018                  ");
  gc.addPdf(62, new PDF_yCP_minus_yCP_KP("WA-2015", mix_param),                      "yCP-yCP(KP)  WA       2015                  ");
  gc.addPdf(63, new PDF_yCP_plus_yCP_RS("Belle", mix_param),                         "yCP+yCP(RS)  Belle    2019                  ");
  gc.addPdf(64, new PDF_yCP_minus_yCP_RS("LHCb-R2", mix_param),                      "yCP-yCP(RS)  LHCb     2022                  ");

  if (dy_fsc_hypo == dy_fsc::none) {
    gc.addPdf(70, new PDF_DY("WA2019", dy_fsc_hypo, acp_param, mix_param),           "DY           WA       2019                  ");
    gc.addPdf(73, new PDF_DY("Belle&BaBar", dy_fsc_hypo, acp_param, mix_param),      "DY           B-factories                    ");
  }
  gc.addPdf(71, new PDF_DY("WA2020", dy_fsc_hypo, acp_param, mix_param),             "DY           WA       2020                  ");
  gc.addPdf(72, new PDF_DY("WA2021", dy_fsc_hypo, acp_param, mix_param),             "DY           WA       2021                  ");

  gc.addPdf(80, new PDF_DY_RS("LHCb2021", mix_param),                                "DY(RS)       LHCb     2021                  ");

  if (dy_fsc_hypo == dy_fsc::none) {
    gc.addPdf(85, new PDF_DY_pipipi0("LHCb-R2", mix_param),                          "DY(pipipi0)  LHCb     Run2                  ");
  }

  gc.addPdf(90, new PDF_AcpHH_LHCb_Run12(dy_fsc_hypo, acp_param, mix_param),                   "ACP(KK/PP)   LHCb     Run1+2                ");
  gc.addSubsetPdf(93, new PDF_AcpHH_LHCb_Run12(dy_fsc_hypo, acp_param, mix_param), 0, 1, 4, 5, "ACP(KK/PP)   LHCb     Run1                  ");

  gc.addPdf(100, new PDF_scan_DY_RS(mix_param),                                      "ScanDYRS     This is just a nuisance parameter");

  gc.addPdf(110, new PDF_yCP("WA-biased-2019", mix_param),                           "yCP          WA       2019     [biased]     ");
  gc.addPdf(111, new PDF_yCP("WA-biased-2022", mix_param),                           "yCP-yCP(RS)  WA       2022     [biased]     ");
  // clang-format on

  ///////////////////////////////////////////////////
  //
  // Define combinations
  //
  ///////////////////////////////////////////////////

  gc.newCombiner(0, "empty", "empty");

  // WA 2020
  gc.newCombiner(1, "WA-2020", "World average (Dec 2020)", 1, 2, 3, 4, 10, 11, 20, 21, 30, 31, 32, 35, 37, 50, 51);
  gc.getCombiner(1)->addPdf(gc[60]);
  gc.getCombiner(1)->addPdf(gc[61]);
  gc.getCombiner(1)->addPdf(gc[62]);
  gc.getCombiner(1)->addPdf(gc[63]);
  gc.getCombiner(1)->addPdf(gc[71]);

  // WA June 2021
  gc.cloneCombiner(20, 1, "WA-2021", "World average (June 2021)");
  gc.getCombiner(20)->addPdf(gc[22]);  // bin-flip run 2
  gc.getCombiner(20)->delPdf(gc[71]);  // DY WA 2020
  gc.getCombiner(20)->addPdf(gc[72]);  // DY WA 2021
  gc.getCombiner(20)->delPdf(gc[11]);  // LHCb K3pi (x2 + y2)/4
  gc.getCombiner(20)->addPdf(gc[5]);   // LHCb K3pi full
  gc.getCombiner(20)->addPdf(gc[54]);  // BES3 + CLEO K3pi, Kpipi0

  // WA after LHCb 2022 yCP measurement
  gc.cloneCombiner(30, 20, "WA-2022-02", "World average (Feb 2022)");
  gc.getCombiner(30)->addPdf(gc[64]);  // yCP LHCb 2022

  // WA after LHCb 2022 yCP measurement - biased
  gc.cloneCombiner(31, 30, "WA-2022-02-biased",
                   "World average (Feb 2022) #minus no #it{y}_{#it{CP}}^{#it{K^{#minus}#pi^{+}}} correction");
  gc.getCombiner(31)->delPdf(gc[60]);
  gc.getCombiner(31)->delPdf(gc[61]);
  gc.getCombiner(31)->delPdf(gc[62]);
  gc.getCombiner(31)->delPdf(gc[63]);
  gc.getCombiner(31)->delPdf(gc[64]);
  gc.getCombiner(31)->addPdf(gc[111]);

  // WA September 2022
  gc.cloneCombiner(40, 30, "WA-2022-09", "World average (Sept 2022)");
  gc.getCombiner(40)->delPdf(gc[51]);  // old BESIII measurement of delta_Kpi
  gc.getCombiner(40)->addPdf(gc[52]);  // new BESIII measurement of delta_Kpi
  gc.getCombiner(40)->addPdf(gc[53]);  // F+_pipipi0
  gc.getCombiner(40)->delPdf(gc[22]);  // bin-flip LHCb Run 2 prompt
  gc.getCombiner(40)->addPdf(gc[24]);  // bin-flip LHCb Run 2
  gc.getCombiner(40)->addPdf(gc[90]);  // ACP(KK) + DeltaACP LHCb Run 1+2

  // WA March 2024 March before WS/RS
  gc.cloneCombiner(49, 40, "WA-2024-02", "World average (Feb 2024)");
  if (dy_fsc_hypo == dy_fsc::none) gc.getCombiner(49)->addPdf(gc[85]);  // DY(pi+ pi- pi0) from LHCb Run 2

  // WA March 2024 - no FSC
  gc.cloneCombiner(50, 40, "WA-2024-03", "World average (March 2024)");
  gc.getCombiner(50)->delPdf(gc[37]);                                   // WS/RS in D0 -> Kpi from LHCb 2011-2016
  gc.getCombiner(50)->addPdf(gc[39]);                                   // WS/RS in D0 -> Kpi from LHCb Run 1+2
  if (dy_fsc_hypo == dy_fsc::none) gc.getCombiner(50)->addPdf(gc[85]);  // DY(pi+ pi- pi0) from LHCb Run 2

  // WA March 2024 with parametrisation of prompt LHCb WS/RS decays from Sec. 9 - no FSC
  gc.cloneCombiner(51, 50, "WA-2024-03-WSsec9", "World average (March 2024, prompt WS/RS from Sec. 9)");
  gc.getCombiner(51)->delPdf(gc[39]);  // WS/RS in D0 -> Kpi from LHCb Run 1+2
  gc.getCombiner(51)->addPdf(gc[38]);  // WS/RS in D0 -> Kpi from LHCb Run 1+2

  // WA Sept 2024 (new BESIII F+(pi+pi-pi0))
  gc.cloneCombiner(53, 50, "WA-2024-09", "World average (Sep 2024)");
  gc.getCombiner(53)->addPdf(gc[55]);  // BESIII measurement of Fp_pipipi0

  // WA October 2024 (new WS/RS with DT Run 2 data)
  gc.cloneCombiner(54, 53, "WA-2024-10", "World average (October 2024)");
  gc.getCombiner(54)->delPdf(gc[35]);  // WS/RS in D0 -> Kpi from LHCb Run 1 DT
  gc.getCombiner(54)->addPdf(gc[41]);  // WS/RS in D0 -> Kpi from LHCb Run 1+2 DT

  // WA October 2025 (new BinFlip from Belle + Belle 2, new BESIII Delta_Kpi, no new LHCb D0 -> K3pi Run 2) ------------
  gc.cloneCombiner(55, 54, "WA-2025-10", "World average (October 2025)");
  gc.getCombiner(55)->delPdf(gc[20]);  // D0 -> KS hh from Belle
  gc.getCombiner(55)->addPdf(gc[6]);   // D0 -> KS pi pi BinFlip Belle + Belle 2
  gc.getCombiner(55)->delPdf(gc[52]);  // D0 -> Kpi BESIII 3   fb
  gc.getCombiner(55)->addPdf(gc[56]);  // D0 -> Kpi BESIII 3+7 fb

  // LHCb-only averages ------------------------------------------------------------------------------------------------

  gc.newCombiner(300, "LHCb-2024-05", "LHCb average (May 2024)");
  for (const auto imeas : get_lhcb_pdfs("run12", dy_fsc_hypo)) gc.getCombiner(300)->addPdf(gc[imeas]);

  // LHCb-only + charm factories averages ------------------------------------------------------------------------------

  gc.cloneCombiner(400, 300, "LHCb-CF-2024-05", "LHCb + Charm factories average (May 2024)");
  for (auto imeas : {50, 52, 53}) gc.getCombiner(400)->addPdf(gc[imeas]);

  // Impact of LHCb upgrades -------------------------------------------------------------------------------------------

  // WA before LHCb Run 2
  if (dy_fsc_hypo == dy_fsc::none) {
    gc.newCombiner(500, "LHCb-Run1", "World average before LHCb Run 2");
    for (auto imeas : {1, 2, 3, 4, 10, 11, 20, 21, 30, 31, 32, 35, 36, 50, 52, 60, 61, 62, 63, 70, 93})
      gc.getCombiner(500)->addPdf(gc[imeas]);
  }

  // WA after LHCb Run 2
  gc.cloneCombiner(501, 50, "LHCb-Run2", "World average after LHCb Run 2");

  ///////////////////////////////////////////////////
  //
  // Run
  //
  ///////////////////////////////////////////////////

  gc.run();

  return 0;
}
