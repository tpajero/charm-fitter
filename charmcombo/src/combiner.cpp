#include "combiner.h"
#include "ParametersCharmCombo.h"
#include "PDF_AcpHH_LHCb_Run12.h"
#include "PDF_BES_Kpi.h"
#include "PDF_BES_Kpi_1d.h"
#include "PDF_BinFlip.h"
#include "PDF_Cleo.h"
#include "PDF_Fp_pipipi0.h"
#include "PDF_DY.h"
#include "PDF_DY_RS.h"
#include "PDF_Kpipi0.h"
#include "PDF_Kshh.h"
#include "PDF_RM.h"
#include "PDF_scan_DY_RS.h"
#include "PDF_XY.h"
#include "PDF_yCP.h"
#include "PDF_yCP_minus_yCP_RS.h"
#include "PDF_yCP_minus_yCP_KP.h"
#include "PDF_yCP_plus_yCP_RS.h"
#include "PDF_WS.h"
#include "PDF_WS_NoCPV.h"

#include <string>

#include "GammaComboEngine.h"


void run_combiner(std::string name, theory_config th_cfg, int argc, char* argv[]) {

    GammaComboEngine gc(name, argc, argv);

    ///////////////////////////////////////////////////
    //
    // define PDFs
    //
    ///////////////////////////////////////////////////
"YCP        LHCb  3fb  D0  -> hh         [B+  -> D0mu+nu] "
    gc.addPdf(1,  new PDF_XY      ("BaBar_Kshh",            th_cfg), "XY KShh     BaBar                      ");
    gc.addPdf(2,  new PDF_XY      ("BaBar_pipipi0",         th_cfg), "XY pipipi0  BaBar                      ");
    gc.addPdf(3,  new PDF_XY      ("LHCb_KSpipi",           th_cfg), "XY KSpipi   LHCb   2011   [D* -> D0 pi]");
    gc.addPdf(4,  new PDF_Kpipi0  ("BaBar",                 th_cfg), "Kpipi0      BaBar                      ");

    gc.addPdf(10, new PDF_RM      ("HFLAV2016",             th_cfg), "R_M         HFLAV  2016                ");
    gc.addPdf(11, new PDF_RM      ("LHCb_K3pi_Run1",        th_cfg), "R_M K3pi    LHCb                       ");

    gc.addPdf(20, new PDF_Kshh    ("Belle",                 th_cfg), "KShh        Belle                      ");
    gc.addPdf(21, new PDF_BinFlip ("LHCb_Run1",             th_cfg), "Bin-flip    LHCb   Run 1               ");
    gc.addPdf(22, new PDF_BinFlip ("LHCb_Run2_prompt",      th_cfg), "Bin-flip    LHCb   Run 2  [D* -> D0 pi]");
    gc.addPdf(23, new PDF_BinFlip ("LHCb_Run2_sl",          th_cfg), "Bin-flip    LHCb   Run 2  [B -> D0 mu] ");
    gc.addPdf(24, new PDF_BinFlip ("LHCb_Run2",             th_cfg), "Bin-flip    LHCb   Run 2               ");

    gc.addPdf(30, new PDF_WS_NoCPV("CDF",                   th_cfg), "WS/RS       CDF                        ");
    gc.addPdf(31, new PDF_WS_NoCPV("BaBar",                 th_cfg), "WS/RS       BaBar                      ");
    gc.addPdf(32, new PDF_WS_NoCPV("Belle",                 th_cfg), "WS/RS       Belle                      ");
    gc.addPdf(33, new PDF_WS      ("LHCb_Prompt_2011_2016", th_cfg), "WS/RS       LHCb   2011-6 [B -> D* mu] ");
    gc.addPdf(34, new PDF_WS      ("LHCb_DT_Run1",          th_cfg), "WS/RS       LHCb   Run 1  [B -> D* mu] ");

    gc.addPdf(50,  new PDF_Cleo      ("Cleo-c",             th_cfg), "Delta_Kpi   Cleo-c                     ");
    gc.addPdf(51,  new PDF_BES_Kpi_1d(                      th_cfg), "Delta_Kpi   BES    3fb                 ");
    gc.addPdf(52,  new PDF_BES_Kpi   (                      th_cfg), "Delta_Kpi   BES    3fb                 ");
    gc.addPdf(53,  new PDF_Fp_pipipi0(                            ), "Fpipipi0    Cleo-c                     ");

    gc.addPdf(60, new PDF_yCP             ("WA2020",        th_cfg), "yCP         WA     2020                ");
    gc.addPdf(61, new PDF_yCP_minus_yCP_RS("WA2020",        th_cfg), "yCP-yCP(RS) WA     2020                ");
    gc.addPdf(62, new PDF_yCP_minus_yCP_KP("WA2020",        th_cfg), "yCP-yCP(KP) WA     2020                ");
    gc.addPdf(63, new PDF_yCP_plus_yCP_RS ("WA2020",        th_cfg), "yCP+yCP(RS) WA     2020                ");
    gc.addPdf(64, new PDF_yCP_minus_yCP_RS("LHCb2022",      th_cfg), "yCP-yCP(RS) LHCb   2022                ");

    gc.addPdf(70, new PDF_DY      ("WA2020",                th_cfg), "DY          WA     2020                ");
    gc.addPdf(71, new PDF_DY      ("WA2021",                th_cfg), "DY          WA     2021                ");

    gc.addPdf(80, new PDF_DY_RS   ("LHCb2021",              th_cfg), "DY(RS)      LHCb   2021                ");

    gc.addPdf(100, new PDF_scan_DY_RS(                      th_cfg), "ScanDYRS    This is just a nuisance parameter");

    gc.addPdf(110, new PDF_yCP    ("WA2020_biased",         th_cfg), "yCP         WA     2020   [biased]     ");
    gc.addPdf(111, new PDF_yCP    ("LHCb2022_biased",       th_cfg), "yCP-yCP(RS) LHCb   2022   [biased]     ");

    gc.addPdf(120, new PDF_AcpHH_LHCb_Run12(                th_cfg), "ACP(KK/PP)  LHCb   Run1+2              ");

    ///////////////////////////////////////////////////
    //
    // Define combinations
    //
    ///////////////////////////////////////////////////

    gc.newCombiner(0, "emtpy", "empty");

    // WA 2020
    gc.newCombiner(1, "WA2020", "World average (Dec 2020)", 1, 2, 3, 4, 10, 11, 20, 21, 30, 31, 32, 33, 34, 50, 51);
    gc.getCombiner(1)->addPdf(gc[60]);
    gc.getCombiner(1)->addPdf(gc[61]);
    gc.getCombiner(1)->addPdf(gc[62]);
    gc.getCombiner(1)->addPdf(gc[63]);
    gc.getCombiner(1)->addPdf(gc[70]);

    // WA June 2021
    gc.cloneCombiner(20, 1, "WA2021", "World average (June 2021)");
    gc.getCombiner(20)->addPdf(gc[22]);  // bin-flip run 2
    gc.getCombiner(20)->delPdf(gc[70]);  // DY WA 2020
    gc.getCombiner(20)->addPdf(gc[71]);  // DY WA 2021

    // WA after LHCb 2022 yCP measurement
    gc.cloneCombiner(30, 20, "WAFeb2022", "World average (Feb 2022)");
    gc.getCombiner(30)->addPdf(gc[64]);  // yCP LHCb 2022

    // WA after LHCb 2022 yCP measurement - biased
    gc.cloneCombiner(31, 30, "WAFeb2022_biased", "World average (Feb 2022) #minus no #it{y}_{#it{CP}}^{#it{K^{#minus}#pi^{+}}} correction");
    gc.getCombiner(31)->delPdf(gc[60]);
    gc.getCombiner(31)->delPdf(gc[61]);
    gc.getCombiner(31)->delPdf(gc[62]);
    gc.getCombiner(31)->delPdf(gc[63]);
    gc.getCombiner(31)->delPdf(gc[64]);
    gc.getCombiner(31)->addPdf(gc[110]);
    gc.getCombiner(31)->addPdf(gc[111]);

    // WA September 2022
    gc.cloneCombiner(40, 30, "WASept2022", "World average (Sept 2022)");
    gc.getCombiner(40)->delPdf(gc[51]);  // old BESIII measurement of delta_Kpi
    gc.getCombiner(40)->addPdf(gc[52]);  // new BESIII measurement of delta_Kpi
    gc.getCombiner(40)->addPdf(gc[53]);  // F+_pipipi0
    gc.getCombiner(40)->delPdf(gc[22]);  // bin-flip LHCb Run 2 prompt
    gc.getCombiner(40)->addPdf(gc[24]);  // bin-flip LHCb Run 2
    gc.getCombiner(40)->addPdf(gc[120]); // ACP(KK) + DeltaACP LHCb Run 1+2

    ///////////////////////////////////////////////////
    //
    // Run
    //
    ///////////////////////////////////////////////////

    gc.run();
}
