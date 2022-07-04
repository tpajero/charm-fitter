#include "combiner.h"


void run_combiner(std::string name, theory_config th_cfg, int argc, char* argv[])
{
    GammaComboEngine gc(name, argc, argv);

    ///////////////////////////////////////////////////
    //
    // define PDFs
    //
    ///////////////////////////////////////////////////

    gc.addPdf(1,  new PDF_XY      ("BaBar_Kshh", "BaBar_Kshh", "BaBar_Kshh", th_cfg), "BaBar Kshh");
    gc.addPdf(2,  new PDF_XY      ("BaBar_pipipi0", "BaBar_pipipi0", "BaBar_pipipi0", th_cfg), "BaBar pipipi0");
    gc.addPdf(3,  new PDF_XY      ("LHCb_KSpipi", "LHCb_KSpipi", "LHCb_KSpipi", th_cfg), "LHCb KSpipi 2011");
    gc.addPdf(4,  new PDF_Kpipi0  ("BaBar", "BaBar", "BaBar", th_cfg), "BaBar Kpipi0");

    gc.addPdf(10, new PDF_RM      ("HFLAV2016", "HFLAV2016", "HFLAV2016", th_cfg), "R_M HFLAV 2016");
    gc.addPdf(11, new PDF_RM      ("LHCb_K3pi_Run1", "LHCb_K3pi_Run1", "LHCb_K3pi_Run1", th_cfg), "R_M K3pi LHCb");

    gc.addPdf(20, new PDF_Kshh    ("Belle", "Belle", "Belle", th_cfg), "Belle Kshh");
    gc.addPdf(21, new PDF_BinFlip ("LHCb_Run1", "LHCb_Run1", "LHCb_Run1", th_cfg), "Bin-flip LHCb Run 1");
    gc.addPdf(22, new PDF_BinFlip ("LHCb_Run2_prompt", "LHCb_Run2_prompt", "LHCb_Run2_prompt", th_cfg), "Bin-flip LHCb Run 2 prompt");

    gc.addPdf(30, new PDF_WS_NoCPV("CDF", "CDF", "CDF", th_cfg), "WS/RS CDF");
    gc.addPdf(31, new PDF_WS_NoCPV("BaBar", "BaBar", "BaBar", th_cfg), "WS/RS BaBar");
    gc.addPdf(32, new PDF_WS_NoCPV("Belle", "Belle", "Belle", th_cfg), "WS/RS Belle");
    gc.addPdf(33, new PDF_WS      ("LHCb_Prompt_2011_2016", "LHCb_Prompt_2011_2016", "LHCb_Prompt_2011_2016", th_cfg), "WS/RS LHCb");
    gc.addPdf(34, new PDF_WS      ("LHCb_DT_Run1", "LHCb_DT_Run1", "LHCb_DT_Run1", th_cfg), "WS/RS LHCb DT");

    gc.addPdf(50,  new PDF_Cleo   ("Cleo-c", "Cleo-c", "Cleo-c", th_cfg), "Cleo-c");
    gc.addPdf(51,  new PDF_BES    ("BES", "BES", "BES", th_cfg), "BES");

    gc.addPdf(60, new PDF_yCP             ("WA2020", "WA2020", "WA2020", th_cfg), "y_CP (WA 2020)");
    gc.addPdf(61, new PDF_yCP_minus_yCP_RS("WA2020", "WA2020", "WA2020", th_cfg), "y_CP-yCP(RS) (WA 2020)");
    gc.addPdf(62, new PDF_yCP_minus_yCP_KP("WA2020", "WA2020", "WA2020", th_cfg), "y_CP-yCP(KP) (WA 2020)");
    gc.addPdf(63, new PDF_yCP_plus_yCP_RS ("WA2020", "WA2020", "WA2020", th_cfg), "y_CP+yCP(RS) (WA 2020)");
    gc.addPdf(64, new PDF_yCP_minus_yCP_RS("LHCb2022", "LHCb2022", "LHCb2022", th_cfg), "y_CP-yCP(RS) (LHCb 2022)");

    gc.addPdf(70, new PDF_DY      ("WA2020", "WA2020", "WA2020", th_cfg), "DY (2020 WA)");
    gc.addPdf(71, new PDF_DY      ("WA2021", "WA2021", "WA2021", th_cfg), "DY (2021 WA)");

    gc.addPdf(80, new PDF_DY_RS   ("LHCb2021", "LHCb2021", "LHCb2021", th_cfg), "DYRS (LHCb 2021)");

    gc.addPdf(100, new PDF_scan_DY_RS(th_cfg), "ScanDYRS");

    gc.addPdf(110, new PDF_yCP    ("WA2020_biased", "WA2020_biased", "WA2020_biased", th_cfg), "y_CP (WA 2020, biased)");
    gc.addPdf(111, new PDF_yCP    ("LHCb2022_biased", "LHCb2022_biased", "LHCb2022_biased", th_cfg), "y_CP-yCP(RS) (LHCb 2022, biased)");

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
    gc.cloneCombiner(2, 1, "WA2021", "World average (June 2021)");
    gc.getCombiner(2)->addPdf(gc[22]);  // bin-flip run 2
    gc.getCombiner(2)->delPdf(gc[70]);  // DY WA 2020
    gc.getCombiner(2)->addPdf(gc[71]);  // DY WA 2021

    // WA after LHCb 2022 yCP measurement
    gc.cloneCombiner(3, 2, "WAFeb2022", "World average (Feb 2022)");
    gc.getCombiner(3)->addPdf(gc[64]);  // yCP LHCb 2022

    // WA after LHCb 2022 yCP measurement - biased
    gc.cloneCombiner(100, 2, "WAFeb2022_biased", "World average (Feb 2022) #minus no #it{y}_{#it{CP}}^{#it{K^{#minus}#pi^{+}}} correction");
    gc.getCombiner(100)->delPdf(gc[60]);
    gc.getCombiner(100)->delPdf(gc[61]);
    gc.getCombiner(100)->delPdf(gc[62]);
    gc.getCombiner(100)->delPdf(gc[63]);
    gc.getCombiner(100)->addPdf(gc[110]);
    gc.getCombiner(100)->addPdf(gc[111]);

    ///////////////////////////////////////////////////
    //
    // Run
    //
    ///////////////////////////////////////////////////

    gc.run();
}
