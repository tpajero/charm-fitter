/**
 * Combination of all WS/RS D0 -> K pi measurements.
 *
 * Needed to plot the y' vs. x'2 contours plot.
 */

// Core
#include <Combiner.h>
#include <GammaComboEngine.h>

// CharmFitter
#include <CharmUtils.h>
#include <PDF_WS.h>
#include <PDF_WS_NoCPV.h>

#include <format>
#include <stdexcept>
#include <string>
#include <vector>

/**
 * Main function to combine WS/RS D0 -> K pi measurements using the (y', x'2) parametrisation.
 */
int main(int argc, char* argv[]) {

  GammaComboEngine gc("ws-combo", argc, &argv[0]);

  // Define the PDFs
  const auto mix_param = parametrisations::mix::d0_to_kpi;
  using parametrisations::kpi;

  // clang-format off
  gc.addPdf(0, new PDF_WS_NoCPV("CDF",   mix_param),                           "WS/RS  CDF    (no CPV)              ");
  gc.addPdf(1, new PDF_WS_NoCPV("BaBar", mix_param),                           "WS/RS  BaBar  no CPV                ");
  gc.addPdf(2, new PDF_WS_NoCPV("Belle", mix_param),                           "WS/RS  Belle  no CPV                ");

  gc.addPdf(10, new PDF_WS("BaBar", mix_param),                                "WS/RS  BaBar                        ");
  gc.addPdf(11, new PDF_WS("Belle", mix_param),                                "WS/RS  Belle                        ");

  gc.addPdf(20, new PDF_WS("LHCb_DT_Run1",           mix_param),               "WS/RS  LHCb   Run 1    [B -> D* mu] ");
  gc.addPdf(21, new PDF_WS("LHCb_Run1",              mix_param),               "WS/RS  LHCb   Run 1                 ");
  gc.addPdf(22, new PDF_WS("LHCb_Prompt_2011_2016",  mix_param),               "WS/RS  LHCb   2011-6   [D* -> D0 pi]");
  gc.addPdf(23, new PDF_WS("LHCb_Prompt_Run12_sec9", mix_param, kpi::ccprime), "WS/RS  LHCb   Run 1-2  [D* -> D0 pi]");
  gc.addPdf(24, new PDF_WS("LHCb_DT_Run2",           mix_param),               "WS/RS  LHCb   Run 2    [B -> D* mu] ");
  gc.addPdf(25, new PDF_WS("LHCb_DT_Run12",          mix_param),               "WS/RS  LHCb   Run 1-2  [B -> D* mu] ");
  // clang-format on

  // Define the combinations
  gc.newCombiner(0, "empty", "empty");

  gc.newCombiner(1, "NonLHCb2025", "CDF + BaBar + Belle", {0, 10, 11});

  gc.newCombiner(10, "LHCb2024", "LHCb 2024", {20, 23});
  gc.newCombiner(11, "LHCb2025", "LHCb Run 1+2", {23, 25});

  gc.newCombiner(20, "WA2025", "World average 2025", {0, 10, 11, 23, 25});

  // Run the combination
  gc.run();

  return 0;
}
