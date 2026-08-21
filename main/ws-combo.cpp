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
  const auto mix = parametrisations::mix::d0_to_kpi;
  using parametrisations::kpi;

  // clang-format off
  gc.addPdf(0, new PDF_WS_NoCPV("CDF",   mix),                           "WS/RS  CDF    No CPV                      ");
  gc.addPdf(1, new PDF_WS_NoCPV("BaBar", mix),                           "WS/RS  BaBar  No CPV                      ");  // Alternative to 10
  gc.addPdf(2, new PDF_WS_NoCPV("Belle", mix),                           "WS/RS  Belle  No CPV                      ");  // Incompatible with 11 and based on a larger sample

  gc.addPdf(10, new PDF_WS("BaBar", mix),                                "WS/RS  BaBar                              ");  // Alternative to 1
  gc.addPdf(11, new PDF_WS("Belle", mix),                                "WS/RS  Belle                              ");  // Incompatible with 2 and based on a smaller sample

  gc.addPdf(30, new PDF_WS("LHCb-R1-prompt-biased",  mix),               "WS/RS  LHCb   Run 1    [D -> D* pi] biased");  // Ghost pions neglected. Superseded by 31
  gc.addPdf(31, new PDF_WS("LHCb-R1-prompt",         mix),               "WS/RS  LHCb   Run 1    [D -> D* pi]       ");  // Supersedes 30; superseded by 34
  gc.addPdf(32, new PDF_WS("LHCb-R1-DT",             mix),               "WS/RS  LHCb   Run 1    [B -> D* mu]       ");  // Superseded by 37
  gc.addPdf(33, new PDF_WS("LHCb-R1-biased",         mix),               "WS/RS  LHCb   Run 1                       ");
  gc.addPdf(34, new PDF_WS("LHCb-2011-2016-prompt",  mix),               "WS/RS  LHCb   2011-6   [D* -> D0 pi]      ");  // Supersedes 31; superseded by 35
  gc.addPdf(35, new PDF_WS("LHCb-R12-prompt-sec9",   mix, kpi::ccprime), "WS/RS  LHCb   Run 1+2  [D* -> D0 pi]      ");
  gc.addPdf(36, new PDF_WS("LHCb-R2-DT",             mix),               "WS/RS  LHCb   Run 2    [B -> D* mu]       ");  // Superseded by 37
  gc.addPdf(37, new PDF_WS("LHCb-R12-DT",            mix),               "WS/RS  LHCb   Run 1+2  [B -> D* mu]       ");
  // clang-format on

  // Define the combinations
  gc.newCombiner(0, "empty", "empty");

  gc.newCombiner(1, "NonLHCb-2013", "CDF + BaBar + Belle", {0, 2, 10});
  gc.newCombiner(2, "NonLHCb-2014-ForCPV", "CDF + BaBar + Belle", {0, 10, 11});

  gc.newCombiner(100, "LHCb-2018", "LHCb 2018", {32, 34});

  gc.newCombiner(200, "LHCb-R1", "LHCb Run 1", {31, 32});
  gc.newCombiner(201, "LHCb-R12", "LHCb Run 1+2", {35, 37});

  gc.newCombiner(300, "WA-2024", "World average 2024", {0, 2, 10, 35, 37});

  // Run the combination
  gc.run();

  return 0;
}
