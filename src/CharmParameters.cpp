/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include <CharmParameters.h>

#include <Utils.h>

using Utils::DegToRad;

CharmParameters::CharmParameters() { defineParameters(); }

/**
 * Define all (nuisance) parameters.
 *
 * scan:  defines scan range (for Prob and Plugin methods)
 * phys:  physically allowed range (needs to be set!)
 */
void CharmParameters::defineParameters() {
  Parameter* p = nullptr;

  // CP asymmetry in the decay (``direct'' CP asymmetry) of D0 -> K+ K-
  p = newParameter("Acp_KK");
  p->title = "#it{a}_{#it{K}^{+}#it{K}^{#minus}}^{d}";
  p->startvalue = 8e-4;
  p->unit = "";
  p->scan = range(-2e-2, 2e-2);
  p->phys = range(-1., 1.);

  // CP asymmetry in the decay (``direct'' CP asymmetry) of D0 -> pi+ pi-
  p = newParameter("Acp_PP");
  p->title = "#it{a}_{#it{#pi}^{+}#it{#pi}^{#minus}}^{d}";
  p->startvalue = 2.4e-3;
  p->unit = "";
  p->scan = range(-2e-2, 2e-2);
  p->phys = range(-1., 1.);

  // CP asymmetry in the decay (``direct'' CP asymmetry) of D0 -> K+ pi-
  p = newParameter("Acp_KP");
  p->title = "#it{a}_{#it{K}^{+}#it{#pi}^{#minus}}^{d}";
  p->startvalue = -6e-3;
  p->unit = "";
  p->scan = range(-2.5e-2, 2e-2);
  p->phys = range(-1., 1.);

  // Ratio of the magnitudes of the decay amplitudes of D0 -> K+ pi- to D0 -> K- pi+
  p = newParameter("r_Kpi");
  p->title = "#it{r_{K#pi}}";
  p->startvalue = 5.857e-2;
  p->unit = "";
  p->scan = range(5.74e-2, 6e-2);
  p->phys = range(0., 10.);

  // Ratio of the magnitudes of the decay amplitudes of D0 -> K+ pi- pi- pi+ to D0 -> K- pi- pi+ pi+
  p = newParameter("r_K3pi");
  p->title = "#it{r_{K3#pi}}";
  p->startvalue = 5.5e-2;
  p->unit = "";
  p->scan = range(4e-2, 6e-2);
  p->phys = range(0., 1e2);

  // Ratio of the magnitudes of the decay amplitudes of D0 -> K+ pi- pi0 to D0 -> K- pi+ pi0
  p = newParameter("r_Kpipi0");
  p->title = "#it{r_{K#pi#pi^{0}}}";
  p->startvalue = 4.4e-2;
  p->unit = "";
  p->scan = range(4e-2, 6e-2);
  p->phys = range(0., 1e2);

  // Coherence parameters for the decay amplitudes of D0 -> K+ pi- pi- pi+ to D0 -> K- pi- pi+ pi+
  p = newParameter("k_K3pi");
  p->title = "#it{#kappa_{K3#pi}}";
  p->startvalue = 0.44;
  p->unit = "";
  p->scan = range(0., 1.);
  p->phys = range(0, 1e4);

  // Coherence parameters for the decay amplitudes of D0 -> K+ pi- pi0 to D0 -> K- pi+ pi0
  p = newParameter("k_Kpipi0");
  p->title = "#it{#kappa_{K#pi#pi^{0}}}";
  p->startvalue = 0.79;
  p->unit = "";
  p->scan = range(0., 1.);
  p->phys = range(0, 1e4);

  // Strong-phase difference between the decay amplitudes of D0 -> K+ pi- to D0 -> K- pi+
  // Follows the convention of:
  //   A. L. Kagan and L. Silvestrini, Dispersive and absorptive CP violation in D0/anti-D0 mixing
  //   Phys. Rev. D 103 053008, 2021 (https://inspirehep.net/literature/1776611)
  p = newParameter("Delta_Kpi");
  p->title = "#Delta_{#it{K#pi}} [rad]";
  p->startvalue = -0.25;
  p->unit = "";
  p->scan = range(DegToRad(-180), DegToRad(180));
  p->phys = range(DegToRad(-180), DegToRad(180));

  // Strong-phase difference between the decay amplitudes of D0 -> K+ pi- pi0 to D0 -> K- pi+ pi0
  p = newParameter("Delta_Kpipi0");
  p->title = "#Delta_{#it{K#pi#pi^{0}}} [rad]";
  p->startvalue = -0.24;
  p->unit = "";
  p->scan = range(DegToRad(-180), DegToRad(180));
  p->phys = range(DegToRad(-180), DegToRad(180));

  // Strong-phase difference between the decay amplitudes of D0 -> K+ pi- pi- pi+ and D0 -> K- pi- pi+ pi+
  // Follows the convention of:
  //   A. L. Kagan and L. Silvestrini, Dispersive and absorptive CP violation in D0/anti-D0 mixing
  //   Phys. Rev. D 103 053008, 2021 (https://inspirehep.net/literature/1776611)
  p = newParameter("Delta_K3pi");
  p->title = "#Delta_{#it{K3#pi}} [rad]";
  p->startvalue = 0.4;
  p->unit = "";
  p->scan = range(DegToRad(-180), DegToRad(180));
  p->phys = range(DegToRad(-180), DegToRad(180));

  // Parametrising the decay amplitude of D0 -> K+ K- decays as
  //
  //   A = A_sd (lambda_s - lambda_d) / 2 + A_b lambda_b / 2,
  //
  // delta_KK is defined as arg(A_b / A_sd).
  p = newParameter("cot_delta_KK");
  p->title = "cot(#it{#delta_{#it{K^{+}K^{#minus}}}})";
  p->startvalue = 0.;
  p->unit = "";
  p->scan = range(-1e4, 1e4);
  p->phys = range(-1e4, 1e4);

  // Parametrising the decay amplitude of D0 -> pi+ pi- decays as
  //
  //   A = A_sd (lambda_s - lambda_d) / 2 + A_b lambda_b / 2,
  //
  // delta_PP is defined as arg(A_b / A_sd).
  p = newParameter("cot_delta_PP");
  p->title = "cot(#it{#delta_{#it{#pi^{+}#pi^{#minus}}}})";
  p->startvalue = 0.;
  p->unit = "";
  p->scan = range(-1e4, 1e4);
  p->phys = range(-1e4, 1e4);

  // CP-even fraction of D0 -> pi+ pi- pi0
  p = newParameter("F_pipipi0");
  p->title = "#it{F^{#pi#pi#pi^{0}}_{+}}";
  p->startvalue = 0.943;
  p->unit = "";
  p->scan = range(0, 1);
  p->phys = range(-1, 2);
  p->force = range(0, 1);

  /**
   * Notes on D0 mixing:
   *
   * By default, we rely on the theoretical parametrisation (see https://arxiv.org/abs/2001.07207).
   *
   * The phenomenological parametrisation is also implemented, with the definition
   *
   *     |D_1,2> = p |D0> \pm q |anti-D0>
   *
   * and the following conventions:
   *
   *   1. CP|D0> = - |anti-D0>
   *   2. |D_1(2)> is the approximately CP-odd (CP-even) mass eigenstate
   *   3. x = (m_2 - m_1) / Gamma; y = (Gamma_2 - Gamma_1) / (2 Gamma)
   */

  constexpr double x12_start = 4.1e-3;
  constexpr double y12_start = 6.6e-3;

  // Phenomenological parametrisation ----------------------------------------------------------------------------------

  p = newParameter("x");
  p->title = "#it{x}";
  p->startvalue = x12_start;
  p->unit = "";
  p->scan = range(-4e-3, 1e-2);
  p->phys = range(-1e2, 1e2);

  p = newParameter("y");
  p->title = "#it{y}";
  p->startvalue = y12_start;
  p->unit = "";
  p->scan = range(1e-3, 1e-2);
  p->phys = range(-1e2, 1e2);

  p = newParameter("qop");
  p->title = "|#it{q/p}|";
  p->startvalue = 0.98;
  p->unit = "";
  p->scan = range(0.75, 1.25);
  p->phys = range(0., 1e4);

  p = newParameter("phi");
  p->title = "#it{#phi}_{2} [rad]";
  p->startvalue = -0.03;
  p->unit = "";
  p->scan = range(-0.6, 0.6);
  p->phys = range(DegToRad(-180.), DegToRad(180.));

  // Theoretical parametrisation ---------------------------------------------------------------------------------------

  p = newParameter("x12");
  p->title = "#it{x}_{12}";
  p->startvalue = x12_start;
  p->unit = "";
  p->scan = range(0., 1e-2);
  p->phys = range(0., 1e2);

  p = newParameter("y12");
  p->title = "#it{y}_{12}";
  p->startvalue = y12_start;
  p->unit = "";
  p->scan = range(0., 1.1e-2);
  p->phys = range(0., 1e2);

  p = newParameter("phiM");
  p->title = "#it{#phi}^{#it{M}}_{2} [rad]";
  p->startvalue = 0.02;
  p->unit = "";
  p->scan = range(-0.5, 0.5);
  p->phys = range(-1.5, -1.5 + DegToRad(360.));

  p = newParameter("phiG");
  p->title = "#it{#phi}^{#Gamma}_{2} [rad]";
  p->startvalue = 0.02;
  p->unit = "";
  p->scan = range(-0.3, 0.3);
  p->phys = range(DegToRad(-180.), DegToRad(180.));

  // Other parameters --------------------------------------------------------------------------------------------------

  // Nuisance parameter to get predictions for CP violation in RS decays
  p = newParameter("DY_RS");
  p->title = "#||{#Delta#it{Y}_{#it{K}^{#minus}#pi^{+}}}";
  p->startvalue = 5e-6;
  p->unit = "";
  p->scan = range(0., 1e-4);
  p->phys = range(0., 1.);

  // Parameters for (D0 -> Kpi)-only combination -----------------------------------------------------------------------

  p = newParameter("yp");
  p->title = "#it{y'}";
  p->startvalue = 5.3e-3;
  p->unit = "";
  p->scan = range(0.1, 1);
  p->phys = range(-1e4, 1e4);

  p = newParameter("dyp");
  p->title = "#it{#Delta y'}";
  p->startvalue = 1e-5;
  p->unit = "";
  p->scan = range(0.1, 1);
  p->phys = range(-1e4, 1e4);

  p = newParameter("xp2");
  p->title = "#it{x'}^{2}";
  p->startvalue = 1.8e-5;
  p->unit = "";
  p->scan = range(-4, 4);
  p->phys = range(-1e4, 1e4);

  p = newParameter("dxp2");
  p->title = "#it{#Delta x'}^{2}";
  p->startvalue = 1e-6;
  p->unit = "";
  p->scan = range(-0.4, 1.);
  p->phys = range(-1e4, 1e4);
}
