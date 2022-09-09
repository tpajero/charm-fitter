/**
 * Charm Combination
 * Author: Tommaso Pajero, tommaso.pajero@cern.ch
 * Date: October 2021
 **/

#include "ParametersCharmCombo.h"


ParametersCharmCombo::ParametersCharmCombo(
        theory_config th_cfg, bool predict_dy_rs)
{
    defineParameters(th_cfg, predict_dy_rs);
}


/**
 * Define all (nuisance) parameters.

 * scan:  defines scan range (for Prob and Plugin methods)
 * phys:  physically allowed range (needs to be set!)
 * free:  range applied when no boundary is required - set this far away from any possible value
 **/
void ParametersCharmCombo::defineParameters(
        theory_config th_cfg, bool predict_dy_rs)
{
    Parameter *p = nullptr;

    p = newParameter("AcpKK");
    p->title = "#it{a}_{K^{+}K^{#minus}}^{d} [%]";
    p->startvalue = -0.05;
    p->unit = "";
    p->scan = range(-2, 2);
    p->phys = range(-100, 100);
    p->free = range(-1e4, 1e4);

    p = newParameter("AcpPP");
    p->title = "#it{a}_{#pi^{+}#pi^{#minus}}^{d} [%]";
    p->startvalue = -0.05;
    p->unit = "";
    p->scan = range(-2, 2);
    p->phys = range(-100, 100);
    p->free = range(-1e4, 1e4);

    p = newParameter("R_Kpi");
    p->title = "#it{R_{K#pi}} [%]";
    p->startvalue = 0.345;
    p->unit = "";
    p->scan = range(0.33, 0.36);
    p->phys = range(0, 1e4);
    p->free = range(0, 1e4);

    if (th_cfg != superweak) {
        p = newParameter("A_Kpi");
        p->title = "#it{A_{K#pi}} [%]";
        p->startvalue = -0.6;
        p->unit = "";
        p->scan = range(-2.5, 2);
        p->phys = range(-100, 100);
        p->free = range(-1e4, 1e4);
    }

    p = newParameter("Delta_Kpi");
    p->title = "#Delta_{#it{K#pi}} [rad]";
    p->startvalue = -0.2;
    p->unit = "";
    p->scan = range(DegToRad(-180), DegToRad(180));
    p->phys = range(DegToRad(-180), DegToRad(180));
    p->free = range(DegToRad(-180), DegToRad(180));

    p = newParameter("Delta_Kpipi");
    p->title = "#Delta_{#it{K#pi#pi}} [rad]";
    p->startvalue = -0.4;
    p->unit = "";
    p->scan = range(DegToRad(-180), DegToRad(180));
    p->phys = range(DegToRad(-180), DegToRad(180));
    p->free = range(DegToRad(-180), DegToRad(180));

    switch (th_cfg) {
        case phenomenological:
            p = newParameter("x");
            p->title = "#it{x} [%]";
            p->startvalue = 0.4;
            p->unit = "";
            p->scan = range(-0.4, 1.);
            p->phys = range(-1e4, 1e4);
            p->free = range(-1e4, 1e4);

            p = newParameter("y");
            p->title = "#it{y} [%]";
            p->startvalue = 0.65;
            p->unit = "";
            p->scan = range(0.1,  1);
            p->phys = range(-1e4, 1e4);
            p->free = range(-1e4, 1e4);

            p = newParameter("qop");
            p->title = "|#it{q/p}|#minus1";
            p->startvalue = -0.02;
            p->unit = "";
            p->scan = range(-0.25, 0.25);
            p->phys = range(-1., 1e4);
            p->free = range(-1., 1e4);

            p = newParameter("phi");
            p->title = "#it{#phi}_{2} [rad]";
            p->startvalue = 0.02;
            p->unit = "";
            p->scan = range(-0.6, 0.6);
            p->phys = range(DegToRad(-180.), DegToRad(180.));
            p->free = range(DegToRad(-180.), DegToRad(180.));
            break;

        case theoretical:
            p = newParameter("x12");
            p->title = "#it{x}_{12} [%]";
            p->startvalue = 0.4;
            p->unit = "";
            p->scan = range(0, 1);
            p->phys = range(0, 1e4);
            p->free = range(0, 1e4);

            p = newParameter("y12");
            p->title = "#it{y}_{12} [%]";
            p->startvalue = 0.65;
            p->unit = "";
            p->scan = range(0., 1.1);
            p->phys = range(0, 1e4);
            p->free = range(0, 1e4);

            p = newParameter("phiM");
            p->title = "#it{#phi}^{#it{M}}_{2} [rad]";
            p->startvalue = 0.02;
            p->unit = "";
            p->scan = range(-0.5, 0.5);
            p->phys = range(-1.5, -1.5 + DegToRad(360.));
            p->free = range(-1.5, -1.5 + DegToRad(360.));

            p = newParameter("phiG");
            p->title = "#it{#phi}^{#Gamma}_{2} [rad]";
            p->startvalue = 0.02;
            p->unit = "";
            p->scan = range(-0.3, 0.3);
            p->phys = range(DegToRad(-180.), DegToRad(180.));
            p->free = range(DegToRad(-180.), DegToRad(180.));
            break;

        case superweak:
            p = newParameter("x12");
            p->title = "#it{x}_{12} [%]";
            p->startvalue = 0.4;
            p->unit = "";
            p->scan = range(0, 1);
            p->phys = range(0, 1e4);
            p->free = range(0, 1e4);

            p = newParameter("y12");
            p->title = "#it{y}_{12} [%]";
            p->startvalue = 0.65;
            p->unit = "";
            p->scan = range(0, 1.1);
            p->phys = range(0, 1e4);
            p->free = range(0, 1e4);

            p = newParameter("phiM");
            p->title = "#it{#phi}^{#it{M}}_{2} [rad]";
            p->startvalue = 0.02;
            p->unit = "";
            p->scan = range(-0.5, 0.5);
            p->phys = range(-1.5, -1.5 + DegToRad(360.));
            p->free = range(-1.5, -1.5 + DegToRad(360.));
            break;

        default:
            cout << "ParametersCharmCombo::defineParameters : ERROR : "
                    "theory_config " + to_string(th_cfg) + " not found." << endl;
            exit(1);
    }
    if (predict_dy_rs) {
        p = newParameter("DY_RS");
        p->title = "#||{#Delta#it{Y}_{#it{K}^{#minus}#pi^{+}}}  [10^{#minus4}]";
        p->startvalue = 0.05;
        p->unit = "";
        p->scan = range(0., 1.);
        p->phys = range(0., 1e4);
        p->free = range(0., 1e4);
    }
}
