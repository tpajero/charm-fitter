# !/usr/bin/env python3
"""Extrapolate the precision of a measurement to a future luminosity, taking into account possible improvements in
trigger and selection efficiency.
"""

import argparse

INVFB = 1.0

LUMI = {
    "LHCb-run2": 6 * INVFB,
    "LHCb-UI": 50 * INVFB,
    "LHCb-UII": 300 * INVFB,
}


def extrapolate_precision(precision: float, lumi_old: float, lumi_new: float, trigger_improvement: float) -> float:
    return precision * (lumi_old / (lumi_new * trigger_improvement)) ** 0.5


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_value", type=float, help="Input precision to extrapolate")
    parser.add_argument("--lumi-old", type=str, default="LHCb-run2", help="Old luminosity ID")
    parser.add_argument("--lumi-new", type=str, default="LHCb-UII", help="New luminosity ID")
    parser.add_argument("--lumi-old-val", type=float, default=None, help="Old luminosity")
    parser.add_argument("--lumi-new-val", type=float, default=None, help="New luminosity")
    parser.add_argument(
        "--trigger-improvement",
        type=float,
        default=2.0,
        help="Improvement in yield per invfb (due eg to removal of L0 trigger between LHCb Run 2 and 3)",
    )
    args = parser.parse_args()

    lumi_old = LUMI[args.lumi_old] if args.lumi_old_val is None else args.lumi_old_val
    lumi_new = LUMI[args.lumi_new] if args.lumi_new_val is None else args.lumi_new_val
    precision_new = extrapolate_precision(args.input_value, lumi_old, lumi_new, args.trigger_improvement)

    print(
        f"Extrapolated precision: {args.input_value:.3e} -> {precision_new:.3e} (from {lumi_old / INVFB:.1f} -> {lumi_new / INVFB:.1f} fb-1, trigger improvement of {args.trigger_improvement:.1f})"
    )
