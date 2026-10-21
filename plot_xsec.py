#!/usr/bin/env python3
# Rmatrix — multilevel, multichannel R-matrix code
# Copyright (C) 2026 Alexander Volya
# SPDX-License-Identifier: GPL-3.0-or-later
"""
plot_xsec.py — Plot Rmatrix differential cross section output.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
GETTING MULTIPLE ANGLES FROM Rmatrix
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Edit the &RMATRIX namelist in your .inp file:

    Astart = 20.0,  Afinal = 180.0,  Astep = 20.0,

This produces angles 20°, 40°, 60°, …, 180°.
By default this script creates one figure per angle (excitation function).
Use --overlay to put all angles on one figure.
Use --angular-dist to plot σ(θ) at one or more fixed energies instead.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
EXPERIMENTAL DATA FORMAT (--data)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Two accepted column layouts (whitespace- or comma-separated):

  2-column:  Ecm(MeV)  xsec(mb/sr)
             → plotted on every angle panel, or use with --angle

  3-column:  Ecm(MeV)  angle(deg)  xsec(mb/sr)
             → matched to the correct angle panel automatically

  4-column:  Ecm(MeV)  xsec(mb/sr)  err_lo  err_hi
             → same as 2-col but with asymmetric error bars
             (symmetric errors: use the same value for err_lo and err_hi)

  5-column:  Ecm(MeV)  angle(deg)  xsec(mb/sr)  err_lo  err_hi
             → 3-col with asymmetric error bars

Lines beginning with # or ! are treated as comments.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
USAGE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  python plot_xsec.py <file> [file2 ...] [options]

OPTIONS
  --save <file>        Save figure(s) to file (.png/.pdf/.svg).
                       When multiple figures are produced the angle or
                       energy is embedded in the filename automatically.
  --log                Logarithmic y-axis
  --angle <deg>        Plot only this angle (default: all angles)
  --overlay            Put all angles from each file on one figure
  --angular-dist       Plot σ(θ) at fixed energies instead of σ(E)
  --energy <MeV>       Energy (or energies) for --angular-dist;
                       use the nearest grid point if not exact.
                       Repeat flag or comma-separate: --energy 3.5,4.0
  --labels <l1,l2,...> Legend labels (one per theory file)
  --data <file>        Experimental data file (see format above)
  --data-label <str>   Legend label for experimental data (default: "data")
  --no-legend          Suppress legend

EXAMPLES
  # Single angle, no data
  python plot_xsec.py 18O_dp.output

  # Compare two calculations, log scale
  python plot_xsec.py 18O.output 18O_dp.output --labels "SP,DP" --log

  # All angles on one figure
  python plot_xsec.py 18O_dp.output --overlay

  # Angular distribution at two energies
  python plot_xsec.py 18O_dp.output --angular-dist --energy 3.0,4.0

  # With experimental data, save PDF
  python plot_xsec.py 18O_dp.output --data exp.dat --data-label "Exp. (2020)" --save xsec.pdf

  # Single angle with data
  python plot_xsec.py 18O_dp.output --angle 180 --data exp.dat
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


# ── data loading ────────────────────────────────────────────────────────────

def load_output(path):
    """Load an Rmatrix .output/.out file.
    Returns dict: angle (float) -> (energies array, xsec array)
    """
    data = defaultdict(lambda: ([], []))
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(("!", "#")):
                continue
            parts = line.replace(",", " ").split()
            if len(parts) != 3:
                continue
            try:
                ecm, angle, xsec = float(parts[0]), float(parts[1]), float(parts[2])
            except ValueError:
                continue
            data[angle][0].append(ecm)
            data[angle][1].append(xsec)
    return {angle: (np.array(e), np.array(x)) for angle, (e, x) in data.items()}


def load_data(path):
    """Load experimental data file.
    Returns dict: angle (float or None) -> (ecm, xsec, err_lo, err_hi)
    angle=None means no angle column was present.
    """
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(("!", "#")):
                continue
            parts = line.replace(",", " ").split()
            try:
                vals = [float(v) for v in parts]
            except ValueError:
                continue
            rows.append(vals)

    if not rows:
        return {}

    ncol = len(rows[0])
    result = defaultdict(lambda: ([], [], [], []))

    for vals in rows:
        n = len(vals)
        if n == 2:                          # Ecm  xsec
            ecm, xsec, elo, ehi = vals[0], vals[1], 0.0, 0.0
            angle = None
        elif n == 3:                        # Ecm  angle  xsec
            ecm, angle, xsec, elo, ehi = vals[0], vals[1], vals[2], 0.0, 0.0
        elif n == 4:                        # Ecm  xsec  elo  ehi
            ecm, xsec, elo, ehi = vals
            angle = None
        elif n == 5:                        # Ecm  angle  xsec  elo  ehi
            ecm, angle, xsec, elo, ehi = vals
        else:
            continue
        result[angle][0].append(ecm)
        result[angle][1].append(xsec)
        result[angle][2].append(elo)
        result[angle][3].append(ehi)

    return {
        a: (np.array(e), np.array(x), np.array(el), np.array(eh))
        for a, (e, x, el, eh) in result.items()
    }


# ── plot helpers ─────────────────────────────────────────────────────────────

def nearest_energy(energies, target):
    idx = np.argmin(np.abs(energies - target))
    return energies[idx], idx


def add_exp(ax, exp_data, angle, color="black", label="data"):
    """Overlay experimental data on ax for the given angle (or None-keyed data)."""
    keys_to_try = [angle, None]
    plotted = False
    for key in keys_to_try:
        if key not in exp_data:
            continue
        ecm, xsec, elo, ehi = exp_data[key]
        has_errors = np.any(elo > 0) or np.any(ehi > 0)
        kw = dict(fmt="o", color=color, ms=4, lw=1.0, capsize=3,
                  label=label if not plotted else "_nolegend_")
        if has_errors:
            ax.errorbar(ecm, xsec, yerr=[elo, ehi], **kw)
        else:
            ax.plot(ecm, xsec, "o", color=color, ms=4, label=kw["label"])
        plotted = True
        break
    return plotted


def style_ax(ax, xlabel, ylabel, title, log=False):
    ax.set_xlabel(xlabel, fontsize=13)
    ax.set_ylabel(ylabel, fontsize=13)
    ax.set_title(title, fontsize=13)
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.grid(True, linestyle="--", alpha=0.4)
    if log:
        ax.set_yscale("log")


def save_or_show(fig, save_base, suffix, args):
    if args.save:
        p = Path(save_base)
        if suffix:
            p = p.with_stem(f"{p.stem}_{suffix}")
        fig.savefig(p, dpi=150)
        print(f"Saved: {p}")
    else:
        plt.show()
    plt.close(fig)


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Plot Rmatrix differential cross section output.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("files", nargs="+", help="Rmatrix output file(s)")
    parser.add_argument("--save", metavar="FILE")
    parser.add_argument("--log", action="store_true")
    parser.add_argument("--angle", type=float, default=None)
    parser.add_argument("--overlay", action="store_true",
                        help="All angles on one figure (excitation functions)")
    parser.add_argument("--angular-dist", action="store_true", dest="angdist",
                        help="Plot dσ/dΩ vs θ at fixed energies")
    parser.add_argument("--energy", metavar="E1,E2,...",
                        help="Energies for --angular-dist (comma-sep or repeat)")
    parser.add_argument("--labels", metavar="L1,L2,...")
    parser.add_argument("--data", metavar="FILE", help="Experimental data file")
    parser.add_argument("--data-label", metavar="STR", default="data", dest="data_label")
    parser.add_argument("--no-legend", action="store_true", dest="no_legend")
    args = parser.parse_args()

    labels = args.labels.split(",") if args.labels else None
    if labels and len(labels) != len(args.files):
        print("Error: --labels count must match number of files", file=sys.stderr)
        sys.exit(1)

    # Load theory
    datasets = []
    for i, path in enumerate(args.files):
        label = labels[i] if labels else Path(path).stem
        datasets.append((label, load_output(path)))

    # Load experimental data
    exp_data = load_data(args.data) if args.data else {}

    colors = list(cm.tab10.colors)
    exp_color = "black"

    all_angles = sorted({a for _, d in datasets for a in d})

    # ── angular distribution mode ──────────────────────────────────────────
    if args.angdist:
        if not args.energy:
            # pick a few representative energies from the middle of the range
            ref_ecm = next(iter(next(iter(datasets))[1].values()))[0]
            n = len(ref_ecm)
            idxs = [n // 4, n // 2, 3 * n // 4]
            target_energies = ref_ecm[idxs]
        else:
            target_energies = np.array([float(e) for e in args.energy.replace(",", " ").split()])

        fig, ax = plt.subplots(figsize=(8, 5))

        for ei, etgt in enumerate(target_energies):
            for fi, (label, data) in enumerate(datasets):
                if not data:
                    continue
                ref_ecm = next(iter(data.values()))[0]
                eactual, eidx = nearest_energy(ref_ecm, etgt)
                theta = np.array(sorted(data.keys()))
                xsec = np.array([data[a][1][eidx] for a in theta])

                n_curves = len(target_energies) * len(datasets)
                cidx = ei * len(datasets) + fi
                col = colors[cidx % len(colors)]

                lbl = f"{label}, {eactual:.3f} MeV" if len(datasets) > 1 else f"{eactual:.3f} MeV"
                ax.plot(theta, xsec, lw=1.5, color=col, label=lbl)

        # experimental angular distribution data
        if exp_data:
            # group by energy if 3-col, otherwise plot as-is treating col1 as angle
            for key in exp_data:
                ecm_e, xsec_e, elo, ehi = exp_data[key]
                has_err = np.any(elo > 0) or np.any(ehi > 0)
                lbl = args.data_label
                if has_err:
                    ax.errorbar(ecm_e, xsec_e, yerr=[elo, ehi],
                                fmt="o", color=exp_color, ms=4, capsize=3, label=lbl)
                else:
                    ax.plot(ecm_e, xsec_e, "o", color=exp_color, ms=4, label=lbl)

        style_ax(ax,
                 r"$\theta_\mathrm{cm}$ (deg)",
                 r"$d\sigma/d\Omega$ (mb/sr)",
                 "Angular distribution",
                 log=args.log)
        if not args.no_legend:
            ax.legend(fontsize=10)
        fig.tight_layout()
        save_or_show(fig, args.save or "angdist.pdf", "", args)
        return

    # ── excitation function mode ────────────────────────────────────────────
    plot_angles = [args.angle] if args.angle is not None else all_angles

    if args.overlay:
        # All angles on one figure, one line per (file × angle)
        fig, ax = plt.subplots(figsize=(10, 5))
        curve_idx = 0
        for fi, (label, data) in enumerate(datasets):
            for angle in plot_angles:
                if angle not in data:
                    continue
                ecm, xsec = data[angle]
                col = colors[curve_idx % len(colors)]
                lbl = f"{label}, {angle:.0f}°" if len(datasets) > 1 else f"{angle:.0f}°"
                ax.plot(ecm, xsec, lw=1.5, color=col, label=lbl)
                curve_idx += 1

        if exp_data:
            add_exp(ax, exp_data, angle=None, color=exp_color, label=args.data_label)

        style_ax(ax,
                 r"$E_\mathrm{cm}$ (MeV)",
                 r"$d\sigma/d\Omega$ (mb/sr)",
                 "Differential cross sections",
                 log=args.log)
        if not args.no_legend:
            ax.legend(fontsize=9, ncol=max(1, curve_idx // 8))
        fig.tight_layout()
        save_or_show(fig, args.save or "xsec_overlay.pdf", "", args)

    else:
        # One figure per angle
        for angle in plot_angles:
            if angle not in all_angles:
                print(f"Warning: angle {angle}° not found in data", file=sys.stderr)
                continue

            fig, ax = plt.subplots(figsize=(9, 5))

            for fi, (label, data) in enumerate(datasets):
                if angle not in data:
                    continue
                ecm, xsec = data[angle]
                ax.plot(ecm, xsec, lw=1.5, color=colors[fi % len(colors)], label=label)

            if exp_data:
                add_exp(ax, exp_data, angle=angle, color=exp_color, label=args.data_label)

            show_legend = not args.no_legend and (
                len(datasets) > 1 or bool(exp_data) or (labels and labels[0])
            )
            style_ax(ax,
                     r"$E_\mathrm{cm}$ (MeV)",
                     r"$d\sigma/d\Omega$ (mb/sr)",
                     rf"$d\sigma/d\Omega$ at $\theta_\mathrm{{cm}} = {angle:.1f}°$",
                     log=args.log)
            if show_legend:
                ax.legend(fontsize=11)
            fig.tight_layout()

            suffix = f"{int(angle)}deg" if len(plot_angles) > 1 else ""
            save_or_show(fig, args.save or "xsec.pdf", suffix, args)


if __name__ == "__main__":
    main()
