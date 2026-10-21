#!/usr/bin/env python3
# Rmatrix — multilevel, multichannel R-matrix code
# Copyright (C) 2026 Alexander Volya
# SPDX-License-Identifier: GPL-3.0-or-later
"""Regression tests for Rmatrix on the 18O example.

Groups
  1. build          make (skipped with --exe)
  2. 18O single     tests/18O/18O_sa.inp  vs tests/18O/refs/18O_sa.output
                    (180 deg, 1.0-6.0 MeV: sub-threshold energies for the
                    n+17O channels, bound 0+ level)
  3. 18O multiangle tests/18O/18O_multiangle.inp  vs tests/18O/refs/18O_multiangle.output
  4. AZURE2 15O     tests/azure2_15O/: p + 14N elastic scattering on the parameters of
                    AZURE2's own test input, compared with AZURE2's cross sections
                    (see tests/azure2_15O/README.md); also checks that BcMode = 1
                    reproduces AZURE2's boundary conditions
                    (20-180 deg in steps of 20, 1.0-6.0 MeV)

Both use tests/18O/18O.par (the example's parameter file; tests/18O/18O.inp
is the full 10-180 degree example).  Criterion: every point outside
1.5 keV of a level energy agrees with the reference to REL relative.
The references were produced by this code; differences can come only
from compiler and library rounding.

Usage
  python3 tests/run_tests.py                # build, then test
  python3 tests/run_tests.py --exe ./Rmatrix
  python3 tests/run_tests.py --update-refs  # rewrite refs/ from the current build

tests/18O/refs-frozen-2026/ holds the outputs of the as-received 2026-08 code
(18O_sa, 18O_multiangle) and of the 2026-03 version (18O_dp_2026-03,
2.0-6.0 MeV at 180 deg) for the record; they are not compared here.
"""
import argparse
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CASE = HERE / "18O"
REL = 1e-6
REL_AZURE2 = 2e-6
BC_TOL = 2e-5                      # AZURE2 prints its boundary conditions with 6 digits
POLE_MASK_MEV = 0.0015


def level_energies(par: Path) -> np.ndarray:
    out = []
    for line in par.read_text().splitlines()[4:]:
        t = line.split()
        if len(t) == 3 and not re.fullmatch(r"-?\d+", t[0]):
            out.append(float(t[0]))
    return np.asarray(out)


def run(exe: Path, inp: Path, par: Path, out: Path, want_stdout: bool = False):
    """Rows (E, theta, dsigma/dOmega) of one run (and its stdout)."""
    res = subprocess.run([str(exe), "-i", str(inp), "-p", str(par), "-o", str(out)],
                         capture_output=True, text=True, cwd=str(out.parent), timeout=1200)
    if res.returncode != 0:
        sys.exit(f"Rmatrix failed on {inp.name}:\n{res.stdout}\n{res.stderr}")
    rows = np.loadtxt(out, comments="#", ndmin=2)[:, :3]
    return (rows, res.stdout) if want_stdout else rows


def group_azure2(exe: Path, work: Path, case: Path) -> list:
    meta = json.loads((case / "points.json").read_text())
    inp_text = (case / "model.inp").read_text()
    results = []
    bc_checked = False
    for p in meta["points"]:
        ang = p["angle"]
        txt = re.sub(r"(Astart\s*=\s*)[-+0-9.eEdD]+", rf"\g<1>{ang:g}", inp_text)
        txt = re.sub(r"(Afinal\s*=\s*)[-+0-9.eEdD]+", rf"\g<1>{ang:g}", txt)
        inp = work / f"{case.name}_{ang:g}.inp"
        inp.write_text(txt)
        new, log = run(exe, inp, case / "model.par", work / f"{case.name}_{ang:g}.output", True)
        ecm = np.asarray(p["ecm"])
        a2 = np.asarray(p["azure2_b_sr"]) * 1000.0                   # b/sr -> mb/sr
        idx = np.array([int(np.argmin(np.abs(new[:, 0] - e))) for e in ecm])
        if np.any(np.abs(new[idx, 0] - ecm) > 1e-6):
            print(f"    {ang:g} deg: energy grid mismatch")
            results.append(False)
            continue
        dev = np.abs(new[idx, 2] - a2) / a2
        ok = bool(np.all(dev <= REL_AZURE2))
        print(f"    {meta['reaction']} {ang:5g} deg: {len(ecm)} pts, max rel dev vs AZURE2 "
              f"{dev.max():.2e} -> {'PASS' if ok else 'FAIL'}")
        results.append(ok)
        if not bc_checked:                       # BcMode = 1 table vs AZURE2's check file
            bc_checked = True
            want = meta.get("boundary_conditions_by_jpi", {})
            bad = 0
            n = 0
            for m_ in re.finditer(r"J =\s*([0-9.]+)\s+parity\s*(-?1).*?B_c =\s*(-?[0-9.]+)", log):
                J, par, b = float(m_.group(1)), int(m_.group(2)), float(m_.group(3))
                key = f"{J:g}{'+' if par > 0 else '-'}"
                n += 1
                if key in want and min(abs(b - v) for v in want[key]) > BC_TOL:
                    bad += 1
                    print(f"    B_c mismatch for {key}: {b} not in {want[key]}")
            print(f"    boundary conditions B_c = S_c(E_1): {n} rows checked against AZURE2 "
                  f"-> {'PASS' if bad == 0 and n > 0 else 'FAIL'}")
            results.append(bad == 0 and n > 0)
    return results


def compare(new: np.ndarray, ref: np.ndarray, poles: np.ndarray, label: str) -> bool:
    if new.shape != ref.shape or not np.allclose(new[:, :2], ref[:, :2], atol=1e-6):
        print(f"    {label}: grid differs from the reference ({new.shape} vs {ref.shape})")
        return False
    mask = np.ones(len(ref), bool)
    if poles.size:
        mask = np.min(np.abs(ref[:, :1] - poles[None, :]), axis=1) > POLE_MASK_MEV
    r = ref[mask, 2]
    n = new[mask, 2]
    ok = r != 0
    dev = np.abs(n[ok] - r[ok]) / np.abs(r[ok])
    worst = int(np.argmax(dev))
    passed = bool(np.all(dev <= REL))
    print(f"    {label}: {int(ok.sum())} pts ({int((~mask).sum())} pole-adjacent excluded), "
          f"max rel dev {dev.max():.2e} at E={ref[mask][ok][worst, 0]:.4f} "
          f"th={ref[mask][ok][worst, 1]:.1f} -> {'PASS' if passed else 'FAIL'}")
    return passed


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", type=Path, help="existing executable (skips the build)")
    ap.add_argument("--update-refs", action="store_true",
                    help="rewrite tests/refs from this build instead of comparing")
    a = ap.parse_args()

    if a.exe:
        exe = a.exe.resolve()
        print(f"[1/4] build — skipped (--exe {exe})")
    else:
        print("[1/4] build")
        res = subprocess.run(["make"], cwd=str(ROOT), capture_output=True, text=True)
        if res.returncode != 0:
            print(res.stdout[-2000:], res.stderr[-2000:])
            return 1
        exe = ROOT / "Rmatrix"

    work = Path(tempfile.mkdtemp(prefix="rmatrix_tests_"))
    par = CASE / "18O.par"
    poles = level_energies(par)
    results = []
    for k, case in enumerate(("18O_sa", "18O_multiangle"), start=2):
        print(f"[{k}/4] {case}")
        inp = CASE / f"{case}.inp"
        ref = CASE / "refs" / f"{case}.output"
        out = work / f"{case}.output"
        new = run(exe, inp, par, out)
        if a.update_refs:
            (CASE / "refs").mkdir(exist_ok=True)
            ref.write_text(out.read_text())
            print(f"    reference written ({len(new)} points)")
            results.append(True)
        elif not ref.exists():
            print(f"    no reference {ref}")
            results.append(False)
        else:
            results.append(compare(new, np.loadtxt(ref, comments="#", ndmin=2)[:, :3],
                                   poles, case))

    print("[4/4] azure2_15O (p + 14N, AZURE2's own test input)")
    results += group_azure2(exe, work, HERE / "azure2_15O")

    print()
    if not all(results):
        print(f"FAILED: {results.count(False)} check(s)")
        return 1
    print(f"ALL TESTS PASSED ({len(results)} checks)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
