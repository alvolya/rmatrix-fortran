#!/usr/bin/env python3
# Rmatrix — multilevel, multichannel R-matrix code
# Copyright (C) 2026 Alexander Volya
# SPDX-License-Identifier: GPL-3.0-or-later
"""make_case.py — build the p + 14N (15O) benchmark case from AZURE2's own
test input, so that AZURE2 and Rmatrix evaluate the SAME R-matrix parameters.

    python3 make_case.py /path/to/AZURE2-checkout [--out DIR]

Needs the AZURE2 checkout (https://github.com/AZURE2DEV/AZURE2) built with
its Python binding pyazr.  Everything written here is derived from the file
tests/thick_target_resonance/thick_target_resonance.azr of that checkout:

1. a scratch copy of the input is prepared: all level energies and widths
   are released (so that pyazr exposes the internal parameter vector),
   the widths of the photon channels are set to zero (this code has no
   photon channels), the boundary-condition check file is requested, and
   <segmentsTest> is replaced by one differential elastic segment per
   angle on a uniform c.m. energy grid;
2. AZURE2 runs it in extrapolation mode WITHOUT the Brune parametrisation
   (use_brune=False): its internal parameters are then formal R-matrix
   parameters with the boundary condition B_c = S_c(E_1), E_1 the energy
   of the first level of each J^pi group;
3. those internal parameters are written as model.inp / model.par
   (BcMode = 1 reproduces AZURE2's B_c inside Rmatrix), AZURE2's cross
   sections and its boundary-condition table go to points.json, and the
   scratch input is kept as azure2_input_used.azr.
"""
import argparse
import json
import re
import shutil
import subprocess
import sys
import tempfile
from collections import OrderedDict
from pathlib import Path

import numpy as np

INPUT = "tests/thick_target_resonance/thick_target_resonance.azr"
ANGLES = (60.0, 90.0, 140.0)             # c.m. angles
EMIN, EMAX, NPTS = 0.1, 3.4, 12          # c.m. energies above the p + 14N threshold
DIFF_CM = 5                              # AZURE2 segment type: differential, c.m. angle


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("azure2", help="AZURE2 checkout (built with pyazr)")
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent))
    a = ap.parse_args()
    import pyazr
    try:
        from importlib.metadata import version
        pyazr_version = version("pyazr")
    except Exception:                       # noqa: BLE001
        pyazr_version = getattr(pyazr, "__version__", "unknown")

    repo = Path(a.azure2).resolve()
    src = repo / INPUT
    commit = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"],
                            capture_output=True, text=True).stdout.strip()
    out = Path(a.out).resolve()
    out.mkdir(parents=True, exist_ok=True)

    work = Path(tempfile.mkdtemp(prefix="azure2_15O_"))
    case = work / "case"
    shutil.copytree(src.parent, case)
    for d in ("output", "checks"):
        (case / d).mkdir(exist_ok=True)
    azr = case / src.name

    # --- 1. scratch input -------------------------------------------------
    model = pyazr.AzrModel.from_file(str(azr))
    file_levels = [(norm_jpi(lv.jpi), lv.energy) for lv in model.levels]   # file order
    for lv in model.levels:
        lv.set_fixed(False)
        for ch in lv.channels:
            ch.channel_fixed = False
            if ch.is_photon:
                ch.gamma = 0.0
    model.write(str(azr))
    txt = azr.read_text()
    txt, n = re.subn(r"^none(\s+#Boundary Condition Check)", r"file\1", txt, flags=re.M)
    assert n == 1, "Boundary Condition Check line not found"
    lines = txt.splitlines(keepends=True)
    with pyazr.azure2(str(azr), cwd=str(case), data_mode=True, use_brune=False) as s0:
        pr = next(p for p in s0.pairs if p.is_entrance)
    lab = (pr.M1 + pr.M2) / pr.M2                   # E_lab = E_cm (m1+m2)/m2
    step = (EMAX - EMIN) / (NPTS - 1)
    i0, i1 = lines.index("<segmentsTest>\n"), lines.index("</segmentsTest>\n")
    segs = [f"1 {pr.key} {pr.key} {EMIN * lab!r} {EMAX * lab!r} {step * lab!r} "
            f"{ang!r} {ang!r} 0.0 {DIFF_CM}\n" for ang in ANGLES]
    lines = lines[:i0 + 1] + segs + lines[i1:]
    if "<targetInt>\n" in lines:                     # target integration applies by segment index
        t0, t1 = lines.index("<targetInt>\n"), lines.index("</targetInt>\n")
        lines = lines[:t0 + 1] + lines[t1:]
    azr.write_text("".join(lines))

    # --- 2. AZURE2 ----------------------------------------------------------
    with pyazr.azure2(str(azr), cwd=str(case), data_mode=False, use_brune=False) as s:
        pairs = list(s.pairs)
        params = list(s.parameters)
        x = np.asarray(s.params_rwa, float)
        ecm = [np.asarray(v, float) for v in s.calculate_energies(x)]
        th = [np.asarray(v, float) for v in s.calculate_angles(x)]
        sig = [np.asarray(v, float) for v in s.calculate_rwa(x)]
    chk = (case / "checks" / "boundaryconditions.chk").read_text()
    bc_rows = [(int(t[0]), int(t[1]), float(t[2])) for t in
               (l.split() for l in chk.splitlines()) if len(t) == 3 and t[0].isdigit()]

    ent = next(p for p in pairs if p.is_entrance)
    assert ent.ptype == 0
    thr = ent.sep_energy + ent.excitation
    nfile = {}
    for jpi, _e in file_levels:
        nfile[jpi] = nfile.get(jpi, 0) + 1

    # --- 3. levels: (jpi, level) -> energy, {(s, L): gamma} ----------------
    levels: "OrderedDict[tuple, dict]" = OrderedDict()
    for q in params:
        kind = getattr(q, "kind", "")
        if kind not in ("energy", "width") or q.level > nfile.get(q.jpi, 0):
            continue                                  # final-state levels AZURE2 adds for photon pairs
        lev = levels.setdefault((q.jpi, q.level), {"E": None, "rows": {}})
        if kind == "energy":
            # a released level that AZURE2 keeps fixed anyway (bound states
            # held at their excitation energy) has no entry in the vector;
            # its transformed energy is level_energy, NOT the observed value
            if q.free_index is not None:
                e = x[q.free_index]
            elif q.level_energy is not None:
                e = q.level_energy
                print(f"note: {q.jpi} level {q.level}: energy fixed, using level_energy {e:.6f}")
            else:
                e = q.value
            lev["E"] = float(e)
        else:
            if lev["E"] is None and q.level_energy is not None:
                lev["E"] = float(q.level_energy)
            if q.pair != ent.number or q.radiation_type != "P":
                continue
            g = x[q.free_index] if q.free_index is not None else 0.0
            lev["rows"][(float(q.S), int(q.L))] = float(g)
    for key, lev in levels.items():
        if lev["E"] is None:                          # fixed ground state without an energy parameter
            jpi, k = key
            lev["E"] = [e for j, e in file_levels if j == jpi][k - 1]
    spins = sorted({s for lev in levels.values() for (s, _L) in lev["rows"]})
    chn = {s: i + 1 for i, s in enumerate(spins)}
    Lmax = max(L for lev in levels.values() for (_s, L) in lev["rows"])

    # --- .inp -----------------------------------------------------------------
    ecm0 = ecm[0]
    dstep = float(np.median(np.diff(ecm0)))
    cards = ["! Channel cards: Ap (u)   Zp    Sp     At (u)    Zt    St     Q(MeV)  ChSp  ChR(fm)",
             f"! p + 14N from AZURE2 {INPUT} (commit {commit[:12]}), pyazr {pyazr_version}"]
    for i, s_ in enumerate(spins):
        Ap = ent.M1 if i < len(spins) - 1 else -ent.M1
        cards.append(f"  {Ap:10.5f} {ent.Z1:5.1f} {ent.J1 * ent.parity1:6.2f} "
                     f"{ent.M2:10.5f} {ent.Z2:5.1f} {ent.J2 * ent.parity2:6.2f} "
                     f"{0.0:12.6f} {s_:5.1f} {ent.channel_radius:8.4f}")
    inp = "\n".join(cards + [
        "", "&RMATRIX",
        f"  Lmax    = {Lmax},",
        f"  Estart  = {ecm0[0]:.12g},  Efinal = {ecm0[-1] + 1.5 * dstep:.12g},  Estep = {dstep:.12g},",
        f"  Astart  = {ANGLES[0]:g},  Afinal = {ANGLES[0]:g},  Astep = 1.0,",
        "  ExitCh  = 1,", "  NRDef   = 0,", "  DKEY    = 0,", "  QEx     = 0.0,",
        "  ecmlow  = 0.0,  ecmhigh = 0.0,",
        "  imphase = 0,  APA  = -1.0,  APAL = 0.3,  APEP = 3.0,  APAP = 0.3,",
        "  rdchan  = 1,",
        "  BcMode  = 1,          ! B_c = S_c(E_1): AZURE2's boundary condition",
        "/", ""])
    (out / "model.inp").write_text(inp)

    # --- .par -----------------------------------------------------------------
    plines = ["; Resonance parameters: Eres, Jtotal, Parity   (c.m. energy above the p + 14N threshold)",
              ";                       ChN, L, Bc, G          (Bc ignored: BcMode = 1 in model.inp)",
              f"; AZURE2 internal (formal) parameters of {INPUT}, commit {commit[:12]}, photon widths zeroed",
              ";                      -ChN  - last row of the level"]
    order = sorted(levels, key=lambda k: levels_index(file_levels, k))   # file order
    for key in order:
        lev = levels[key]
        jpi = key[0]
        J = float(jpi[:-1])
        par = 1.0 if jpi.endswith("+") else -1.0
        rows = sorted(lev["rows"].items(), key=lambda kv: (chn[kv[0][0]], kv[0][1]))
        plines.append(f"  {lev['E'] - thr:.17g}  {J:6.3f}  {par:6.3f}")
        for j, ((s_, L), g) in enumerate(rows):
            sign = -1 if j == len(rows) - 1 else 1
            plines.append(f"  {sign * chn[s_]:4d} {L:3d} {0.0:.3f} {g:.17g}")
    (out / "model.par").write_text("\n".join(plines) + "\n")

    # --- points.json ------------------------------------------------------------
    jgroups = list(OrderedDict.fromkeys(j for j, _e in file_levels))      # AZURE2 J-group order
    bc = {}
    for jg, ch, v in bc_rows:
        bc.setdefault(jgroups[jg - 1] if jg - 1 < len(jgroups) else str(jg), []).append(v)
    meta = dict(
        system="15O", reaction="p14N->p14N", entrance="p + 14N", threshold_MeV=thr,
        azure2_repository="https://github.com/AZURE2DEV/AZURE2", azure2_commit=commit,
        azure2_input=INPUT, pyazr_version=pyazr_version, use_brune=False,
        photon_widths="set to zero", angles_cm=list(ANGLES),
        boundary_conditions_by_jpi={k: sorted(set(round(v, 6) for v in vs)) for k, vs in bc.items()},
        points=[dict(angle=ang, ecm=e.tolist(), theta=t.tolist(), azure2_b_sr=sg.tolist())
                for ang, e, t, sg in zip(ANGLES, ecm, th, sig)])
    (out / "points.json").write_text(json.dumps(meta, indent=1))
    shutil.copy2(azr, out / "azure2_input_used.azr")
    shutil.rmtree(work, ignore_errors=True)
    print(f"{out}: {len(levels)} levels, Lmax {Lmax}, channel spins {spins}, "
          f"{len(ecm0)} energies x {len(ANGLES)} angles; AZURE2 {commit[:12]}, pyazr {pyazr_version}")
    return 0


def norm_jpi(label: str) -> str:
    """'1/2-' or '0.5-' -> '0.5-' (the engine's parameter labels)."""
    sign = label[-1]
    body = label[:-1]
    if "/" in body:
        n, d = body.split("/")
        J = float(n) / float(d)
    else:
        J = float(body)
    return f"{J:g}{sign}"


def levels_index(file_levels, key):
    """Position in file order of the k-th level of J^pi `key[0]`."""
    jpi, k = key
    n = 0
    for i, (j, _e) in enumerate(file_levels):
        if j == jpi:
            n += 1
            if n == k:
                return i
    raise KeyError(key)


if __name__ == "__main__":
    sys.exit(main())
