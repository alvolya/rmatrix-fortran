# Benchmark against AZURE2: p + 14N elastic scattering (15O)

The same R-matrix parameters evaluated by AZURE2 and by Rmatrix.  Any
difference is a difference between the two codes.  Everything here is
derived from public material only:

| item | value |
|------|-------|
| AZURE2 | https://github.com/AZURE2DEV/AZURE2, commit `9441e3d3c82f49b7bdd02c995f91dca1cc73b18b` (2026-10-03) |
| Python binding | pyazr 2.8.0 (part of that checkout, `USE_API=ON`) |
| input | `tests/thick_target_resonance/thick_target_resonance.azr` of that checkout: p + 14N, 33 levels in seven J^π groups (1/2±, 3/2±, 5/2±, 7/2+), channel spins 1/2 and 3/2, L ≤ 3, channel radius 5.5 fm, six photon pairs |
| generator | `make_case.py <AZURE2 checkout>` (needs pyazr); rewrites every file in this directory |

## What the generator does

1. Copies the input and releases all level energies and widths, so that
   pyazr exposes AZURE2's internal parameter vector; sets the widths of
   the photon channels to zero (this code has no photon channels; the
   test is about particle channels); asks AZURE2 for its
   boundary-condition check file; replaces `<segmentsTest>` by one
   differential elastic segment per angle (60°, 90°, 140° c.m.) on a
   uniform grid of 12 c.m. energies from 0.1 to 3.4 MeV above the
   p + 14N threshold.  The modified input is kept as
   `azure2_input_used.azr`.
2. Runs AZURE2 in extrapolation mode **without** the Brune
   parametrisation.  Its internal parameters are then formal R-matrix
   parameters (E_λ, γ_λc) for the boundary condition B_c = S_c(E_1),
   E_1 the energy of the first level of each J^π group
   (`CNuc::CalcBoundaryConditions`).  A bound level that AZURE2 keeps
   fixed has no vector entry; its transformed energy is read from
   `level_energy`.
3. Writes those parameters as `model.inp` / `model.par` (energies
   relative to the p + 14N threshold, 7.297 MeV; `BcMode = 1` makes
   Rmatrix compute the same B_c internally) and AZURE2's cross sections,
   the energies and angles it actually used, and its boundary-condition
   table into `points.json`.

## What the test checks (`tests/run_tests.py`, group 4)

- dσ/dΩ of Rmatrix at every (E, θ) against AZURE2 to 2e-6 relative
  (observed: 4e-7, AZURE2's own numerical floor);
- the B_c table Rmatrix prints under `BcMode = 1` against AZURE2's
  check file (6 digits).

Rmatrix needs no AZURE2 to run the test; AZURE2 is needed only to
regenerate the files.  A note on AZURE2's two formalisms: run on this
input with and without the Brune parametrisation, AZURE2 gives the same
cross sections to 1e-13, so the choice above is one of convenience.
