# Example: alpha + 14C / n + 17O (18O compound nucleus)

Five channel-spin groups (alpha + 14C; n + 17O with s = 2, 3; n + 17O*
with s = 0, 1), Lmax = 8, 27 levels including a bound 0+ level at
-2.0 MeV.  Energies 1.0-6.0 MeV (c.m., relative to the alpha + 14C
threshold; the neutron channels open at 1.8173 and 2.6880 MeV), angles
10-180 degrees.  The parameter set is illustrative (G.V. Rogachev); the
two triangle-forbidden widths of the bound level in the original file
are set to zero.

```sh
../../Rmatrix -i 18O.inp -p 18O.par -o 18O.out      # full example, 10-180 degrees
python3 ../../plot_xsec.py 18O.out --angle 180
```

Writes `18O.out` (E, angle, dsigma/dOmega in mb/sr), `18O.total`
(angle-integrated nuclear term and S-factor) and `18O.screen`.
`NRDef = 10, rdchan = 2` also runs the resonance finder on level 10 in
the n + 17O channel and prints its observed energy and width.

The regression test (`tests/run_tests.py`) runs the same parameter file
with `18O_sa.inp` (180 degrees) and `18O_multiangle.inp` (20-180 degrees)
against `refs/`; `refs-frozen-2026/` holds the outputs of the as-received
2026-08 code and of the 2026-03 version for the record.
