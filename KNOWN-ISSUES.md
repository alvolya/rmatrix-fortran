# Known limitations and conventions

Items listed below are documented behaviour. 

## Conventions

- **Masses are in atomic mass units.**  k = sqrt(2 m_u mu E)/hbar c and
  eta = alpha Z1 Z2 sqrt(m_u mu / 2E) with mu the reduced mass in u
  (`channel_data.f90`: DWE, DET).  Integer mass numbers are accepted as
  an approximation (0.03–0.1 % in mu).  Before 2026-10 the constants
  implied a mass unit of 938.4 MeV.
- **Partial waves without a level** carry no hard-sphere scattering:
  the T-matrix is built only for the J^pi blocks that appear in the
  parameter file, and within a block only for the (channel, L) rows
  listed there (zero-width rows do get the hard-sphere phase).  AZURE2
  follows the same convention; `rmatrixpy` offers both
  (`partial_waves: levels | all`).
- **Rows are positional.**  All levels of one J^pi block must list the
  same (ChN, L, Bc) rows in the same order; the reader enforces this.
- **Boundary condition** B_c is a free per-row input in the `.par` file
  (dimensionless; often B_c = -L), or, with `BcMode = 1`, B_c = S_c(E_1)
  at the first level of each J^pi block (AZURE2's convention).  Nothing
  converts between formal and observed parameters; level energies
  E_lambda are the formal pole positions for the chosen B_c.
- **Threshold handling.**  Within 0.1 keV of any threshold the channel
  energy is nudged 0.1 keV away from it (E_c = 0 counts as open).
- **Closed channels.**  L_c = (shift from the Whittaker function
  W_{-eta, L+1/2}(2 kappa a)) - B_c, P_c = 0, through GSL's confluent
  hypergeometric U.  The complex Coulomb functions are still evaluated
  for closed channels by WCLBES and then discarded.
- **Energy grid.**  `nsteps = int((Efinal - Estart)/Estep)` points
  starting at Estart; Efinal itself is not computed.

## Limitations

- **No ANC input, no observed-parameter input.**  Bound and
  sub-threshold levels are ordinary poles with a formal reduced-width
  amplitude gamma; the ANC of a bound state must be converted to gamma
  outside the code (Whittaker normalisation, Thomas factor — see
  docs/PHYSICS.md; `rmatrixpy` accepts ANCs and observed widths directly).
- **Entrance channel = first final-state group**; the Coulomb amplitude
  uses eta and k of channel 1.  Any exit group can be selected (ExitCh).
- **Identical particles are not symmetrised** (no Mott interference,
  no factor for identical-particle channels).  The reader warns when
  the entrance pair consists of identical particles.
- **No photon channels, no radiative capture.**
- **Angle-integrated output (`.total`).**  For a charged-particle
  elastic exit it contains the nuclear Blatt–Biedenharn B_0 term only
  (the Rutherford and interference parts diverge); the S-factor column
  is built from it.
- **DKEY = 1 chi-squared mode** reads a hard-coded `data.dat`, divides
  by the number of points, not degrees of freedom, and prints one
  number; it is a remnant of the original fitting workflow.
- **RDef** (resonance finder) works on one level at a time, scans in 1
  keV steps from the threshold of channel `rdchan`, and reports the
  peak of |1 - W_cc|^2 and its FWHM in that channel.
- **Array limits** are compile-time constants in `channel_data.f90`:
  MaxCh = 50 channels, MaxNL = 16 (Lmax), MaxTJL = 20 (2J), MaxRes =
  1000 levels; the reader rejects inputs that exceed them.  The
  T-matrix array TJL is static (about 240 MB).
- **Angles below 1 degree** stop the run (the Rutherford amplitude is
  evaluated for every exit, elastic or not).
- **Absorption phase** (`imphase = 1`, APSRB) is carried over from the
  original code and is untested.

## Numerics

- Coulomb functions: CERNLIB WCLBES (COULCC) in double precision,
  cross-checked against mpmath to 1e-11 for open channels.
- Output precision: 7 significant digits (`ES15.7`).
- Pole term: 1/(E_lambda - E) is used as is; only an exact coincidence
  |E_lambda - E| < 1e-9 MeV is replaced by a sign-preserving 1e-9 MeV.
