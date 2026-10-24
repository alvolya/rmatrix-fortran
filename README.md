# Rmatrix

Multilevel, multichannel R-matrix code for differential cross sections of
low-energy nuclear reactions.  Original code G.V. Rogachev (Texas A&M
University, 1998–2004; update 2026); Fortran 90
modernisation, fixes and tests A. Volya (Florida State University, 2026).
Parts of the 2026 revision were prepared with the help of Claude (Anthropic) 
under the authors' direction and review.
See `KNOWN-ISSUES.md` for conventions and limitations.
License: GNU GPL version 3 or later, SPDX `GPL-3.0-or-later` (`LICENSE`); third-party
components in `NOTICE`; citation in `CITATION.cff`, DOI
[10.5281/zenodo.23617712](https://doi.org/10.5281/zenodo.23617712).

## Physics

The code implements the R-matrix formalism for nuclear reactions:

- **Input**: channel definitions (masses, charges, spins, channel radii), resonance
  parameters (energy, J, parity, reduced widths per channel and L)
- **Output**: differential cross section dσ/dΩ (mb/sr) or phase shifts as a function
  of centre-of-mass energy and angle
- **Angular distributions**: Blatt–Biedenharn decomposition using Legendre polynomials
- **Coulomb wave functions**: CERNLIB `WCLBES` (COULCC, Barnett/Steed); closed
  channels through the Whittaker function (GSL confluent hypergeometric U)
- **Angular momentum coupling**: Clebsch–Gordan and Racah W coefficients (inline)
- **Masses** are entered in atomic mass units; k and η use the CODATA atomic mass unit.

### Formulas and conventions

The formalism is that of Lane and Thomas [1]: formal R-matrix parameters
(level energies E_λ and reduced-width amplitudes γ_λc) with a boundary
condition B_c per channel, the collision matrix from R, the hard-sphere
phases and the logarithmic derivative L_c = S_c − B_c + iP_c with the usual
P_c = ρ/(F_L² + G_L²) and S_c = ρ(F_L F_L′ + G_L G_L′)/(F_L² + G_L²) at ρ = k_c a_c;
angular distributions after Blatt and Biedenharn [2].  Conventions specific
to this code:

- masses in u, k and η from the CODATA constants (`channel_data.f90`);
- B_c from the `.par` file, or B_c = S_c(E_1) at the first level of each J^π
  block with `BcMode = 1` (the convention of AZURE2 [6]); level energies are
  formal, no observed-parameter (Brune [5]) input;
- closed channels: S_c from the logarithmic derivative of the Whittaker
  function, P_c = 0;
- Coulomb functions from COULCC [3,4]; Coulomb phases σ_L by recursion;
- the Legendre expansion of dσ/dΩ runs to 2·Lmax; partial waves enter only
  through the (channel, L) rows listed in the parameter file;
- dσ/dΩ in mb/sr, averaged over the entrance spins; for charged particles the
  Rutherford amplitude and its interference with the nuclear amplitude are
  included.

[1] A. M. Lane and R. G. Thomas, Rev. Mod. Phys. 30, 257 (1958).
[2] J. M. Blatt and L. C. Biedenharn, Rev. Mod. Phys. 24, 258 (1952).
[3] A. R. Barnett, D. H. Feng, J. W. Steed and L. J. B. Goldfarb, Comput. Phys.
    Commun. 8, 377 (1974).
[4] I. J. Thompson and A. R. Barnett, Comput. Phys. Commun. 36, 363 (1985).
[5] C. R. Brune, Phys. Rev. C 66, 044611 (2002).
[6] R. E. Azuma et al., Phys. Rev. C 81, 045805 (2010).
[7] P. Descouvemont and D. Baye, Rep. Prog. Phys. 73, 036301 (2010).

### Example: ¹⁸O test case

`tests/18O/` holds the example: α + ¹⁴C / n + ¹⁷O multichannel system with
5 channel-spin groups and a bound 0⁺ level, energy range 1.0–6.0 MeV (cm),
angles 10°–180° (`18O.inp`, `18O.par`; see its README).  The same files
serve the regression test.

## Build requirements

- **gfortran** — GCC Fortran compiler
  - macOS: `brew install gcc gsl`
  - Ubuntu/Debian: `sudo apt install gfortran make libgsl-dev liblapack-dev`

  Tested on macOS (Apple Silicon, GCC 16, GSL 2.8) and Ubuntu 24.04 (gfortran 13,
  GSL 2.7.1 and 2.8).
- **GSL** (closed-channel Whittaker function) and **LAPACK/BLAS** (complex matrix
  inversion; Accelerate on macOS).
- The CERNLIB routine `WCLBES` (C309, the complex Coulomb-function routine COULCC
  of Barnett and Thompson) is included as source together with the seven helper
  routines and the two gamma-function routines it calls (CERNLIB 2006, GNU GPL).

## Building

```sh
make            # builds ./Rmatrix
make test       # runs tests/run_tests.py on the build
make clean
make debug      # clean rebuild with bounds checking
```

GSL is located with `pkg-config`; without pkg-config the Makefile assumes
`/usr/local`, and `make GSL_PREFIX=/path` overrides both.  On Linux the
LAPACK/BLAS link line is `LAPACK_LIBS` (default `-llapack -lblas`); macOS uses
the Accelerate framework.  `make test` needs Python 3 with numpy.

### Quick start

```sh
make && make test
```

The test script runs the 18O example against stored references and a p + 14N
case against AZURE2's cross sections (see `tests/azure2_15O/README.md`).

## Running

```sh
./Rmatrix -i <input.inp> -p <params.par> -o <output.out>
```

All three flags are optional — if omitted the program prompts interactively:

```sh
./Rmatrix
  Input file name:   tests/18O/18O.inp
  Parameter file name: tests/18O/18O.par
  Output file name:  18O_test.out
```

Or pipe the filenames in the old way:

```sh
echo "tests/18O/18O.inp
tests/18O/18O.par
18O_test.out" | ./Rmatrix
```

Use `-h` or `--help` for a usage summary.

## Output files

The program writes several output files derived from the output name given with `-o`.
For `-o 18O_test.out` the following files are created:

| File | Contents |
|------|----------|
| `18O_test.out` | Main output: `Ecm (MeV)`, `angle (deg)`, `dσ/dΩ (mb/sr)` |
| `18O_test.total` | Total cross section and astrophysical S-factor vs energy |

Both files include a `#`-prefixed header line with column names, compatible with
gnuplot, Python pandas/numpy, and similar tools.

## Input file format (`.inp`)

```
! Comment lines starting with ! are ignored anywhere in the file

! Channel cards: Ap  Zp  Sp  At  Zt  St  Q(MeV)  ChSp  ChR(fm)
  <Ap> <Zp> <Sp> <At> <Zt> <St> <Q> <ChSp> <ChR>
  ...  (one line per channel; last channel has negative Ap)

&RMATRIX
  Lmax    = 8,           ! max orbital angular momentum in partial-wave expansion
  Estart  = 2.00,  Efinal = 6.00,  Estep = 0.005,   ! cm energy range (MeV)
  Astart  = 180.0, Afinal = 180.0, Astep = 20.0,    ! angular range (deg)
  ExitCh  = 1,           ! exit channel index (negative = inverse kinematics)
  NRDef   = 0,           ! resonance finder: 0=off, N=level N, -1=all levels
  DKEY    = 0,           ! 0 = uniform energy grid,  1 = use data.dat
  QEx     = 0.0,         ! threshold energy shift applied to data.dat (MeV)
  ecmlow  = 2.8,  ecmhigh = 4.5,   ! chi-squared fit window (MeV, used if DKEY=1)
  imphase = 0,           ! absorption phase: 0=off, 1=on
  APA     = -1.0,        ! absorption amplitude (radians; negative = absorptive)
  APAL    = 0.3,         ! diffuseness of angular-momentum cutoff
  APEP    = 3.0,         ! centre of rho (=k*R) cutoff
  APAP    = 0.3,         ! diffuseness of rho cutoff
  rdchan  = 1,           ! channel used to report resonance widths (NRDef)
  BcMode  = 0,           ! 0: Bc from the .par file; 1: Bc = S_c(E_1) at the first
                         !    level of each J^pi block (AZURE2's convention)
/
```

**Channel card columns**: `Ap/At` = mass in u (integer mass numbers are an acceptable
approximation), `Zp/Zt` = charge, `Sp/St` = spin (negative value encodes negative
parity), `Q` = Q-value (MeV) relative to the first channel, `ChSp` = channel spin,
`ChR` = channel radius (fm).  The last channel is flagged by a negative `Ap`.
Channels with the same `Q` form one final-state group; the first group is the
entrance, `ExitCh` selects the exit group.  `Lmax` must not exceed 16 (MaxNL);
the Legendre expansion of the cross section runs to 2·Lmax.

**Absorption phase** (`imphase`, `APA`, `APAL`, `APEP`, `APAP`): an optional
optical-model-like term that adds a complex imaginary part to the hard-sphere phase
shifts, simulating flux loss to unobserved open channels.  The phase follows a product
of two Fermi-Dirac steps in $\ell$ and $\rho = k_c a_c$ respectively — see
`PHYSICS.md` §4 for the full formula.  Set `imphase = 0` to disable (default).

## Parameter file format (`.par`)

```
; four comment/header lines (skipped)
; Eres  Jtotal  Parity
  <E(MeV)>  <Jtotal>  <parity(±1)>
    <ChN>  <L>  <Bc>  <G>   ! channel N, orbital L, boundary cond., reduced width
    ...
   -<ChN> <L>  <Bc>  <G>   ! negative ChN = last channel-spin group for this level
  ...                       ! repeat for each resonance level
```

Reduced-width amplitudes `G` have units MeV$^{1/2}$ (formal parameters for the given
`Bc`).  The boundary condition `Bc` is dimensionless, often `Bc = −L`.  Level
energies are c.m. energies relative to the first channel threshold; negative
values are bound levels.  Every row must satisfy the parity rule
π_p π_t (−1)^L = π and the triangle rule |J−s| ≤ L ≤ J+s; a forbidden row with
non-zero `G` is rejected, a forbidden row with `G = 0` is ignored.  All levels of
one J^π must list the same `(ChN, L, Bc)` rows in the same order.

## Source file overview

### Fortran 90 modules

| File | Module | Contents |
|------|--------|---------|
| `io_filenames.f90` | `io_filenames` | Shared storage for the three file paths (`inp_file`, `par_file`, `out_file`) |
| `channel_data.f90` | `channel_data` | Channel definitions, energy/angle ranges, run parameters, constants |
| `resonance_data.f90` | `resonance_data` | Resonance energies, reduced widths, level count |
| `tmat_data.f90` | `tmat_data` | T-matrix array `TJL(2J, ch_in, ch_out, L_in, L_out)` |
| `prep_data.f90` | `prep_data` | Kinematic arrays: reduced masses, wave numbers, Coulomb parameters |

### R-matrix physics routines

| File | Subroutine | Purpose |
|------|-----------|---------|
| `rm.f90` | `MAIN` | Entry point: parses CLI flags, calls PREAD → REREAD → DOCALC |
| `rmpread.f90` | `PREAD` | Read channel input file (`.inp`) |
| `rmreread.f90` | `REREAD` | Read resonance parameter file (`.par`) |
| `rmdocalc.f90` | `DOCALC` | Main energy/angle loop, open output files |
| `pretjl.f90` | `PreTJL` | Build R-matrix, invert → collision matrix W → T-matrix |
| `rdef.f90` | `RDef` | Locate resonance energy and observed width |
| `rmblatt.f90` | `rmblatt` | Blatt–Biedenharn angular distribution → dσ/dΩ |
| `rmbbb.f90` | `BBB` | Blatt–Biedenharn B-coefficients |
| `rmnrphase_complex.f90` | `NResPh` | Hard-sphere (non-resonant) phase shifts |
| `rmcoulph.f90` | `CPC` | Coulomb phase shifts σ_L |
| `apsrb_phase.f90` | `APSRB` | Artificial absorption phase (optical-model attenuation) |
| `intgsl.f90` | `intgsl` | Whittaker function W_{κ,μ}(z) through GSL `gsl_sf_hyperg_U_e` |

### CERNLIB routines (`*.F`, unmodified CERNLIB 2006 sources, GNU GPL)

The code calls only `WCLBES`, for the Coulomb wave functions and their derivatives
at the channel radius; the other files are its internal dependencies.

| File | Subroutine | Purpose |
|------|-----------|---------|
| `wclbes.F` | `WCLBES` | Coulomb wave functions F_L, G_L, F′_L, G′_L (COULCC, Barnett/Steed) |
| `r1dp.F`–`r8dp.F` | internal | Continued-fraction helpers for WCLBES |
| `clogam64.F` | `WLOGAM`/`WLGAMA` | Complex log-gamma |
| `cdigam64.F` | `WDIGAM` | Complex digamma |
| `mtlprt.for` | `MTLPRT` | Error message printer (stub) |

Clebsch–Gordan and Racah W coefficients are implemented inline in `rmbbb.f90`
following the CERNLIB `DCLEBG`/`DRACAW` algorithms; Legendre polynomials by
3-term recursion in `rmblatt.f90`.

## Tests

```sh
python3 tests/run_tests.py                 # builds, then tests (needs numpy)
python3 tests/run_tests.py --exe ./Rmatrix # existing build
```

Each case is self-contained in its directory:

- `tests/18O/refs/`: outputs of the 18O case (180°, and 20°–180°) produced by this
  code; every point outside 1.5 keV of a level must agree to 1e-6.
- `tests/18O/refs-frozen-2026/`: outputs of the as-received 2026-08 code and of the
  2026-03 version, kept for the record.  They differ from the current code
  (mass unit, Legendre order, pole term).
- `tests/azure2_15O/`: p + 14N elastic scattering on the parameters of AZURE2's
  own test input (public repository, pinned commit), compared with AZURE2's cross
  sections at 60°, 90° and 140°; agreement 2e-7.  `make_case.py` regenerates the
  case from an AZURE2 checkout with its Python binding; running the test needs
  no AZURE2.

## Plotting

`plot_xsec.py` is a Python script for visualising the output files.
Requires **numpy** and **matplotlib**.

### Excitation function (σ vs energy at fixed angle)

```sh
# Single file — one figure per angle present in the output
python3 plot_xsec.py 18O_dp.output

# Restrict to one angle
python3 plot_xsec.py 18O_dp.output --angle 180

# All angles overlaid on one figure
python3 plot_xsec.py 18O_dp.output --overlay

# Compare two calculations
python3 plot_xsec.py 18O.output 18O_dp.output --labels "SP,DP"

# Logarithmic y-axis and save to PDF
python3 plot_xsec.py 18O_dp.output --log --save xsec.pdf
```

### Getting multiple angles from Rmatrix

Edit the `&RMATRIX` namelist in your `.inp` file:

```fortran
Astart = 20.0,  Afinal = 180.0,  Astep = 20.0,
```

This produces angles 20°, 40°, …, 180° in a single run.

### Angular distribution (σ vs angle at fixed energy)

```sh
# Auto-select a few representative energies
python3 plot_xsec.py 18O_dp.output --angular-dist

# Specify energies explicitly (MeV, nearest grid point is used)
python3 plot_xsec.py 18O_dp.output --angular-dist --energy 3.0,4.1,5.0
```

### Overlaying experimental data

Provide a plain-text data file with `--data`.  Supported column layouts
(whitespace- or comma-separated; lines beginning with `#` or `!` are comments):

| Columns | Meaning |
|---------|---------|
| `Ecm  xsec` | 2-col: plotted on every angle panel |
| `Ecm  angle  xsec` | 3-col: matched to the correct angle panel |
| `Ecm  xsec  err_lo  err_hi` | 4-col: 2-col with asymmetric error bars |
| `Ecm  angle  xsec  err_lo  err_hi` | 5-col: 3-col with asymmetric error bars |

```sh
python3 plot_xsec.py 18O_dp.output --data exp.dat --data-label "Smith (2020)"

# Single angle with data, save to PDF
python3 plot_xsec.py 18O_dp.output --angle 180 --data exp.dat --save xsec.pdf
```

### All options

| Flag | Description |
|------|-------------|
| `--save FILE` | Save figure(s) to file (`.png` / `.pdf` / `.svg`). When multiple figures are produced the angle or energy is embedded in the filename. |
| `--log` | Logarithmic y-axis |
| `--angle DEG` | Plot only this angle (excitation-function mode) |
| `--overlay` | All angles on one figure (excitation-function mode) |
| `--angular-dist` | Plot σ(θ) at fixed energies instead of σ(E) |
| `--energy E1,E2,...` | Energies for `--angular-dist` (nearest grid point used) |
| `--labels L1,L2,...` | Legend labels, one per theory file |
| `--data FILE` | Experimental data file |
| `--data-label STR` | Legend label for experimental data (default: `data`) |
| `--no-legend` | Suppress legend |

## Acknowledgment

The 2026 modernisation, review and release of this code were carried out within
the DOE Genesis Mission project "Quantum-Enhanced Evaluation of Nuclear Reaction Cross
Sections in R-Matrix Theory" (Florida State University, Texas A&M University,
University of Notre Dame).

This material is based upon work supported by the Department of Energy under
Award Number DE-SC0026615.

This report was prepared as an account of work sponsored by an agency of the
United States Government. Neither the United States Government nor any agency
thereof, nor any of their employees, makes any warranty, express or implied, or
assumes any legal liability or responsibility for the accuracy, completeness, or
usefulness of any information, apparatus, product, or process disclosed, or
represents that its use would not infringe privately owned rights. Reference
herein to any specific commercial product, process, or service by trade name,
trademark, manufacturer, or otherwise does not necessarily constitute or imply
its endorsement, recommendation, or favoring by the United States Government or
any agency thereof. The views and opinions of authors expressed herein do not
necessarily state or reflect those of the United States Government or any agency
thereof.
