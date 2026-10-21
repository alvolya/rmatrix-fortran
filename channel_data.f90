! Rmatrix — multilevel, multichannel R-matrix code
! Copyright (C) 1998-2004, 2026 Grigory V. Rogachev
! Copyright (C) 2026 Alexander Volya
! SPDX-License-Identifier: GPL-3.0-or-later
! This program is free software under the GNU General Public License,
! version 3 or (at your option) any later version; see LICENSE.
! Third-party components: see NOTICE.
      MODULE channel_data
        IMPLICIT NONE

        INTEGER, PARAMETER :: MaxCh  = 50
        INTEGER, PARAMETER :: MaxNL  = 16 
        INTEGER, PARAMETER :: MaxTJL = 20 
        INTEGER, PARAMETER :: MaxRes = 1000
        DOUBLE PRECISION, PARAMETER :: PI  = 3.14159265358979323846d0
        DOUBLE PRECISION, PARAMETER :: d2r = PI/180.d0
        ! Physical constants (CODATA 2018)
        DOUBLE PRECISION, PARAMETER :: HBARC = 197.3269804d0   ! MeV fm
        DOUBLE PRECISION, PARAMETER :: ALPHA = 1.d0/137.035999084d0
        DOUBLE PRECISION, PARAMETER :: AMU   = 931.49410242d0    ! MeV, atomic mass unit
        ! Kinematic constants for masses in u:
        !   k [fm^-1] = DWE sqrt(mu E),  eta = DET sqrt(mu) Z1 Z2 / sqrt(E)
        ! (the 2004 code used DWE = 0.21954, DET = 0.158052, i.e. a mass unit
        !  of 938.4 MeV; with masses in u that made k 0.37 % too large)
        DOUBLE PRECISION, PARAMETER :: DWE = sqrt(2.d0*AMU)/HBARC
        DOUBLE PRECISION, PARAMETER :: DET = ALPHA*sqrt(AMU/2.d0)

        ! Channel particle data
        DOUBLE PRECISION :: Ap(MaxCh), Zp(MaxCh), Sp(MaxCh)
        DOUBLE PRECISION :: At(MaxCh), Zt(MaxCh), St(MaxCh)
        DOUBLE PRECISION :: Q(MaxCh)
        INTEGER          :: Pp(MaxCh), Pt(MaxCh), FS(MaxCh)
        DOUBLE PRECISION :: ChSp(MaxCh), ChR(MaxCh)

        ! Geometry and control
        INTEGER :: Lmax, lastch, nfs, inv_out
        LOGICAL :: Inverse

        ! Energy and angle ranges
        DOUBLE PRECISION :: Estart, Efinal, Estep
        DOUBLE PRECISION :: Astart, Afinal, Astep
        INTEGER          :: ExitCh, NRDef
        ! boundary conditions: 0 = Bc column of the .par file,
        ! 1 = B_c = S_c(E_1) at the first level of each J^pi block (AZURE2)
        INTEGER          :: BcMode = 0

        ! Experimental data
        DOUBLE PRECISION :: expdat(0:10000, 4)
        DOUBLE PRECISION :: QEx
        INTEGER :: numdp, DKEY

        ! Chi-squared energy window
        DOUBLE PRECISION :: ecmlow, ecmhigh

        ! Absorption phase parameters
        INTEGER :: imphase
        DOUBLE PRECISION :: APA, APAL, APEP, APAP

        ! Channel definition for resonance finder
        INTEGER :: rdchan

      END MODULE channel_data
