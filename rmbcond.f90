! Rmatrix — multilevel, multichannel R-matrix code
! Copyright (C) 1998-2004, 2026 Grigory V. Rogachev
! Copyright (C) 2026 Alexander Volya
! SPDX-License-Identifier: GPL-3.0-or-later
! This program is free software under the GNU General Public License,
! version 3 or (at your option) any later version; see LICENSE.
! Third-party components: see NOTICE.
SUBROUTINE SetBoundary
! ----------------------------------------------------------------------
! Purpose: BcMode = 1 — boundary condition B_c = S_c(E_1) for every row
!          of a J^pi block, where E_1 is the energy of the FIRST level of
!          that block in the parameter file (the convention of AZURE2,
!          CNuc::CalcBoundaryConditions).  Open channel:
!          S = Re[rho (F F' + G G')/(F^2 + G^2)]; closed channel: the
!          logarithmic derivative of the Whittaker function (as PreTJL).
!          Overrides the Bc column of the .par file for all levels of the
!          block; the values are printed.
! ----------------------------------------------------------------------
  USE channel_data
  USE resonance_data
  USE prep_data
  IMPLICIT NONE

  INTEGER :: k, m, n, ch, L, J2, PA
  LOGICAL :: first
  COMPLEX*16 :: NRPh(0:MaxNL, MaxCh), SIG(0:MaxNL, MaxCh), p(MaxCh)
  COMPLEX*16 :: F(0:MaxNL, MaxCh), FP(0:MaxNL, MaxCh)
  COMPLEX*16 :: G(0:MaxNL, MaxCh), GP(0:MaxNL, MaxCh)
  COMPLEX*16, PARAMETER :: i = (0.d0, 1.d0)
  DOUBLE PRECISION :: S, KAW, ZZ, MU, W1, W2, E1
  DOUBLE PRECISION, EXTERNAL :: intgsl

  WRITE(*,'(2x,a)') 'Boundary conditions B_c = S_c(E_1), E_1 = first level of each J^pi (BcMode = 1):'
  DO k = 1, lastre
    J2 = nint(2.d0*eres(k,2))
    PA = nint(eres(k,3))
    first = .TRUE.
    DO m = 1, k-1
      if (nint(2.d0*eres(m,2)) == J2 .and. nint(eres(m,3)) == PA) first = .FALSE.
    ENDDO
    if (.not. first) cycle

    E1 = eres(k,1)
    CALL Kinematics(E1)
    CALL NResPh(NRPh, F, G, FP, GP, SIG, p, E1)

    DO n = 1, rdim(k)
      ch = nint(gres(k,n,1))
      L  = nint(gres(k,n,2))
      if ((E1 + Q(ch)).lt.0.d0) then
        KAW = dble((-i)*eta(ch))
        ZZ  = dble(2.d0*(-i)*p(ch))
        MU  = dble(L) + 0.5d0
        W1  = intgsl(KAW, MU, ZZ)
        W2  = intgsl(KAW-1.d0, MU, ZZ)
        S   = ((KAW-0.5d0*ZZ)*W1 - (MU*MU-(KAW-0.5d0)**2)*W2) / W1
      else
        S = dble(p(ch) * (F(L,ch)*FP(L,ch) + G(L,ch)*GP(L,ch)) / &
                 (F(L,ch)*F(L,ch) + G(L,ch)*G(L,ch)))
      endif
      DO m = 1, lastre
        if (nint(2.d0*eres(m,2)) == J2 .and. nint(eres(m,3)) == PA) gres(m,n,3) = S
      ENDDO
      WRITE(*,'(4x,a,f5.1,a,i2,a,i3,a,i2,a,f10.4,a,f13.8)') 'J =', eres(k,2), &
            '  parity', PA, '  channel', ch, '  L =', L, '  E_1 =', E1, '  B_c =', S
    ENDDO
  ENDDO

END SUBROUTINE SetBoundary
