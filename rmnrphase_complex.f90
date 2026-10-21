! Rmatrix — multilevel, multichannel R-matrix code
! Copyright (C) 1998-2004, 2026 Grigory V. Rogachev
! Copyright (C) 2026 Alexander Volya
! SPDX-License-Identifier: GPL-3.0-or-later
! This program is free software under the GNU General Public License,
! version 3 or (at your option) any later version; see LICENSE.
! Third-party components: see NOTICE.
SUBROUTINE NResPh(NonRPh, F, G, FP, GP, SIG, p, ecm)
! ----------------------------------------------------------------------
! Purpose: Hard-sphere (non-resonant) phase shifts and Coulomb
!          wave functions F_L, G_L, F'_L, G'_L, sigma_L.
! Created 7-OCT-1998, Modified 07-MAR-2004 Grigory Rogachev, 2026 A. Volya
! ----------------------------------------------------------------------
  USE channel_data
  USE prep_data
  IMPLICIT NONE

  COMPLEX*16, INTENT(OUT) :: NonRPh(0:MaxNL, MaxCh)
  COMPLEX*16, INTENT(OUT) :: F(0:MaxNL, MaxCh), FP(0:MaxNL, MaxCh)
  COMPLEX*16, INTENT(OUT) :: G(0:MaxNL, MaxCh), GP(0:MaxNL, MaxCh)
  COMPLEX*16, INTENT(OUT) :: SIG(0:MaxNL, MaxCh)
  COMPLEX*16, INTENT(OUT) :: p(MaxCh)
  DOUBLE PRECISION, INTENT(IN) :: ecm

  INTEGER    :: nchan, il, JFAIL
  DOUBLE PRECISION :: Aphase(0:MaxNL, MaxCh)
  COMPLEX*16 :: zlmin, imone
  COMPLEX*16 :: Ftmp(0:MaxNL), FPtmp(0:MaxNL)
  COMPLEX*16 :: Gtmp(0:MaxNL), GPtmp(0:MaxNL), SIGtmp(0:MaxNL)
  DOUBLE PRECISION :: rel

  JFAIL = 0

  DO nchan = 1, lastch
    p(nchan) = ChR(nchan) * wavnum(nchan)

    zlmin = cmplx(0.d0, 0.d0)
    imone = cmplx(0.d0, 1.d0)

    if (ABS(p(nchan)).lt.1.d-6) p(nchan) = p(nchan) + 1.d-6

    CALL WCLBES(p(nchan), eta(nchan), zlmin, MaxNL, &
                Ftmp, Gtmp, FPtmp, GPtmp, SIGtmp, 0, 1, JFAIL, 1)

    IF (JFAIL.NE.0) THEN
      PRINT *, 'Error in WCLBES, code = ', JFAIL
      STOP
    ENDIF

    CALL APSRB(Aphase, p, nchan)

    DO il = 0, MaxNL
      F(il,nchan)   = Ftmp(il)
      FP(il,nchan)  = FPtmp(il)
      G(il,nchan)   = Gtmp(il)
      GP(il,nchan)  = GPtmp(il)
      SIG(il,nchan) = SIGtmp(il)

      if ((ecm + Q(nchan)).gt.1.d-4) then
        rel = dreal(F(il,nchan)) / dreal(G(il,nchan))
        if (dreal(F(il,nchan)).gt.0.d0) then
          if (dreal(G(il,nchan)).gt.0.d0) then
            NonRPh(il,nchan) = DATAN(rel)
          else
            NonRPh(il,nchan) = PI + DATAN(rel)
          endif
        else
          if (dreal(G(il,nchan)).gt.0.d0) then
            NonRPh(il,nchan) = 2.d0*PI + DATAN(rel)
          else
            NonRPh(il,nchan) = PI + DATAN(rel)
          endif
        endif
        if (imphase.EQ.1) then
          NonRPh(il,nchan) = NonRPh(il,nchan) + &
                              imone * Aphase(il,nchan)
        endif
      else
        NonRPh(il,nchan) = 0.d0
      endif
    ENDDO
  ENDDO

END SUBROUTINE NResPh
