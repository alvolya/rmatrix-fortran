! Rmatrix — multilevel, multichannel R-matrix code
! Copyright (C) 1998-2004, 2026 Grigory V. Rogachev
! Copyright (C) 2026 Alexander Volya
! SPDX-License-Identifier: GPL-3.0-or-later
! This program is free software under the GNU General Public License,
! version 3 or (at your option) any later version; see LICENSE.
! Third-party components: see NOTICE.
SUBROUTINE APSRB(Aphase, p, nchan)
! ----------------------------------------------------------------------
! Purpose: Artificial absorption phase (imaginary addition to
!          hard-sphere phase for optical-model-like attenuation).
!
!   delta = APA * (1 - (1+exp((rho-L)/APAL))^-1)
!               * (1 - (1+exp((rho-APEP)/APAP))^-1)
!
!   Parameters from channel_data: APA, APAL, APEP, APAP
! ----------------------------------------------------------------------
  USE channel_data
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(OUT) :: Aphase(0:MaxNL, MaxCh)
  COMPLEX*16, INTENT(IN)  :: p(MaxCh)
  INTEGER,    INTENT(IN)  :: nchan
  INTEGER    :: L
  COMPLEX*16 :: x

  x = p(nchan)
  DO L = 0, MaxNL
    Aphase(L,nchan) = dble(APA * (1.d0 - (1.d0+exp((x-L)/APAL))**(-1)) * &
                                 (1.d0 - (1.d0+exp((x-APEP)/APAP))**(-1)))
  ENDDO

END SUBROUTINE APSRB
