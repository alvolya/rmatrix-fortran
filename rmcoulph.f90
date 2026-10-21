! Rmatrix — multilevel, multichannel R-matrix code
! Copyright (C) 1998-2004, 2026 Grigory V. Rogachev
! Copyright (C) 2026 Alexander Volya
! SPDX-License-Identifier: GPL-3.0-or-later
! This program is free software under the GNU General Public License,
! version 3 or (at your option) any later version; see LICENSE.
! Third-party components: see NOTICE.
SUBROUTINE CPC(Ph)
! ----------------------------------------------------------------------
! Purpose: Coulomb phase shifts sigma_L = sum_{k=1}^{L} atan(|eta|/k)
! Created 7-OCT-1998, Modified 07-MAR-2004 Grigory Rogachev
! ----------------------------------------------------------------------
  USE channel_data
  USE prep_data
  IMPLICIT NONE

  COMPLEX*16, INTENT(OUT) :: Ph(0:MaxNL, MaxCh)
  INTEGER   :: k, nchan
  DOUBLE PRECISION :: reta, cphase

  DO nchan = 1, lastch
    Ph(0,nchan) = 0.d0
    DO k = 1, MaxNL
      reta   = abs(eta(nchan))
      cphase = ATAN(reta / dble(k))
      if (cphase.lt.0.d0) cphase = PI + cphase
      Ph(k,nchan) = Ph(k-1,nchan) + cphase
    ENDDO
  ENDDO

END SUBROUTINE CPC
