! Rmatrix — multilevel, multichannel R-matrix code
! Copyright (C) 1998-2004, 2026 Grigory V. Rogachev
! Copyright (C) 2026 Alexander Volya
! SPDX-License-Identifier: GPL-3.0-or-later
! This program is free software under the GNU General Public License,
! version 3 or (at your option) any later version; see LICENSE.
! Third-party components: see NOTICE.
SUBROUTINE Kinematics(ecm)
! ----------------------------------------------------------------------
! Purpose: wave numbers k_c and Sommerfeld parameters eta_c of all
!          channels at the c.m. energy ecm (relative to the first channel
!          threshold).  Closed channels get imaginary k and eta.  Within
!          0.1 keV of a threshold the channel energy is nudged by
!          SMALLDELTA away from it; E_c = 0 counts as open.
!          Shared by DOCALC, RDef and SetBoundary.
! ----------------------------------------------------------------------
  USE channel_data
  USE prep_data
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: ecm
  DOUBLE PRECISION, PARAMETER  :: SMALLDELTA = 0.0001d0
  COMPLEX*16, PARAMETER        :: mone = (0.d0, 1.d0)
  DOUBLE PRECISION :: CEnergy
  INTEGER          :: nchan

  DO nchan = 1, lastch
    CEnergy = ecm + Q(nchan)
    if (CEnergy.GE.0.0d0) then
      if (CEnergy.LT.SMALLDELTA) CEnergy = CEnergy + SMALLDELTA
      wavnum(nchan) = DWE * dsqrt(RMass(nchan)*CEnergy)
      eta(nchan)    = DET * dsqrt(RMass(nchan)) * &
                      Zp(nchan)*Zt(nchan) / dsqrt(CEnergy)
    else
      if (CEnergy.GT.-SMALLDELTA) CEnergy = CEnergy - SMALLDELTA
      wavnum(nchan) = mone * DWE * dsqrt(RMass(nchan)*(-1.d0)*CEnergy)
      eta(nchan)    = (-1.d0)*mone * DET * dsqrt(RMass(nchan)) * &
                      Zp(nchan)*Zt(nchan) / dsqrt((-1.d0)*CEnergy)
    endif
  ENDDO

END SUBROUTINE Kinematics
