! Rmatrix — multilevel, multichannel R-matrix code
! Copyright (C) 1998-2004, 2026 Grigory V. Rogachev
! Copyright (C) 2026 Alexander Volya
! SPDX-License-Identifier: GPL-3.0-or-later
! This program is free software under the GNU General Public License,
! version 3 or (at your option) any later version; see LICENSE.
! Third-party components: see NOTICE.
SUBROUTINE DOCALC
! ----------------------------------------------------------------------
! Purpose: Main calculation loop over energy and angle
! ----------------------------------------------------------------------
  USE channel_data
  USE resonance_data
  USE tmat_data
  USE prep_data
  USE io_filenames
  IMPLICIT NONE

  DOUBLE PRECISION :: ecm, acm, cshi, xi2
  DOUBLE PRECISION :: csexp, cserror
  INTEGER          :: nchan, istep, nsteps, nxist, ios, dot_pos
  CHARACTER(LEN=256) :: fname, base, errmsg

  IF (LEN_TRIM(out_file) == 0) THEN
    WRITE(6,'(2x,a)', advance='no') 'Output file name: '
    READ(5,*) out_file
  END IF
  fname = out_file

  ! Derive base name by stripping the last extension (e.g. "18O_test.out" -> "18O_test")
  base    = TRIM(fname)
  dot_pos = INDEX(base, '.', BACK=.TRUE.)
  IF (dot_pos > 1) base = base(1:dot_pos-1)

  OPEN(1,  file=TRIM(fname),                status='unknown', iostat=ios, iomsg=errmsg)
  IF (ios /= 0) THEN
    WRITE(*,'(a,a)') 'ERROR: cannot open output file: ', TRIM(errmsg)
    STOP 1
  END IF
  OPEN(2,  file=TRIM(base)//'.total',  status='unknown')
  OPEN(4,  file=TRIM(base)//'.screen', status='unknown')

  ! Column headers (# prefix is recognised by gnuplot, Python pandas, etc.)
  WRITE(1,'(a)') '# Ecm(MeV)         angle(deg)     dsigma/dOmega(mb/sr)'
  WRITE(2,'(a)') '# Ecm(MeV)      TotalCS(mb)   S-factor(MeV.mb)'
  WRITE(4,'(a)') '# E+Q(MeV)      Ue(MeV)       escreen'

  DO nchan = 1, lastch
    RMass(nchan) = Ap(nchan)*At(nchan) / (Ap(nchan)+At(nchan))
  ENDDO

  IF (BcMode.EQ.1) CALL SetBoundary

  IF ((NRDef.gt.0).and.(NRDef.le.lastre)) THEN
    PRINT *
    CALL RDef
  ENDIF

  IF (NRDef.eq.-1) THEN
    PRINT *
    DO NRDef = 1, lastre
      CALL RDef
    ENDDO
  ENDIF

  IF (DKEY.NE.1) THEN
    nsteps = (Efinal - Estart) / Estep
  ELSE
    nsteps = numdp
    xi2    = 0.0d0
    nxist  = 0
  ENDIF

  DO istep = 1, nsteps

    IF (DKEY.NE.1) THEN
      ecm = Estart + (istep-1)*Estep
      acm = 0.0d0
    ELSE
      ecm     = expdat(istep,1)
      csexp   = expdat(istep,2)
      cserror = expdat(istep,3)
      acm     = expdat(istep,4)
    ENDIF

    if ((ecm.GE.Estart).AND.(ecm.LE.Efinal)) then
    
      CALL Kinematics(ecm)

      CALL PreTJL(ecm)
      CALL rmblatt(ecm, acm, cshi)

      IF (DKEY.EQ.1) THEN
        if ((ecm.GT.ecmlow).and.(ecm.LT.ecmhigh)) then
          xi2   = xi2 + (csexp-cshi)**2.d0 / cserror**2.d0
          nxist = nxist + 1
        endif
      ENDIF

    endif
  ENDDO

  IF (DKEY.EQ.1) THEN
    PRINT *, 'xi^2 = ', xi2/nxist
  ENDIF

  CLOSE (UNIT=1)
  CLOSE (UNIT=2)

END SUBROUTINE DOCALC
