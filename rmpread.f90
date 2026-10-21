! Rmatrix — multilevel, multichannel R-matrix code
! Copyright (C) 1998-2004, 2026 Grigory V. Rogachev
! Copyright (C) 2026 Alexander Volya
! SPDX-License-Identifier: GPL-3.0-or-later
! This program is free software under the GNU General Public License,
! version 3 or (at your option) any later version; see LICENSE.
! Third-party components: see NOTICE.
SUBROUTINE PREAD
! ----------------------------------------------------------------------
! Purpose: Read channel input file (.inp)
!
! File format:
!   - Lines beginning with ! or blank are ignored (anywhere in file)
!   - One line per channel: Ap Zp Sp At Zt St Q ChSp ChR
!     (last channel marked with negative Ap)
!   - &RMATRIX namelist block with all scalar parameters
! ----------------------------------------------------------------------
  USE channel_data
  USE io_filenames
  IMPLICIT NONE

  INTEGER            :: i, n, ios
  LOGICAL            :: LCH
  CHARACTER(LEN=256) :: fname, errmsg, line

  NAMELIST /RMATRIX/ Lmax,                          &
                     Estart, Efinal, Estep,          &
                     Astart, Afinal, Astep,          &
                     ExitCh, NRDef, DKEY, QEx,       &
                     ecmlow, ecmhigh,                &
                     imphase, APA, APAL, APEP, APAP, &
                     rdchan, BcMode

  IF (LEN_TRIM(inp_file) == 0) THEN
    WRITE(6,'(2x,a)', advance='no') 'Input file name: '
    READ(5,*) inp_file
  END IF
  fname = inp_file

  OPEN(1, file=TRIM(fname), status='old', iostat=ios, iomsg=errmsg)
  IF (ios /= 0) THEN
    WRITE(*,'(a,a)') 'ERROR: cannot open input file: ', TRIM(errmsg)
    STOP 1
  END IF

  ! Read channel cards, skipping blank lines and ! comments.
  ! Each data line is read into a buffer first, then parsed internally.
  LCH = .TRUE.
  i   = 1
  nfs = 1
  DO WHILE (LCH)
    ! Skip blank/comment lines to reach the next data line
    DO
      READ(1, '(a)', iostat=ios) line
      IF (ios /= 0) THEN
        WRITE(*,'(a)') 'ERROR: unexpected end of file reading channel cards'
        STOP 1
      END IF
      line = ADJUSTL(line)
      IF (LEN_TRIM(line) > 0 .AND. line(1:1) /= '!') EXIT
    END DO
    IF (i > MaxCh) THEN
      WRITE(*,'(a,i0,a)') 'ERROR: more than MaxCh = ', MaxCh, ' channel cards'
      STOP 1
    END IF
    READ(line, *, iostat=ios) Ap(i), Zp(i), Sp(i), At(i), Zt(i), &
                  St(i), Q(i), ChSp(i), ChR(i)
    IF (ios /= 0) THEN
      WRITE(*,'(a,i0)') 'ERROR: cannot parse channel card ', i
      STOP 1
    END IF

    if (i > 1) then
      if ((Q(i) - Q(i-1)) /= 0) nfs = nfs + 1
    endif

    if (Sp(i) < 0.) then
      Sp(i) = ABS(Sp(i))
      Pp(i) = -1
    else
      Pp(i) = +1
    endif

    if (St(i) < 0.) then
      St(i) = ABS(St(i))
      Pt(i) = -1
    else
      Pt(i) = +1
    endif

    FS(i) = nfs

    IF (Ap(i) < 0.) THEN
      Ap(i)  = ABS(Ap(i))
      LCH    = .FALSE.
      lastch = i
    ENDIF
    i = i + 1
  ENDDO

  ! Read scalar parameters from &RMATRIX namelist block.
  ! NAMELIST natively supports ! comments inside the block.
  READ(1, NML=RMATRIX, iostat=ios, iomsg=errmsg)
  IF (ios /= 0) THEN
    WRITE(*,'(a,a)') 'ERROR: reading &RMATRIX namelist: ', TRIM(errmsg)
    STOP 1
  END IF

  CLOSE(UNIT=1)

  IF (BcMode < 0 .or. BcMode > 1) THEN
    WRITE(*,'(a,i0)') 'ERROR: BcMode must be 0 or 1, got ', BcMode
    STOP 1
  END IF
  IF (Lmax < 0 .or. Lmax > MaxNL) THEN
    WRITE(*,'(a,i0,a,i0)') 'ERROR: Lmax = ', Lmax, ' must be between 0 and MaxNL = ', MaxNL
    STOP 1
  END IF
  IF (abs(ExitCh) < 1 .or. abs(ExitCh) > nfs) THEN
    WRITE(*,'(a,i0,a,i0,a)') 'ERROR: ExitCh = ', ExitCh, ' but there are ', nfs, &
                             ' final-state groups'
    STOP 1
  END IF
  DO i = 1, lastch
    IF (ChSp(i) < abs(Sp(i)-St(i)) - 1.d-6 .or. ChSp(i) > Sp(i)+St(i) + 1.d-6 .or. &
        abs(2.d0*ChSp(i) - nint(2.d0*ChSp(i))) > 1.d-6 .or. &
        mod(nint(2.d0*(ChSp(i)-Sp(i)-St(i))), 2) /= 0) THEN
      WRITE(*,'(a,i0,a,f5.1,a,f4.1,a,f4.1)') 'ERROR: channel ', i, ': channel spin ', &
            ChSp(i), ' not allowed for spins ', Sp(i), ' and ', St(i)
      STOP 1
    END IF
  ENDDO
  IF (abs(Ap(1)-At(1)) < 1.d-6 .and. abs(Zp(1)-Zt(1)) < 1.d-6 .and. &
      abs(Sp(1)-St(1)) < 1.d-6 .and. Pp(1) == Pt(1)) THEN
    WRITE(*,'(a)') 'WARNING: identical particles in the entrance channel; this code'
    WRITE(*,'(a)') '         does not symmetrise the amplitude (no Mott interference).'
  END IF

  ! Resolve inverse-reaction flag from sign of ExitCh
  if (ExitCh < 0) then
    ExitCh  = ABS(ExitCh)
    Inverse = .TRUE.
    do i = 1, lastch
      if (FS(i) == ExitCh) inv_out = i
    enddo
  else
    Inverse = .FALSE.
  endif

  IF (DKEY == 1) THEN
    OPEN(3, file='data.dat', status='old', iostat=ios, iomsg=errmsg)
    IF (ios /= 0) THEN
      WRITE(*,'(a,a)') 'ERROR: cannot open data.dat: ', TRIM(errmsg)
      STOP 1
    END IF
    expdat(1,1) = 0.
    i = 0
    DO WHILE (expdat(i,1) < Efinal)
      i = i + 1
      READ(3,*,END=50) (expdat(i,n), n=1,4)
      expdat(i,1) = expdat(i,1) - QEx
    ENDDO
50  numdp = i
  ENDIF

END SUBROUTINE PREAD
