! Rmatrix — multilevel, multichannel R-matrix code
! Copyright (C) 1998-2004, 2026 Grigory V. Rogachev
! Copyright (C) 2026 Alexander Volya
! SPDX-License-Identifier: GPL-3.0-or-later
! This program is free software under the GNU General Public License,
! version 3 or (at your option) any later version; see LICENSE.
! Third-party components: see NOTICE.
SUBROUTINE REREAD
! ----------------------------------------------------------------------
! Purpose: Read the resonance parameter file (.par) and validate it.
!
! File format: four header lines, then per level
!     E_lambda  J  parity
!     ChN  L  Bc  gamma        (one row per channel-spin/L combination;
!    -ChN  L  Bc  gamma         a negative ChN marks the last row)
! E in MeV (c.m., relative to the first channel threshold), gamma in
! MeV^1/2 (reduced-width amplitude), Bc dimensionless.
!
! Validation (fatal unless stated): level count <= MaxRes, rows per
! level <= MaxCh, 2J <= MaxTJL, parity +-1, channel number in range,
! L <= Lmax, and the selection rules of every row — parity
! pi_p pi_t (-1)^L = parity and the triangle |J-s| <= L <= J+s with
! the channel spin s of the channel card.  A row violating them with a
! non-zero width is an error; with zero width it is ignored (rowok =
! .FALSE.; the 2004 code gave such rows a hard-sphere phase).  All
! levels of one J^pi block must list the same (ChN, L) rows in the
! same order with the same Bc: the engine assembles the block
! positionally.
! ----------------------------------------------------------------------
  USE channel_data
  USE resonance_data
  USE io_filenames
  IMPLICIT NONE

  INTEGER            :: k, n, m, ccount, ios, ch, L, J2, S2, nbad, nzero
  CHARACTER(LEN=256) :: fname, errmsg
  DOUBLE PRECISION   :: ChN, Bc, G, Lr, e1, e2, e3
  LOGICAL            :: ok

  IF (LEN_TRIM(par_file) == 0) THEN
    WRITE(6,'(2x,a)', advance='no') 'Parameter file name: '
    READ(5,*) par_file
  END IF
  fname = par_file

  OPEN(1, file=TRIM(fname), status='old', iostat=ios, iomsg=errmsg)
  IF (ios /= 0) THEN
    WRITE(*,'(a,a)') 'ERROR: cannot open parameter file: ', TRIM(errmsg)
    STOP 1
  END IF

  READ(1,*)
  READ(1,*)
  READ(1,*)
  READ(1,*)
  lastre = 0
  DO
    READ(1,*,iostat=ios) e1, e2, e3
    IF (ios /= 0) EXIT
    IF (lastre >= MaxRes) THEN
      WRITE(*,'(a,i0,a)') 'ERROR: more than MaxRes = ', MaxRes, ' levels'
      STOP 1
    END IF
    lastre = lastre + 1
    eres(lastre,1) = e1
    eres(lastre,2) = e2
    eres(lastre,3) = e3
    ccount = 0
    ChN    = 0
    DO WHILE (ChN >= 0)
      READ(1,*,iostat=ios) ChN, Lr, Bc, G
      IF (ios /= 0) THEN
        WRITE(*,'(a,i0)') 'ERROR: unexpected end of parameter file in level ', lastre
        STOP 1
      END IF
      IF (ccount >= MaxCh) THEN
        WRITE(*,'(a,i0,a,i0)') 'ERROR: more than MaxCh = ', MaxCh, &
                               ' rows in level ', lastre
        STOP 1
      END IF
      ccount                   = ccount + 1
      gres(lastre,ccount,1)    = ABS(ChN)
      gres(lastre,ccount,2)    = Lr
      gres(lastre,ccount,3)    = Bc
      gres(lastre,ccount,4)    = G
    ENDDO
    rdim(lastre) = ccount
  ENDDO
  CLOSE(UNIT=1)

  ! ---- validation ----------------------------------------------------
  rowok = .FALSE.
  nbad  = 0
  nzero = 0
  DO k = 1, lastre
    J2 = nint(2.d0*eres(k,2))
    IF (abs(2.d0*eres(k,2) - J2) > 1.d-6 .or. J2 < 0) THEN
      WRITE(*,'(a,i0,a,f8.3)') 'ERROR: level ', k, ' has J = ', eres(k,2)
      STOP 1
    END IF
    IF (J2 > MaxTJL) THEN
      WRITE(*,'(a,i0,a,f6.1,a,i0)') 'ERROR: level ', k, ' has J = ', &
            eres(k,2), ' but 2J must not exceed MaxTJL = ', MaxTJL
      STOP 1
    END IF
    IF (abs(nint(eres(k,3))) /= 1 .or. abs(eres(k,3) - nint(eres(k,3))) > 1.d-6) THEN
      WRITE(*,'(a,i0,a,f8.3)') 'ERROR: level ', k, ' has parity ', eres(k,3)
      STOP 1
    END IF
    DO n = 1, rdim(k)
      ch = nint(gres(k,n,1))
      L  = nint(gres(k,n,2))
      IF (ch < 1 .or. ch > lastch) THEN
        WRITE(*,'(a,i0,a,i0,a,i0)') 'ERROR: level ', k, ' row ', n, &
              ': channel number out of range: ', ch
        STOP 1
      END IF
      IF (L < 0 .or. L > Lmax) THEN
        WRITE(*,'(a,i0,a,i0,a,i0,a,i0)') 'ERROR: level ', k, ' row ', n, &
              ': L = ', L, ' exceeds Lmax = ', Lmax
        STOP 1
      END IF
      S2 = nint(2.d0*ChSp(ch))
      ok = (Pp(ch)*Pt(ch)*(-1)**L == nint(eres(k,3)))        ! parity
      ok = ok .and. (mod(J2+S2, 2) == 0)                      ! J, s both (half-)integer
      ok = ok .and. (2*L >= abs(J2-S2)) .and. (2*L <= J2+S2)  ! triangle
      rowok(k,n) = ok
      IF (.not. ok) THEN
        IF (gres(k,n,4) /= 0.d0) THEN
          nbad = nbad + 1
          WRITE(*,'(a,i0,a,i0,a,f5.1,a,i2,a,i0,a,i0,a,f5.1,a,es10.3)') &
            'ERROR: level ', k, ' row ', n, ' (J = ', eres(k,2), ', parity ', &
            nint(eres(k,3)), '): channel ', ch, ' L = ', L, ' s = ', ChSp(ch), &
            ' violates the selection rules but has gamma = ', gres(k,n,4)
        ELSE
          nzero = nzero + 1
        END IF
      END IF
    ENDDO
  ENDDO
  IF (nbad > 0) THEN
    WRITE(*,'(a)') 'ERROR: parity/triangle-forbidden widths in the parameter file'
    STOP 1
  END IF
  IF (nzero > 0) WRITE(*,'(2x,a,i0,a)') 'Note: ', nzero, &
      ' zero-width rows violate the selection rules and are ignored'

  ! all levels of one J^pi block: identical (ChN, L, Bc) rows
  DO k = 1, lastre
    DO m = 1, k-1
      IF (nint(2.d0*eres(m,2)) == nint(2.d0*eres(k,2)) .and. &
          nint(eres(m,3)) == nint(eres(k,3))) THEN
        ok = (rdim(m) == rdim(k))
        IF (ok) THEN
          DO n = 1, rdim(k)
            IF (nint(gres(m,n,1)) /= nint(gres(k,n,1)) .or. &
                nint(gres(m,n,2)) /= nint(gres(k,n,2)) .or. &
                (BcMode == 0 .and. abs(gres(m,n,3) - gres(k,n,3)) > 1.d-9)) ok = .FALSE.
          ENDDO
        END IF
        IF (.not. ok) THEN
          WRITE(*,'(a,i0,a,i0,a)') 'ERROR: levels ', m, ' and ', k, &
            ' have the same J and parity but different (ChN, L, Bc) rows;'
          WRITE(*,'(a)') '       all levels of one J^pi block must list the same rows'
          STOP 1
        END IF
      END IF
    ENDDO
  ENDDO

END SUBROUTINE REREAD
