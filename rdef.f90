! Rmatrix — multilevel, multichannel R-matrix code
! Copyright (C) 1998-2004, 2026 Grigory V. Rogachev
! Copyright (C) 2026 Alexander Volya
! SPDX-License-Identifier: GPL-3.0-or-later
! This program is free software under the GNU General Public License,
! version 3 or (at your option) any later version; see LICENSE.
! Third-party components: see NOTICE.
SUBROUTINE RDef
! ----------------------------------------------------------------------
! Purpose: Resonance finder. Scans energy, computes S-matrix for
!          level NRDef, finds peak in |1-W_cc|^2 and its FWHM.
! ----------------------------------------------------------------------
  USE channel_data
  USE resonance_data
  USE prep_data
  IMPLICIT NONE

  INTEGER :: nres, incount, outcount, inCh, outCh, nchan, irow
  DOUBLE PRECISION :: delta, ECM, Emax
  DOUBLE PRECISION :: KAW, ZZ, MU, Whitt1, Whitt2
  DOUBLE PRECISION, EXTERNAL :: intgsl
  COMPLEX*16 :: p(MaxCh)

  INTEGER :: IPIV(MaxCh), INFO
  COMPLEX*16 :: WORK(MaxCh)
  INTEGER :: J, L, LPr, PA, MDIM, ichan, zflag

  COMPLEX*16 :: NRPh(0:MaxNL, MaxCh)
  COMPLEX*16 :: CoulPh(0:MaxNL, MaxCh), i
  COMPLEX*16 :: F(0:MaxNL, MaxCh), FP(0:MaxNL, MaxCh)
  COMPLEX*16 :: G(0:MaxNL, MaxCh), GP(0:MaxNL, MaxCh)
  COMPLEX*16 :: SIG(0:MaxNL, MaxCh)
  COMPLEX*16 :: Lmat(MaxCh,MaxCh), Pmat(MaxCh,MaxCh), Wmat(MaxCh,MaxCh)
  COMPLEX*16 :: Omat(MaxCh,MaxCh), twoimat(MaxCh,MaxCh)
  COMPLEX*16 :: RL(MaxCh,MaxCh)
  COMPLEX*16 :: Rmatrix(MaxCh,MaxCh), Unitmat(MaxCh,MaxCh)
  COMPLEX*16 :: UmRL(MaxCh,MaxCh), PUmRL(MaxCh,MaxCh)
  COMPLEX*16 :: PIR(MaxCh,MaxCh), PIRP(MaxCh,MaxCh), PIRPW(MaxCh,MaxCh)
  COMPLEX*16 :: Smat(MaxCh,MaxCh)
  LOGICAL    :: cond
  DOUBLE PRECISION :: fste, Tnew, Told
  DOUBLE PRECISION :: Elevel, csmax, cshalf, Tcomp
  DOUBLE PRECISION :: leftb, rightb

  fste   = 0.001d0
  zflag  = 0
  ! rdchan is a CHANNEL number: the scan starts at its threshold and the
  ! peak is searched in the diagonal element of the first row of level
  ! NRDef that belongs to this channel (the 2004 code used rdchan as a
  ! row index as well)
  ichan  = rdchan
  IF (ichan < 1 .or. ichan > lastch) THEN
    WRITE(*,'(a,i0)') 'RDef: rdchan out of range: ', ichan
    RETURN
  END IF
  irow = 0
  DO incount = 1, rdim(NRDef)
    if (irow == 0 .and. nint(gres(NRDef,incount,1)) == ichan .and. &
        gres(NRDef,incount,4) /= 0.d0 .and. rowok(NRDef,incount)) irow = incount
  ENDDO
  IF (irow == 0) THEN
    DO incount = 1, rdim(NRDef)
      if (irow == 0 .and. gres(NRDef,incount,4) /= 0.d0 .and. &
          rowok(NRDef,incount)) irow = incount
    ENDDO
  END IF
  IF (irow == 0) THEN
    WRITE(*,'(a,i0,a)') 'RDef: level ', NRDef, ' has no active width; skipped'
    RETURN
  END IF
  Emax   = max(Efinal, eres(NRDef,1)) + 5.d0
  Tcomp  = 0.01d0
  Told   = -1.d0
  ecm    = fste - Q(ichan)
  cond   = .TRUE.
  i      = (0.d0, 1.d0)
  cshalf = 10000.d0
  Elevel = 10000.d0
  leftb  = 0.d0
  rightb = 0.d0

  DO WHILE (cond)

    CALL Kinematics(ecm)

    CALL NResPh(NRPh, F, G, FP, GP, SIG, p, ecm)
    CALL CPC(CoulPh)
    i = (0.d0, 1.d0)

    J    = nint(2.d0 * eres(NRDef,2))
    PA   = nint(eres(NRDef,3))
    MDIM = rdim(NRDef)

    ! Zero all working matrices
    Rmatrix = (0.d0,0.d0)
    Omat    = (0.d0,0.d0)
    twoimat = (0.d0,0.d0)
    Lmat    = (0.d0,0.d0)
    PUmRL   = (0.d0,0.d0)
    Pmat    = (0.d0,0.d0)
    RL      = (0.d0,0.d0)
    Smat    = (0.d0,0.d0)
    Unitmat = (0.d0,0.d0)
    UmRL    = (0.d0,0.d0)

    ! Unit and 2i diagonals
    DO incount = 1, MaxCh
      Unitmat(incount,incount) = (1.d0,0.d0)
      twoimat(incount,incount) = (0.d0,2.d0)
    ENDDO

    ! Build matrices for this resonance
    DO incount = 1, MDIM
      DO outcount = incount, MDIM
        inCh  = int(gres(NRDef,incount,1))
        outCh = int(gres(NRDef,outcount,1))
        L     = int(gres(NRDef,incount,2))
        LPr   = int(gres(NRDef,outcount,2))

        delta = eres(NRDef,1) - ecm
        if (abs(delta).lt.1.d-9) delta = sign(1.d-9, delta)

        if (rowok(NRDef,incount) .and. rowok(NRDef,outcount)) then
        Rmatrix(incount,outcount) = &
          Rmatrix(incount,outcount) + &
          gres(NRDef,incount,4)*gres(NRDef,outcount,4) / delta

        if (incount.eq.outcount) then

          Pmat(incount,outcount) = sqrt(p(outCh) / &
            (F(LPr,outCh)*F(LPr,outCh) + &
             G(LPr,outCh)*G(LPr,outCh)))

          Smat(incount,outcount) = &
            p(outCh) * &
            (F(LPr,outCh)*FP(LPr,outCh) + &
             G(LPr,outCh)*GP(LPr,outCh)) / &
            (F(LPr,outCh)*F(LPr,outCh) + &
             G(LPr,outCh)*G(LPr,outCh)) - &
            gres(NRDef,outcount,3)

          Lmat(incount,outcount) = Smat(incount,outcount) + &
            i * Pmat(incount,outcount) * Pmat(incount,outcount)

          Omat(incount,outcount) = &
            exp(i*(CoulPh(LPr,outCh) - NRPh(LPr,outCh)))

          ! closed channel: logarithmic derivative of the Whittaker
          ! function, P = 0 (as in PreTJL)
          if ((ECM+Q(outCh)).lt.0.d0) then
            KAW = dble((-i)*eta(outCh))
            ZZ  = dble(2.d0*(-i)*p(outCh))
            MU  = dble(LPr) + 0.5d0
            Whitt1 = intgsl(KAW, MU, ZZ)
            Whitt2 = intgsl(KAW-1.d0, MU, ZZ)
            Lmat(incount,outcount) = &
              ((KAW-0.5d0*ZZ)*Whitt1 - &
               (MU*MU-(KAW-0.5d0)**2)*Whitt2) / Whitt1 - &
              gres(NRDef,outcount,3)
            Pmat(incount,outcount)    = (0.d0,0.d0)
            twoimat(incount,outcount) = (0.d0,0.d0)
            Omat(incount,outcount)    = (0.d0,0.d0)
          endif

        else
          Pmat(incount,outcount)    = (0.d0,0.d0)
          Lmat(incount,outcount)    = (0.d0,0.d0)
          twoimat(incount,outcount) = (0.d0,0.d0)
          Omat(incount,outcount)    = (0.d0,0.d0)
          Unitmat(incount,outcount) = (0.d0,0.d0)
        endif
        endif

        ! Symmetrize
        Rmatrix(outcount,incount) = Rmatrix(incount,outcount)
        Pmat(outcount,incount)    = Pmat(incount,outcount)
        Lmat(outcount,incount)    = Lmat(incount,outcount)
        twoimat(outcount,incount) = twoimat(incount,outcount)
        Omat(outcount,incount)    = Omat(incount,outcount)
        Unitmat(outcount,incount) = Unitmat(incount,outcount)

      ENDDO
    ENDDO

    ! All matrices ready. Main algebra:

    IF (MDIM.GE.1) THEN

      ! RL = R * L
      CALL ZGEMM('N','N', MDIM, MDIM, MDIM, &
                 (1.d0,0.d0), Rmatrix, MaxCh, &
                 Lmat, MaxCh, (0.d0,0.d0), RL, MaxCh)

      ! UmRL = I - RL
      UmRL(1:MDIM,1:MDIM) = Unitmat(1:MDIM,1:MDIM) - RL(1:MDIM,1:MDIM)

      ! UmRL = (I - RL)^-1
      CALL ZGETRF(MDIM, MDIM, UmRL, MaxCh, IPIV, INFO)
      if (INFO.ne.0) then
        PRINT *, 'ZGETRF failed for 1-RL, INFO=', INFO
        STOP
      endif
      CALL ZGETRI(MDIM, UmRL, MaxCh, IPIV, WORK, MaxCh, INFO)
      if (INFO.ne.0) then
        PRINT *, 'ZGETRI failed for 1-RL, INFO=', INFO
        STOP
      endif

      ! PUmRL = P * (I-RL)^-1
      CALL ZGEMM('N','N', MDIM, MDIM, MDIM, &
                 (1.d0,0.d0), Pmat, MaxCh, &
                 UmRL, MaxCh, (0.d0,0.d0), PUmRL, MaxCh)

      ! PIR = P*(I-RL)^-1 * R
      CALL ZGEMM('N','N', MDIM, MDIM, MDIM, &
                 (1.d0,0.d0), PUmRL, MaxCh, &
                 Rmatrix, MaxCh, (0.d0,0.d0), PIR, MaxCh)

      ! PIRP = P*(I-RL)^-1*R * P
      CALL ZGEMM('N','N', MDIM, MDIM, MDIM, &
                 (1.d0,0.d0), PIR, MaxCh, &
                 Pmat, MaxCh, (0.d0,0.d0), PIRP, MaxCh)

      ! PIRPW = P*(I-RL)^-1*R*P * 2i
      CALL ZGEMM('N','N', MDIM, MDIM, MDIM, &
                 (1.d0,0.d0), PIRP, MaxCh, &
                 twoimat, MaxCh, (0.d0,0.d0), PIRPW, MaxCh)

      ! Wmat = I + PIRPW
      Wmat(1:MDIM,1:MDIM) = Unitmat(1:MDIM,1:MDIM) + PIRPW(1:MDIM,1:MDIM)

    ENDIF

    Tnew = (1-Wmat(irow,irow)) * (1-CONJG(Wmat(irow,irow)))

    if ((Tnew.LE.Told).AND.(Told.GT.Tcomp)) then
      Elevel = ecm - fste
      csmax  = Told
      cshalf = csmax / 2.d0
      Tcomp  = 10000000.d0
      ecm    = 0.d0 - Q(ichan)
    endif

    if ((ecm.lt.Elevel).and.(Tnew.le.cshalf)) leftb = ecm
    if ((ecm.gt.Elevel).and.(Tnew.le.cshalf)) then
      rightb = ecm
      if (Tcomp.eq.10000000.d0) cond = .FALSE.
    endif

    ecm   = ecm + fste
    Told  = Tnew
    if (ecm > Emax) then
      WRITE(*,'(2x,a,i0,a,f8.3,a)') 'RDef: no resonance peak found for level #', &
            NRDef, ' below ', Emax, ' MeV'
      RETURN
    endif

  ENDDO

  write(6,'(2x,a,I2,2x,F4.1,a,F6.3,a,F6.3,a)') &
    'Energy and Width of level #', NRDef, eres(NRDef,2)*eres(NRDef,3), &
    ' are ', Elevel, ' MeV and ', rightb-leftb, ' MeV'

END SUBROUTINE RDef
