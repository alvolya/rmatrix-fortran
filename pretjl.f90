! Rmatrix — multilevel, multichannel R-matrix code
! Copyright (C) 1998-2004, 2026 Grigory V. Rogachev
! Copyright (C) 2026 Alexander Volya
! SPDX-License-Identifier: GPL-3.0-or-later
! This program is free software under the GNU General Public License,
! version 3 or (at your option) any later version; see LICENSE.
! Third-party components: see NOTICE.
SUBROUTINE PreTJL(ECM)
! ----------------------------------------------------------------------
! Purpose: Core R-matrix engine. Builds T-matrix TJL for all J^pi.
! For each (J, parity) block: constructs R, L, P, O matrices from
! resonance data and Coulomb wave functions, then computes
! S = O * (I + P*(I-RL)^-1*R*P*2i) * O  and stores T = e^{2i*sigma} - S
! ----------------------------------------------------------------------
  USE channel_data
  USE resonance_data
  USE tmat_data
  USE prep_data
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: ECM

  INTEGER :: nres, incount, outcount, inCh, outCh, refres
  DOUBLE PRECISION :: delta, KAW, ZZ, MU, Whitt1, Whitt2
  DOUBLE PRECISION, EXTERNAL :: intgsl
  COMPLEX*16 :: p(MaxCh)

  INTEGER :: IPIV(MaxCh), INFO
  COMPLEX*16 :: WORK(MaxCh)
  INTEGER :: J, L, LPr, Jfirst, Jlast, S2max, PA, MDIM

  COMPLEX*16 :: NRPh(0:MaxNL, MaxCh)
  COMPLEX*16 :: CoulPh(0:MaxNL, MaxCh), i
  COMPLEX*16 :: F(0:MaxNL, MaxCh), FP(0:MaxNL, MaxCh)
  COMPLEX*16 :: G(0:MaxNL, MaxCh), GP(0:MaxNL, MaxCh)
  COMPLEX*16 :: SIG(0:MaxNL, MaxCh)
  COMPLEX*16 :: Lmat(MaxCh,MaxCh), Pmat(MaxCh,MaxCh), Wmat(MaxCh,MaxCh)
  COMPLEX*16 :: Omat(MaxCh,MaxCh), twoimat(MaxCh,MaxCh)
  COMPLEX*16 :: RL(MaxCh,MaxCh)
  COMPLEX*16 :: PIR(MaxCh,MaxCh), PIRP(MaxCh,MaxCh), PIRPW(MaxCh,MaxCh)
  COMPLEX*16 :: Rmatrix(MaxCh,MaxCh), Unitmat(MaxCh,MaxCh)
  COMPLEX*16 :: UmRL(MaxCh,MaxCh), PUmRL(MaxCh,MaxCh)
  COMPLEX*16 :: OW(MaxCh,MaxCh), OWO(MaxCh,MaxCh)
  COMPLEX*16 :: OWOPH(MaxCh,MaxCh), PHmat(MaxCh,MaxCh)
  COMPLEX*16 :: OmatPH(MaxCh,MaxCh)
  COMPLEX*16 :: Smat



  CALL NResPh(NRPh, F, G, FP, GP, SIG, p, ecm)
  CALL CPC(CoulPh)

  i = (0.d0, 1.d0)

  ! Zero-initialise full T-matrix
  DO J = 0, MaxTJL
    DO inCh = 1, lastch
      DO outCh = 1, lastch
        DO L = 0, Lmax
          DO LPr = 0, Lmax
            TJL(J,inCh,outCh,L,LPr) = (0.d0,0.d0)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  ! J range: parity of 2J from the channel spins, maximum from the
  ! largest channel spin of ANY channel (the 2004 code used channel 1
  ! only and silently skipped levels beyond its range)
  S2max = 0
  DO inCh = 1, lastch
    S2max = max(S2max, nint(2.d0*ChSp(inCh)))
  ENDDO
  Jfirst = mod(S2max, 2)
  Jlast  = min(2*Lmax + S2max, MaxTJL)

  DO J = Jfirst, Jlast, 2
    DO PA = -1, 1, 2

      ! Zero all working matrices
      Rmatrix = (0.d0,0.d0)
      Omat    = (0.d0,0.d0)
      OmatPH  = (0.d0,0.d0)
      OWOPH   = (0.d0,0.d0)
      twoimat = (0.d0,0.d0)
      Lmat    = (0.d0,0.d0)
      Pmat    = (0.d0,0.d0)
      PHmat   = (0.d0,0.d0)
      RL      = (0.d0,0.d0)
      UmRL    = (0.d0,0.d0)
      PUmRL   = (0.d0,0.d0)
      Unitmat = (0.d0,0.d0)
      PIR     = (0.d0,0.d0)
      PIRP    = (0.d0,0.d0)
      PIRPW   = (0.d0,0.d0)

      ! Unit and 2i diagonal matrices
      DO incount = 1, MaxCh
        Unitmat(incount,incount) = (1.d0,0.d0)
        twoimat(incount,incount) = (0.d0,2.d0)
      ENDDO

      MDIM = 0

      ! Build R-matrix and channel matrices for matching J, PA
      DO nres = 1, lastre
        if ((nint(2*eres(nres,2)).eq.J) .and. &
             (PA.eq.nint(eres(nres,3)))) then

          MDIM   = rdim(nres)
          refres = nres

          DO incount = 1, MDIM
            DO outcount = incount, MDIM
              inCh  = nint(gres(nres,incount,1))
              outCh = nint(gres(nres,outcount,1))
              L     = nint(gres(nres,incount,2))
              LPr   = nint(gres(nres,outcount,2))

              ! pole term; guard only against an exact coincidence E = E_lambda
              ! (the 2004 code added 1 keV whenever |E_lambda - E| < 1 keV, which
              ! moved the pole onto the grid point 1 keV above it)
              delta = eres(nres,1) - ecm
              if (abs(delta).lt.1.d-9) delta = sign(1.d-9, delta)

              ! rows failing the parity or triangle rule (REREAD) do not
              ! take part: R, L, P and Omega stay zero, Unitmat keeps 1
              if (rowok(nres,incount) .and. rowok(nres,outcount)) then

                Rmatrix(incount,outcount) = &
                  Rmatrix(incount,outcount) + &
                  gres(nres,incount,4)*gres(nres,outcount,4) / delta

                if (incount.eq.outcount) then

                  Pmat(incount,outcount) = sqrt(p(outCh) / &
                    (F(LPr,outCh)*F(LPr,outCh) + &
                     G(LPr,outCh)*G(LPr,outCh)))

                  Omat(incount,outcount) = &
                    exp(i*(CoulPh(L,outCh) - NRPh(LPr,outCh)))

                  OmatPH(incount,outcount) = &
                    exp(i*(-NRPh(LPr,outCh)))


                  Smat = p(outCh) * &
                    (F(LPr,outCh)*FP(LPr,outCh) + &
                     G(LPr,outCh)*GP(LPr,outCh)) / &
                    (F(LPr,outCh)*F(LPr,outCh) + &
                     G(LPr,outCh)*G(LPr,outCh)) - &
                    gres(nres,outcount,3)

                  Lmat(incount,outcount) = Smat + i*p(outCh) / &
                    (F(LPr,outCh)*F(LPr,outCh) + &
                     G(LPr,outCh)*G(LPr,outCh))

                  ! Closed-channel logarithmic derivative from the
                  ! Whittaker W function.
                  if ((ECM+Q(outCh)).lt.0.d0) then
                    KAW = dble((-i)*eta(outCh))
                    ZZ  = dble(2.d0*(-i)*p(outCh))
                    MU  = dble(LPr) + 0.5d0
                    Whitt1 = intgsl(KAW, MU, ZZ)
                    Whitt2 = intgsl(KAW-1.d0, MU, ZZ)
                    Lmat(incount,outcount) = &
                      ((KAW-0.5d0*ZZ)*Whitt1 - &
                       (MU*MU-(KAW-0.5d0)**2)*Whitt2) / Whitt1 - &
                      gres(nres,outcount,3)
                    Pmat(incount,outcount)    = (0.d0,0.d0)
                    twoimat(incount,outcount) = (0.d0,0.d0)
                    Omat(incount,outcount)    = (0.d0,0.d0)
                  endif

                else
                  Pmat(incount,outcount)    = (0.d0,0.d0)
                  Lmat(incount,outcount)    = (0.d0,0.d0)
                  twoimat(incount,outcount) = (0.d0,0.d0)
                  Omat(incount,outcount)    = (0.d0,0.d0)
                  OmatPH(incount,outcount)  = (0.d0,0.d0)
                  Unitmat(incount,outcount) = (0.d0,0.d0)
                endif
              endif

              ! Symmetrize
              Rmatrix(outcount,incount) = Rmatrix(incount,outcount)
              Pmat(outcount,incount)    = Pmat(incount,outcount)
              Lmat(outcount,incount)    = Lmat(incount,outcount)
              twoimat(outcount,incount) = twoimat(incount,outcount)
              Omat(outcount,incount)    = Omat(incount,outcount)
              OmatPH(outcount,incount)  = OmatPH(incount,outcount)
              Unitmat(outcount,incount) = Unitmat(incount,outcount)

            ENDDO
          ENDDO

        endif
      ENDDO

      ! All matrices are ready. Main algebra:

      IF (MDIM.GE.1) THEN

        ! RL = R * L
        CALL ZGEMM('N','N', MDIM, MDIM, MDIM, &
                   (1.d0,0.d0), Rmatrix, MaxCh, &
                   Lmat, MaxCh, (0.d0,0.d0), RL, MaxCh)

        ! UmRL = I - RL
        UmRL(1:MDIM,1:MDIM) = Unitmat(1:MDIM,1:MDIM) - RL(1:MDIM,1:MDIM)

        ! UmRL = (I - RL)^-1  (in-place via LU factorisation)
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

        ! Wmat = I + P*(I-RL)^-1*R*P*2i
        Wmat(1:MDIM,1:MDIM) = Unitmat(1:MDIM,1:MDIM) + PIRPW(1:MDIM,1:MDIM)

        ! OW = O * Wmat
        CALL ZGEMM('N','N', MDIM, MDIM, MDIM, &
                   (1.d0,0.d0), Omat, MaxCh, &
                   Wmat, MaxCh, (0.d0,0.d0), OW, MaxCh)

        ! OWO = O * Wmat * O  (S-matrix)
        CALL ZGEMM('N','N', MDIM, MDIM, MDIM, &
                   (1.d0,0.d0), OW, MaxCh, &
                   Omat, MaxCh, (0.d0,0.d0), OWO, MaxCh)

        ! OWOPH = OmatPH * Wmat  (nuclear phase only)
        CALL ZGEMM('N','N', MDIM, MDIM, MDIM, &
                   (1.d0,0.d0), OmatPH, MaxCh, &
                   Wmat, MaxCh, (0.d0,0.d0), OWOPH, MaxCh)

        ! PHmat = OmatPH * Wmat * OmatPH  (nuclear phase matrix)
        CALL ZGEMM('N','N', MDIM, MDIM, MDIM, &
                   (1.d0,0.d0), OWOPH, MaxCh, &
                   OmatPH, MaxCh, (0.d0,0.d0), PHmat, MaxCh)

        ! Store T-matrix and write phase shifts for selected wave
        DO incount = 1, MDIM
          DO outcount = 1, MDIM
            inCh  = int(gres(refres,incount,1))
            outCh = int(gres(refres,outcount,1))
            L     = int(gres(refres,incount,2))
            LPr   = int(gres(refres,outcount,2))
            if (.not.(rowok(refres,incount) .and. rowok(refres,outcount))) cycle
            if ((inCh.eq.outCh).and.(L.eq.LPr)) then
              TJL(J,inCh,outCh,L,LPr) = &
                exp(2*i*CoulPh(L,inCh)) - OWO(incount,outcount)
            else
              TJL(J,inCh,outCh,L,LPr) = -OWO(incount,outcount)
            endif


          ENDDO
        ENDDO

      ENDIF

    ENDDO  ! PA
  ENDDO    ! J

END SUBROUTINE PreTJL
