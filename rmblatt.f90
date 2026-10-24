! Rmatrix — multilevel, multichannel R-matrix code
! Copyright (C) 1998-2004, 2026 Grigory V. Rogachev
! Copyright (C) 2026 Alexander Volya
! SPDX-License-Identifier: GPL-3.0-or-later
! This program is free software under the GNU General Public License,
! version 3 or (at your option) any later version; see LICENSE.
! Third-party components: see NOTICE.
SUBROUTINE rmblatt(ecm, acm, cshi)
! ----------------------------------------------------------------------
! Purpose: Differential cross section via Blatt-Biedenharn
!          decomposition. dσ/dΩ in mb/sr.
!          ecm - cm energy (MeV), acm - cm angle (degrees)
! ----------------------------------------------------------------------
  USE channel_data
  USE resonance_data
  USE tmat_data
  USE prep_data
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN)    :: ecm
  DOUBLE PRECISION, INTENT(INOUT) :: acm
  DOUBLE PRECISION, INTENT(OUT)   :: cshi

  INTEGER          :: chindex, chx, incount, outcount, lang, Jtot, iang, nang
  DOUBLE PRECISION :: sigma, sigI, sigII, sigIII, BBB
  DOUBLE PRECISION :: cscms(MaxCh), LP(0:2*MaxNL), cshi_loc
  DOUBLE PRECISION :: Bcoef(1:MaxCh,1:MaxCh,0:2*MaxNL), Bvalue, TCS, Btot
  DOUBLE PRECISION :: MainCl, AnCos, ClSin, Sfactor
  INTEGER          :: nchan, ChSpin, ChSpPr, nch2pr
  COMPLEX*16       :: ClExp, ClAmpl, i

  i = cmplx(0.d0, 1.d0)

  ! Pre-compute Blatt-Biedenharn B-coefficients. The Legendre expansion
  ! of dsigma/dOmega runs to 2*Lmax (the 2004 code truncated it at Lmax,
  ! which distorts angular distributions once partial waves above Lmax/2
  ! contribute and can make the cross section negative).
  DO incount = 1, lastch
    DO outcount = incount, lastch
      if ((FS(incount).EQ.1).and.(FS(outcount).EQ.ExitCh)) then
        ChSpin = nint(2.d0*ChSp(incount))
        ChSpPr = nint(2.d0*ChSp(outcount))
        DO lang = 0, 2*Lmax
          Bcoef(outcount,incount,lang) = &
              BBB(outcount, ChSpPr, incount, ChSpin, lang)
          Bcoef(incount,outcount,lang) = Bcoef(outcount,incount,lang)
        ENDDO
      endif
    ENDDO
  ENDDO

  IF (acm.GT.0.0d0) THEN
    Astart = acm
    Afinal = acm
    Astep  = acm
  ENDIF

  ! angle loop (integer form of the former real DO loop: same iteration
  ! count MAX(INT((Afinal-Astart+Astep)/Astep),0), same incremented angles)
  nang = max(int((Afinal - Astart + Astep)/Astep), 0)
  acm  = Astart
  DO iang = 1, nang

    AnCos = DCOS(acm*d2r)

    IF (acm.LT.1.0d0) THEN
      PRINT *, 'Too small c.m.s. scattering angle!'
      PRINT *, 'Elastic scattering cannot be calculated.'
      STOP
    ENDIF

    ClSin  = DSIN(0.5d0*acm*d2r)
    MainCl = (ETA(1)/dsqrt(4.d0*PI)) * (1.0d0/ClSin)**2.
    ClExp  = CDEXP(-2.d0*i*ETA(1)*DLOG(ClSin))
    ClAmpl = -MainCl * ClExp

    sigma = 0.d0
    Btot  = 0.d0

    ! Legendre polynomials P_l(cos θ) via 3-term recursion
    LP(0) = 1.d0
    LP(1) = AnCos
    do lang = 2, 2*Lmax
      LP(lang) = ((2*lang-1)*AnCos*LP(lang-1) - (lang-1)*LP(lang-2)) / lang
    enddo

    chindex = ExitCh
    ! first channel of the exit group (ExitCh numbers final-state groups,
    ! while Q, eta, Zp, Zt are per channel; the 2004 code mixed the two)
    chx = 1
    DO incount = lastch, 1, -1
      if (FS(incount).EQ.chindex) chx = incount
    ENDDO

    DO incount = 1, lastch
      DO outcount = 1, lastch
        sigI   = 0.d0
        sigII  = 0.d0
        sigIII = 0.d0
        if ((FS(incount).EQ.1).and.(FS(outcount).EQ.chindex)) then
          ChSpin = nint(2.d0*ChSp(incount))
          ChSpPr = nint(2.d0*ChSp(outcount))
          Btot   = Btot + Bcoef(outcount,incount,0)

          DO lang = 0, 2*Lmax
            sigII = sigII + Bcoef(outcount,incount,lang)*LP(lang)
          ENDDO
          sigII = sigII / (WAVNUM(incount)*WAVNUM(incount))

          IF ((chindex.EQ.1).and.(incount.EQ.outcount)) THEN
            DO lang = 0, Lmax
              DO Jtot = ABS(2*lang-ChSpin),(2*lang+ChSpin),2
                if (Jtot > MaxTJL) cycle
                sigIII = sigIII + (Jtot+1)*2 * &
                  dreal(i * &
                    TJL(Jtot,incount,outcount,lang,lang) * &
                    CONJG(ClAmpl) * LP(lang))
              ENDDO
            ENDDO
            sigIII = sigIII * sqrt(PI) / &
                     (2.d0*WAVNUM(1)*WAVNUM(1))
            sigI   = (ChSpin+1)*PI*(CDABS(ClAmpl)**2.) / &
                     (WAVNUM(incount)*WAVNUM(incount))
          ENDIF
        endif

        sigma = sigma + sigI + sigII + sigIII

      ENDDO
    ENDDO

    if (Inverse) then
      cscms(chindex) = sigma*10.d0/(2.*Sp(inv_out)+1) / &
                       (2.*St(inv_out)+1) * &
                       wavnum(1)*wavnum(1) / &
                       wavnum(inv_out)/wavnum(inv_out)
    else
      cscms(chindex) = sigma*10.d0/(2.*Sp(1)+1)/(2.*St(1)+1)
    endif
    cshi  = cscms(chindex)
    sigma = 0.d0

    write(1,'(2x,3ES15.7)') ecm, acm, cscms(chindex)

    if (Inverse) then
      Bvalue = 10.d0/wavnum(inv_out)/wavnum(inv_out) / &
               (2.*Sp(inv_out)+1)/(2.*St(inv_out)+1)*Btot
    else
      Bvalue = 10.d0/wavnum(1)/wavnum(1) / &
               (2.*Sp(1)+1)/(2.*St(1)+1)*Btot
    endif

    if (Inverse) then
      Sfactor = 4.d0*PI*Bvalue*(ecm+Q(chx)) * &
                exp(2.d0*PI*eta(chx))*0.001d0
    else
      Sfactor = 4.d0*PI*Bvalue*(ecm+Q(1)) * &
                exp(2.d0*PI*eta(1))*0.001d0
    endif

    write(2,*) ecm, 4.d0*PI*Bvalue*0.001d0, Sfactor

    acm = acm + Astep
  ENDDO

END SUBROUTINE rmblatt
