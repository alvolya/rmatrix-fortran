! Rmatrix — multilevel, multichannel R-matrix code
! Copyright (C) 1998-2004, 2026 Grigory V. Rogachev
! Copyright (C) 2026 Alexander Volya
! SPDX-License-Identifier: GPL-3.0-or-later
! This program is free software under the GNU General Public License,
! version 3 or (at your option) any later version; see LICENSE.
! Third-party components: see NOTICE.
DOUBLE PRECISION FUNCTION BBB(ChNPr, ChSpPr, ChN, ChSpin, lang)
! ----------------------------------------------------------------------
! Purpose: Blatt-Biedenharn B-coefficient.
! B_lang = (1/4) * sum_{L1,J1,L2,J2} Z(L1,J1,L2,J2,S,lang) *
!                                      Z(L1Pr,J1,L2Pr,J2,S',lang) *
!                                      Re[ T(J1,ch,ch',L1,L1Pr) *
!                                          T*(J2,ch,ch',L2,L2Pr) ]
! lang runs to 2*Lmax (the Legendre order of a product of two partial
! waves up to Lmax; the 2004 code stopped at Lmax).  Z-coefficients
! are cached on first call (energy-independent) in an array sized by
! the actual Lmax and the largest channel spin.
! ----------------------------------------------------------------------
  USE channel_data
  USE tmat_data
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ChNPr, ChSpPr, ChN, ChSpin, lang

  INTEGER          :: J1, J2, L1, L2, L1Pr, L2Pr, sdet, SSpin, LLL, n
  DOUBLE PRECISION :: ZZZ, ZZZPr, thesign
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: zcoef(:,:,:,:,:,:)
  INTEGER, SAVE    :: S2max = 0
  COMPLEX*16       :: Btmp

  DOUBLE PRECISION, EXTERNAL :: ZBLAT

  Btmp = (0.d0, 0.d0)

  ! Build Z-coefficient cache on first call (geometry only)
  if (.not. allocated(zcoef)) then
    S2max = 0
    do n = 1, lastch
      S2max = max(S2max, nint(2.d0*ChSp(n)))
    enddo
    allocate(zcoef(0:Lmax, 0:MaxTJL, 0:Lmax, 0:MaxTJL, 0:S2max, 0:2*Lmax))
    do L1 = 0, Lmax
      do J1 = 0, MaxTJL
        do L2 = 0, Lmax
          do J2 = 0, MaxTJL
            do SSpin = 0, S2max
              do LLL = 0, 2*Lmax
                zcoef(L1,J1,L2,J2,SSpin,LLL) = &
                    ZBLAT(L1,J1,L2,J2,SSpin,LLL)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  endif

  DO L1 = 0, Lmax
    DO J1 = ABS(2*L1-ChSpin), (2*L1+ChSpin), 2
      if (J1 > MaxTJL) cycle
      DO L2 = ABS(L1-lang), min(L1+lang, Lmax)
        DO J2 = ABS(2*L2-ChSpin), (2*L2+ChSpin), 2
          if (J2 > MaxTJL) cycle
          ZZZ = zcoef(L1,J1,L2,J2,ChSpin,lang)
          if (ZZZ == 0.d0) cycle
          DO L1Pr = ABS(J1-ChSpPr)/2, min((J1+ChSpPr)/2, Lmax)
            DO L2Pr = ABS(J2-ChSpPr)/2, min((J2+ChSpPr)/2, Lmax)
              ZZZPr = zcoef(L1Pr,J1,L2Pr,J2,ChSpPr,lang)
              Btmp = Btmp + ZZZ*ZZZPr * &
                dble(TJL(J1,ChN,ChNPr,L1,L1Pr) * &
                CONJG(TJL(J2,ChN,ChNPr,L2,L2Pr)))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  sdet    = ABS(ChSpin - ChSpPr) / 2
  thesign = (-1.d0)**sdet
  BBB     = 0.25d0 * dble(Btmp) * thesign

END FUNCTION BBB


DOUBLE PRECISION FUNCTION ZBLAT(L1, J1, L2, J2, ChSpin, lang)
! ----------------------------------------------------------------------
! Purpose: Z-coefficient for Blatt-Biedenharn decomposition.
! Z = sqrt[(2L1+1)(2L2+1)(J1+1)(J2+1)] * CG(L1,0,L2,0;lang,0)
!        * W(L1, J1/2, L2, J2/2; ChSpin/2, lang)
! ----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: L1, J1, L2, J2, ChSpin, lang
  DOUBLE PRECISION :: Ztmp, CG, RACAH
  DOUBLE PRECISION, EXTERNAL :: CG_m0, RACAHW

  Ztmp  = SQRT((2.d0*L1+1.d0)*(2.d0*L2+1.d0)*(J1+1.d0)*(J2+1.d0))
  CG    = CG_m0(L1, L2, lang)
  ! RACAHW takes arguments as 2*j integers: W(L1, J1/2, L2, J2/2; ChSpin/2, lang)
  RACAH = RACAHW(2*L1, J1, 2*L2, J2, ChSpin, 2*lang)

  ZBLAT = Ztmp * CG * RACAH

END FUNCTION ZBLAT


! ----------------------------------------------------------------------
! Log-factorial: ln(n!) using Fortran 2008 log_gamma intrinsic
! ----------------------------------------------------------------------
DOUBLE PRECISION FUNCTION lnfact(n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  if (n <= 0) then
    lnfact = 0.d0
  else
    lnfact = log_gamma(dble(n+1))
  endif
END FUNCTION lnfact


! ----------------------------------------------------------------------
! CG_m0: Clebsch-Gordan coefficient <j1,0;j2,0|J,0> (all integer j).
! Implements the CERNLIB DCLEBG algorithm for m1=m2=M=0.
! Returns 0 if triangle inequality fails or j1+j2+J is odd.
! ----------------------------------------------------------------------
DOUBLE PRECISION FUNCTION CG_m0(j1, j2, J)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: j1, j2, J
  INTEGER :: ia, ib, ic, k0, k1, k2, k3, k4, k5, k6, k7, k8, k9
  INTEGER :: k10, k11, ka, kz, k
  DOUBLE PRECISION :: w, s, q
  DOUBLE PRECISION, EXTERNAL :: lnfact

  CG_m0 = 0.d0
  ia = 2*j1;  ib = 2*j2;  ic = 2*J

  ! Triangle inequality and parity
  if (ia+ib < ic .or. ia+ic < ib .or. ib+ic < ia) return
  if (mod(ia+ib+ic, 2) /= 0) return

  ! K-values (CERNLIB notation, ix=iy=iz=0 since m1=m2=M=0)
  k0  = ia+ib+ic+2
  k1  = ia+ib-ic;  k2  = ia-ib+ic;  k3  = ib+ic-ia
  k4  = ia;        k5  = ib;         k6  = ic
  k7  = ia;        k8  = ib;         k9  = ic
  k10 = ib-ic;     k11 = ia-ic
  ka  = max(0, k10, k11)
  kz  = min(k1, k5, k7)
  if (ka > kz) return

  w = 0.5d0*(lnfact(k1/2)+lnfact(k2/2)+lnfact(k3/2) &
            +lnfact(k4/2)+lnfact(k5/2)+lnfact(k6/2) &
            +lnfact(k7/2)+lnfact(k8/2)+lnfact(k9/2) &
            -lnfact(k0/2))
  s = 0.d0
  q = dble((-1)**(ka/2))
  do k = ka, kz, 2
    s = s + q * exp(w - (lnfact(k/2) + lnfact((k1-k)/2) &
                        + lnfact((k5-k)/2) + lnfact((k7-k)/2) &
                        + lnfact((k-k10)/2) + lnfact((k-k11)/2)))
    q = -q
  enddo
  CG_m0 = sqrt(dble(ic+1)) * s

END FUNCTION CG_m0


! ----------------------------------------------------------------------
! RACAHW: Racah W coefficient W(A,B,C,X;Y,Z) with half-integer args.
! Arguments passed as integers = 2*j:
!   ia1=2A, ib1=2B, ic1=2C, ix1=2X, iy1=2Y, iz1=2Z
! Implements the CERNLIB DRACAW algorithm exactly, including its
! argument permutation: A=A1, B=B1, C=Y1, X=X1, Y=C1, Z=Z1.
! ----------------------------------------------------------------------
DOUBLE PRECISION FUNCTION RACAHW(ia1, ib1, ic1, ix1, iy1, iz1)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ia1, ib1, ic1, ix1, iy1, iz1
  INTEGER :: ia, ib, ic, ix, iy, iz
  INTEGER :: k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12
  INTEGER :: iabc, iayz, ixbz, ixyc, iabxy, ibcyz, icazx, ka, kz, j1, k
  DOUBLE PRECISION :: w, s, q
  DOUBLE PRECISION, EXTERNAL :: lnfact

  RACAHW = 0.d0

  ! CERNLIB DRACAW argument permutation: C gets Y1, Y gets C1
  ia = ia1;  ib = ib1;  ic = iy1;  ix = ix1;  iy = ic1;  iz = iz1

  if (ia<0 .or. ib<0 .or. ic<0 .or. ix<0 .or. iy<0 .or. iz<0) return

  iabc = ia+ib+ic;  iayz = ia+iy+iz
  ixbz = ix+ib+iz;  ixyc = ix+iy+ic
  if (mod(iabc,2)/=0 .or. mod(iayz,2)/=0) return
  if (mod(ixbz,2)/=0 .or. mod(ixyc,2)/=0) return

  ! All four triangles
  k1=ia+ib-ic;  k2=ia-ib+ic;  k3=ib+ic-ia
  if (k1<0 .or. k2<0 .or. k3<0) return
  k4=ia+iy-iz;  k5=ia-iy+iz;  k6=iy+iz-ia
  if (k4<0 .or. k5<0 .or. k6<0) return
  k7=ix+ib-iz;  k8=ix-ib+iz;  k9=ib+iz-ix
  if (k7<0 .or. k8<0 .or. k9<0) return
  k10=ix+iy-ic; k11=ix-iy+ic; k12=iy+ic-ix
  if (k10<0 .or. k11<0 .or. k12<0) return

  iabxy = ia+ib+ix+iy;  ibcyz = ib+ic+iy+iz;  icazx = ic+ia+iz+ix
  ka = max(iabc, iayz, ixbz, ixyc)
  kz = min(iabxy, ibcyz, icazx)
  j1 = ka + iabxy    ! Extra Racah phase vs 6j
  if (ka > kz) return

  w = 0.5d0*(lnfact(k1/2)+lnfact(k2/2)+lnfact(k3/2)-lnfact(iabc/2+1) &
            +lnfact(k4/2)+lnfact(k5/2)+lnfact(k6/2)-lnfact(iayz/2+1) &
            +lnfact(k7/2)+lnfact(k8/2)+lnfact(k9/2)-lnfact(ixbz/2+1) &
            +lnfact(k10/2)+lnfact(k11/2)+lnfact(k12/2)-lnfact(ixyc/2+1))
  s = 0.d0
  q = dble((-1)**(j1/2))
  do k = ka, kz, 2
    s = s + q * exp(w + lnfact(k/2+1) &
                   - (lnfact((k-iabc)/2) + lnfact((k-iayz)/2) &
                    + lnfact((k-ixbz)/2) + lnfact((k-ixyc)/2) &
                    + lnfact((iabxy-k)/2) + lnfact((ibcyz-k)/2) &
                    + lnfact((icazx-k)/2)))
    q = -q
  enddo
  RACAHW = s

END FUNCTION RACAHW
