! Rmatrix — multilevel, multichannel R-matrix code
! Copyright (C) 1998-2004, 2026 Grigory V. Rogachev
! Copyright (C) 2026 Alexander Volya
! SPDX-License-Identifier: GPL-3.0-or-later
! This program is free software under the GNU General Public License,
! version 3 or (at your option) any later version; see LICENSE.
! Third-party components: see NOTICE.
DOUBLE PRECISION FUNCTION intgsl(kappa, mu, z)
! ----------------------------------------------------------------------
! Whittaker W_{kappa,mu}(z), evaluated through GSL's confluent
! hypergeometric U function:
!   W_{kappa,mu}(z) = exp(-z/2) z^(mu+1/2)
!                     U(mu-kappa+1/2, 1+2*mu, z)
! This interface is for the real, positive-z closed-channel case.
! GSL's default error handler (which aborts the program) is switched
! off; a non-zero status stops with a message naming the arguments.
! ----------------------------------------------------------------------
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_DOUBLE, C_INT, C_PTR
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN) :: kappa, mu, z
  DOUBLE PRECISION :: a, b

  TYPE, BIND(C) :: gsl_sf_result
    REAL(C_DOUBLE) :: val
    REAL(C_DOUBLE) :: err
  END TYPE gsl_sf_result

  INTERFACE
    FUNCTION gsl_sf_hyperg_U_e(a_arg, b_arg, x_arg, res) &
        BIND(C, NAME='gsl_sf_hyperg_U_e') RESULT(status)
      IMPORT :: C_DOUBLE, C_INT, gsl_sf_result
      REAL(C_DOUBLE), VALUE :: a_arg, b_arg, x_arg
      TYPE(gsl_sf_result) :: res
      INTEGER(C_INT) :: status
    END FUNCTION gsl_sf_hyperg_U_e
    FUNCTION gsl_set_error_handler_off() &
        BIND(C, NAME='gsl_set_error_handler_off') RESULT(old)
      IMPORT :: C_PTR
      TYPE(C_PTR) :: old
    END FUNCTION gsl_set_error_handler_off
  END INTERFACE

  TYPE(gsl_sf_result) :: res
  TYPE(C_PTR)         :: old
  INTEGER(C_INT)      :: status
  LOGICAL, SAVE       :: first = .TRUE.

  if (first) then
    old   = gsl_set_error_handler_off()
    first = .FALSE.
  endif

  if (z.le.0.d0) then
    PRINT *, 'INTGSL error: Whittaker argument must be positive, z=', z
    STOP 1
  endif

  a = mu - kappa + 0.5d0
  b = 1.d0 + 2.d0*mu
  status = gsl_sf_hyperg_U_e(a, b, z, res)
  if (status /= 0) then
    WRITE(*,'(a,i0,a,3es14.6)') 'INTGSL error: gsl_sf_hyperg_U status ', status, &
          ' for kappa, mu, z =', kappa, mu, z
    STOP 1
  endif
  intgsl = exp(-0.5d0*z) * z**(mu+0.5d0) * res%val

END FUNCTION intgsl
