! Rmatrix — multilevel, multichannel R-matrix code
! Copyright (C) 1998-2004, 2026 Grigory V. Rogachev
! Copyright (C) 2026 Alexander Volya
! SPDX-License-Identifier: GPL-3.0-or-later
! This program is free software under the GNU General Public License,
! version 3 or (at your option) any later version; see LICENSE.
! Third-party components: see NOTICE.
MODULE io_filenames
! ----------------------------------------------------------------------
! Shared storage for the three user-specified file paths.
! Populated by MAIN from command-line flags (-i/-p/-o); if left empty,
! each reader subroutine falls back to an interactive prompt.
! ----------------------------------------------------------------------
  IMPLICIT NONE
  CHARACTER(LEN=256) :: inp_file = ''   ! channel input    (-i)
  CHARACTER(LEN=256) :: par_file = ''   ! resonance params (-p)
  CHARACTER(LEN=256) :: out_file = ''   ! output           (-o)
END MODULE io_filenames
