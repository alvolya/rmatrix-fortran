! Rmatrix — multilevel, multichannel R-matrix code
! Copyright (C) 1998-2004, 2026 Grigory V. Rogachev
! Copyright (C) 2026 Alexander Volya
! SPDX-License-Identifier: GPL-3.0-or-later
! This program is free software under the GNU General Public License,
! version 3 or (at your option) any later version; see LICENSE.
! Third-party components: see NOTICE.
! -------------------------------------------------------------------
! Multilevel-Multichannel R-matrix calculation
! Any number of channels and levels are accepted
! Created by G.V. Rogachev (2004), revised by A. Volya (2026)
! -------------------------------------------------------------------
PROGRAM MAIN
  USE io_filenames
  IMPLICIT NONE

  INTEGER            :: nargs, iarg
  CHARACTER(LEN=256) :: arg
  LOGICAL            :: show_help

  show_help = .FALSE.
  nargs     = COMMAND_ARGUMENT_COUNT()
  iarg      = 1

  DO WHILE (iarg <= nargs)
    CALL GET_COMMAND_ARGUMENT(iarg, arg)
    SELECT CASE (TRIM(arg))
      CASE ('-h', '--help')
        show_help = .TRUE.
        EXIT
      CASE ('-i')
        iarg = iarg + 1
        IF (iarg > nargs) THEN
          WRITE(*,'(a)') 'ERROR: -i requires a filename argument'
          STOP 1
        END IF
        CALL GET_COMMAND_ARGUMENT(iarg, inp_file)
      CASE ('-p')
        iarg = iarg + 1
        IF (iarg > nargs) THEN
          WRITE(*,'(a)') 'ERROR: -p requires a filename argument'
          STOP 1
        END IF
        CALL GET_COMMAND_ARGUMENT(iarg, par_file)
      CASE ('-o')
        iarg = iarg + 1
        IF (iarg > nargs) THEN
          WRITE(*,'(a)') 'ERROR: -o requires a filename argument'
          STOP 1
        END IF
        CALL GET_COMMAND_ARGUMENT(iarg, out_file)
      CASE DEFAULT
        WRITE(*,'(a,a)') 'ERROR: unknown option: ', TRIM(arg)
        show_help = .TRUE.
        EXIT
    END SELECT
    iarg = iarg + 1
  END DO

  IF (show_help) THEN
    WRITE(*,'(a)') 'Usage: Rmatrix [-i <input.inp>] [-p <params.par>] [-o <output.out>]'
    WRITE(*,'(a)') ''
    WRITE(*,'(a)') 'Options:'
    WRITE(*,'(a)') '  -i <file>   Channel input file (.inp)'
    WRITE(*,'(a)') '  -p <file>   Resonance parameter file (.par)'
    WRITE(*,'(a)') '  -o <file>   Output file'
    WRITE(*,'(a)') '  -h, --help  Show this help message'
    WRITE(*,'(a)') ''
    WRITE(*,'(a)') 'If options are omitted, the program prompts interactively.'
    STOP
  END IF

  CALL PREAD
  CALL REREAD
  CALL DOCALC

END PROGRAM MAIN
