* Rmatrix — multilevel, multichannel R-matrix code
* Copyright (C) 1998-2004, 2026 Grigory V. Rogachev
* Copyright (C) 2026 Alexander Volya
* SPDX-License-Identifier: GPL-3.0-or-later
* This program is free software under the GNU General Public License,
* version 3 or (at your option) any later version; see LICENSE.
* Third-party components: see NOTICE.
C Stub for CERNLIB MTLPRT error message printer.
C Prints the routine name, error code, and message to stderr.
      SUBROUTINE MTLPRT(NAME,CODE,TEXT)
      CHARACTER*(*) NAME,CODE,TEXT
      WRITE(0,*) '***** CERNLIB ',NAME,' [',CODE,'] ',TEXT
      RETURN
      END
