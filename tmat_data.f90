! Rmatrix — multilevel, multichannel R-matrix code
! Copyright (C) 1998-2004, 2026 Grigory V. Rogachev
! Copyright (C) 2026 Alexander Volya
! SPDX-License-Identifier: GPL-3.0-or-later
! This program is free software under the GNU General Public License,
! version 3 or (at your option) any later version; see LICENSE.
! Third-party components: see NOTICE.
      MODULE tmat_data
        USE channel_data, ONLY: MaxTJL, MaxCh, MaxNL
        IMPLICIT NONE

        ! T-matrix indexed by (2J, ch_in, ch_out, L_in, L_out)
        COMPLEX*16 :: TJL(0:MaxTJL, 1:MaxCh, 1:MaxCh, 0:MaxNL, 0:MaxNL)

      END MODULE tmat_data
