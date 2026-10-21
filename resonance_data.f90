! Rmatrix — multilevel, multichannel R-matrix code
! Copyright (C) 1998-2004, 2026 Grigory V. Rogachev
! Copyright (C) 2026 Alexander Volya
! SPDX-License-Identifier: GPL-3.0-or-later
! This program is free software under the GNU General Public License,
! version 3 or (at your option) any later version; see LICENSE.
! Third-party components: see NOTICE.
      MODULE resonance_data
        USE channel_data, ONLY: MaxRes, MaxCh
        IMPLICIT NONE

        DOUBLE PRECISION :: eres(MaxRes, 3)
        DOUBLE PRECISION :: gres(MaxRes, MaxCh, 4)
        INTEGER :: lastre
        INTEGER :: rdim(MaxRes)
        ! row passes the parity and triangle selection rules (set by REREAD)
        LOGICAL :: rowok(MaxRes, MaxCh)

      END MODULE resonance_data
