! Rmatrix — multilevel, multichannel R-matrix code
! Copyright (C) 1998-2004, 2026 Grigory V. Rogachev
! Copyright (C) 2026 Alexander Volya
! SPDX-License-Identifier: GPL-3.0-or-later
! This program is free software under the GNU General Public License,
! version 3 or (at your option) any later version; see LICENSE.
! Third-party components: see NOTICE.
      MODULE prep_data
        USE channel_data, ONLY: MaxCh
        IMPLICIT NONE

        ! Kinematic arrays computed once per energy in DOCALC,
        ! shared by PreTJL, RDef, NResPh, CPC, rmblatt
        DOUBLE PRECISION :: RMass(MaxCh)
        COMPLEX*16       :: wavnum(MaxCh)
        COMPLEX*16       :: eta(MaxCh)

      END MODULE prep_data
