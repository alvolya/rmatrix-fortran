* Rmatrix — multilevel, multichannel R-matrix code
* Copyright (C) 1998-2004, 2026 Grigory V. Rogachev
* Copyright (C) 2026 Alexander Volya
* SPDX-License-Identifier: GPL-3.0-or-later
* This program is free software under the GNU General Public License,
* version 3 or (at your option) any later version; see LICENSE.
* Third-party components: see NOTICE.
*
* c309prec.inc - machine precision parameters for WCLBES (REAL*8 / double precision)
*
      PARAMETER(ACC8   = 0.5D0**52)
      PARAMETER(ACC16  = 0.5D0**105)
      PARAMETER(ACCU   = 1.0D-14)
      PARAMETER(FPMAX  = 1.7976931348623157D+308 * ACC8)
      PARAMETER(FPMIN  = 2.2250738585072014D-308 / ACC8)
