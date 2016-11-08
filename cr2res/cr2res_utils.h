/*
 * This file is part of the CR2RES Pipeline
 * Copyright (C) 2002,2003 European Southern Observatory
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02111-1307  USA
 */

#ifndef CR2RES_UTILS_H
#define CR2RES_UTILS_H

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>

/*-----------------------------------------------------------------------------
                                   Define
 -----------------------------------------------------------------------------*/

#define CR2RES_NB_DETECTORS     4

/*-----------------------------------------------------------------------------
                                       Prototypes
 -----------------------------------------------------------------------------*/

const char * cr2res_get_license(void) ;

cpl_frameset * cr2res_extract_frameset(const cpl_frameset *, const char *) ;

cpl_error_code cr2res_detector_shotnoise_model(
        const cpl_image *   ima_data,
        const double        gain,
        const double        ron,
        cpl_image       **  ima_errs) ;

cpl_table * cr2res_orders_fit(
        cpl_table           *   table);

cpl_mask * cr2res_signal_detect(
        const cpl_image     *   image,
        int                     ordersep,
        double                  smoothfactor,
        double                  thresh);

#endif