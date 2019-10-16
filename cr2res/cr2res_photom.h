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

#ifndef CR2RES_PHOTOM_H
#define CR2RES_PHOTOM_H

/*-----------------------------------------------------------------------------
   								Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>

/*-----------------------------------------------------------------------------
   								Define
 -----------------------------------------------------------------------------*/

#define CR2RES_PHOTOM_STARSDIST_ARCSEC  30.0

/* TODO : Gain values have to be set  - these are the old CRIRES ones */
#define CRIRES_GAIN_CHIP1       		7.737 * 0.9
#define CRIRES_GAIN_CHIP2       		7.664 * 0.9
#define CRIRES_GAIN_CHIP3       		7.689 * 0.9

/*-----------------------------------------------------------------------------
   							        Prototypes
 -----------------------------------------------------------------------------*/

int cr2res_photom_engine(   
        const cpl_table     *   extr,
        const char          *   std_star_file,
        double                  ra,
        double                  dec,
        double                  gain,
        double                  exptime,
        int                     display_order,
        int                     display_trace,
        cpl_table           **  throughput) ;

cpl_bivector * cr2res_photom_conv_get_star(
        const cpl_table     *   tab,
        double                  ra,
        double                  dec) ;


#endif
