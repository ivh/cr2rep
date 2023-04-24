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

#ifndef CR2RES_IDP_H
#define CR2RES_IDP_H

/*-----------------------------------------------------------------------------
                                   Includes
-----------------------------------------------------------------------------*/

#include <cpl.h>

/*-----------------------------------------------------------------------------
                                    Define
-----------------------------------------------------------------------------*/

#define SPEC_RESOL_SLIT02       86000.0
#define SPEC_RESOL_SLIT04       43000.0


/*-----------------------------------------------------------------------------
                                   Functions prototypes
-----------------------------------------------------------------------------*/
int cr2res_wl_is_ghost(const char * setting, double wl);
cpl_bivector * cr2res_get_ghosts(const char * setting);

int cr2res_idp_save(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   rawframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        cpl_propertylist        *  ext_plist[3],
        const char              *   recipe) ;

cpl_table * cr2res_idp_create_table(
        cpl_table               **  tables,
        const char              *   recipe,
        const char              *   setting) ;

int cr2res_idp_compute_mjd(
        cpl_frameset        *   fset,
        double              *   mjd_start,
        double              *   mjd_end) ;

#endif

