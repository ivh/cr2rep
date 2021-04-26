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

#ifndef CR2RES_POL_H
#define CR2RES_POL_H

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>
#include <hdrl.h>

/*-----------------------------------------------------------------------------
                                   Define
 -----------------------------------------------------------------------------*/

#define     CR2RES_POLARIMETRY_GROUP_SIZE   4

/*-----------------------------------------------------------------------------
                                       Prototypes
 -----------------------------------------------------------------------------*/

cpl_bivector * cr2res_pol_demod_stokes(
        cpl_vector  **  intens,
        cpl_vector  **  wl,
        cpl_vector  **  errors,
        int             n) ;
cpl_bivector * cr2res_pol_demod_null(
        cpl_vector  **  intens,
        cpl_vector  **  wl,
        cpl_vector  **  errors,
        int             n) ;
cpl_bivector * cr2res_pol_demod_intens(
        cpl_vector  **  intens,
        cpl_vector  **  wl,
        cpl_vector  **  errors,
        int             n) ;
cpl_table * cr2res_pol_POL_SPEC_create(
        int             *   orders,
        cpl_vector      **  wl,
        cpl_bivector    **  stokes,
        cpl_bivector    **  null,
        cpl_bivector    **  intens,
        int                 norders) ;
int * cr2res_pol_sort_frames(
        const cpl_frame *   frame1,
        const cpl_frame *   frame2,
        const cpl_frame *   frame3,
        const cpl_frame *   frame4) ;

cpl_table * cr2res_pol_spec_pol_merge(
        const cpl_table **  pol_spec_list,
        int                 pol_spec_nb) ;

cpl_error_code cr2res_pol_subtract_background(
          cpl_frameset * rawframes_a,
          cpl_frameset * rawframes_b);


#endif
