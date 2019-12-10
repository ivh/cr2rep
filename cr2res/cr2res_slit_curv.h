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

#ifndef CR2RES_SLIT_CURV_H
#define CR2RES_SLIT_CURV_H

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>
#include "hdrl.h"
#include "cr2res_utils.h"

/*-----------------------------------------------------------------------------
                                       Prototypes
 -----------------------------------------------------------------------------*/

cpl_polynomial ** cr2res_slit_curv_compute_order_trace(
        cpl_table           *   trace_wave,
        int                     order,
        int                     trace_id,
        int                     max_curv_degree,
        int                     display) ;

int cr2res_slit_curv_fit_coefficients(
        cpl_polynomial  **  curvatures,
        int                 nb_polys,
        cpl_polynomial  **  slit_polya,
        cpl_polynomial  **  slit_polyb,
        cpl_polynomial  **  slit_polyc) ;

hdrl_image * cr2res_slit_curv_gen_map(
        const cpl_table *   trace_wave,
        int                 order,
        int                 trace_id,
        int                 spacing_pixels,
        int                 full_trace) ;

cpl_polynomial * cr2res_slit_curv_build_poly(
        cpl_polynomial  *   slit_poly_a,
        cpl_polynomial  *   slit_poly_b,
        cpl_polynomial  *   slit_poly_c,
        cpl_size            x) ;


int cr2res_slit_curv_from_image(
        const hdrl_image * img,
        const cpl_table * trace_wave,
        int order,
        int trace,
        int height,
        int window,
        cpl_size degree,
        int fit_second_order,
        cpl_polynomial  **   slit_poly_a,
        cpl_polynomial  **   slit_poly_b,
        cpl_polynomial  **   slit_poly_c);

#endif
