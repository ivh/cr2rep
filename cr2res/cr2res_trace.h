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

#ifndef CR2RES_TRACE_H
#define CR2RES_TRACE_H

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>
#include "cr2res_utils.h"

/*-----------------------------------------------------------------------------
                                       Prototypes
 -----------------------------------------------------------------------------*/

cpl_table * cr2res_trace(
        cpl_image       *   ima,
        int                 smooth_x,
        int                 smooth_y,
        double              threshold,
        int                 opening,
        int                 degree,
        int                 min_cluster) ;

cpl_mask * cr2res_trace_clean(
        cpl_mask    *   mask,
        int             opening,
        int             min_cluster) ;

cpl_image * cr2res_trace_gen_image(
        cpl_table   *   trace,
        int             nx,
        int             ny) ;

int * cr2res_trace_get_order_numbers(
        const cpl_table   *   trace, 
        int         *   nb_orders) ;

cpl_table * cr2res_trace_merge(
        const cpl_table     *   trace_wave1,
        const cpl_table     *   trace_wave2) ;

cpl_polynomial * cr2res_get_trace_wave_poly(
        const cpl_table     *   trace_wave,
        const char          *   poly_column,
        int                     order,
        int                     trace_nb) ;

cpl_vector * cr2res_trace_get_wl(
        const cpl_table *   trace_wave,
        int                 order,
        int                 trace_nb,
        int                 size) ;

cpl_size cr2res_get_trace_table_index(
        const cpl_table     *   trace_wave,
        int                     order,
        int                     trace_nb) ;

cpl_size cr2res_get_nb_traces(
        const cpl_table     *   trace_wave,
        int                     order) ;

int * cr2res_get_trace_numbers(
        const cpl_table     *   trace_wave,
        int                     order,
        int *                   nb_traces) ;

cpl_size cr2res_get_nb_traces_with_wavelength(
        const cpl_table     *   trace_wave,
        int                     order) ;

cpl_vector * cr2res_trace_get_ycen(
        const cpl_table *   trace,
        cpl_size            order_nb,
        cpl_size            trace_nb,
        int                 size) ;

int cr2res_trace_get_height(
        const cpl_table *   trace,
        cpl_size            order_nb,
        cpl_size            trace_nb) ;

cpl_vector * cr2res_trace_compute_middle(
        cpl_polynomial  *   trace1,
        cpl_polynomial  *   trace2,
        int                 vector_size) ;

int cr2res_trace_compute_height(
        cpl_polynomial  *   trace1,
        cpl_polynomial  *   trace2,
        int                 vector_size) ;

double cr2res_trace_get_trace_ypos(
        const cpl_table *   traces,
        int                 idx) ;

int cr2res_trace_add_extra_columns(
        cpl_table           *   traces,
        const char          *   infile,
        int                     det_nr) ;

cpl_table * cr2res_trace_new_slit_fraction(
        const cpl_table     *   traces,
        const cpl_array     *   new_slit_fraction) ;

cpl_table * cr2res_trace_adjust(
        const cpl_table     *   trace_wave,
        const cpl_frameset  *   flat_raw,
        int                     det_nr) ;

cr2res_decker cr2res_trace_slit_fraction_info(
        const cpl_array *   slit_frac,
        int             *   up_or_down) ;

cpl_array * cr2res_trace_slit_fraction_create(
        cr2res_decker   decker_position,
        int             up_or_down) ;

cpl_table * cr2res_trace_split(
        cpl_table   *   trace_wave,
        int             order,
        int             nb_subtraces) ;

#endif
