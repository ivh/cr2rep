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
        double              smoothfactor,
        int                 opening,
        int                 degree,
        int                 min_cluster,
        int                 split_single_trace_orders) ;

cpl_mask * cr2res_trace_clean(
        cpl_mask    *   mask,
        int             opening,
        int             min_cluster) ;

cpl_table * cr2res_trace_fit(
        cpl_image   *   labels,
        int             degree, 
        int             split_single_trace_orders) ;

cpl_image * cr2res_trace_gen_image(
        cpl_table   *   trace,
        int             nx,
        int             ny) ;

int * cr2res_trace_get_order_numbers(
        cpl_table   *   trace, 
        int         *   nb_orders) ;

cpl_polynomial * cr2res_get_trace_wave_poly(
        const cpl_table     *   trace_wave,
        const char          *   poly_column,
        int                     order,
        int                     trace_nb) ;

cpl_size cr2res_get_trace_table_index(
        const cpl_table     *   trace_wave,
        int                     order,
        int                     trace_nb) ;

cpl_size cr2res_get_nb_traces(
        const cpl_table     *   trace_wave,
        int                     order) ;

cpl_size cr2res_get_nb_traces_with_wavelength(
        const cpl_table     *   trace_wave,
        int                     order) ;

cpl_vector * cr2res_trace_get_ycen(
            cpl_table   *   trace,
            cpl_size        order_nb,
            cpl_size        trace_nb,
            int             size);

int cr2res_trace_get_height(
            cpl_table   *   trace,
            cpl_size        order_nb,
            cpl_size        trace_nb);

cpl_vector * cr2res_trace_compute_middle(
        cpl_polynomial  *   trace1,
        cpl_polynomial  *   trace2,
        int                 vector_size) ;

int cr2res_trace_compute_height(
        cpl_polynomial  *   trace1,
        cpl_polynomial  *   trace2,
        int                 vector_size) ;

double cr2res_trace_get_trace_ypos(
        cpl_table   *   traces,
        int             idx) ;

int cr2res_trace_add_ord_tra_wav_curv_columns(
        cpl_table           *   traces,
        const char          *   file_for_wl,
        int                     det_nr) ;

#endif
