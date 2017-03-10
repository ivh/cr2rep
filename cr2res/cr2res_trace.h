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

cpl_table * cr2res_trace_cpl(
        cpl_image       *   ima,
        cr2res_decker       decker,
        double              smoothfactor,
        int                 opening,
        int                 degree,
        int                 min_cluster) ;

cpl_table * cr2res_trace_nocpl(
        cpl_image       *   ima,
        cr2res_decker       decker,
        double              smoothfactor,
        int                 opening,
        int                 degree,
        int                 min_cluster) ;

cpl_mask * cr2res_trace_detect(
        cpl_image   *   ima,
        double          smoothfactor,
        int             opening,
        int             min_cluster) ;

cpl_image * cr2res_trace_labelize(cpl_mask * mask) ;

cpl_table * cr2res_trace_fit(
        cpl_image   *   labels,
        int             degree) ;

cpl_vector * cr2res_trace_compare(
        cpl_table   *   trace1,
        cpl_table   *   trace2) ;

cpl_table * cr2res_trace_combine(
        cpl_table   *   td_13,
        cpl_table   *   td_24) ;

cpl_image * cr2res_trace_gen_image(
        cpl_table   *   trace,
        int             nx,
        int             ny) ;

int * cr2res_trace_get_order_numbers(
        cpl_table   *   trace, 
        int         *   nb_orders) ;

cpl_polynomial ** cr2res_trace_open_get_polynomials(
        cpl_table   *   trace,
        cpl_size        order_nb) ;

cpl_vector * cr2res_trace_compute_middle(
        cpl_polynomial  *   trace1,
        cpl_polynomial  *   trace2,
        int                 vector_size) ;

int cr2res_trace_compute_height(
        cpl_polynomial  *   trace1,
        cpl_polynomial  *   trace2,
        int                 vector_size) ;

#endif
