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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <math.h>
#include <string.h>

#include <cpl.h>

#include "cr2res_dfs.h"
#include "cr2res_io.h"
#include "cr2res_pfits.h"
#include "cr2res_utils.h"

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_tilt    Slit Tilt
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Compute the slit tilt polynomial for all x along the order
  @param    trace_wave  The trace wave table
  @param    order       The order that needѕ to be computed
  @param    display     For Display
  @return   A list of CR2RES_DETECTOR_SIZE polynomials or NULL in error case




 */
/*----------------------------------------------------------------------------*/
cpl_polynomial ** cr2res_tilt(
        cpl_table           *   trace_wave,
        int                     order,
        int                     display)
{
    cpl_polynomial  **  out_polys ;
    cpl_polynomial  *   wave1 ;
    cpl_polynomial  *   wave2 ;
    int                 trace_nb, trace_nb1, trace_nb2, nrows, cur_order ;
    int                 i ;

    /* Check Entries */
    if (trace_wave == NULL) return NULL ;

    /* Initialise */
    nrows = cpl_table_get_nrow(trace_wave) ;
    trace_nb1 = trace_nb2 = -1 ;

    /* Check that there are exactly 2 traces for this order */
    for (i=0 ; i<nrows ; i++) {
        cur_order = cpl_table_get(trace_wave, CR2RES_COL_ORDER, i, NULL) ;
        trace_nb = cpl_table_get(trace_wave, CR2RES_COL_TRACENB, i, NULL) ;
        if (order == cur_order) {
            /* Fill trace #1 */
            if (trace_nb1 < 0) {
                trace_nb1 = trace_nb ;
                continue ;
            }
            /* Fill trace #2 */
            if (trace_nb2 < 0) {
                trace_nb2 = trace_nb ;
                continue ;
            }
            /* We encounter here the 3rd trace - Exit for the moment */
            cpl_msg_error(__func__, "Third trace found for this Order"); 
            return NULL ;
        }
    }

    /* We need 2 traces to compute the tilt */
    if (trace_nb1 < 0 || trace_nb2 < 0) {
        cpl_msg_error(__func__, "Two traces are needed to compute the tilt"); 
        return NULL ;
    }

    /* Check that the wavelength polynomials are available */
    wave1 = cr2res_get_trace_wave_poly(trace_wave, CR2RES_COL_WAVELENGTH, order,
            trace_nb1) ;
    wave2 = cr2res_get_trace_wave_poly(trace_wave, CR2RES_COL_WAVELENGTH, order,
            trace_nb2) ;
    if (wave1 == NULL || wave2 == NULL) {
        cpl_msg_error(__func__, "The wavelength solution is missing"); 
        if (wave1 != NULL) cpl_polynomial_delete(wave1) ;
        if (wave2 != NULL) cpl_polynomial_delete(wave2) ;
        return NULL ;
    }

    /* Create Output polynomials array */
    out_polys = cpl_calloc(CR2RES_DETECTOR_SIZE, sizeof(cpl_polynomial *)) ;

    out_polys[10] = cpl_polynomial_duplicate(wave1) ; 





    cpl_polynomial_delete(wave1) ;
    cpl_polynomial_delete(wave2) ;

    return out_polys ;
}

/**@}*/
