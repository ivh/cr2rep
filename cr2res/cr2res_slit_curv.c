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

#include "cr2res_slit_curv.h"
#include "cr2res_dfs.h"
#include "cr2res_io.h"
#include "cr2res_pfits.h"
#include "cr2res_utils.h"
#include "cr2res_trace.h"

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_slit_curv   Slit Curvature
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief 
  @param 
  @return 

 */
/*----------------------------------------------------------------------------*/
cpl_polynomial ** cr2res_slit_curv_compute_order_trace(
        cpl_table           *   trace_wave,
        int                     order,
        int                     trace_nb,
        int                     display)
{
    cpl_polynomial  **  out_polys ;
    cpl_size            ntraces, i ;
    cpl_polynomial  *   current_trace ;
    cpl_polynomial  **  order_other_traces ;

    /* Check Entries */
    if (trace_wave == NULL) return NULL ;

    /* Allocate Output array */
    out_polys = cpl_malloc(CR2RES_DETECTOR_SIZE * sizeof(cpl_polynomial*)) ;

    /* Get the number of traces */
    ntraces = cr2res_get_traces_number(trace_wave, order) ;






    /* Loop on all x positiona */
    for (i=0 ; i<CR2RES_DETECTOR_SIZE ; i++) {
        out_polys[i] = cpl_polynomial_new(1) ;














    }
    return out_polys ;
}

/**@}*/

