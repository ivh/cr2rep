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

#ifndef CR2RES_EXTRACT_H
#define CR2RES_EXTRACT_H

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>
#include "hdrl.h"

#include "cr2res_utils.h"

/*-----------------------------------------------------------------------------
                                       Prototypes
 -----------------------------------------------------------------------------*/

int cr2res_extract_slitdec_vert(
        cpl_image   *   img_hdrl,
        cpl_table   *   trace_tab,
        int             order,
        int             trace_id,
        int             height,
        int             swath,
        int             oversample,
        double          smooth_slit,
        cpl_vector  **  slit_func,
        cpl_vector  **  spec,
        hdrl_image  **  model) ;

int cr2res_extract_slitdec_curved(
        hdrl_image  *   img_hdrl,
        cpl_table   *   trace_tab,
        int             order,
        int             trace_id,
        int             height,
        int             swath,
        int             oversample,
        double          smooth_slit,
        cpl_vector  **  slit_func,
        cpl_vector  **  spec,
        hdrl_image  **  model);

int cr2res_extract_sum_vert(
        cpl_image   *   img_in,
        cpl_table   *   trace_tab,
        int             order,
        int             trace_id,
        int             height,
        cpl_vector  **  slit_func,
        cpl_vector  **  spec,
        hdrl_image  **  model) ;

cpl_table * cr2res_extract_SLITFUNC_create(
        cpl_vector      **  slit_func,
        cpl_table       *   trace_table) ;

cpl_table * cr2res_extract_EXTRACT1D_create(
        cpl_vector      **  spectrum,
        cpl_table       *   trace_table) ;

cpl_vector * cr2res_extract_EXTRACT1D_get_spectrum(
        cpl_table   *   tab,
        int             order,
        int             trace_nb) ;

#endif
