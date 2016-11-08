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

cpl_polynomial ** cr2res_trace(cpl_image * ima, cr2res_decker decker, int * npolys);

cpl_image * cr2res_trace_detect(cpl_image * ima, int min_cluster_size);

cpl_image * cr2res_trace_labelize(cpl_image * bin_ima);

cpl_polynomial ** cr2res_trace_fit(cpl_image * ima, int degree);

cpl_vector * cr2res_trace_compare(cpl_table * trace1, cpl_table * trace2);

cpl_table * cr2res_trace_combine(cpl_table * td_13, cpl_table * td_24);

cpl_image * cr2res_trace_gen_image(cpl_table * trace_open,
                                   cpl_table * trace_decker_1_3,
                                   cpl_table * trace_decker_2_4
                               );

#endif
