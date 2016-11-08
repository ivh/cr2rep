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

#include <cpl.h>
#include "cr2res_trace.h"
#include "cr2res_utils.h"

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------
                                Trace Functions
 -----------------------------------------------------------------------------*/
 /**@{*/

 /*----------------------------------------------------------------------------*/
 /**
   @brief Main function for running all parts of the trace algorithm on an img.
   @param ima input image
   @param decker slit layout
   @param npolys [out] the number of trace polynomials determined
   @return set of trace polynomials that describe the orders

  */
 /*----------------------------------------------------------------------------*/
cpl_polynomial ** cr2res_trace(cpl_image * ima, cr2res_decker decker, int * npolys){
    return NULL;
}


 /*----------------------------------------------------------------------------*/
 /**
   @brief Determine which pixels belong to an order and which not
   @param ima input image
   @param min_cluster_size smaller clusters of connected pixels with signal are disregarded
   @return binary image with 1s for pixels in orders

  */
 /*----------------------------------------------------------------------------*/
cpl_image * cr2res_trace_detect(cpl_image * ima, int min_cluster_size){
    return NULL;
}



 /*----------------------------------------------------------------------------*/
 /**
   @brief Find out which pixels belong to the same order, label orders
   @param bin_ima input binary image
   @return image with labelled pixels.

  */
 /*----------------------------------------------------------------------------*/
cpl_image * cr2res_trace_labelize(cpl_image * bin_ima){
    return NULL;
}



 /*----------------------------------------------------------------------------*/
 /**
   @brief Fit polynomials to pixel coordinates in each order
   @param ima input image with labels
   @param degree polynomial degree
   @return fit result polynomials

  */
 /*----------------------------------------------------------------------------*/
cpl_polynomial ** cr2res_trace_fit(cpl_image * ima, int degree){
    return NULL;
}



 /*----------------------------------------------------------------------------*/
 /**
   @brief compare two traces
   @param trace1 first trace
   @param trace2 second trace
   @return sqrt(sum(distances along x-axis ^2))

  */
 /*----------------------------------------------------------------------------*/
cpl_vector * cr2res_trace_compare(cpl_table * trace1, cpl_table * trace2){
    return NULL;
}



 /*----------------------------------------------------------------------------*/
 /**
   @brief Combine two decker traces into an open trace
   @param td_13 1-3-decker trace
   @param td_24 2-4-decker trace
   @return open trace

  */
 /*----------------------------------------------------------------------------*/
cpl_table * cr2res_trace_combine(cpl_table * td_13, cpl_table * td_24){
    return NULL;
}



 /*----------------------------------------------------------------------------*/
 /**
   @brief Make an image out of the trace solution
   @param trace_open (or NULL)
   @param trace_decker_1_3 (or NULL)
   @param trace_decker_2_4 (or NULL)
   @return image for trace visualization

  */
 /*----------------------------------------------------------------------------*/
cpl_image * cr2res_trace_gen_image(cpl_table * trace_open,
                                   cpl_table * trace_decker_1_3,
                                   cpl_table * trace_decker_2_4
                               ){
    return NULL;
}


/**@}*/
