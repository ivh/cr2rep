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

#include <string.h>
#include <math.h>

#include <cpl.h>

#include "cr2res_qc.h"
#include "cr2res_trace.h"
#include "cr2res_dfs.h"

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_qc  QC related functions
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    The Read Out Noise computation
  @param    ima1        the first input image
  @param    ima2        the second input image
  @param    hsize
  @param    nsamples
  @param    ndit        the NDIT
  @return   the RON or -1 in error case
 */
/*----------------------------------------------------------------------------*/
double cr2res_qc_ron(
        const cpl_image     *   ima1,
        const cpl_image     *   ima2,
        int                     hsize,
        int                     nsamples,
        int                     ndit)
{
    cpl_image       *   ima ;
    double              norm, ron ;

    /* Test entries */
    if (ima1 == NULL || ima2 == NULL || ndit < 1)   return -1.0 ;

    /* Compute norm */
    norm = 0.5 * ndit ;
    norm = sqrt(norm) ;

    /* Subtraction */
    if ((ima = cpl_image_subtract_create(ima2, ima1)) == NULL) return -1.0 ;

    /* RON measurement */
    cpl_flux_get_noise_window(ima, NULL, hsize, nsamples, &ron, NULL) ;
    cpl_image_delete(ima) ;
    return norm*ron ;
}

/*---------------------------------------------------------------------------*/
/**
* @brief    get the mean dark current of the image
* @param    dark   dark frame image, i.e. with no light exposure
* @return   mean   mean dark current
*/
/*---------------------------------------------------------------------------*/
double cr2res_qc_dark_mean(cpl_image * dark)
{
    if (dark == NULL) return -1;

    double mean = cpl_image_get_mean(dark);
    return mean;
}

/*---------------------------------------------------------------------------*/
/**
* @brief    get the median of the dark current in the image
* @param    dark   dark frame image, i.e. with no light exposure
* @return   median median dark current
*/
/*---------------------------------------------------------------------------*/
double cr2res_qc_dark_median(cpl_image * dark)
{
    if (dark == NULL) return -1;

    double median = cpl_image_get_median(dark);
    return median;
}

/*---------------------------------------------------------------------------*/
/**
* @brief    get the standard deviation of the dark current
* @param    dark    dark frame image, i.e. with no light exposure
* @return   stddev  standard deviation of dark current
*/
/*---------------------------------------------------------------------------*/
double cr2res_qc_dark_stddev(cpl_image * dark)
{
    if (dark == NULL) return -1;

    double stddev = cpl_image_get_stdev(dark);
    return stddev;
}


/*---------------------------------------------------------------------------*/
/**
* @brief    count the trace orders in the table
* @param    tracewave    table with traces as polynomials
* @return   nb_orders    number of individual orders
*/
/*---------------------------------------------------------------------------*/
int cr2res_qc_trace_count_orders(cpl_table * tracewave)
{
    int nb_orders = -1;
    int * res = cr2res_trace_get_order_numbers(tracewave, &nb_orders);
    cpl_free(res);
    return nb_orders; 
}

/*---------------------------------------------------------------------------*/
/**
* @brief    count the number of traces in the table
* @param    tracewave    table with traces as polynomials
* @return   ntraces      number of traces
*/
/*---------------------------------------------------------------------------*/
int cr2res_qc_trace_count_traces(cpl_table * tracewave)
{
    if (tracewave == NULL) return -1;

    int ntraces = cpl_table_get_nrow(tracewave);
    return ntraces;
}

/*---------------------------------------------------------------------------*/
/**
* @brief    get the central y position of a given trace and order
* @param    tracewave    table with traces as polynomials
* @param    order        order to get values for
* @param    trace        trace of that order
* @return   ycen         y value of central pixel of trace and order
*/
/*---------------------------------------------------------------------------*/
int cr2res_qc_trace_get_ypos(cpl_table * tracewave, int order, int trace)
{
    cpl_vector *center = cr2res_trace_get_ycen(tracewave, order, trace, CR2RES_DETECTOR_SIZE);
    
    if (center == NULL) return -1;

    int ycen = cpl_vector_get(center, CR2RES_DETECTOR_SIZE/2-1);
    cpl_vector_delete(center);
    return ycen;
}

/*---------------------------------------------------------------------------*/
/**
* @brief    get the zeropoint (i.e. y(x=0)) for a given order and trace 
* @param    tracewave    table with traces as polynomials
* @param    order        order to get values for
* @param    trace        trace of that order
* @return   y0           y position of the center leftmost pixel of the trace and order
*/
/*---------------------------------------------------------------------------*/
double cr2res_qc_wave_zeropoint(cpl_table * tracewave, int order, int trace)
{
    cpl_vector *center = cr2res_trace_get_ycen(tracewave, order, trace, 1);
 
    if (center == NULL) return -1;

    int y0 = cpl_vector_get(center, 0);
    cpl_vector_delete(center);
    return y0;
}

/**@}*/

