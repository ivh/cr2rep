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

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

int cr2res_qc_count_badpix(cpl_image * bpm, int type);

double cr2res_qc_read_out_noise(cpl_image * im1, cpl_image * im2);
double cr2res_qc_dark_mean(cpl_image * dark);
double cr2res_qc_dark_median(cpl_image * dark);
double cr2res_qc_dark_stddev(cpl_image * dark);

int cr2res_qc_trace_count_orders(cpl_table * tracewave);
int cr2res_qc_trace_count_traces(cpl_table * tracewave);
int cr2res_qc_trace_get_ypos(cpl_table * tracewave, int order, int trace);

double cr2res_qc_wave_zeropoint(cpl_table * tracewave, int order, int trace);

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_qc  QC related functions
 *
 * TBD
 */
/*----------------------------------------------------------------------------*/

/**@{*/


/*---------------------------------------------------------------------------*/
/**
* @brief    count number of bad pixels in image
* @param    bpm    image with bad pixels
* @param    type   value to count as bad pixels
*/
/*---------------------------------------------------------------------------*/
int cr2res_qc_count_badpix(cpl_image * bpm, int type)
{
    cpl_image * tmp = cpl_image_copy(bpm);
    cpl_image_reject_value(tmp, type);
    int nbp = cpl_image_count_rejected(tmp);
    return nbp;
}

/*---------------------------------------------------------------------------*/
/**
* @brief    get the read out noise of image1 and 2
* @param    im1    first image
* @param    im2    second image
*/
/*---------------------------------------------------------------------------*/
double cr2res_qc_read_out_noise(cpl_image * im1, cpl_image * im2)
{
    cpl_image *tmp = cpl_image_copy(im1);
    cpl_image_subtract(tmp, im2);
    double rdnoise = cpl_image_get_stdev(tmp);
}

/*---------------------------------------------------------------------------*/
/**
* @brief    get the mean dark current of the image
* @param    dark   dark frame image, i.e. with no light exposure
*
*/
/*---------------------------------------------------------------------------*/
double cr2res_qc_dark_mean(cpl_image * dark)
{
    double result = cpl_image_get_mean(dark);
    return result;
}

/*---------------------------------------------------------------------------*/
/**
* @brief    get the median of the dark current in the image
* @param    dark   dark frame image, i.e. with no light exposure
*
*/
/*---------------------------------------------------------------------------*/
double cr2res_qc_dark_median(cpl_image * dark)
{
    double result = cpl_image_get_median(dark);
    return result;
}

/*---------------------------------------------------------------------------*/
/**
* @brief    get the standard deviation of the dark current
* @param    dark   dark frame image, i.e. with no light exposure
*
*/
/*---------------------------------------------------------------------------*/
double cr2res_qc_dark_stddev(cpl_image * dark)
{
    double result = cpl_image_get_stdev(dark);
    return result;
}


/*---------------------------------------------------------------------------*/
/**
* @brief    count the trace orders in the table
* @param    tracewave    table with traces as polynomials
*
*/
/*---------------------------------------------------------------------------*/
int cr2res_qc_trace_count_orders(cpl_table * tracewave)
{
    int nb_orders = 0;
    int * res = cr2res_trace_get_order_numbers(tracewave, &nb_orders);
    cpl_free(res);
    return nb_orders; 
}

/*---------------------------------------------------------------------------*/
/**
* @brief    count the number of traces in the table
* @param    tracewave    table with traces as polynomials
*
*/
/*---------------------------------------------------------------------------*/
int cr2res_qc_trace_count_traces(cpl_table * tracewave)
{
    int ntraces = cpl_table_get_nrow(tracewave);
    return ntraces;
}

/*---------------------------------------------------------------------------*/
/**
* @brief    get the central y position of a given trace and order
* @param    tracewave    table with traces as polynomials
* @param    order        order to get values for
* @param    trace        trace of that order
*/
/*---------------------------------------------------------------------------*/
int cr2res_qc_trace_get_ypos(cpl_table * tracewave, int order, int trace)
{
    cpl_vector *center = cr2res_trace_get_ycen(tracewave, order, trace, CR2RES_DETECTOR_SIZE);
    int result = cpl_vector_get(center, CR2RES_DETECTOR_SIZE/2-1);
    cpl_vector_delete(center);
    return result;
}

/*---------------------------------------------------------------------------*/
/**
* @brief    get the zeropoint (i.e. y(x=0)) for a given order and trace 
* @param    tracewave    table with traces as polynomials
* @param    order        order to get values for
* @param    trace        trace of that order
*/
/*---------------------------------------------------------------------------*/
double cr2res_qc_wave_zeropoint(cpl_table * tracewave, int order, int trace)
{
    cpl_vector *center = cr2res_trace_get_ycen(tracewave, order, trace, 1);
    int result = cpl_vector_get(center, 0);
    cpl_vector_delete(center);
    return result;
}

/**@}*/

