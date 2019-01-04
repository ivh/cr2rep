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
double cr2res_dark_qc_ron(
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

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the detlin median non linearity 
  @param    coeffs  The detector non linearity coefficients
  @return   The computed median or -1.0 in error case
 */
/*----------------------------------------------------------------------------*/
double cr2res_qc_detlin_median(
        const cpl_imagelist     *   coeffs)
{
    double      qc_detlin_median ;

    /* Check Entries */
    if (coeffs == NULL) return -1.0 ;

    /* Initialise */
    qc_detlin_median = -1.0 ;

    /* TODO */
    
    return qc_detlin_median ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the detlin gain
  @param    coeffs  The detector non linearity coefficients
  @return   The computed gain or -1.0 in error case
 */
/*----------------------------------------------------------------------------*/
double cr2res_qc_detlin_gain(
        const cpl_imagelist     *   coeffs)
{
    double      qc_detlin_gain ;

    /* Check Entries */
    if (coeffs == NULL) return -1.0 ;

    /* Initialise */
    qc_detlin_gain = -1.0 ;

    /* TODO */
    
    return qc_detlin_gain ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the detlin min and max level
  @param    ima         input image
  @param    min_level   [out] The computed min level
  @param    max_level   [out] The computed max level
  @return   0 if ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
int cr2res_qc_detlin_min_max_level(
        const cpl_image     *   ima,
        double              *   min_level,
        double              *   max_level)
{
    double      min_level_loc, max_level_loc ;

    /* Check Entries */
    if (ima == NULL || min_level == NULL || max_level == NULL) return -1 ;

    /* Initialise */
    *min_level = -1.0 ;
    *max_level = -1.0 ;

    /* TODO */
    
    *min_level = min_level_loc ;
    *max_level = max_level_loc ;
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the total lamp intensity over a limited spectral region
  @param    ima     input image
  @return   The computed intensity or -1.0 in error case
 */
/*----------------------------------------------------------------------------*/
double cr2res_qc_flat_lamp_ints(
        const cpl_image     *   ima)
{
    double      qc_lamp_ints ;

    /* Check Entries */
    if (ima == NULL) return -1.0 ;

    /* Initialise */
    qc_lamp_ints = -1.0 ;

    /* TODO */
    
    return qc_lamp_ints ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the mean relative background level
  @param    ima     input image
  @return   The computed level or -1.0 in error case
 */
/*----------------------------------------------------------------------------*/
double cr2res_qc_flat_mean_level(
        const cpl_image     *   ima)
{
    double      qc_mean_level ;

    /* Check Entries */
    if (ima == NULL) return -1.0 ;

    /* Initialise */
    qc_mean_level = -1.0 ;

    /* TODO */
    
    return qc_mean_level ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the median SNR per pixel
  @param    ima     input image
  @return   The computed median SNR or -1.0 in error case
 */
/*----------------------------------------------------------------------------*/
double cr2res_qc_flat_med_snr(
        const cpl_image     *   ima)
{
    double      qc_med_snr ;

    /* Check Entries */
    if (ima == NULL) return -1.0 ;

    /* Initialise */
    qc_med_snr = -1.0 ;

    /* TODO */
    
    return qc_med_snr ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the mean and median flux
  @param    ima     input image
  @param    mean_flux   [out] The computed mean flux
  @param    med_flux    [out] The computed median flux
  @return   0 if ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
int cr2res_qc_flat_mean_med_flux(
        const cpl_image     *   ima,
        double              *   mean_flux,
        double              *   med_flux)
{
    double      mean_flux_loc, med_flux_loc ;

    /* Check Entries */
    if (ima == NULL || mean_flux == NULL || med_flux == NULL) return -1 ;

    /* Initialise */
    *mean_flux = -1.0 ;
    *med_flux = -1.0 ;

    /* TODO */
    
    *mean_flux = mean_flux_loc ;
    *med_flux = med_flux_loc ;
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the mean Y coord of the central order 
  @param    trace   the trace table
  @return   The computed Y coordinate or -1.0 in error case
 */
/*----------------------------------------------------------------------------*/
double cr2res_qc_flat_trace_center_y(
        const cpl_table     *   trace)
{
    double      qc_trace_center_y ;

    /* Check Entries */
    if (trace == NULL) return -1.0 ;

    /* Initialise */
    qc_trace_center_y = -1.0 ;

    /* TODO */
    
    return qc_trace_center_y ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the number of overexposed pixels in the first raw frame
  @param    ima     the first raw image
  @return   the computed number or -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_qc_flat_nb_overexposed(
        const cpl_image     *   ima)
{
    int     qc_overexposed ;

    /* Check Entries */
    if (ima == NULL) return -1 ;

    /* Initialise */
    qc_overexposed = -1 ;
    
    return qc_overexposed ;
}

/**@}*/

