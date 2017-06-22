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
#include "cr2res_wave.h"
#include "cr2res_io.h"
#include "cr2res_pfits.h"

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_wave        Wavelength Calibration
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Top level function that takes spectrum, returns solution.
  @param    spectrum        Input spectrum: arc lamp, etalon etc.
  @param    initial_guess   Starting wavelength solution
  @param    catalog         Line catalog, wavelengths, strengths
  @param    template        Template spectrum (flux, wavelengths)
  @return  Wavelength solution, i.e. polynomial that translates pixel
            values to wavelength.
 */
/*----------------------------------------------------------------------------*/
cpl_polynomial * cr2res_wave(
        cpl_vector      *   spectrum,
        cpl_polynomial  *   initial_guess,
        cpl_table       *   catalog,
        cpl_bivector    *   template)
{


    cpl_polynomial_dump(initial_guess, stdout) ;
    cpl_plot_vector("", "w lines", "", spectrum) ;

    return NULL ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief  Find solution by cross-correlating template spectrum
  @param    spectrum        Input spectrum
  @param    initial_guess   Starting wavelength solution
  @param    template        Template spectrum (flux, wavelengths)
  @return  Wavelength solution, i.e. polynomial that translates pixel
            values to wavelength.

    TODO: Summarize method
 */
/*----------------------------------------------------------------------------*/
cpl_polynomial * cr2res_wave_xcorr(
        cpl_vector      *   spectrum,
        cpl_polynomial  *   initial_guess,
        cpl_bivector    *   template)
{
    return NULL ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief   Find solution by finding lines and fitting
  @param    spectrum        Input spectrum
  @param    initial_guess   Starting wavelength solution
  @param    catalog         Line catalog, wavelengths, strengths
  @param    template        Template spectrum (flux, wavelengths)
  @return  Wavelength solution, i.e. polynomial that translates pixel
            values to wavelength.
 */
/*----------------------------------------------------------------------------*/
cpl_polynomial * cr2res_wave_line_fitting(
        cpl_vector      *   spectrum,
        cpl_polynomial  *   initial_guess,
        cpl_table       *   catalog)
{
    return NULL ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief   Find solution from etalon
  @param    spectrum        Input spectrum: etalon
  @param    initial_guess   Starting wavelength solution
  @return  Wavelength solution, i.e. polynomial that translates pixel
            values to wavelength.

    This function uses the intrinsic property of the etalon spectrum,
    that the lines are *supposed* to be equi-spaced, to refine the wavelength
    solution. The input solution needs to be good enough for the zero-point
    since the etalon spectrum carries no information about the absolute
    wavelength scale.

    The method involves these steps:
    * Identify lines (thresholding)
    * Determine line centers (ceter of gravity or gauss fit, needs testing)
    * Subtract x-coods from subsequent lines, i.e. measure d many times
        d is the distance between fringes in pixels, assumed to be constant
    * Determine mis-counts in fringes, by looking at outliers in d-distribution,
        re-count with e.g. half-distance between two fringes, if one is missing.
    * Fit the d-distibution to determine "true d"
    * Use d to calculate new x-coodinates, use input-solution to translate to
        wavelength
    * Fit measured x-coords to new lamdas to get new solution.
 */
/*----------------------------------------------------------------------------*/
cpl_polynomial * cr2res_wave_etalon(
        cpl_vector      *   spectrum,
        cpl_polynomial  *   initial_guess)
{
    cpl_vector  *   Ds;

    Ds = cr2res_wave_etalon_measure_Ds(spectrum);

    cpl_vector_delete(Ds);
    return NULL ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief Measure the distances between etalon lines
  @param
  @return
 */
/*----------------------------------------------------------------------------*/

cpl_vector * cr2res_wave_etalon_measure_Ds(
            cpl_vector * spectrum)
{
    cpl_vector  *   D_vec;
    cpl_vector  *   spec_thresh;
    cpl_vector  *   cur_peak;
    cpl_vector  *   X;
    int             numD = 0 ;
    int             smooth = 70 ;   // TODO: make free parameter?
                                    // interfringe ~30 in Y, ~70 in K
    int             thresh = 10 ;
    int             i, j, k ;
    cpl_size        nx ;
    double          spec_i;
    double      *   x0, *sigma, *area, *offset;

    nx = cpl_vector_get_size(spectrum) ;
    spec_thresh = cr2res_threshold_spec(spectrum, smooth, thresh) ;

    if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
        cpl_vector_save(spec_thresh, "debug_thresh.fits", CPL_TYPE_DOUBLE, NULL,
                CPL_IO_CREATE);
    }

    cur_peak = cpl_vector_new(256);
    X = cpl_vector_new(256);
    for (i=0;i<256;i++) cpl_vector_set(X, i, (double)i ) ;

    i = 0;
    while (i < nx){
        j = 0;
        while ( (spec_i = cpl_vector_get(spec_thresh, i)) > 0 ) {
            i++;
            j++;
        }
        cur_peak = cpl_vector_new(j);
        X = cpl_vector_new(j);
        for (k=0; k<i-j; k++){
            cpl_vector_set(cur_peak, k, cpl_vector_get(spec_thresh, k+(i-j) ));
            cpl_vector_set(X, k, (double)k );
        }
        cpl_vector_fit_gaussian(X, NULL, cur_peak, NULL, CPL_FIT_ALL,
                                x0, sigma, area, offset,
                                NULL,NULL,NULL);
        cpl_msg_debug(__func__,"Fit: %.2f, %.2f, %.2f, %.2f",x0, sigma, area, offset);
        cpl_vector_delete(cur_peak);
        cpl_vector_delete(X);
    }

    cpl_vector_delete(spec_thresh);
    D_vec = cpl_vector_new(1);
    return D_vec;
}

/*----------------------------------------------------------------------------*/
/**
  @brief
  @param
  @return
 */
/*----------------------------------------------------------------------------*/
cpl_vector * cr2res_wave_line_detection(
        cpl_vector      *   spectrum)
{
    return NULL ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief
  @param
  @return
 */
/*----------------------------------------------------------------------------*/
cpl_vector * cr2res_wave_gen_spectrum(
        cpl_table       *   catalog,
        cpl_polynomial  *   initial_guess)
{
    return NULL ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Compute the wavelength polynomial from boundaries
  @param    wmin    First pixel wavelength
  @param    wmax    Last pixel wavelength
  @return   the polynomial or NULL in error case

  The returned polynomial must be deallocated with cpl_polynomial_delete()
 */
/*----------------------------------------------------------------------------*/
cpl_array * cr2res_wave_get_estimate(
        const char  *   filename,
        int             detector,
        int             order)
{
    double                  wmin, wmax ;
    cpl_array           *   wl ;
    cpl_propertylist    *   plist ;
    int                     wished_ext_nb ;

    /* Check Entries */
    if (filename == NULL) return NULL ;
    if (order < 0 || detector <= 0 || detector > CR2RES_NB_DETECTORS)
        return NULL ;

    /* Load the propertylist */
    wished_ext_nb = cr2res_io_get_ext_idx(filename, detector) ;
    plist = cpl_propertylist_load(filename, wished_ext_nb) ;

    /* Get the values for this order */
    wmin = cr2res_pfits_get_wmin(plist, order, detector) ;
    wmax = cr2res_pfits_get_wmax(plist, order, detector) ;
    cpl_propertylist_delete(plist) ;
    if (cpl_error_get_code() != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot Get the WMIN/WMAX from the header") ;
        return NULL ;
    }

    /* Create the array */
    wl = cpl_array_new(2, CPL_TYPE_DOUBLE) ;
    cpl_array_set(wl, 0, wmin) ;
    cpl_array_set(wl, 1, wmax) ;

    return wl ;
}

/**@}*/
