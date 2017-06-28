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
#include "cr2res_utils.h"

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
        cpl_vector          *   spectrum,
        cpl_polynomial      *   initial_guess,
        cr2res_wavecal_type     wavecal_type,
        int                     line_fitting,
        const char          *   static_file)
{
    cpl_polynomial  *   solution ;

    /* Initialise */
    solution = NULL ;

    cpl_plot_vector("", "w lines", "", spectrum) ;
    cpl_polynomial_dump(initial_guess, stdout) ;

    /* Switch on the possible methods */
    if (wavecal_type == CR2RES_LAMP) {
        if (line_fitting) {
            solution = cr2res_wave_line_fitting(spectrum, initial_guess,
                    NULL) ;
        } else {
            solution = cr2res_wave_xcorr(spectrum, initial_guess, NULL) ;
        }

    } else if (wavecal_type == CR2RES_GAS) {
        solution = cr2res_wave_xcorr(spectrum, initial_guess, NULL) ;
    } else if (wavecal_type == CR2RES_ETALON) {
        solution = cr2res_wave_etalon(spectrum, initial_guess) ;
    }
    if (solution != NULL) cpl_polynomial_dump(solution, stdout) ;
    return solution ;
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
    cpl_array   *   peaks;
    cpl_vector  *   Ds;
    cpl_vector  *   spec_thresh;
    cpl_vector  *   cur_peak;
    cpl_vector  *   X_all, *X_peak;
    int             i, j, k ;
    int             numD = 0 ;
    int             smooth = 70 ;   // TODO: make free parameter?
                                    // interfringe ~30 in Y, ~70 in K
    int             thresh = 10 ;   // TODO: derive from read-out noise
    int             max_num_peaks = 256 ;
    int             max_len_peak = 256 ;
    int             min_len_peak = 5 ; //TODO: tweak or make parameter?;
    cpl_size        nx ;
    double          spec_i;
    double          prev_peak;
    double      *   x0, *sigma, *area, *offset;

    nx = cpl_vector_get_size(spectrum) ;
    spec_thresh = cr2res_threshold_spec(spectrum, smooth, thresh) ;

    if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
        cpl_vector_save(spec_thresh, "debug_thresh.fits", CPL_TYPE_DOUBLE,
                NULL, CPL_IO_CREATE);
    }

    /*Output array, values are invalid until set.*/
    peaks = cpl_array_new(max_num_peaks, CPL_TYPE_DOUBLE);

    /* X-axis to cut out from for each peak */
    X_all = cpl_vector_new(max_len_peak);
    for (i=0; i<max_len_peak; i++) cpl_vector_set(X_all, i, (double)i) ;

    i = 0;
    while (i < nx){
        j = 0;
        while ( (spec_i = cpl_vector_get(spec_thresh, i)) > 0 ) {
            i++;
            j++;
        }
        if (j < min_len_peak) continue;
        cur_peak = cpl_vector_extract(spec_thresh, i-j, i, 1) ;
        X_peak = cpl_vector_extract(X_all, 0, j, 1) ;
        cur_peak = cpl_vector_new(j);
        cpl_vector_fit_gaussian(X_peak, NULL, cur_peak, NULL, CPL_FIT_ALL,
                                x0, sigma, area, offset,
                                NULL,NULL,NULL);
        cpl_msg_debug(__func__,"Fit: %.2f, %.2f, %.2f, %.2f",
                                    x0, sigma, area, offset);

        if ((k = cpl_array_count_invalid(peaks)) <1)
            cpl_msg_error(__func__,"Output array overflow!");
        cpl_array_set_double(peaks, max_num_peaks - k, *x0);

        cpl_vector_delete(cur_peak);
        cpl_vector_delete(X_peak);
    }

    /* Copy into output array, treating first element outside the loop */
    k = max_num_peaks - cpl_array_count_invalid(peaks) ;
    Ds = cpl_vector_new(k) ;
    prev_peak = cpl_array_get(peaks, 0, NULL) ;
    cpl_vector_set(Ds, 0, prev_peak) ;
    for (i=1; i<k; i++)
        cpl_vector_set(Ds, i,
            cpl_array_get(peaks, 0, NULL) - prev_peak ) ;

    cpl_vector_delete(spec_thresh) ;
    cpl_vector_delete(X_all) ;
    cpl_array_delete(peaks) ;
    return Ds;
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
  @return   the array with two polynomial coeffs, or NULL in error case

  The returned array must be deallocated with cpl_array_delete()
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
    wmin = cr2res_pfits_get_wmin(plist, order) ;
    wmax = cr2res_pfits_get_wmax(plist, order) ;
    cpl_propertylist_delete(plist) ;
    if (cpl_error_get_code() != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot get WMIN/WMAX from header, order "
                            "%d, detector %d", order, detector) ;
        return NULL ;
    }

    /* Create the array */
    wl = cpl_array_new(2, CPL_TYPE_DOUBLE) ;
    cpl_array_set(wl, 0, wmin) ;
    cpl_array_set(wl, 1, (wmax-wmin) / CR2RES_DETECTOR_SIZE ) ; // linear slope
    return wl ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Compute the wavelength map from the trace_wave table
  @param    trace_wave      The trace wave table
  @return   the wave_map image or NULL in error case

  The returned image must be deallocated with hdrl_image_delete()
 */
/*----------------------------------------------------------------------------*/
hdrl_image * cr2res_wave_gen_wave_map(
        const cpl_table *   trace_wave)
{

    /* Check Entries */
    if (trace_wave == NULL) return NULL ;

    /* TODO */
    return NULL ;
}
/**@}*/
