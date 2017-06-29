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
  @param    display         Flag to display results
  @return  Wavelength solution, i.e. polynomial that translates pixel
            values to wavelength.
 */
/*----------------------------------------------------------------------------*/
cpl_polynomial * cr2res_wave(
        cpl_vector          *   spectrum,
        cpl_polynomial      *   initial_guess,
        cr2res_wavecal_type     wavecal_type,
        int                     line_fitting,
        const char          *   static_file,
        int                     display)
{
    cpl_polynomial  *   solution ;
    cpl_bivector    *   ref_spectrum ;

    /* Initialise */
    solution = NULL ;

    /* Switch on the possible methods */
    if (wavecal_type == CR2RES_LAMP) {
        if (line_fitting) {
            solution = cr2res_wave_line_fitting(spectrum, initial_guess, NULL) ;
        } else {
            /* Create the lines spectrum from the lines list */
            ref_spectrum = cr2res_wave_gen_lines_spectrum(static_file,
                    initial_guess) ;
            solution = cr2res_wave_xcorr(spectrum, initial_guess,
                    ref_spectrum, display) ;
            cpl_bivector_delete(ref_spectrum) ;
        }
    } else if (wavecal_type == CR2RES_GAS) {
        solution = cr2res_wave_xcorr(spectrum, initial_guess, NULL, display) ;
    } else if (wavecal_type == CR2RES_ETALON) {
        solution = cr2res_wave_etalon(spectrum, initial_guess) ;
    }
    return solution ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief  Find solution by cross-correlating template spectrum
  @param    spectrum        Input spectrum
  @param    initial_guess   Starting wavelength solution
  @param    lines_list      Lines List (flux, wavelengths)
  @param    display         Flag to display results
  @return  Wavelength solution, i.e. polynomial that translates pixel
            values to wavelength.

    TODO: Summarize method
 */
/*----------------------------------------------------------------------------*/
cpl_polynomial * cr2res_wave_xcorr(
        cpl_vector      *   spectrum,
        cpl_polynomial  *   initial_guess,
        cpl_bivector    *   lines_list,
        int                 display)
{
    /* Check Entries */
    if (spectrum == NULL || initial_guess == NULL || lines_list == NULL)
        return NULL ;


    if (display) {
        /* Plot Reference Spectrum */

        /* Plot Extracted Spectrum */
    }

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
    cpl_vector  *   peaks;
	double			trueD;

    peaks = cr2res_wave_etalon_measure_fringes(spectrum);
	trueD = cr2res_wave_etalon_fringe_stats(peaks, initial_guess);
    cpl_vector_delete(peaks);
    return NULL ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief Find the true D from fringe statistics
  @param
  @return
 */
/*----------------------------------------------------------------------------*/

double cr2res_wave_etalon_fringe_stats(
            cpl_vector      * peaks,
            cpl_polynomial  * initial_guess)
{
	int				i;
	cpl_size		num_peaks;
    double      	trueD=-1.0;
	cpl_vector	*	diffs;
	cpl_vector	*	waves;


    /* cpl_plot_vector("", "w lines", "", peaks) ; */


	num_peaks = cpl_vector_get_size(peaks);
    waves = cr2res_polynomial_eval_vector(initial_guess, peaks);
	diffs = cpl_vector_new(num_peaks-1);
	for (i=1; i<num_peaks; i++){
		cpl_vector_set(diffs,i-1,
			cpl_vector_get(waves,i) - cpl_vector_get(waves,i-1) );
	}

    if (cpl_msg_get_level() == CPL_MSG_DEBUG){
        cpl_table   *   tab;
        tab = cpl_table_new(num_peaks-1);
        cpl_table_new_column(tab, "wavediff", CPL_TYPE_DOUBLE) ;
        for(i=0; i<num_peaks-1; i++) {
            cpl_table_set_double(tab, "wavediff", i, cpl_vector_get(diffs, i));
        }

        if ( cpl_table_save(tab, NULL, NULL, "debug_wavediffs.fits",
                CPL_IO_EXTEND) == CPL_ERROR_FILE_NOT_FOUND) {
            cpl_error_reset() ;
            cpl_table_save(tab, NULL, NULL, "debug_wavediffs.fits",
                                CPL_IO_CREATE);
        }
        cpl_table_delete(tab);
    }
	cpl_vector_delete(diffs);
	cpl_vector_delete(waves);
    return trueD;
}

/*----------------------------------------------------------------------------*/
/**
  @brief Identify and fit etalon lines
  @param spectrum The input spectrum vector
  @return Vector with the fitted peak positions, in pixels.

  The peak positions start with 1 for the first pixel !!
    TODO : Add unit test with an atificial vector
 */
/*----------------------------------------------------------------------------*/

cpl_vector * cr2res_wave_etalon_measure_fringes(
            cpl_vector * spectrum)
{
    cpl_array   *   peaks;
    cpl_vector  *   peak_vec;
    cpl_vector  *   spec_thresh;
    cpl_vector  *   cur_peak;
    cpl_vector  *   X_all, *X_peak;
    int             i, j, k ;
    int             numD = 0 ;
    int             smooth = 70 ;   // TODO: make free parameter?
                                    // interfringe ~30 in Y, ~70 in K
    int             thresh = 1 ;   // TODO: derive from read-out noise
    int             max_num_peaks = 256 ;
    int             max_len_peak = 256 ;
    int             min_len_peak = 5 ; //TODO: tweak or make parameter?;
    cpl_size        nx ;
    double          spec_i;
    double          prev_peak;
    double          x0, sigma, area, offset;

    nx = cpl_vector_get_size(spectrum) ;
    spec_thresh = cr2res_threshold_spec(spectrum, smooth, thresh) ;
    //cpl_plot_vector("", "w lines", "", spec_thresh) ;

    if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
        cpl_vector_save(spec_thresh, "debug_thresh.fits", CPL_TYPE_DOUBLE,
                NULL, CPL_IO_CREATE);
    }

    /*Output array, values are invalid until set.*/
    peaks = cpl_array_new(max_num_peaks, CPL_TYPE_DOUBLE);

    /* X-axis to cut out from for each peak */
    X_all = cpl_vector_new(max_len_peak);
    for (i=0; i<max_len_peak; i++) cpl_vector_set(X_all, i, (double)i+1) ;

    for (i=0; i < nx; i++){
        j = 0;
        while ( (spec_i = cpl_vector_get(spec_thresh, i)) > 0 ) {
            j++;
            i++;
        }
        if (j < min_len_peak) continue;
        cpl_msg_debug(__func__, "Peak length j=%d at i=%d",j,i);
        cur_peak = cpl_vector_extract(spec_thresh, i-j, i, 1) ;
        X_peak = cpl_vector_extract(X_all, 0, j, 1) ;

        if (cpl_vector_fit_gaussian(X_peak, NULL, cur_peak, NULL, CPL_FIT_ALL,
                                &x0, &sigma, &area, &offset,
                                NULL,NULL,NULL) != CPL_ERROR_NONE ) {
            cpl_msg_warning(__func__, "Fit at j=%d i=%d failed",j,i);
            cpl_vector_delete(cur_peak);
            cpl_vector_delete(X_peak);
            continue;
        }

        // Shift x0 to absolute position
        x0 += i-j ;
        cpl_msg_debug(__func__,"Fit: %.2f, %.2f, %.2f, %.2f",
                                    x0, sigma, area, offset);
        if ((k = cpl_array_count_invalid(peaks)) <1)
            cpl_msg_error(__func__,"Output array overflow!");
        cpl_msg_debug(__func__,"k=%d, x0=%d",k,x0);
        cpl_array_set_double(peaks, max_num_peaks - k, x0);

        cpl_vector_delete(cur_peak);
        cpl_vector_delete(X_peak);
    }

    /* Copy into output array, treating first element outside the loop */
    k = max_num_peaks - cpl_array_count_invalid(peaks);
    peak_vec = cpl_vector_new(k) ;
    for (i=1; i<k; i++)
        cpl_vector_set(peak_vec, i, cpl_array_get(peaks, i, NULL) );

    cpl_vector_delete(spec_thresh) ;
    cpl_vector_delete(X_all) ;
    cpl_array_delete(peaks) ;
    return peak_vec;
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
cpl_bivector * cr2res_wave_gen_lines_spectrum(
        const char      *   catalog,
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

  wmin = poly(1)
  wmax = poly((CR2RES_DETECTOR_SIZE)

  The returned array must be deallocated with cpl_array_delete()
 */
/*----------------------------------------------------------------------------*/
cpl_array * cr2res_wave_get_estimate(
        const char  *   filename,
        int             detector,
        int             order)
{
    double                  wmin, wmax, a, b ;
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

    /* Compute polynomial coefficients */
    b = (wmax - wmin) / (CR2RES_DETECTOR_SIZE-1) ;
    a = wmin - b ;

    /* Create the array */
    wl = cpl_array_new(2, CPL_TYPE_DOUBLE) ;
    cpl_array_set(wl, 0, a) ;
    cpl_array_set(wl, 1, b) ;
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
