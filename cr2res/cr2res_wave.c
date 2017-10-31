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

#include "irplib_wlxcorr.h"

#include "cr2res_dfs.h"
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
  @param    degree          The polynomial degree of the solution
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
        int                     degree,
        int                     display)
{
    cpl_polynomial  *   solution ;
    cpl_bivector    *   ref_spectrum ;
    int                 wl_error ;

    /* Initialise */
    solution = NULL ;
    wl_error = 50 ;

    /* Switch on the possible methods */
    if (wavecal_type == CR2RES_LAMP) {
        if (line_fitting) {
            solution = cr2res_wave_line_fitting(spectrum, initial_guess, NULL) ;
        } else {
            /* Create the lines spectrum from the lines list */
            if ((ref_spectrum = cr2res_wave_gen_lines_spectrum(static_file,
                    initial_guess, wl_error)) == NULL) {
                cpl_msg_error(__func__, "Cannot generate catalog spectrum");
                return NULL ;
            }
            solution = cr2res_wave_xcorr(spectrum, initial_guess,
                    wl_error, ref_spectrum, degree, display) ;
            cpl_bivector_delete(ref_spectrum) ;
        }
    } else if (wavecal_type == CR2RES_GAS) {
        solution = cr2res_wave_xcorr(spectrum, initial_guess, wl_error, NULL,
                degree, display) ;
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
  @param    wl_error        Max error in pixels of the initial guess
  @param    lines_list      Lines List (flux, wavelengths)
  @param    degree          The polynomial degree
  @param    display         Flag to display results
  @return  Wavelength solution, i.e. polynomial that translates pixel
            values to wavelength.

    TODO: Summarize method
 */
/*----------------------------------------------------------------------------*/
cpl_polynomial * cr2res_wave_xcorr(
        cpl_vector      *   spectrum,
        cpl_polynomial  *   initial_guess,
        int                 wl_error,
        cpl_bivector    *   lines_list,
        int                 degree,
        int                 display)
{
    cpl_vector          *   wl_errors ;
    cpl_polynomial      *   sol ;
    cpl_polynomial      *   sol_guess ;
    cpl_vector          *   spec_clean ;
    double              *   pspec_clean ;
    cpl_vector          *   filtered ;
    cpl_vector          *   xcorrs ;
    cpl_table           *   spc_table ;
    double                  wl_min, wl_max, wl_error_wl, wl_error_pix ;
    double                  xc, slit_width, fwhm ;
    int                     i, nsamples, clean_spec, filt_size, degree_loc ;

    /* Check Entries */
    if (spectrum == NULL || initial_guess == NULL || lines_list == NULL)
        return NULL ;

    /* Initialise */
    clean_spec = 1 ;
    slit_width = 2.0 ;
    fwhm = 2.0 ;
    filt_size = 9 ;

    /* Compute wl boundaries */
    wl_min = cpl_polynomial_eval_1d(initial_guess, 1, NULL);
    wl_max = cpl_polynomial_eval_1d(initial_guess, CR2RES_DETECTOR_SIZE, NULL);

    /* Clean the spectrum from the low frequency signal if requested */
    if (clean_spec) {
        cpl_msg_info(__func__, "Low Frequency removal from spectrum") ;
        cpl_msg_indent_more() ;
        /* Subrtract the low frequency part */
        if ((filtered=cpl_vector_filter_median_create(spectrum,
                        filt_size))==NULL){
            cpl_msg_error(__func__, "Cannot filter the spectrum") ;
            spec_clean = cpl_vector_duplicate(spectrum) ;
        } else {
            spec_clean = cpl_vector_duplicate(spectrum) ;
            cpl_vector_subtract(spec_clean, filtered) ;
            cpl_vector_delete(filtered) ;
        }
        cpl_msg_indent_less() ;
    } else {
        spec_clean = cpl_vector_duplicate(spectrum) ;
    }

    /* Remove Negative values */
    pspec_clean = cpl_vector_get_data(spec_clean) ;
    for (i=0 ; i<cpl_vector_get_size(spec_clean) ; i++)
        if (pspec_clean[i] < 0.0) pspec_clean[i] = 0 ;

    /* Display */
    if (display) {
        /* Plot Catalog Spectrum */
        cpl_plot_bivector(
                "set grid;set xlabel 'Wavelength (nm)';set ylabel 'Emission';",
                "t 'Catalog Spectrum' w impulses", "", lines_list);
        /* Plot Extracted Spectrum */
        /*
        cpl_plot_vector(
    "set grid;set xlabel 'Position (Pixel)';set ylabel 'Intensity (ADU/sec)';",
                "t 'Extracted spectrum' w lines", "", spectrum);
        */
        /* Plot Extracted Spectrum */
        cpl_plot_vector(
    "set grid;set xlabel 'Position (Pixel)';set ylabel 'Intensity (ADU/sec)';",
                "t 'Cleaned Extracted spectrum' w lines", "", spec_clean);
    }

    /* Pass #1 */
    degree_loc = 1 ;
    nsamples = 100 ;
    wl_error_pix = wl_error ;
    sol_guess = initial_guess ;

    wl_error_wl = (wl_max-wl_min)*wl_error_pix/CR2RES_DETECTOR_SIZE ;
    wl_errors = cpl_vector_new(degree_loc+1) ;
    cpl_vector_fill(wl_errors, wl_error_wl) ;
    cpl_msg_info(__func__,
    "Pass #1 : Degree %d / Error %g nm (%g pix) / %d samples -> %g Polys",
            degree_loc,wl_error_wl,wl_error_pix,nsamples,pow(nsamples,
                degree_loc+1)) ;
    cpl_msg_indent_more() ;
    if ((sol=irplib_wlxcorr_best_poly(spec_clean, lines_list, degree_loc,
                    sol_guess, wl_errors, nsamples, slit_width, fwhm, &xc, NULL,
                    &xcorrs)) == NULL) {
        cpl_msg_error(__func__, "Cannot get the best polynomial") ;
        cpl_msg_indent_less() ;
        cpl_vector_delete(wl_errors) ;
        cpl_vector_delete(spec_clean) ;
        if (xcorrs != NULL) cpl_vector_delete(xcorrs) ;
        cpl_error_reset() ;
        return NULL ;
    }
    cpl_vector_delete(wl_errors) ;
    cpl_msg_info(__func__, "Cross-Correlation factor: %g", xc) ;

    /* Plot the correlation values */
    if (display) {
        cpl_plot_vector("set grid;", "t 'Correlation values (Pass #1)' w lines",
                "", xcorrs) ;
    }
    if (xcorrs != NULL) cpl_vector_delete(xcorrs) ;
    cpl_msg_indent_less() ;

    /* Pass #2 */
    degree_loc = degree ;
    nsamples = 10 ;
    wl_error_pix = wl_error/2.0 ;
    sol_guess = sol ;

    wl_error_wl = (wl_max-wl_min)*wl_error_pix/CR2RES_DETECTOR_SIZE ;
    wl_errors = cpl_vector_new(degree_loc+1) ;
    cpl_vector_fill(wl_errors, wl_error_wl) ;
    cpl_msg_info(__func__,
    "Pass #2 : Degree %d / Error %g nm (%g pix) / %d samples -> %g Polys",
            degree_loc,wl_error_wl,wl_error_pix,nsamples,
            pow(nsamples,degree_loc+1)) ;
    cpl_msg_indent_more() ;
    if ((sol=irplib_wlxcorr_best_poly(spec_clean, lines_list, degree_loc,
                    sol_guess, wl_errors, nsamples, slit_width, fwhm, &xc, NULL,
                    &xcorrs)) == NULL) {
        cpl_msg_error(__func__, "Cannot get the best polynomial") ;
        cpl_msg_indent_less() ;
        cpl_vector_delete(wl_errors) ;
        cpl_vector_delete(spec_clean) ;
        if (xcorrs != NULL) cpl_vector_delete(xcorrs) ;
        cpl_polynomial_delete(sol_guess) ;
        cpl_error_reset() ;
        return NULL ;
    }
    cpl_vector_delete(wl_errors) ;
    cpl_msg_info(__func__, "Cross-Correlation factor: %g", xc) ;

    /* Plot the correlation values */
    if (display) {
        cpl_plot_vector("set grid;", "t 'Correlation values (Pass #2)' w lines",
                "", xcorrs) ;
    }
    if (xcorrs != NULL) cpl_vector_delete(xcorrs) ;
    cpl_msg_indent_less() ;

    /* Plot results table */
    if (display) {
        /* Create the spc_table  */
        if ((spc_table = irplib_wlxcorr_gen_spc_table(spec_clean, lines_list,
                        slit_width, fwhm, sol_guess, sol)) == NULL) {
            cpl_msg_error(cpl_func, "Cannot generate infos table") ;
        } else {
			/* Plot */
			irplib_wlxcorr_plot_spc_table(spc_table, "XC", 0, 0) ;
			cpl_table_delete(spc_table) ;
        }
    }
    cpl_vector_delete(spec_clean) ;
    cpl_polynomial_delete(sol_guess) ;

    return sol ;
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
    cpl_vector  *   peaks_shouldbe;
	double			trueD, x_0;
    int             npeaks, i;

    peaks = cr2res_wave_etalon_measure_fringes(spectrum);
    npeaks=cpl_vector_get_size(peaks);
    x_0 = cpl_vector_get(peaks,0);

	trueD = cr2res_wave_etalon_fringe_stats(peaks, initial_guess);

    peaks_shouldbe = cpl_vector_new(npeaks);
    for (i=0; i<npeaks; i++) {
        cpl_vector_set(peaks_shouldbe, i, x_0 + (trueD*i));
    }

    if (cpl_msg_get_level() == CPL_MSG_DEBUG){
        cpl_vector_dump(peaks, stdout);
        cpl_vector_dump(peaks_shouldbe, stdout);
    }

    cpl_vector_delete(peaks);
    cpl_vector_delete(peaks_shouldbe);
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

    trueD = cpl_vector_get_median(diffs);
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
        cpl_msg_debug(__func__,"k=%d, x0=%g",k,x0);
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
  @brief    Load the emission lines in a bivector
  @param    catalog         The catalog file
  @param    initial_guess   The wavelength polynomial
  @param    wl_error        Max error in pixels of the initial guess
  @return   The lines spectrum
 */
/*----------------------------------------------------------------------------*/
cpl_bivector * cr2res_wave_gen_lines_spectrum(
        const char      *   catalog,
        cpl_polynomial  *   initial_guess,
        int                 wl_error)
{
    cpl_bivector    *   lines ;
    cpl_bivector    *   lines_sub ;
    double          *   lines_sub_wl ;
    double          *   lines_sub_intens ;
    double              wl_error_wl, wl_min, wl_max ;
    int                 i ;

    /* Check Entries */
    if (catalog == NULL || initial_guess == NULL) return NULL ;

    /* Load the lines */
    lines = cr2res_io_load_EMISSION_LINES(catalog) ;

    /* Extract the needed spectrum */
    wl_min = cpl_polynomial_eval_1d(initial_guess, 1, NULL);
    wl_max = cpl_polynomial_eval_1d(initial_guess, CR2RES_DETECTOR_SIZE, NULL);
    wl_error_wl = (wl_max-wl_min)*wl_error/CR2RES_DETECTOR_SIZE ;
    lines_sub = irplib_wlxcorr_cat_extract(lines, wl_min-wl_error_wl,
            wl_max+wl_error_wl) ;

	/* Zero the beginning and the end */
    lines_sub_wl = cpl_bivector_get_x_data(lines_sub) ;
    lines_sub_intens = cpl_bivector_get_y_data(lines_sub) ;
    for (i=0 ; i<cpl_bivector_get_size(lines_sub) ; i++) {
        if (wl_min > 0)
            if (lines_sub_wl[i] < wl_min) lines_sub_intens[i] = 0.0 ;
        if (wl_max > 0)
            if (lines_sub_wl[i] > wl_max) lines_sub_intens[i] = 0.0 ;
    }

    /* Free and return */
    cpl_bivector_delete(lines) ;
    return lines_sub ;
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
    wmin = cr2res_pfits_get_wstrt(plist, order) ;
    wmax = cr2res_pfits_get_wend(plist, order) ;
    cpl_propertylist_delete(plist) ;
    if (cpl_error_get_code() != CPL_ERROR_NONE) {
        cpl_msg_error(__func__,
                "Cannot get WSTRT/WEND from header for Detector %d / Order %d",
                detector, order) ;
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
    hdrl_image      *   out ;
    cpl_image       *   out_ima ;
    double          *   pout_ima ;
    const cpl_array *   tmp_array ;
    cpl_polynomial  *   wave_poly ;
    cpl_polynomial  *   upper_poly ;
    cpl_polynomial  *   lower_poly ;
    double              upper_pos, lower_pos, wavelength ;
    cpl_size            i, j, k, nrows, nx, ny ;

    /* Check Entries */
    if (trace_wave == NULL) return NULL ;

    /* Initialise */
    nrows = cpl_table_get_nrow(trace_wave) ;

    /* Create the image */
    out = hdrl_image_new(CR2RES_DETECTOR_SIZE, CR2RES_DETECTOR_SIZE) ;
    out_ima = hdrl_image_get_image(out) ;
    nx = cpl_image_get_size_x(out_ima) ;
    ny = cpl_image_get_size_y(out_ima) ;
    pout_ima = cpl_image_get_data_double(out_ima) ;

    /* Loop on the traces */
    for (k=0 ; k<nrows ; k++) {
        /* Check if there is a Wavelength Polynomial available */
        tmp_array = cpl_table_get_array(trace_wave, CR2RES_COL_WAVELENGTH, k) ;
        wave_poly = cr2res_convert_array_to_poly(tmp_array) ;
        if (wave_poly != NULL) {
            /* Get the Upper Polynomial */
            tmp_array = cpl_table_get_array(trace_wave, CR2RES_COL_UPPER, k) ;
            upper_poly = cr2res_convert_array_to_poly(tmp_array) ;

            /* Get the Lower Polynomial */
            tmp_array = cpl_table_get_array(trace_wave, CR2RES_COL_LOWER, k) ;
            lower_poly = cr2res_convert_array_to_poly(tmp_array) ;

            /* Check if all Polynomials are available */
            if (upper_poly == NULL || lower_poly == NULL) {
                if (upper_poly != NULL) cpl_polynomial_delete(upper_poly) ;
                if (lower_poly != NULL) cpl_polynomial_delete(lower_poly) ;
                cpl_msg_warning(__func__, "Cannot get UPPER/LOWER information");
                cpl_polynomial_delete(wave_poly) ;
                continue ;
            }

            /* Set the Pixels in the trace */
            for (i=0 ; i<nx ; i++) {
                upper_pos = cpl_polynomial_eval_1d(upper_poly, i+1, NULL) ;
                lower_pos = cpl_polynomial_eval_1d(lower_poly, i+1, NULL) ;
                wavelength = cpl_polynomial_eval_1d(wave_poly, i+1, NULL) ;
                for (j=0 ; j<ny ; j++) {
                    if (j+1 >= lower_pos && j+1 <= upper_pos)
                        pout_ima[i+j*nx] = wavelength ;
                }
            }
            cpl_polynomial_delete(wave_poly) ;
            cpl_polynomial_delete(upper_poly) ;
            cpl_polynomial_delete(lower_poly) ;
        }
    }
    return out ;
}
/**@}*/
