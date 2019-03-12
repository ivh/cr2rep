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

static int poly(const double x[], const double a[], double * result) ;
static int deriv_poly(const double x[], const double a[], double * result) ;
static int gauss(const double x[], const double a[], double * result) ;
static int gauss_derivative(const double x[], const double a[], 
        double * result) ;
cpl_polynomial  * polyfit_1d(
    cpl_matrix  * px, 
    cpl_vector  * py,
    cpl_vector  * sigma_py,
    int degree, 
    const cpl_polynomial * solution_init,
    cpl_array   ** wavelength_error,
    cpl_vector  ** sigma_fit,
    cpl_matrix  ** cov);
int cr2res_wave_extract_lines(
    cpl_bivector    *   spectrum,
    cpl_bivector    *   spectrum_err,
    cpl_polynomial  *   wavesol_init,
    const cpl_array *   wave_error_init,
    cpl_bivector    *   lines_list,
    int                 display,
    cpl_matrix      **  px,
    cpl_vector      **  py,
    cpl_vector      **  sigma_py,
    cpl_vector      **  heights);

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_wave        Wavelength Calibration
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    1D Wavelength Calibration
  @param    spectrum        Extracted spectrum
  @param    spectrum        Extracted spectrum error
  @param    wavesol_init    Initial wavelength solution
  @param    wave_error_init Initial wavelength error (can be NULL)
  @param    catalog         Line catalog or template spectrum or NULL
  @param    degree          The polynomial degree of the solution
  @param    display         Flag to display results
  @param    wavelength_error    [out] array of wave_mean_error, wave_max_error
  @return   Wavelength solution, i.e. polynomial that translates pixel
            values to wavelength.
 */
/*----------------------------------------------------------------------------*/
cpl_polynomial * cr2res_wave_1d(
        cpl_bivector        *   spectrum,
        cpl_bivector        *   spectrum_err,
        cpl_polynomial      *   wavesol_init,
        const cpl_array     *   wave_error_init,
        cr2res_wavecal_type     wavecal_type,
        const char          *   static_file,
        int                     degree,
        int                     display,
        cpl_array           **  wavelength_error)
{
    cpl_polynomial  *   solution ;
    cpl_bivector    *   ref_spectrum ;
    cpl_bivector    *   simple_ref_spectrum ;
    const cpl_bivector    **  plot;
    int                 wl_error ;

    /* Check Inputs */
    if (spectrum == NULL || spectrum_err == NULL || wavesol_init == NULL)
        return NULL ;
    if ((wavecal_type == CR2RES_XCORR || wavecal_type == CR2RES_LINE1D) &&
            static_file == NULL) return NULL ;

    /* Initialise */
    solution = NULL ;
    wl_error = 100 ;
    *wavelength_error = NULL ;

    /* Create the lines spectrum from the lines list */
    ref_spectrum = cr2res_wave_gen_lines_spectrum(static_file, wavesol_init,
            wl_error) ;
            
    /* Just Extract the lines from the catalog */
    simple_ref_spectrum = cr2res_io_load_EMISSION_LINES(static_file) ;

    /* Switch on the possible methods */
    if (wavecal_type == CR2RES_XCORR) {
        solution = cr2res_wave_xcorr(spectrum, wavesol_init, wl_error,
                ref_spectrum, degree, display) ;
    } else if (wavecal_type == CR2RES_LINE1D) {
        solution = cr2res_wave_line_fitting(spectrum, spectrum_err,
                wavesol_init, wave_error_init, simple_ref_spectrum, degree,
                display, NULL, wavelength_error) ;
    } else if (wavecal_type == CR2RES_ETALON) {
        solution = cr2res_wave_etalon(spectrum, spectrum_err, wavesol_init, 
                degree, wavelength_error);
    }

    if (ref_spectrum != NULL) cpl_bivector_delete(ref_spectrum) ;
    if (simple_ref_spectrum != NULL) cpl_bivector_delete(simple_ref_spectrum) ;

    return solution ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Compute the 2D wavelength polynomial based on a line spectrum
            and a reference catalog by finding lines and fitting
  @param    spectra         List of extracted spectra 
  @param    spectra_err     List of extracted spectra errors
  @param    wavesol_init    List of Initial wavelength solutions
  @param    wavesol_init_err List of Initial wavelength error (can be NULL)
  @param    orders          List of orders of the various spectra
  @param    ninputs         Number of entries in the previous parameters
  @param    catalog_spec    Catalog spectrum
  @param    degree_x        The polynomial degree in x
  @param    degree_y        The polynomial degree in y
  @param    display         Flag to display results
  @param    wavesol_error   [out] array of wave_mean_error, wave_max_error
  @return   Wavelength solution, i.e. polynomial that translates pixel
            values to wavelength.

    NOTES
        - Some input spectra might be NULL - be robust against this.
        - Avoid hardcoded  values (2)
        - Limit the line length to 80 characters
        - Keep the interface stable
        - Add Documentation in the function code
        - Remove dead code 
        - All declarations at the beginning of a function
        - Make sure to deallocate everything
 */
/*----------------------------------------------------------------------------*/
cpl_polynomial * cr2res_wave_2d(
        cpl_bivector        **  spectra,
        cpl_bivector        **  spectra_err,
        cpl_polynomial      **  wavesol_init,
        const cpl_array     **  wavesol_init_err,
        int                 *   orders,
        int                     ninputs,
        cpl_bivector        *   catalog_spec,
        cpl_size                degree_x,
        cpl_size                degree_y,
        int                     display,
        cpl_array           **  wavesol_error)
{
    cpl_vector      *   diff;
    cpl_size            old, new, i, j, k, spec_size ;
    cpl_polynomial  *   result ;
    cpl_error_code      error;
    cpl_size            degree_2d[2];
    cpl_matrix      *   tmp_x;
    cpl_vector      *   tmp_y;
    cpl_vector      *   tmp_sigma;
    cpl_vector      *   pos;
    cpl_matrix      *   px ;
    cpl_matrix      *   sigma_px ;
    cpl_vector      *   py ;
    cpl_vector      *   sigma_py ;
    int                 n ;

    /* Check Inputs */
    if (spectra==NULL || spectra_err==NULL || wavesol_init==NULL || 
            orders==NULL || catalog_spec==NULL) 
        return NULL ;

    /* Initialise */
    n = cpl_bivector_get_size(catalog_spec);
    px = sigma_px = NULL ;
    py = sigma_py = NULL ;

    result = cpl_polynomial_new(2);

    for (i = 0; i < ninputs; i++){
        // extract line data in 1 spectrum
        cr2res_wave_extract_lines(spectra[i], spectra_err[i], wavesol_init[i],
            wavesol_init_err[i], catalog_spec, display, &tmp_x, &tmp_y, 
            &tmp_sigma, NULL);

        // append new data onto existing vectors/matrices
        new = cpl_vector_get_size(tmp_y);

        if (px == NULL){
            // First order to run
            px = cpl_matrix_new(2, new);
            py = cpl_vector_new(new);
            sigma_py = cpl_vector_new(new);
            old = 0;
        } else {
            old = cpl_vector_get_size(py);
            cpl_vector_set_size(py, old + new);
            cpl_vector_set_size(sigma_py, old + new);
            cpl_matrix_set_size(px, 2, old + new);
        }

        for (j = 0; j < new; j++){
            cpl_vector_set(py, old + j, cpl_vector_get(tmp_y, j));
            cpl_vector_set(sigma_py, old + j, cpl_vector_get(tmp_sigma, j));
            cpl_matrix_set(px, 0, old + j, cpl_matrix_get(tmp_x, j, 0));
            cpl_matrix_set(px, 1, old + j, orders[i]);
        }
        cpl_matrix_delete(tmp_x);
        cpl_vector_delete(tmp_y);
        cpl_vector_delete(tmp_sigma);
    }

    degree_2d[0] = degree_x ;
    degree_2d[1] = degree_y ;
    error = cpl_polynomial_fit(result, px, NULL, py, NULL, TRUE, NULL,
            degree_2d);

    if (error == CPL_ERROR_NONE){
        if (wavesol_error != NULL){
            // Calculate absolute difference between polynomial and 
            // catalog value for each line
            // use px and py, so that only good lines are used
            diff = cpl_vector_new(cpl_vector_get_size(py));
            pos = cpl_vector_new(2);
            for (i = 0; i < cpl_vector_get_size(py); i++){
                cpl_vector_set(pos, 0, cpl_matrix_get(px, i, 0));
                cpl_vector_set(pos, 1, cpl_matrix_get(px, i, 1));
                cpl_vector_set(diff, i, abs(
                    cpl_polynomial_eval(result, pos)
                    - cpl_vector_get(py, i)));
            }
            // Set wavesol_error to mean and max difference
            cpl_array_set_double(*wavesol_error, 0, 
                    cpl_vector_get_mean(diff));
            cpl_array_set_double(*wavesol_error, 1, 
                    cpl_vector_get_max(diff));
            cpl_vector_delete(diff);
            cpl_vector_delete(pos);
        }
    }
    cpl_matrix_delete(px);
    cpl_vector_delete(py);
    cpl_vector_delete(sigma_py);

    // in case something went wrong during fitting
    if (error != CPL_ERROR_NONE){
        cpl_polynomial_delete(result);
        cpl_error_reset();
        return NULL;
    }
    return result;
}

/*----------------------------------------------------------------------------*/
/**
  @brief  Find solution by cross-correlating template spectrum
  @param    spectrum        Input spectrum
  @param    wavesol_init    Starting wavelength solution
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
        cpl_bivector    *   spectrum,
        cpl_polynomial  *   wavesol_init,
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
    if (spectrum == NULL || wavesol_init == NULL || lines_list == NULL)
        return NULL ;

    /* Initialise */
    clean_spec = 1 ;
    slit_width = 2.0 ;
    fwhm = 2.0 ;
    filt_size = 9 ;

    /* Compute wl boundaries */
    wl_min = cpl_polynomial_eval_1d(wavesol_init, 1, NULL);
    wl_max = cpl_polynomial_eval_1d(wavesol_init, CR2RES_DETECTOR_SIZE, NULL);

    cpl_msg_info(__func__, "Wl Range Input : %g - %g", wl_min, wl_max) ;

    /* Clean the spectrum from the low frequency signal if requested */
    if (clean_spec) {
        cpl_msg_info(__func__, "Low Frequency removal from spectrum") ;
        cpl_msg_indent_more() ;
        /* Subrtract the low frequency part */
        if ((filtered=cpl_vector_filter_median_create(
                        cpl_bivector_get_y(spectrum),
                        filt_size))==NULL){
            cpl_msg_error(__func__, "Cannot filter the spectrum") ;
            spec_clean = cpl_vector_duplicate(cpl_bivector_get_y(spectrum)) ;
        } else {
            spec_clean = cpl_vector_duplicate(cpl_bivector_get_y(spectrum)) ;
            cpl_vector_subtract(spec_clean, filtered) ;
            cpl_vector_delete(filtered) ;
        }
        cpl_msg_indent_less() ;
    } else {
        spec_clean = cpl_vector_duplicate(cpl_bivector_get_y(spectrum)) ;
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
        cpl_plot_vector(
    "set grid;set xlabel 'Position (Pixel)';set ylabel 'Intensity (ADU/sec)';",
                "t 'Cleaned Extracted spectrum' w lines", "", spec_clean);
    }

    /* Pass #1 */
    degree_loc = 1 ;
    nsamples = 100 ;
    wl_error_pix = wl_error ;
    sol_guess = wavesol_init ;

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
    wl_error_pix = wl_error ;
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

            /* TODO : Return the table as a product CR2RES_PROTYPE_XCORR */
			cpl_table_delete(spc_table) ;
        }
    }
    cpl_vector_delete(spec_clean) ;
    cpl_polynomial_delete(sol_guess) ;

    return sol ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief  Extract line positions in spectrum
  @param    spectrum        Input spectrum
  @param    spectrum_err    Input spectrum error
  @param    wavesol_init    Starting wavelength solution
  @param    wl_error        Max error in pixels of the initial guess
  @param    lines_list      Lines List (flux, wavelengths)
  @param    degree          The polynomial degree
  @param    display         Flag to display results
  @return  Wavelength solution, i.e. polynomial that translates pixel
            values to wavelength.

    TODO: Summarize method
 */
/*----------------------------------------------------------------------------*/
int cr2res_wave_extract_lines(
    cpl_bivector    *   spectrum,
    cpl_bivector    *   spectrum_err,
    cpl_polynomial  *   wavesol_init,
    const cpl_array *   wave_error_init,
    cpl_bivector    *   lines_list,
    int                 display,
    cpl_matrix      **  px,
    cpl_vector      **  py,
    cpl_vector      **  sigma_py,
    cpl_vector      **  heights)
{

    /* Check Entries */
    if (spectrum == NULL || spectrum_err == NULL || wavesol_init == NULL ||
            lines_list == NULL){
        return -1;
    }

    cpl_size power = 1;
    int window_size = 500;
    /* set window_size using the wave_error_init, scaled by the initial guess */

    if (wave_error_init != NULL &&0)
        if (cpl_array_get_double(wave_error_init, 1, NULL) > 0){
            window_size = 10 * ceil(cpl_array_get_double(wave_error_init, 1, NULL) /
                            fabs(cpl_polynomial_get_coeff(wavesol_init, &power)));
        }
    
    cpl_size i, j, k, ngood, spec_size, npossible;
    double pixel_pos, pixel_new, red_chisq, dbl, res;
    int n = cpl_bivector_get_size(lines_list);
    cpl_error_code error;
    cpl_vector * wave_vec, * pixel_vec, *width_vec, *flag_vec, *height_vec;
    const cpl_vector *spec, *unc;
    double * wave, width;
    const double *height;
    double value, diff;
    double max_wl, min_wl;
    cpl_vector * fit, *fit_x;
    const cpl_vector ** plot;
    if (heights != NULL) *heights = cpl_vector_new(n);

    // For gaussian fit of each line
    // gauss = A * exp((x-mu)^2/(2*sig^2)) + cont
    cpl_matrix * x = cpl_matrix_new(window_size, 1);
    cpl_matrix * sigma_x = NULL;
    cpl_vector * y = cpl_vector_new(window_size);
    cpl_vector * sigma_y = cpl_vector_new(window_size);
    cpl_vector * a = cpl_vector_new(4);
    int ia[] = {1, 1, 1, 1};
    double x0, sigma, area, offset;

    spec = cpl_bivector_get_y_const(spectrum);
    unc = cpl_bivector_get_y_const(spectrum_err);

    // Prepare fit data vectors
    pixel_vec = cpl_vector_new(n);
    height_vec = cpl_vector_new(n);
    width_vec = cpl_vector_new(n);
    flag_vec = cpl_vector_new(n);
    cpl_vector_fill(flag_vec, 1);

    // The number of good lines, start with all
    ngood = 0;
    // The number of possible lines to fit, for debugging only
    npossible = 0;

    // evaluate the initial wavelength solution for all pixels
    // so that we can find the closest pixel position of each line
    spec_size = cpl_vector_get_size(spec);
    wave_vec = cpl_vector_new(spec_size);
    for (i = 0; i< spec_size; i++){
        cpl_vector_set(wave_vec, i, cpl_polynomial_eval_1d(wavesol_init, i, NULL));
    }
    max_wl = cpl_vector_get_max(wave_vec);
    min_wl = cpl_vector_get_min(wave_vec);


    // get line data
    wave = cpl_bivector_get_x_data(lines_list);
    height = cpl_bivector_get_y_data_const(lines_list);
    // TODO width is not provided in the catalog at the moment,
    // use half window size instead?
    width = 1;

    // for each line fit a gaussian around guessed position
    // and find actual pixel position
    for (i = 0; i < n; i++){

        // skip lines that are outside this wavelength region
        if ((wave[i] < min_wl) | (wave[i] > max_wl)){
            cpl_vector_set(flag_vec, i, 0);
            continue;
        }
        // The number of possible lines to fit, for debugging
        npossible++;

        // cut out a part of the spectrum around each line
        // assumes that the wavelength vector is ascending !!!
        pixel_pos = cpl_vector_find(wave_vec, wave[i]);
        value = 0;
        for (j = 0; j < window_size; j++){
            k = pixel_pos - window_size / 2 + j;
            if (k < 0 | k >= spec_size){
                // if the window reaches outside the spectrum
                // don't use the line
                cpl_vector_set(flag_vec, i, 0);
                break;
            }

            value = cpl_vector_get(spec, k);
            if (value < 0) value = 0;
            cpl_matrix_set(x, j, 0, k);
            cpl_vector_set(y, j, value);
            cpl_vector_set(sigma_y, j, cpl_vector_get(unc, k));
        }

        if (cpl_vector_get(flag_vec, i) == 0){
            // if the line was flagged as bad, skip the fit
            continue;
        }

        // Filter out bad pixels
        for (j = 1; j < window_size-1; j++){
            diff = 2 * cpl_vector_get(y, j) - cpl_vector_get(y, j-1) - cpl_vector_get(y, j+1);
            diff = fabs(diff);
            if (diff > 300){
                value = (cpl_vector_get(y, j-1) + cpl_vector_get(y, j+1)) / 2.;
                cpl_vector_set(y, j, value);
            }
        }
        cpl_vector_set(y, 0, 0);
        cpl_vector_set(y, window_size-1, 0);


        // get initial guess for gaussian fit
        value = pixel_pos - window_size / 2 + cpl_vector_get_maxpos(y);
        cpl_vector_set(a, 0, value);
        cpl_vector_set(a, 1, width);
        cpl_vector_set(a, 2, cpl_vector_get_max(y) - cpl_vector_get_min(y));
        cpl_vector_set(a, 3, cpl_vector_get_min(y));

        error = cpl_fit_lvmq(x, sigma_x, y, sigma_y, a, ia, &gauss, &gauss_derivative,
                        CPL_FIT_LVMQ_TOLERANCE, CPL_FIT_LVMQ_COUNT,
                        CPL_FIT_LVMQ_MAXITER, NULL, &red_chisq, NULL);

        // Set new pixel pos based on gaussian fit
        pixel_new = cpl_vector_get(a, 0);
        cpl_vector_set(pixel_vec, i, pixel_new);
        // width == uncertainty of wavelength position?
        cpl_vector_set(width_vec, i, cpl_vector_get(a, 1));
        cpl_vector_set(height_vec, i, cpl_vector_get(a, 2));

        // if fit to bad set flag to 0(False)
        // fit is bad, when 
        // 1) it caused an error
        // 2) its chi square is large
        // 3) the gaussian is centered outside the window   
        // 4) the fitted height is negative
        // 5) Peak is smaller than the noise level (SNR > 1) 

        if (error != CPL_ERROR_NONE 
            // | red_chisq > 100
            | fabs(pixel_new - pixel_pos) > window_size
            | cpl_vector_get(a, 2) < 0
            | cpl_vector_get(a, 2) < cpl_vector_get(a, 3) * 5.
        ){
            cpl_vector_set(flag_vec, i, 0);
            cpl_error_reset();
            continue;
        }
        ngood++;

        /* Display */
        if (display) {
            /* Plot Observation and Fit */
            // TODO: Currently this will create a plot window for each line
            // afaik, to get them all in one plot, one needs to save all the fits
            // and then run one big plot_vectors (maybe bivectors) command with set multiplot
            // alternatively create a file for each line, instead of plot window
            plot = cpl_malloc(3 * sizeof(cpl_vector*));
            fit = cpl_vector_new(window_size);
            fit_x = cpl_vector_new(window_size);

            for (j = 0; j < window_size; j++){
                dbl = cpl_matrix_get(x, j, 0);
                gauss(&dbl, cpl_vector_get_data_const(a), &res);
                cpl_vector_set(fit, j, res);
                // dbl = cpl_polynomial_eval_1d(wavesol_init, dbl, NULL);
                cpl_vector_set(fit_x, j, dbl);
            }
            plot[0] = fit_x;
            plot[1] = y;
            plot[2] = fit;
            cpl_plot_vectors(
            "set grid;set xlabel 'Position (Pixel)';set ylabel 'Intensity (ADU/sec)';",
                "title 'Observed' w lines", "q", plot, 3);
            cpl_vector_delete(fit);
            cpl_vector_delete(fit_x);
            cpl_free(plot);
            cpl_error_reset();
        }
    }

    cpl_msg_debug(__func__, "Using %lli out of %lli lines", ngood, npossible);

    // Set vectors/matrices for polyfit
    // only need space for good lines, ignoring bad ones
    *px = cpl_matrix_new(ngood, 1);
    *py = cpl_vector_new(ngood);
    *sigma_py = cpl_vector_new(ngood);

    k = 0;
    for (i = 0; i < n; i++){
        // Skip bad lines
        if (cpl_vector_get(flag_vec, i) == 1){
            cpl_matrix_set(*px, k, 0, cpl_vector_get(pixel_vec, i));
            cpl_vector_set(*py, k, wave[i]);
            cpl_vector_set(*sigma_py, k, fabs(cpl_vector_get(width_vec, i)));
            if (heights != NULL) cpl_vector_set(*heights, k, cpl_vector_get(height_vec, i));
            k++;
        }
    }
    cpl_matrix_delete(x);
    cpl_vector_delete(y);
    cpl_vector_delete(sigma_y);
    cpl_vector_delete(a);

    cpl_vector_delete(wave_vec);
    cpl_vector_delete(pixel_vec);
    cpl_vector_delete(width_vec);
    cpl_vector_delete(flag_vec);
    cpl_vector_delete(height_vec);

    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Compute the wavelength polynomial based on a line spectrum
            and a reference catalog by finding lines and fitting
  @param    spectrum        Observed spectrum
  @param    spectrum_err    Observed spectrum error
  @param    wavesol_init    Initial wavelength solution
  @param    wave_error_init Initial wavelength error (can be NULL)
  @param    lines_list      Lines List (flux, wavelengths)
  @param    degree          The polynomial degree
  @param    display         Flag to display results
  @param    sigma_fit       [out] uncertainties of the polynomial fit
                            parameters (may be NULL)
  @param    wavelength_error [out] array of wave_mean_error, wave_max_error (may be NULL),
                            if pointer to NULL, will create cpl_array which needs to be deleted
  @return  Wavelength solution, i.e. polynomial that translates pixel
            values to wavelength.
  The returned polynomial must be deallocated with cpl_polynomial_delete()
 */
/*----------------------------------------------------------------------------*/
cpl_polynomial * cr2res_wave_line_fitting(
        cpl_bivector    *   spectrum,
        cpl_bivector    *   spectrum_err,
        cpl_polynomial  *   wavesol_init,
        const cpl_array *   wave_error_init,
        cpl_bivector    *   lines_list,
        int                 degree,
        int                 display,
        cpl_vector      **  sigma_fit,
        cpl_array       **  wavelength_error)
{
    cpl_polynomial  *   result;
    cpl_matrix      *   cov;
    cpl_matrix      *   px;
    cpl_vector      *   py;
    cpl_vector      *   sigma_py;
    cpl_vector      *   heights;

    /* Check Entries */
    if (spectrum == NULL || spectrum_err == NULL || wavesol_init == NULL ||
            lines_list == NULL)
        return NULL;

    // extract line data in 1 spectrum
    if (cr2res_wave_extract_lines(spectrum, spectrum_err, wavesol_init, 
                wave_error_init, lines_list, display, &px, &py, &sigma_py, 
                &heights) != 0) {
        cpl_msg_error(__func__, "Cannot extract lines") ;
        return NULL ;
    }

    // fit polynomial to data points
    result = polyfit_1d(px, py, sigma_py, degree, wavesol_init, 
            wavelength_error, sigma_fit, &cov);

    if (display){
        cpl_size n = cpl_bivector_get_size(spectrum);
        cpl_size nlines = cpl_vector_get_size(py);
        const cpl_bivector ** plot = cpl_malloc(2 * sizeof(cpl_bivector*));

        cpl_vector * wave = cpl_bivector_get_x(spectrum);
        for (cpl_size i = 0; i < n; i++) 
            cpl_vector_set(wave, i, cpl_polynomial_eval_1d(result, i, NULL));

        cpl_bivector * lines = cpl_bivector_new(nlines);
        cpl_vector * pos = cpl_bivector_get_x(lines);
        cpl_vector * val = cpl_bivector_get_y(lines);

        for (cpl_size i = 0; i < nlines; i++){
            cpl_vector_set(pos, i, cpl_vector_get(py, i)); // Wavelength
            cpl_vector_set(val, i, cpl_vector_get(heights, i));        
        }
        plot[0] = spectrum;
        plot[1] = lines;

        const char* options[] = 
            {"title 'Observed' w lines", "title 'Lines' w points"};

        cpl_plot_bivectors(
"set grid;set xlabel 'Position (Wavelength)';set ylabel 'Intensity (ADU/sec)';",
            options , "", plot, 2);
        cpl_bivector_delete(lines);
        cpl_free(plot);
    }
    cpl_vector_delete(heights);
    cpl_matrix_delete(px);
    cpl_vector_delete(py);
    cpl_vector_delete(sigma_py);

    if (result != NULL) {
        cpl_matrix_delete(cov);
    }
    return result;
}
/*----------------------------------------------------------------------------*/
/**
  @brief   Find solution from etalon
  @param    spectrum        Input spectrum: etalon
  @param    wavesol_init   Starting wavelength solution
  @param    wavelength_error    [out] array of wave_mean_error, wave_max_error
  @return  Wavelength solution, i.e. polynomial that translates pixel
            values to wavelength.

    This function uses the intrinsic property of the etalon spectrum,
    that the lines are *supposed* to be equi-spaced, to refine the wavelength
    solution. The input solution needs to be good enough for the zero-point
    since the etalon spectrum carries no information about the absolute
    wavelength scale.

    The method involves these steps:
    * Identify lines (thresholding)
    * Determine line centers (gauss fit)
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
        cpl_bivector    *   spectrum,
        cpl_bivector    *   spectrum_err,
        cpl_polynomial  *   wavesol_init,
        int                 degree,
        cpl_array       **  wavelength_error)
{
    cpl_bivector *  is_should;
    cpl_vector  *   xi;
    cpl_vector  *   li;
    // cpl_vector  *   mi;
    cpl_vector  *   li_true;
    cpl_matrix  *   px;
    cpl_polynomial * result;
	double			l0, trueD;
    int             nxi, i, npeaks;

    if (spectrum == NULL | spectrum_err == NULL | 
        wavesol_init == NULL | degree < 0 | wavelength_error == NULL) return NULL;

    // Find etalon peaks xi in spectrum
    // TODO: Use Spectrum Error
    xi = cr2res_wave_etalon_measure_fringes(
            cpl_bivector_get_y(spectrum)); 
    nxi=cpl_vector_get_size(xi);

    /* apply initial solution to get wavelength li at each point xi*/
    li = cr2res_polynomial_eval_vector(wavesol_init, xi);

    /* Calculate delta lambda between peaks */
	trueD = cr2res_wave_etalon_get_D(li);
    cpl_msg_debug(__func__,"trueD: %e", trueD);

    /* Set vector with correct wavelength values */
    // expected number of peaks (+2 just to be sure)
    npeaks = (cpl_vector_get(li, nxi-1) - l0) / trueD + 2;
    li_true = cpl_vector_new(npeaks);
    for (i=0; i<npeaks; i++) {
        cpl_vector_set(li_true, i,l0 + (trueD*i));
    }

    // For each peak find the closest expected wavelength value
    is_should = cr2res_wave_etalon_assign_fringes(li, li_true);

    // polynomial fit to points, xi, li
    px = cpl_matrix_wrap(nxi, 1, cpl_vector_get_data(xi));
    result = polyfit_1d(px, cpl_bivector_get_y(is_should), NULL, 
                degree, wavesol_init, wavelength_error, NULL, NULL);

    cpl_matrix_unwrap(px);
    cpl_vector_delete(li_true);

    cpl_bivector_delete(is_should);
    cpl_vector_delete(xi);
    cpl_vector_delete(li);

    return result;
}

/*----------------------------------------------------------------------------*/
/**
  @brief Associate found fringes with the best "should"-value
  @param
  @return
 */
/*----------------------------------------------------------------------------*/

cpl_bivector * cr2res_wave_etalon_assign_fringes(
            const cpl_vector      * li,
            const cpl_vector      * li_true)
{
    int i, j, n;
    int * best_idx;
    double x,y;
    cpl_bivector * is_should;
    cpl_vector * is;
    cpl_vector * should;

    n = cpl_vector_get_size(li);
    is_should = cpl_bivector_new(n);
    is = cpl_bivector_get_x(is_should);
    should = cpl_bivector_get_y(is_should);
    for (i=0; i<n; i++) {
        x = cpl_vector_get(li, i);
        j = cpl_vector_find(li_true, x);
        cpl_vector_set(is, i, x);
        cpl_vector_set(should, i, cpl_vector_get(li_true, j));
    }
    return is_should;
}
/*----------------------------------------------------------------------------*/
/**
  @brief Find x0
  @param
  @return
 */
/*----------------------------------------------------------------------------*/

double cr2res_wave_etalon_get_x0(
            cpl_vector      * xi,
            cpl_polynomial  * wavesol_init)
{
    double x0, D;
    cpl_vector * li;
    cpl_vector * xs;
    int i, n;

    li = cr2res_polynomial_eval_vector(wavesol_init, xi);
    D = cr2res_wave_etalon_get_D(li);

    n = cpl_vector_get_size(xi);
    for (i=0; i<n; i++) {
        cpl_vector_set(xs, i, cpl_vector_get(li,i)-(i*D) );
    }

    x0 = cpl_vector_get_median(xs);

    cpl_vector_delete(li);
    cpl_vector_delete(xs);
    return x0;
}
/*----------------------------------------------------------------------------*/
/**
  @brief Find the true D from fringe statistics
  @param
  @return
 */
/*----------------------------------------------------------------------------*/

double cr2res_wave_etalon_get_D(
            cpl_vector      * li)
{
	int				i;
	cpl_size		nxi;
    double      	trueD=-1.0;
	cpl_vector	*	diffs;

	nxi = cpl_vector_get_size(li);
	diffs = cpl_vector_new(nxi-1);
	for (i=1; i<nxi; i++){
		cpl_vector_set(diffs,i-1,
			cpl_vector_get(li,i) - cpl_vector_get(li,i-1) );
	}

    if (cpl_msg_get_level() == CPL_MSG_DEBUG){
        cpl_table   *   tab;
        tab = cpl_table_new(nxi-1);
        cpl_table_new_column(tab, "wavediff", CPL_TYPE_DOUBLE) ;
        for(i=0; i<nxi-1; i++) {
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
    int             smooth = 35 ;   // TODO: make free parameter?
                                    // interfringe ~30 in Y, ~70 in K
    double          thresh = 1.0 ;   // TODO: derive from read-out noise
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
        cpl_vector_save(spectrum, "debug_spectrum.fits", CPL_TYPE_DOUBLE,
                NULL, CPL_IO_CREATE);
    }

    /*Output array, values are invalid until set.*/
    peaks = cpl_array_new(max_num_peaks, CPL_TYPE_DOUBLE);

    /* X-axis to cut out from for each peak */
    X_all = cpl_vector_new(nx);
    for (i=0; i<nx; i++) cpl_vector_set(X_all, i, (double)i+1) ;

    for (i=0; i < nx; i++){
        j = 0;
        while ( (spec_i = cpl_vector_get(spec_thresh, i)) > -1 ) {
            j++;
            i++;
        }
        if (j < min_len_peak) continue;
        // cpl_msg_debug(__func__, "Peak length j=%d at i=%d",j,i);
        cur_peak = cpl_vector_extract(spec_thresh, i-j, i, 1) ;
        X_peak = cpl_vector_extract(X_all, i-j, i, 1) ;

        if (cpl_vector_fit_gaussian(X_peak, NULL, cur_peak, NULL, CPL_FIT_ALL,
                                &x0, &sigma, &area, &offset,
                                NULL,NULL,NULL) != CPL_ERROR_NONE ) {
            cpl_msg_warning(__func__, "Fit at j=%d i=%d failed",j,i);
            cpl_vector_delete(cur_peak);
            cpl_vector_delete(X_peak);
            cpl_error_reset();
            continue;
        }

        //cpl_msg_debug(__func__,"Fit: %.2f, %.2f, %.2f, %.2f",
        //                            x0, sigma, area, offset);
        if ((k = cpl_array_count_invalid(peaks)) <1)
            cpl_msg_error(__func__,"Output array overflow!");
        //cpl_msg_debug(__func__,"k=%d, x0=%g",k,x0);
        cpl_array_set_double(peaks, max_num_peaks - k, x0);

        cpl_vector_delete(cur_peak);
        cpl_vector_delete(X_peak);
    }

    /* Copy into output array */
    k = max_num_peaks - cpl_array_count_invalid(peaks);
    peak_vec = cpl_vector_new(k) ;
    for (i=0; i<k; i++)
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
  @param    wavesol_init   The wavelength polynomial
  @param    wl_error        Max error in pixels of the initial guess
  @return   The lines spectrum
 */
/*----------------------------------------------------------------------------*/
cpl_bivector * cr2res_wave_gen_lines_spectrum(
        const char      *   catalog,
        cpl_polynomial  *   wavesol_init,
        int                 wl_error)
{
    cpl_bivector    *   lines ;
    cpl_bivector    *   lines_sub ;
    double          *   lines_sub_wl ;
    double          *   lines_sub_intens ;
    double              wl_error_wl, wl_min, wl_max ;
    int                 i ;

    /* Check Entries */
    if (catalog == NULL || wavesol_init == NULL) return NULL ;

    /* Load the lines */
    lines = cr2res_io_load_EMISSION_LINES(catalog) ;

    /* Extract the needed spectrum */
    wl_min = cpl_polynomial_eval_1d(wavesol_init, 1, NULL);
    wl_max = cpl_polynomial_eval_1d(wavesol_init, CR2RES_DETECTOR_SIZE, NULL);
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
        /* TODO */
        if (lines_sub_intens[i] > 5000) lines_sub_intens[i] = 0.0 ;

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
    wished_ext_nb = cr2res_io_get_ext_idx(filename, detector, 1) ;
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


/*----------------------------------------------------------------------------*/
/**
  @brief    Fit a 1D polynomial to data px, py, sigma_py
  @param    px              X data, i.e. pixel positions
  @param    py              Y data, i.e. wavelengths
  @param    sigma_py        Y error, i.e. uncertainty on each wavelength
  @param    degree          Degree of the polynomial fit
  @param    solution_init   Initial solution of the polynomial fit
  @param    wavelength_error [out] array of wave_mean_error, wave_max_error (may be NULL)
  @param    sigma_fit       [out] uncertainties of the polynomial fit
                            parameters (may be NULL)
  @param    cov             Covariance matrix of the polynomial fit (may be NULL if sigma_fit is NULL)

  @return  Wavelength solution, i.e. polynomial that translates pixel
            values to wavelength.
  The returned polynomial must be deallocated with cpl_polynomial_delete()
 */
/*----------------------------------------------------------------------------*/
cpl_polynomial * polyfit_1d(
        cpl_matrix * px, 
        cpl_vector * py,
        cpl_vector *sigma_py,
        int degree,
        const cpl_polynomial * solution_init,
        cpl_array ** wavelength_error,
        cpl_vector ** sigma_fit,
        cpl_matrix ** cov)
{

    // For polynomial fit
    cpl_vector * diff;
    cpl_polynomial * result = cpl_polynomial_new(1);
    cpl_matrix * sigma_px = NULL;
    cpl_vector * pa = cpl_vector_new(degree + 1 + 1);
    cpl_error_code error;
    cpl_size i, j;
    int *pia;
    
    // first parameter of polynomial fit is the number of degrees (the value is fixed though)
    // the number of parameters is then polynomial degree + 1 (constant term) + 1 (number of degrees)
    // i.e. pa = degrees, a0, a1, ...
    cpl_vector_set(pa, 0, degree);
    pia = cpl_malloc((degree + 1 + 1) * sizeof(int));
    pia[0] = 0;
    for (i = 1; i < degree + 1 + 1; i++){
        pia[i] = 1;
    }

    // initial guess for polynomial fit, based on passed initial guess
    for (j = 0; j < degree+1; j++){
        cpl_vector_set(pa, j+1, cpl_polynomial_get_coeff(solution_init, &j));
    }

    // polynomial fit of px, py
    // with px: line pixel, py: line wavelength
    // I would use cpl_polynomial_fit, but that does not support error estimation
    error = cpl_fit_lvmq(px, sigma_px, py, sigma_py, pa, pia, &poly, &deriv_poly,
                        CPL_FIT_LVMQ_TOLERANCE, CPL_FIT_LVMQ_COUNT,
                        CPL_FIT_LVMQ_MAXITER, NULL, NULL, cov);


    if (error == CPL_ERROR_NONE){
        // Everything is fine

        // if there are not enough data points left for the polynomial fit
        // the lvmq fit will complain about a singular matrix

        // errors of the fit are the square root of the diagonal of the covariance matrix
        // assuming parameters are uncorrelated, better use the whole matrix
        // errors on the wavelength are then given using standard error propagation
        // s_wl**2 = s0**2 + s1**2 * x**2 + s2**2 * x**4 + ... si**2 * x**(2*i)
        for (i = 0; i < degree + 1; i++){
            cpl_polynomial_set_coeff(result, &i, cpl_vector_get(pa, i+1));
            if (sigma_fit != NULL)
                cpl_vector_set(*sigma_fit, i, sqrt(cpl_matrix_get(*cov, i+1, i+1)));
        }

        if (wavelength_error != NULL){
            // Calculate absolute difference between polynomial and catalog value for each line
            // use px and py, so that only good lines are used
            diff = cpl_vector_new(cpl_vector_get_size(py));
            if (*wavelength_error == NULL)
                *wavelength_error = cpl_array_new(2, CPL_TYPE_DOUBLE);


            for (i = 0; i < cpl_vector_get_size(py); i++){
                cpl_vector_set(diff, i, fabs(
                    cpl_polynomial_eval_1d(result, cpl_matrix_get(px, i, 0), NULL)
                    - cpl_vector_get(py, i)));
            }
            // Set wavelength_error to mean and max difference
            cpl_array_set_double(*wavelength_error, 0, cpl_vector_get_mean(diff));
            cpl_array_set_double(*wavelength_error, 1, cpl_vector_get_max(diff));

            cpl_msg_debug(__func__, "Wave Error Mean: %g", cpl_array_get(*wavelength_error, 0, NULL));
            cpl_msg_debug(__func__, "Wave Error Max : %g", cpl_array_get(*wavelength_error, 1, NULL));

            cpl_vector_delete(diff);
        }
    }


    cpl_free(pia);
    cpl_vector_delete(pa);

    if (error != CPL_ERROR_NONE){
        cpl_msg_info(__func__, "Can't fit Polynomial");
        cpl_polynomial_delete(result);
        cpl_error_reset();
        return NULL;
    }
    return result;
}

// Functions for polynomial fit in _wave_catalog()
static int poly(const double x[], const double a[], double * result)
{
    int j;
    int degree = a[0];

    *result = 0;
    for (j = degree; j > 0; j--){
        *result = x[0] * (a[1 + j] + *result);
    }
    *result = *result + a[1];

    return 0;
}

static int deriv_poly(const double x[], const double a[], double * result)
{
    int i, j;
    const double degree = a[0];

    result[0] = 0;
    for (j = 0; j < degree + 1; j++){
        result[j+1] = 1;
        for (i = 0; i < j; i++){
            result[j+1] *= x[0];
        }
    }

    return 0;
}

// following two are shamelessly taken from cpl_vector_fit_gauss

/*----------------------------------------------------------------------------*/
/**
   @internal
   @brief   Evaluate a gaussian
   @param   x             The evaluation point
   @param   a             The parameters defining the gaussian
   @param   result        The function value

   @return  0 iff okay.

   This function computes

   @em a3 +  @em a2 *
   exp( -(@em x0 - @em a0)^2/(2 @em a1^2)).

   where @em a0, ..., @em a3 are the first four elements of @em a, and @em
   x0 is the first element of @em x .

   The function fails iff @em a1 is zero and @em x0 is equal to @em a0.

*/
/*----------------------------------------------------------------------------*/
static int gauss(const double x[], const double a[], double *result)
{
    const double my    = a[0];
    const double sigma = a[1];

    if (sigma != 0.0) {

        const double A = a[2];
        const double B = a[3];

        *result = B + A * exp(- (x[0] - my)*(x[0] - my)
                / (2*sigma*sigma));

    } else {

        /* Dirac's delta function */
        *result = x[0] != my ? 0.0 : DBL_MAX;
    }

    return 0;
}

/*----------------------------------------------------------------------------*/
/**
   @internal
   @brief   Evaluate the derivatives of a gaussian
   @param   x             The evaluation point
   @param   a             The parameters defining the gaussian
   @param   result        The derivatives wrt to parameters

   @return  0 iff okay.

   This function computes the partial derivatives of
   @em f(@em x0,@em a) =
   @em a3 +  @em a2  *
   exp( -(@em x0 - @em a0)^2/(2 @em a1^2))
   with respect to @em a0, ..., @em a3.
   On successful evaluation, the i'th element of the @em result vector
   contains df/da_i.

   The function never returns failure.

*/
/*----------------------------------------------------------------------------*/
static int gauss_derivative(const double x[], const double a[], double result[])
{

    if (a[1] != 0.0) {

        const double my    = a[0];
        const double sigma = a[1];
        const double A     = a[2];
        /* a[3] not used */

        /* f(x) = B + A/sqrt(2 pi s^2) exp(-(x-my)^2/2s^2)
         *
         * df/d(my) = A/sqrt(2 pi s^2) exp(-(x-my)^2/2s^2) * (x-my)  / s^2
         *          = A * fac. * (x-my)  / s^2
         * df/ds    = A/sqrt(2 pi s^2) exp(-(x-my)^2/2s^2) * ((x-my)^2/s^3 - 1/s)
         *          = A * fac. * ((x-my)^2 / s^2 - 1) / s
         * df/dA    = 1/sqrt(2 pi s^2) exp(-(x-my)^2/2s^2)
         *          = fac.
         * df/dB    = 1
         */


        const double factor = exp( -(x[0] - my)*(x[0] - my)/(2.0*sigma*sigma) );

        result[0] = A * factor * (x[0]-my) / (sigma*sigma);
        result[1] = A * factor * ((x[0]-my)*(x[0]-my) / (sigma*sigma*sigma));
        result[2] = factor;
        result[3] = 1.0;

    } else {
        /* Derivative of Dirac's delta function */
        result[0] = 0.0;
        result[1] = 0.0;
        result[2] = 0.0;
        result[3] = 0.0;
    }

    return 0;
}

/**@}*/
