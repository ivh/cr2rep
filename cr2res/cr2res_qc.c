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
#include "hdrl.h"

#include "cr2res_qc.h"
#include "cr2res_qc_lines.h"
#include "cr2res_trace.h"
#include "cr2res_dfs.h"
#include "cr2res_extract.h"
#include "cr2res_calib.h"
#include "cr2res_detlin.h"
#include "cr2res_wave.h"

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
  @param    min_level   [out] min level
  @param    max_level   [out] max level
  @return   The computed median or -1.0 in error case
 */
/*----------------------------------------------------------------------------*/
double cr2res_qc_detlin(
        const hdrl_imagelist    *   hdrl_coeffs,
        double                      bpm_thresh,
        cpl_mask               **   outmask,
        double                  *   min_level,
        double                  *   max_level) 
{
    double      qc_detlin_median, qc_detlin_min, qc_detlin_max ;
    double      level;
    int         nimgs = 3, width, height;
    hdrl_image * img;
    hdrl_value value = {CR2RES_NONLIN_LEVEL, 0};
    cpl_mask * tmpmask;


    /* Check Entries */
    if (hdrl_coeffs==NULL || min_level==NULL || max_level==NULL) return -1.0 ;

    /* Initialise */
    qc_detlin_median = -1.0 ;
    qc_detlin_min = -1.0 ;
    qc_detlin_max = -1.0 ;

    // Apply detlin correction on an image with constant value
    width = hdrl_image_get_size_x(hdrl_imagelist_get(hdrl_coeffs, 0));
    height = hdrl_image_get_size_y(hdrl_imagelist_get(hdrl_coeffs, 0));
    img = hdrl_image_new(width, height); 
    hdrl_image_add_scalar(img, value);
    cr2res_detlin_correct(img, hdrl_coeffs);
    // mask the Nan values, that exist in the out of order pixels
    // as well as bad pixels and wherever detlin failed
    hdrl_image_reject_value(img, CPL_VALUE_NAN);

    if (cpl_msg_get_level() == CPL_MSG_DEBUG){ 
        cpl_image_save(
            hdrl_image_get_image(img), "debug_qccorrection.fits",
        CPL_TYPE_DOUBLE, NULL, CPL_IO_CREATE);
    }

    tmpmask = cpl_mask_new(width, height);
    cpl_mask_threshold_image(tmpmask, hdrl_image_get_image(img),
                (double)CR2RES_NONLIN_LEVEL,
                (double)CR2RES_NONLIN_LEVEL * (1.0 + (bpm_thresh/100)),
                CPL_BINARY_0);
    hdrl_image_reject_from_mask(img, tmpmask);
    cpl_mask_delete(tmpmask);

    // Then determine the median of that corrected image
    qc_detlin_median = cpl_image_get_median(hdrl_image_get_image(img));
    qc_detlin_min = cpl_image_get_min(hdrl_image_get_image(img)) ;
    qc_detlin_max = cpl_image_get_max(hdrl_image_get_image(img)) ;

    tmpmask = cpl_mask_duplicate(hdrl_image_get_mask(img));
    cpl_msg_info(__func__, "BPM has %d pix rejected",
                                (int)cpl_mask_count(tmpmask));
    hdrl_image_delete(img);

    if (isnan(qc_detlin_median) || 
            isnan(qc_detlin_min) ||
            isnan(qc_detlin_max)) return -1.0 ;

    *min_level = qc_detlin_min ;
    *max_level = qc_detlin_max ;
    *outmask = tmpmask;
    return qc_detlin_median ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Compute the central orders positions
  @param    tw          trave Wave table
  @param    order_nb    [out] Array of orders numbers
  @param    order_pos   [out] Array of orders positions
  @param    nbvals      [out] Size of the arrays
  @return   0 if ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
int cr2res_qc_flat_order_positions(
        const cpl_table *   tw,
        int             **  order_nb,
        double          **  order_pos,
        int             *   nbvals)
{
    int             *   order_idx_vals ;
    const cpl_array *   slit_frac ;
    int             *   order_nb_loc ;
    double          *   order_pos_loc ;
    int                 i, j, cur_order, nb_order_idx_vals, n_full ;
    cpl_size            nrows ;

    /* Check Entries */
    if (tw == NULL || order_nb==NULL || order_pos==NULL || nbvals==NULL) 
        return -1 ;

    /* Initialise */
    *nbvals = 0 ;
    *order_nb = NULL ;
    *order_pos = NULL ;
    nrows = cpl_table_get_nrow(tw) ;
    if (nrows <= 0) return 0 ;
    
    /* Get the list of different orders */
    order_idx_vals = cr2res_trace_get_order_idx_values(tw, &nb_order_idx_vals) ;

    /* Count the number of full slit orders in tw */
    n_full = 0 ;
    /* Loop on the different orders */
    for (i=0 ; i<nb_order_idx_vals ; i++) {
        cur_order = order_idx_vals[i] ;
        /* Search an open slit of this order */
        for (j=0 ; j<nrows ; j++) {
            if (cpl_table_get(tw, CR2RES_COL_ORDER, j, NULL) == cur_order) {
                slit_frac = cpl_table_get_array(tw, CR2RES_COL_SLIT_FRACTION,j);
                if (cr2res_trace_slit_fraction_info(slit_frac, NULL) == 
                        CR2RES_DECKER_NONE) {
                    n_full++ ;
                    /* Go to next order */
                    break ;
                }
            }
        }
    }

    /* Allocate output arrays */
    order_nb_loc = cpl_malloc(n_full * sizeof(int)) ;
    order_pos_loc = cpl_malloc(n_full * sizeof(double)) ;

    /* Fill the arrays */
    n_full = 0 ;

    /* Loop on the different orders */
    for (i=0 ; i<nb_order_idx_vals ; i++) {
        cur_order = order_idx_vals[i] ;

        /* Search an open slit of this order */
        for (j=0 ; j<nrows ; j++) {
            if (cpl_table_get(tw, CR2RES_COL_ORDER, j, NULL) == cur_order) {
                slit_frac = cpl_table_get_array(tw, CR2RES_COL_SLIT_FRACTION,j);
                if (cr2res_trace_slit_fraction_info(slit_frac, NULL) == 
                        CR2RES_DECKER_NONE) {
                    order_nb_loc[n_full] = cur_order ;
                    order_pos_loc[n_full] = cr2res_trace_get_trace_ypos(tw, j);
                    cpl_msg_debug(__func__, "Order %d Pos : %g",
                            order_nb_loc[n_full], order_pos_loc[n_full]) ;
                    n_full++ ;
                    /* Go to next order */
                    break ;
                }
            }
        }
    }
    cpl_free(order_idx_vals); 

    /* Return */
    *nbvals = n_full ;
    *order_nb = order_nb_loc ;
    *order_pos = order_pos_loc ;
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
    cpl_vector * vector;
    cpl_array * array;
    int * order_idx_values, nb_order_idx_values, central_order_idx, i;
    int * traces, nb_traces;
    double      qc_trace_center_y ;

    /* Check Entries */
    if (trace == NULL) return -1.0 ;

    /* Initialise */
    qc_trace_center_y = 0;
    // Step 1: find central order
    order_idx_values = cr2res_trace_get_order_idx_values((cpl_table*) trace, 
            &nb_order_idx_values);
    array = cpl_array_wrap_int(order_idx_values, nb_order_idx_values);

/* TODO : Is the median really the CENTRAL order ?? */
    central_order_idx = cpl_array_get_median(array);
    cpl_array_unwrap(array);

    // Step 2: Sum all traces together
    traces = cr2res_get_trace_numbers(trace, central_order_idx, &nb_traces);
    for(cpl_size i = 0; i < nb_traces; i++) {
      vector = cr2res_trace_get_ycen(trace, central_order_idx, traces[i], 
                CR2RES_DETECTOR_SIZE);
      qc_trace_center_y += cpl_vector_get_mean(vector);
      cpl_vector_delete(vector);
    }
    
    // Step 3: take the mean
    qc_trace_center_y /= nb_traces;

    cpl_free(order_idx_values);
    cpl_free(traces);

    return qc_trace_center_y ;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the S2N on the flat
  @param    master_flat     The master flat
  @return   
 */
/*----------------------------------------------------------------------------*/
double cr2res_qc_flat_s2n(
        const cpl_image     *   master_flat)
{
    return -1.0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the central WLEN of a given order
  @param    tw      the TW table
  @param    order   the order index
  @return   the computed number or -1.0 in error case
 */
/*----------------------------------------------------------------------------*/
double cr2res_qc_wave_central(
        const cpl_table *   tw,
        int                 order_idx)
{
    cpl_vector  *   wls ;
    double          wl_central ;

    /* Get the wavelengths */
    if ((wls = cr2res_trace_get_wl(tw, order_idx, CR2RES_QC_TRACE,
                    CR2RES_DETECTOR_SIZE)) == NULL) {
        cpl_msg_warning(__func__, 
                "QC CENTWL - cannot find wl solution for order/trace %d/%d",
                order_idx, CR2RES_QC_TRACE) ;
        return -1.0 ;
    }
    wl_central = cpl_vector_get(wls, (int)(CR2RES_DETECTOR_SIZE/2)) ;
    cpl_vector_delete(wls) ;
    return wl_central ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the dispersion of a given order
  @param    tw      the TW table
  @param    order   the order index
  @return   the computed number or -1.0 in error case
 */
/*----------------------------------------------------------------------------*/
double cr2res_qc_wave_disp(
        const cpl_table *   tw,
        int                 order_idx)
{
    cpl_vector  *   wls ;
    double          wl_disp ;
    int             nbins ;

    if (tw == NULL) return -1.0 ;

    /* Get the wavelengths */
    if ((wls = cr2res_trace_get_wl(tw, order_idx, CR2RES_QC_TRACE,
                    CR2RES_DETECTOR_SIZE)) == NULL) {
        cpl_msg_warning(__func__, 
                "QC DISPWL - cannot find wl solution for order/trace %d/%d",
                order_idx, CR2RES_QC_TRACE) ;
        return -1.0 ;
    }
    nbins = cpl_vector_get_size(wls) ;
    wl_disp = (cpl_vector_get(wls, nbins-1) - cpl_vector_get(wls, 0)) / nbins;
    cpl_vector_delete(wls) ;
    return wl_disp ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Collect lines between 2 wavelengths
  @param    wstart
  @param    wstop
  @return   the vector with the collected lines
 */
/*----------------------------------------------------------------------------*/
cpl_vector * cr2res_qc_lines_collect(double wmin, double wmax)
{
    cpl_vector  *   lines ;
    int             i, nb_found ;

    /* Initialise */
    nb_found = 0 ;

    /* Find the Interesting lines */ 
    for (i=0 ; i< QC_LINES_NUMBER ; i++) 
        if (qc_lines[i] > wmin && qc_lines[i] < wmax) nb_found++ ;

    if (nb_found == 0) return NULL ;

    /* Allocate the container */
    lines = cpl_vector_new(nb_found); 

    /* Fill with the good lines */
    nb_found = 0 ;
    for (i=0 ; i< QC_LINES_NUMBER ; i++) 
        if (qc_lines[i] > wmin && qc_lines[i] < wmax) 
            cpl_vector_set(lines, nb_found++, qc_lines[i]) ;
    return lines; 
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes one line Intensity
  @param    spec        spectrum
  @param    wl          line position
  @return   the computed intensity
 */
/*----------------------------------------------------------------------------*/
double cr2res_qc_wave_line_intens(
        const cpl_bivector  *   spec,
        double                  wl)
{
    // TODO Thomas / Ansgar
    //cpl_plot_bivector("set grid;set xlabel 'Wavelength (nm)';
    // set ylabel 'Spec';", "t 'Spectrum' w lines", "",spec) ;
    const cpl_vector * wave;
    const cpl_vector * flux;
    cpl_vector * tmp;
    cpl_size pixel_pos;
    cpl_size window_size;
    cpl_size k, n_inner, n_outer;
    double sum_inner, sum_outer;
    double value;

    if (spec == NULL) return -1.0;

    // TODO: How large should this window be?
    window_size = CR2RES_QC_WINDOW / 2;

    wave = cpl_bivector_get_x_const(spec);
    flux = cpl_bivector_get_y_const(spec);

    // Determine pixel pos
    // TODO: If we are sure that wave is sorted, we can also use cpl_vector_find
    tmp = cpl_vector_duplicate(wave);
    cpl_vector_subtract_scalar(tmp, wl);
    cpl_vector_multiply(tmp, tmp);
    pixel_pos = cpl_vector_get_minpos(tmp);
    cpl_vector_delete(tmp);

    // If the wavelength value is outside the spectrum
    if (pixel_pos == 0 || pixel_pos == cpl_vector_get_size(wave)) return -1.0;

    // Sum up the values of the spectrum
    // inside the window and outside the window
    sum_inner = 0;
    sum_outer = 0;
    n_inner = 0;
    n_outer = 0;
    for (cpl_size i = -window_size * 2; i < 2 * window_size; i++)
    {
        k = pixel_pos - i;
        if (k < 0 | k > cpl_vector_get_size(flux)){
            continue;
        }
        value = cpl_vector_get(flux, k);
        if (isnan(value)){
            continue;
        }
        if (fabs(i) < window_size){
            // Inner sum
            n_inner++;
            sum_inner += value;
        } else {
            // Outer sum
            n_outer++;
            sum_outer += value;
        }
    }
    if (n_inner + n_outer == 0){
        // No valid points
        return -1.0;
    }
    // Take the mean
    if (n_inner != 0) sum_inner /= n_inner;
    if (n_outer != 0) sum_outer /= n_outer;
    // return the difference
    value = sum_inner - sum_outer;

    return value;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes one line Fwhm 
  @param    spec        spectrum
  @param    wl          line position
  @param    peak_height [out] Fitted Peak height
  @return   the computed fwhm
 */
/*----------------------------------------------------------------------------*/
double cr2res_qc_wave_line_fwhm(
        const cpl_bivector  *   spec,
        double                  wl,
        double              *   peak_height)
{
    // TODO Thomas / Ansgar
    //cpl_plot_bivector("set grid;set xlabel 'Wavelength (nm)';
    // set ylabel 'Spec';", "t 'Spectrum' w lines", "",spec) ;
    const cpl_vector * wave;
    const cpl_vector * flux;
    cpl_vector * unc;
    cpl_vector * result;
    cpl_vector * tmp;
    cpl_size pixel_pos;
    cpl_size window_width;
    double fwhm ;

    /* Check Entries */
    if (spec == NULL || peak_height == NULL) return -1.0;

    /* Initialise */
    *peak_height = -1.0 ;

    wave = cpl_bivector_get_x_const(spec);
    flux = cpl_bivector_get_y_const(spec);
    unc = NULL;
    result = NULL;
    // TODO: what should the size of the window be?
    window_width = CR2RES_QC_WINDOW;

    // Determine pixel pos
    tmp = cpl_vector_duplicate(wave);
    cpl_vector_subtract_scalar(tmp, wl);
    cpl_vector_multiply(tmp, tmp);
    pixel_pos = cpl_vector_get_minpos(tmp);
    cpl_vector_delete(tmp);

    // If the wavelength value is outside the spectrum
    if (pixel_pos == 0 || pixel_pos == cpl_vector_get_size(wave)) return -1.0;

    // Fit the line with a gaussian
    if (cr2res_wave_fit_single_line(flux, unc, pixel_pos, window_width, 1, 0, 
                &result)){
        // Could not determine the line fit
        cpl_error_reset();
        return -1.0 ;
    }

    // Get the FWHM from the results
    // and multiply with the conversion factor from sigma to FWHM
    fwhm = cpl_vector_get(result, 1);
    fwhm = fwhm * 2.355;

    *peak_height = cpl_vector_get(result, 2) ;
    cpl_vector_delete(result);
    return fwhm ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the lines intensities
  @param    extracted   extracted spectrum
  @return   the computed resolution
 */
/*----------------------------------------------------------------------------*/
double cr2res_qc_wave_lamp_effic(
        const cpl_bivector  *   spec)
{
    cpl_vector  *   ref_lines ;
    cpl_vector  *   ref_lines_intens ;
    double          wmin, wmax, intens_med, intens;
    int             i, n, nall ;

    /* Get the reference lines */
    wmin = cpl_vector_get(cpl_bivector_get_x_const(spec), 0);
    wmax = cpl_vector_get(cpl_bivector_get_x_const(spec),
            cpl_bivector_get_size(spec)-1) ;
    ref_lines = cr2res_qc_lines_collect(wmin, wmax) ;
    if (ref_lines == NULL) return -1.0 ;

    //cpl_vector_dump(ref_lines, stdout) ;

    /* Loop on the lines */
    n = 0;
    nall = cpl_vector_get_size(ref_lines);
    ref_lines_intens = cpl_vector_new(nall) ;
    for (i=0 ; i < nall ; i++) {
        intens = cr2res_qc_wave_line_intens(spec, cpl_vector_get(ref_lines, i));
        if (intens > 0.0) {
            cpl_vector_set(ref_lines_intens, n, intens); 
            n++;
        }
    }
    cpl_vector_delete(ref_lines) ;
    if (n == 0){
        cpl_vector_delete(ref_lines_intens) ;
        return -1.0 ;
    }
    cpl_vector_set_size(ref_lines_intens, n);
    cpl_msg_info(__func__, 
        "Using %i of %i lines to estimate Intensity", n, nall);


    /* Compute the median */
    intens_med = cpl_vector_get_median(ref_lines_intens) ;
    cpl_vector_delete(ref_lines_intens) ;

    return intens_med ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the lines Fwhm and return the smallest
  @param    extracted   extracted spectrum
  @param    wl          [out] the wavelength of the thinest line
  @return   the thinest fwhm
 */
/*----------------------------------------------------------------------------*/
double cr2res_qc_wave_resol_fwhm(
        const cpl_bivector  *   spec,
        double              *   wl)
{
    cpl_vector  *   ref_lines ;
    cpl_vector  *   ref_lines_fwhm ;
    cpl_vector  *   ref_lines_pos ;
    double          wmin, wmax, fwhm_med, fwhm, peak_height, min_fwhm_val, 
                    min_fwhm_pos;
    cpl_size        idx ;
    int             i, n, nall;

    /* Check Entries */
    if (spec == NULL || wl == NULL) return -1.0 ;

    /* Initialise */
    *wl = -1.0 ;

    /* Get the reference lines */
    wmin = cpl_vector_get(cpl_bivector_get_x_const(spec), 0);
    wmax = cpl_vector_get(cpl_bivector_get_x_const(spec),
            cpl_bivector_get_size(spec)-1) ;
    ref_lines = cr2res_qc_lines_collect(wmin, wmax) ;
    if (ref_lines == NULL) return -1.0 ;

    //cpl_vector_dump(ref_lines, stdout) ;

    /* Loop on the lines */
    n = 0;
    nall = cpl_vector_get_size(ref_lines);
    ref_lines_fwhm = cpl_vector_new(nall) ;
    ref_lines_pos = cpl_vector_new(nall) ;
    for (i=0 ; i<cpl_vector_get_size(ref_lines) ; i++) {
        fwhm = cr2res_qc_wave_line_fwhm(spec, cpl_vector_get(ref_lines, i), 
                &peak_height);
        if (fwhm > 0.1 && fwhm < CR2RES_QC_WINDOW && peak_height > 0 
                && peak_height<CR2RES_DETECTOR_OVEREXP_THRESH) {
            cpl_vector_set(ref_lines_fwhm, n, fwhm); 
            cpl_vector_set(ref_lines_pos, n, cpl_vector_get(ref_lines, i)); 
            n++;
        }
    }
    cpl_vector_delete(ref_lines) ;
    if (n == 0){
        cpl_vector_delete(ref_lines_fwhm) ;
        cpl_vector_delete(ref_lines_pos) ;
        return -1.0;
    }
    cpl_vector_set_size(ref_lines_fwhm, n);
    cpl_vector_set_size(ref_lines_pos, n);

    /* Return the smallest FWHM and its position */
    idx = cpl_vector_get_minpos(ref_lines_fwhm) ;
    if (idx < 0) {
        cpl_vector_delete(ref_lines_fwhm) ;
        cpl_vector_delete(ref_lines_pos) ;
        return -1.0 ;
    } 
    min_fwhm_val = cpl_vector_get(ref_lines_fwhm, idx) ;
    min_fwhm_pos = cpl_vector_get(ref_lines_pos, idx) ;
    cpl_vector_delete(ref_lines_fwhm) ;
    cpl_vector_delete(ref_lines_pos) ;

    cpl_msg_info(__func__, "Using the smallest FWHM %g of line at %g nm",
            min_fwhm_val, min_fwhm_pos) ;

    /* Return */
    *wl = min_fwhm_pos ;
    return min_fwhm_val ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the Overexposed fraction
  @param    ima         the image (with BPM)
  @param    tw          the TW table
  @param    order_idx   the order index
  @return   the computed overexposed value
 */
/*----------------------------------------------------------------------------*/
double cr2res_qc_overexposed(
        const cpl_image *   ima,
        const cpl_table *   tw,
        int                 order_idx)
{
    cpl_image       *   labels ;
    const double    *   pima ;
    int             *   plabel ;
    int                 i, j, nb_over, nb_total ;

    /* Check Entries */
    if (ima == NULL || tw == NULL) return -1.0 ;
    if (cpl_image_get_size_x(ima) != CR2RES_DETECTOR_SIZE ||
            cpl_image_get_size_y(ima) != CR2RES_DETECTOR_SIZE) {
        return -1.0 ;
    }
    pima = cpl_image_get_data_double_const(ima) ;

    /* Create Label image */
    labels = cr2res_trace_gen_image((cpl_table *)tw, CR2RES_DETECTOR_SIZE, 
            CR2RES_DETECTOR_SIZE);
    plabel = cpl_image_get_data_int(labels) ;

    nb_over = nb_total = 0 ;
    /* Loop on th pixels */
    for (j=0 ; j<CR2RES_DETECTOR_SIZE ; j++) {
        for (i=0 ; i<CR2RES_DETECTOR_SIZE ; i++) {
            if (plabel[i+j*CR2RES_DETECTOR_SIZE] == order_idx &&
                    !cpl_image_is_rejected(ima, i+1, j+1)) {
                nb_total++ ;
                if (pima[i+j*CR2RES_DETECTOR_SIZE] > 
                        CR2RES_DETECTOR_OVEREXP_THRESH) {
                    nb_over ++ ;
                }
            }
        }
    }
    cpl_image_delete(labels) ;
    if (nb_total == 0) return -1.0 ;
    return (double)nb_over/(double)nb_total ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the integrated flux over part of the spectrum
  @param    extracted   Extracted spectrum table
  @return   the computed signal
 */
/*----------------------------------------------------------------------------*/
double cr2res_qc_obs_nodding_signal(
        const cpl_table     *   extracted)
{
    double          qc_signal ;
    char        *   colname;
    double      *   data;
    cpl_vector  *   vector, *vector2;
    int             nrows, size;

    /* Check Entries */
    if (extracted == NULL) return -1.0 ;

    /* Initialise */
    qc_signal = -1.0 ;
    colname = cr2res_dfs_SPEC_colname(CR2RES_QC_ORDER, CR2RES_QC_TRACE);
    data = cpl_table_get_data_double((cpl_table*) extracted, colname);
    cpl_free(colname);

    if (cpl_error_get_code() == CPL_ERROR_DATA_NOT_FOUND) {
      // QC order or trace not found in the table
      cpl_error_reset();
      return -1;
    }
    
    nrows = cpl_table_get_nrow(extracted);
    size = CR2RES_QC_SIZE;
    if (nrows < size) size = nrows;

    vector = cpl_vector_wrap(nrows, data);
    vector2 = cpl_vector_extract(vector, (nrows - size)/ 2, (nrows + size) 
                / 2 - 1, 1);
    qc_signal = cpl_vector_get_median(vector2);

    cpl_vector_unwrap(vector);
    cpl_vector_delete(vector2);

    return qc_signal ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the standard flux over part of the spectrum
  @param    extracted   Extracted spectrum table
  @param    setting     The setting
  @return   the computed standard flux
 */
/*----------------------------------------------------------------------------*/
double cr2res_qc_obs_nodding_standard_flux(
        const cpl_table     *   extracted,
        char                *   setting)
{
    double          qc_flux ;
    char        *   colname;
    double      *   ext_data;
    double      *   wl_data;
    double          wl_start, wl_stop ;
    cpl_vector  *   vector, *vector2;
    int             nrows, wl_start_index, wl_stop_index, i ;

    /* Check Entries */
    if (extracted == NULL) return -1.0 ;

    /* Initialise */
    qc_flux = -1.0 ;
    wl_start = -1.0 ;
    wl_stop = -1 ;
    nrows = cpl_table_get_nrow(extracted);

    /* Get the wl boundaries */
    if (!strcmp(setting, "Y1029")) {
        wl_start = 1087.3 ;
        wl_stop = -1 ;
    } else if (!strcmp(setting, "H1559")) {
        wl_start = 1504.8 ;
        wl_stop = -1 ;
    } else if (!strcmp(setting, "K2217")) {
        wl_start = 2205.0 ;
        wl_stop = -1 ;
    } else if (!strcmp(setting, "L3377")) {
        wl_start = 3828.70 ;
        wl_stop = 3829.30 ;
    } else if (!strcmp(setting, "M4266")) {
        wl_start = 3953.00 ;
        wl_stop = 3953.90 ;
    } else if (!strcmp(setting, "Y1028")) {
        wl_start = 1045.70 ;
        wl_stop =  1047.00 ;
    } else if (!strcmp(setting, "J1228")) {
        wl_start = 1245.50 ;
        wl_stop =  1246.40 ;
    } else if (!strcmp(setting, "H1567")) {
        wl_start = 1695.20 ;
        wl_stop =  1696.20 ;
    } else if (!strcmp(setting, "K2148")) {
        wl_start = 2312.60 ;
        wl_stop = 2313.30 ;
    }

    if (wl_start < 0.0) {
        cpl_msg_info(__func__, 
                "QC Standard Flux : No WL specified for setting %s", setting) ;
        return -1.0 ;
    }

    /* Get the Data */
    colname = cr2res_dfs_SPEC_colname(CR2RES_QC_ORDER, CR2RES_QC_TRACE);
    ext_data = cpl_table_get_data_double((cpl_table*) extracted, colname);
    cpl_free(colname);
    colname = cr2res_dfs_WAVELENGTH_colname(CR2RES_QC_ORDER, CR2RES_QC_TRACE);
    wl_data = cpl_table_get_data_double((cpl_table*) extracted, colname);
    cpl_free(colname);
    if (cpl_error_get_code() == CPL_ERROR_DATA_NOT_FOUND) {
        cpl_error_reset();
        return -1;
    }

    wl_start_index = -1 ;
    for (i=0 ; i<nrows ; i++) {
        if (wl_data[i] > wl_start && wl_start_index < 0) {
            wl_start_index = i ;
            break ;
        }
    }
    wl_stop_index = -1;
    if (wl_stop > 0.0) {
        for (i=0 ; i<nrows ; i++) {
            if (wl_data[i] <  wl_stop && wl_stop_index < 0) {
                wl_stop_index = i ;
                break ;
            }
        }
    } else {
        wl_stop_index = wl_start_index + CR2RES_QC_SIZE ;
    }

    if (wl_start_index < 0 || wl_stop_index >= cpl_table_get_nrow(extracted)) {
        cpl_msg_warning(__func__, 
                "QC Standard Flux : Specified WL out of range %g-%g",
                wl_start, wl_stop) ;
        return -1.0 ;
    }
    cpl_msg_info(__func__, 
            "QC STD FLUX : Trace %d Order %d WL start %g, idx range: [%d-%d]", 
            CR2RES_QC_TRACE, CR2RES_QC_ORDER, wl_start, wl_start_index,
            wl_stop_index) ;

    vector = cpl_vector_wrap(nrows, ext_data);
    vector2 = cpl_vector_extract(vector, wl_start_index, wl_stop_index, 1) ;
    qc_flux = cpl_vector_get_mean(vector2);
    cpl_vector_unwrap(vector);
    cpl_vector_delete(vector2);

    return qc_flux ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the SNR of several spectra
  @param    tw          The TW table
  @param    extracted   The Extracted table
  @param    out_order_idx_values    [out] The order values
  @param    out_nb_order_idx_values [out] The number of order values
  @return   An array of size out_nb_order_idx_values with the SNRs values
 */
/*----------------------------------------------------------------------------*/
double * cr2res_qc_snr(
        const cpl_table *   tw,
        const cpl_table *   extracted,
        int             **  out_order_idx_values,
        int             *   out_nb_order_idx_values)
{
    int             *   order_idx_values ;
    int                 nb_order_idx_values ;
    cpl_bivector    *   my_spec,
                    *   my_spec_err ;
    double          *   pmy_spec_err ;
    cpl_vector      *   my_snr_spec ;
    double          *   snrs ;
    int                 i, j ;

    /* Check Entries */
    if (tw==NULL || extracted==NULL || out_order_idx_values==NULL ||
            out_nb_order_idx_values==NULL) return NULL ;

    /* Get the number of orders from TW table */
    order_idx_values = cr2res_trace_get_order_idx_values(tw,
            &nb_order_idx_values) ;

    /* Allocate the output arrays */
    snrs = cpl_malloc(nb_order_idx_values * sizeof(double)) ;

    /* Loop on the orders */
    for (i=0 ; i<nb_order_idx_values ; i++) {
        /* Compute the SNR */
        if (cr2res_extract_EXTRACT1D_get_spectrum(extracted,
                order_idx_values[i], 1, &my_spec, &my_spec_err) == 0) {
            my_snr_spec = cpl_vector_duplicate(cpl_bivector_get_y(my_spec)) ;
            cpl_bivector_delete(my_spec) ;
            /* Clean the error to avoid division by 0.0 */
            pmy_spec_err = cpl_bivector_get_y_data(my_spec_err) ;
            for (j=0 ; j<cpl_bivector_get_size(my_spec_err) ; j++) {
                if (pmy_spec_err[j] == 0) pmy_spec_err[j] = 1.0 ;
            }
            cpl_vector_divide(my_snr_spec, cpl_bivector_get_y(my_spec_err)) ;
            cpl_bivector_delete(my_spec_err) ;
            snrs[i] = cpl_vector_get_median(my_snr_spec) ;
            if (isnan(snrs[i])) snrs[i] = -1.0 ;
            cpl_vector_delete(my_snr_spec) ; 
        } else {
            snrs[i] = -1.0 ;
        }
    }

    *out_order_idx_values = order_idx_values ; 
    *out_nb_order_idx_values = nb_order_idx_values ;
    return snrs ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the FWHM of the PSF along the slit for a given order
  @param    slitfu      Slit func table
  @param    order_idxp  Order index
  @return   The computed FWHM for this order
 */
/*----------------------------------------------------------------------------*/
double cr2res_qc_obs_nodding_slit_psf(
        const cpl_table     *   slitfu,
        int                     order_idxp)
{
    cpl_vector      *   x ;
    cpl_vector      *   y ;
    char            *   colname ;
    const double    *   data ;
    cpl_size            i, j ;
    cpl_fit_mode        fit_pars ;
    cpl_error_code      err;
    int                 nrow ;
    double              qc_fwhm, x0, sigma, area, offset;

    /* Check Entries */
    if (slitfu == NULL) return -1 ;

    /* Get the table column name */
    colname = cr2res_dfs_SLIT_FUNC_colname(order_idxp, 1);
    data = cpl_table_get_data_double_const(slitfu, colname);
    cpl_free(colname);
    if (data == NULL) {
        cpl_msg_warning(__func__, "No slitfunc data: code %d",
                    cpl_error_get_code());
        cpl_error_reset();
        return -1;
    }

    /* Initialise */
    qc_fwhm = -1.0 ;
    fit_pars = CPL_FIT_CENTROID + CPL_FIT_STDEV + CPL_FIT_AREA;
    offset = 0.;

    nrow = cpl_table_get_nrow(slitfu);

    /* x is just the element number */
    x = cpl_vector_new(nrow);
    y = cpl_vector_new(nrow);
    for (i = 0; i < nrow; i++){
        cpl_vector_set(x, i, i);
        cpl_vector_set(y, i, 0);
    }

    /* remove "strange" values, i.e. nan and unreasonably large values */
    /* otherwise the fit will not work */
    /* Slitfunc should be normalised to the oversampling rate(?) */
    for (j = 0; j < nrow; j++) {
        if (isnan(data[j]) | (data[j] > 1)){
            cpl_vector_set(y, j, 0);
        } else {
            cpl_vector_set(y, j, data[j]);
        }
    }
    err = cpl_vector_fit_gaussian(x, NULL, y, NULL, fit_pars, &x0, &sigma,
            &area, &offset, NULL, NULL, NULL);
    if (cpl_error_get_code()) {
        cpl_msg_warning(__func__, "Failed Fit for the slit PSF") ;
        cpl_error_reset() ;
        sigma = 0.0 ;
    }
    qc_fwhm = 2.355 * sigma; // 2.355 = 2 * sqrt(2 * ln(2))
 
    /* Free Memory */
    cpl_vector_delete(x);
    cpl_vector_delete(y);
    return qc_fwhm ;
}
/**@}*/

