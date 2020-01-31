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

#include "cr2res_slit_curv.h"
#include "cr2res_dfs.h"
#include "cr2res_io.h"
#include "cr2res_pfits.h"
#include "cr2res_utils.h"
#include "cr2res_trace.h"
#include "cr2res_extract.h"
#include "cr2res_etalon.h"

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static cpl_polynomial * cr2res_slit_curv_build_poly(
        cpl_polynomial  *   slit_poly_a,
        cpl_polynomial  *   slit_poly_b,
        cpl_polynomial  *   slit_poly_c,
        cpl_size            x) ;

static int cr2res_slit_curv_get_position(
        cpl_polynomial  *   trace,
        cpl_polynomial  *   wave,
        double              ref_wl,
        double          *   xpos,
        double          *   ypos) ;

static int fmodel(const double x[], const double a[], double *result) CPL_ATTR_NONNULL;
static int dmodel_da(const double x[], const double a[], double *result) CPL_ATTR_NONNULL;
static int cr2res_slit_curv_remove_peaks_at_edge(
    cpl_vector ** peaks,
    int width,
    int ncols) ;
static int cr2res_slit_curv_single_peak(
    cpl_image  * img_peak,
    cpl_vector * ycen,
    double       peak,
    cpl_matrix * x,
    cpl_vector * y,
    cpl_vector * a,
    int        * ia,
    double     * value_b,
    double     * value_c) ;
static int cr2res_slit_curv_all_peaks(
    cpl_vector * peaks,
    cpl_image  * img_rect,
    cpl_vector * ycen,
    int window,
    int fit_second_order,
    cpl_vector ** vec_b,
    cpl_vector ** vec_c) ;

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_slit_curv   Slit Curvature
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief Get the slit curvature directly from the image
  @param img   The input image
  @param trace_wave  The trace wave, with the order tracing
  @param order  The order to determine the curvature for
  @param trace  The trace to determine the curvature for
  @return   0 if ok, -1 in error case

  To determine the curvature we perform several steps. First, find all peaks in
  the spectrum. Next for each peak, fit a Gaussian with a shifted peak to the
  image. The shift of the peak determines the curvature.
  Finally fit the curvature to all 

 */
/*----------------------------------------------------------------------------*/
int cr2res_slit_curv_compute_order_trace(
        const hdrl_image    *   img,
        const cpl_table     *   trace_wave,
        int                     order,
        int                     trace,
        int                     height,
        int                     window,
        cpl_size                degree,
        int                     fit_second_order,
        cpl_polynomial      **  slit_poly_a,
        cpl_polynomial      **  slit_poly_b,
        cpl_polynomial      **  slit_poly_c)
{
    const cpl_image     *   img_in;
    cpl_vector          *   sfunc;
    cpl_vector          *   ycen;
    cpl_bivector        *   spec_bi;
    hdrl_image          *   model;
    cpl_vector          *   peaks;
    cpl_image           *   img_rect;
    cpl_vector          *   vec_b;
    cpl_vector          *   vec_c;
    int                     ncols;
    cpl_size                power;
    cpl_matrix          *   samppos;
    cpl_mask            *   kernel;
    cpl_image           *   img_other;
    hdrl_image          *   hdrl_other;
    double                  mad, median;

    if (img == NULL || trace_wave == NULL || slit_poly_a == NULL ||
        slit_poly_b == NULL || slit_poly_c == NULL) return -1;

    img_in = hdrl_image_get_image_const(img);
    ncols = cpl_image_get_size_x(img_in);

    // Median filter the image
    img_other = cpl_image_new(cpl_image_get_size_x(img_in),
        cpl_image_get_size_y(img_in), cpl_image_get_type(img_in));
    kernel = cpl_mask_new(3, 3);
    cpl_mask_not(kernel);
    cpl_image_filter_mask(img_other, img_in, kernel,
            CPL_FILTER_MEDIAN, CPL_BORDER_ZERO);
    cpl_mask_delete(kernel);

    // Discard outliers
    // cpl_image_subtract(img_other, img_in);
    // cpl_image_abs(img_other);
    // cpl_image_get_mad(img_in, &mad);
    // kernel = cpl_mask_threshold_image_create(img_other, -1, 5 * mad);
    // cpl_mask_not(kernel);
    // cpl_image_reject_from_mask(img_in, kernel);
    // cpl_mask_delete(kernel);
    // cpl_image_fill_rejected(img_in, 0);

    // No hdrl_image_wrap ?
    hdrl_other = hdrl_image_create(img_other, NULL);

    // Determine the peaks and remove peaks at the edges of the order
    if (cr2res_extract_sum_vert(hdrl_other, trace_wave, order, trace,
            height, &sfunc, &spec_bi, &model) != 0){
        return -1;
    }

    peaks = cr2res_etalon_get_maxpos(cpl_bivector_get_x(spec_bi));
    cr2res_slit_curv_remove_peaks_at_edge(&peaks, window, ncols);
    cpl_bivector_delete(spec_bi);
    cpl_vector_delete(sfunc);
    hdrl_image_delete(model);
    hdrl_image_delete(hdrl_other);

    // Rectify the image
    if ((ycen = cr2res_trace_get_ycen(trace_wave, order,
            trace, ncols)) == NULL){
        cpl_vector_delete(peaks);
        return -1 ;
    }
    if ((img_rect = cr2res_image_cut_rectify(img_in, ycen, height)) == NULL){
        cpl_vector_delete(peaks);
        cpl_vector_delete(ycen);
        return -1;
    }
    if (cpl_msg_get_level() == CPL_MSG_DEBUG){
        cpl_image_save(img_rect, "debug_image_rect.fits",
            CPL_TYPE_DOUBLE, NULL, CPL_IO_CREATE);
    }
    cpl_image_delete(img_other);


    // Determine the curvature of all peaks
    if (cr2res_slit_curv_all_peaks(peaks, img_rect, ycen, window,
            fit_second_order, &vec_b, &vec_c) != 0){
        cpl_vector_delete(peaks);
        cpl_vector_delete(ycen);
        cpl_image_delete(img_rect);
        return -1;
    }

    // Discard outliers
    // It would be better to make the fit itself more robust,
    // however this is a lot easier and works as long as there is no
    // strong variation within the order
    int npeaks = cpl_vector_get_size(vec_b);
    img_other = cpl_image_wrap_double(1, npeaks, cpl_vector_get_data(vec_b));
    median = cpl_image_get_mad(img_other, &mad);
    cpl_image_unwrap(img_other);
    for (int i = 0; i < npeaks; i++)
    {
        if (fabs(cpl_vector_get(vec_b, i) - median) > mad * 100)
            cpl_vector_set(vec_b, i, 0);
    }
    if (fit_second_order){
        img_other = cpl_image_wrap_double(1, npeaks, cpl_vector_get_data(vec_c));
        median = cpl_image_get_mad(img_other, &mad);
        cpl_image_unwrap(img_other);
        for (int i = 0; i < npeaks; i++)
        {
            if (fabs(cpl_vector_get(vec_c, i) - median) > mad * 100)
                cpl_vector_set(vec_c, i, 0);
        }
    }

    if (cpl_msg_get_level() == CPL_MSG_DEBUG){
        cpl_vector_save(vec_b, "debug_vector_b.fits",
            CPL_TYPE_DOUBLE, NULL, CPL_IO_CREATE);
        cpl_vector_save(vec_c, "debug_vector_c.fits",
            CPL_TYPE_DOUBLE, NULL, CPL_IO_CREATE);
    }
    cpl_image_delete(img_rect);
    cpl_vector_delete(ycen);

    // Fit the curvature to the whole order
    samppos = cpl_matrix_wrap(1, cpl_vector_get_size(peaks), 
                    cpl_vector_get_data(peaks));

    *slit_poly_a = cpl_polynomial_new(1);
    *slit_poly_b = cpl_polynomial_new(1);
    *slit_poly_c = cpl_polynomial_new(1);

    power = 1;
    cpl_polynomial_set_coeff(*slit_poly_a, &power, 1);
    cpl_polynomial_fit(*slit_poly_b, samppos, NULL, vec_b, NULL,
        CPL_TRUE, NULL, &degree);
    cpl_polynomial_fit(*slit_poly_c, samppos, NULL, vec_c, NULL,
        CPL_TRUE, NULL, &degree);
    cpl_matrix_unwrap(samppos);

    // Clean up memory
    cpl_vector_delete(vec_b);
    cpl_vector_delete(vec_c);
    cpl_vector_delete(peaks);


    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Compute the slit_curv map from the trace_wave table
  @param trace_wave     The trace wave table
  @param order          The order number or -1 for all
  @param trace_id       The trace_id number or -1 for all
  @param spacing_pixels The space in pixels between the traces
  @param full_trace     Draw on the full detector or only in the trace
  @return   the slit_curv_map image or NULL in error case

  The returned image must be deallocated with hdrl_image_delete()
 */
/*----------------------------------------------------------------------------*/
hdrl_image * cr2res_slit_curv_gen_map(
        const cpl_table *   trace_wave,
        int                 order,
        int                 trace_id,
        int                 spacing_pixels,
        int                 full_trace)
{
    hdrl_image      *   out ;
    cpl_image       *   out_ima ;
    double          *   pout_ima ;
    const cpl_array *   tmp_array ;
    cpl_polynomial  *   slit_poly_a ;
    cpl_polynomial  *   slit_poly_b ;
    cpl_polynomial  *   slit_poly_c ;
    cpl_polynomial  *   upper_poly ;
    cpl_polynomial  *   lower_poly ;
    cpl_polynomial  *   slit_curv_poly ;
    int                 cur_order, cur_trace_id ;
    double              upper_pos, lower_pos, x_slit_pos, value, val1, val2 ;
    cpl_size            i, j, k, nrows, nx, ny, x1, x2, ref_x ;

    /* Check Entries */
    if (trace_wave == NULL) return NULL ;

    /* Initialise */
    nrows = cpl_table_get_nrow(trace_wave) ;
    value = 10000. ;

    /* Create the image */
    out = hdrl_image_new(CR2RES_DETECTOR_SIZE, CR2RES_DETECTOR_SIZE) ;
    out_ima = hdrl_image_get_image(out) ;
    nx = cpl_image_get_size_x(out_ima) ;
    ny = cpl_image_get_size_y(out_ima) ;
    pout_ima = cpl_image_get_data_double(out_ima) ;

    /* Loop on the traces */
    for (k=0 ; k<nrows ; k++) {
        /* Only specified order / trace */
        cur_order = cpl_table_get(trace_wave, CR2RES_COL_ORDER, k, NULL) ;
        cur_trace_id = cpl_table_get(trace_wave, CR2RES_COL_TRACENB,k,NULL);
        if (order > -1 && order != cur_order) continue ;
        if (trace_id > -1 && trace_id != cur_trace_id) continue ;

        /* Get the Slit curvature for the trace */
        tmp_array = cpl_table_get_array(trace_wave, CR2RES_COL_SLIT_CURV_A, k) ;
        slit_poly_a = cr2res_convert_array_to_poly(tmp_array) ;
        tmp_array = cpl_table_get_array(trace_wave, CR2RES_COL_SLIT_CURV_B, k) ;
        slit_poly_b = cr2res_convert_array_to_poly(tmp_array) ;
        tmp_array = cpl_table_get_array(trace_wave, CR2RES_COL_SLIT_CURV_C, k) ;
        slit_poly_c = cr2res_convert_array_to_poly(tmp_array) ;

        if (slit_poly_a != NULL && slit_poly_b != NULL && slit_poly_c != NULL) {
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
                cpl_polynomial_delete(slit_poly_a) ;
                cpl_polynomial_delete(slit_poly_b) ;
                cpl_polynomial_delete(slit_poly_c) ;
                continue ;
            }

            /* Set the Pixels in the trace */
            for (i=0 ; i<nx ; i++) {
                ref_x = i+1 ;
                if ((ref_x % spacing_pixels) != 0) continue ;

                cpl_msg_debug(__func__, 
                        "Process Order/Trace: %d/%d - ref_x=%g", 
                        cur_order, cur_trace_id, (double)ref_x) ;

                /* Create the slit curvature polynomial at position i+1 */
                slit_curv_poly = cr2res_slit_curv_build_poly(slit_poly_a, 
                        slit_poly_b, slit_poly_c, ref_x) ;

                upper_pos = cpl_polynomial_eval_1d(upper_poly, ref_x, NULL) ;
                lower_pos = cpl_polynomial_eval_1d(lower_poly, ref_x, NULL) ;
                for (j=0 ; j<ny ; j++) {
                    if ((j+1 >= lower_pos && j+1 <= upper_pos) || full_trace) {
                        x_slit_pos = cpl_polynomial_eval_1d(slit_curv_poly, 
                                j+1, NULL) ;
                        x1 = (cpl_size)x_slit_pos ;
                        x2 = x1 + 1 ;
                        val1 = value * (x2-x_slit_pos) ;
                        val2 = value - val1 ;
                        /*cpl_msg_debug(__func__, 
                                "%"CPL_SIZE_FORMAT" %"CPL_SIZE_FORMAT" (%g) / %"CPL_SIZE_FORMAT" %"CPL_SIZE_FORMAT" (%g)", 
                                x1-1, j, val1, x2-1, j, val2) ; */
                        if (x1>=1 && x1<=nx) pout_ima[(x1-1)+j*nx] = val1 ;
                        if (x2>=1 && x2<=nx) pout_ima[(x2-1)+j*nx] = val2 ;
                    }
                }
                cpl_polynomial_delete(slit_curv_poly) ;
            }
            cpl_polynomial_delete(slit_poly_a) ;
            cpl_polynomial_delete(slit_poly_b) ;
            cpl_polynomial_delete(slit_poly_c) ;
            cpl_polynomial_delete(upper_poly) ;
            cpl_polynomial_delete(lower_poly) ;
        } else {
            if (slit_poly_a != NULL) cpl_polynomial_delete(slit_poly_a) ; 
            if (slit_poly_b != NULL) cpl_polynomial_delete(slit_poly_b) ; 
            if (slit_poly_c != NULL) cpl_polynomial_delete(slit_poly_c) ; 
        }
    }
    return out ;
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Create the slit curvature polynomial for x position
  @param slit_poly_a	Polynomial for the a coefficient
  @param slit_poly_b	Polynomial for the b coefficient
  @param slit_poly_c	Polynomial for the c coefficient
  @param x              The x position (1->2048)
  @return   the slit curvture polynomial or NULL in error case
 */
/*----------------------------------------------------------------------------*/
static cpl_polynomial * cr2res_slit_curv_build_poly(
        cpl_polynomial  *   slit_poly_a,
        cpl_polynomial  *   slit_poly_b,
        cpl_polynomial  *   slit_poly_c,
        cpl_size            x)
{
    cpl_polynomial  *   out ;
    cpl_size            power ;

    /* Generate the slit curvature polynomial */
   	out = cpl_polynomial_new(1) ;
    power = 0 ; cpl_polynomial_set_coeff(out, &power, 
            cpl_polynomial_eval_1d(slit_poly_a, x, NULL)) ;
    power = 1 ; cpl_polynomial_set_coeff(out, &power, 
            cpl_polynomial_eval_1d(slit_poly_b, x, NULL)) ;
    power = 2 ; cpl_polynomial_set_coeff(out, &power, 
            cpl_polynomial_eval_1d(slit_poly_c, x, NULL)) ;
    return out ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief Get the (X,Y) image cordinate of a Wavelength on a trace
  @param trace  The trace center position polynomial on which we search
  @param wave   The wavelength polynomial of the trace
  @param ref_wl The wavelength value
  @param xpos   [out] The X position of the trace point with the ref_wl
  @param ypos   [out] The Y position of the trace point with the ref_wl
  @return   0 if ok, -1 in error case

  The exact X position on the trace for which the wl is ref_wl is found
  Y is computed using X and the trace polynomial.
 */
/*----------------------------------------------------------------------------*/
static int cr2res_slit_curv_get_position(
        cpl_polynomial  *   trace,
        cpl_polynomial  *   wave,
        double              ref_wl,
        double          *   xpos,
        double          *   ypos)
{
    double          cur_wl, tmp_wl ;
    cpl_size        x ;

    /* Check entries */
    if (trace == NULL || wave == NULL || xpos == NULL || ypos == NULL) 
        return -1 ;

    /* Loop on the x positions  */
    for (x=0 ; x<=CR2RES_DETECTOR_SIZE+1 ; x++) {
        cur_wl = cpl_polynomial_eval_1d(wave, (double)x, NULL) ;
        /* As soon as the WL is bigger than ref_wl, we keep X */
        if (cur_wl > ref_wl) break ;
    }

    tmp_wl = cpl_polynomial_eval_1d(wave, (double)(x-1), NULL) ;

    /* Linear interpolation */
    if (fabs(cur_wl-tmp_wl) < 1e-5)
        *xpos = (double)x ;
    else 
        *xpos = (double)(x -((cur_wl-ref_wl) / (cur_wl-tmp_wl))) ;

    *ypos = cpl_polynomial_eval_1d(trace, *xpos, NULL) ;
    /* Check results */
    if (*xpos < 1 || *xpos > CR2RES_DETECTOR_SIZE ||
            *ypos < 1 || *ypos > CR2RES_DETECTOR_SIZE) {
        *xpos = *ypos = -1 ;
        return -1 ;
    }
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief Build a model of a shifted Gaussian peak
  @param x       The x and y position to evaluate
  @param a       The coefficients of the model
  @param result  The evaluated value of the model
  @return   0 if ok, -1 in error case

  Evaluates the model at position x using the coefficients a.
  Is used for the levenberg marquard least squares fit
 */
/*----------------------------------------------------------------------------*/
static int fmodel(const double x[], const double a[], double *result){
    // result = B + A * exp(-(x-mu)**2 / (2 * sigma**2))
    // with mu = center + tilt * y + shear * y**2
    double height = a[0];
    double bottom = a[1];
    double sigma = a[2];
    double center = a[3];
    double tilt = a[4];
    double shear = a[5];
    double nrows = a[6];
    double * slitfunc = (double*)(size_t)a[7];
    double * ycen = (double*)(size_t)a[8];

    double pos_x = x[0];
    double pos_y = x[1] - nrows;

    // Shift the center of the gaussian
    double curvature = center + (tilt + shear * pos_y) * pos_y;
    // Calculate the Gaussian
    double b1 = pos_x - curvature;
    *result = bottom + height * exp(-(b1 * b1) / (2 * sigma * sigma));
    *result *= slitfunc[(int)(x[1] - ycen[(int)pos_x] + nrows/2)];
    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief Calculate the derivative of model of a shifted Gaussian peak
  @param x       The x and y position to evaluate
  @param a       The coefficients of the model
  @param result  The evaluated value of the model
  @return   0 if ok, -1 in error case

  Evaluates the derivative of the model at position x using the coefficients a.
  Is used for the levenberg marquard least squares fit
 */
/*----------------------------------------------------------------------------*/
static int dmodel_da(const double x[], const double a[], double *result){
    double height = a[0];
    double bottom = a[1];
    double sigma = a[2];
    double center = a[3];
    double tilt = a[4];
    double shear = a[5];
    double nrows = a[6];

    double pos_x = x[0];
    double pos_y = x[1] - nrows;

    double curvature = center + (tilt + shear * pos_y) * pos_y;
    double b1 = pos_x - curvature;
    double b2 = sigma * sigma;
    double factor = exp(-(b1 * b1) / (2 * sigma * sigma));

    result[0] = factor;
    result[1] = 1.0;
    result[2] = height * factor * (b1 * b1 / b2 - 1) / sigma;
    result[3] = height * factor * b1 / b2;
    result[4] = result[3] * pos_y;
    result[5] = result[4] * pos_y;
    // The others dont matter?
    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief Remove peaks at the edge of the order
  @param peaks  [in/out] The peaks to check
  @param width  The width of the edge of the order at which to remove peaks
  @param ncols  The size of the order
  @return   0 if ok, -1 in error case

  All peaks smaller than width and larger than ncols - width will be removed
  from the peaks vector.
 */
/*----------------------------------------------------------------------------*/
static int cr2res_slit_curv_remove_peaks_at_edge(
    cpl_vector ** peaks,
    int width,
    int ncols)
{
    cpl_vector * peaks_filtered;
    int i, npeaks;
    double peak;

    // Remove peaks at the edges of the order
    peaks_filtered = cpl_vector_new(cpl_vector_get_size(*peaks));
    npeaks = 0;
    for (i = 0; i < cpl_vector_get_size(*peaks); i++){
        peak = cpl_vector_get(*peaks, i);
        if (peak <= width || peak >= ncols - width) continue;
        cpl_vector_set(peaks_filtered, npeaks, peak);
        npeaks++;
    }
    cpl_vector_set_size(peaks_filtered, npeaks);
    cpl_vector_delete(*peaks);
    *peaks = peaks_filtered;
    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief Fit the curvature of a single peak
  @param img_peak The image of the peak (and only the peak) of shape
  @param x        Matrix of shape (width * height, 2) for least-squares fitting
  @param y        Vector of shape (width * height,) for least-squares fitting
  @param a        Vector of shape (7,), for least-squares fitting
  @param ia       array of shape (7,), for least-squares fitting
  @param value_b  [out] fitted first order coefficient of the curvature
  @param value_c  [out] fitted second order coefficient of the curvature
  @return   0 if ok, -1 in error case

  Fits a model of a spectrum with a Gaussian Line Profile, where the peak
  position depends on the curvature and the row.

  The parameters x, y, a do not need to be filled, but have to exist.
  They are passed here purely to avoid declaring memory for them again, when
  the old will be sufficient. ia however needs to be declared as well, and will
  determine which parameters are used for fitting.
  For more details see cpl_fit_lvmq

 */
/*----------------------------------------------------------------------------*/
static int cr2res_slit_curv_single_peak(
    cpl_image  * img_peak,
    cpl_vector * ycen,
    double       peak,
    cpl_matrix * x,
    cpl_vector * y,
    cpl_vector * a,
    int        * ia,
    double     * value_b,
    double     * value_c)
{
    cpl_size n, j, k;
    double posx, posy;
    int width, height, window, badpix;
    cpl_error_code error;
    cpl_image  * img_slitfunc;
    double minimum, maximum;

    width = cpl_image_get_size_x(img_peak);
    height = cpl_image_get_size_y(img_peak);
    window = width / 2;

    // Set x and y matrix/vector
    n = 0;
    for (j = 0; j < width; j++){
        for (k = 0; k < height; k++){
            // We want to use the absolute reference frame
            // as thats what is desired in the output
            posx = (int)peak - window + j;
            posy = cpl_vector_get(ycen, posx) - height / 2 + k;
            cpl_matrix_set(x, n, 0, posx);
            cpl_matrix_set(x, n, 1, posy);
            cpl_vector_set(y, n, 
                cpl_image_get(img_peak, j+1, k+1, &badpix));
            n++;
        }
    }

    // Collapse image along y axis, to get slit illumination
    // normalized to the peak maximum
    img_slitfunc = cpl_image_collapse_median_create(img_peak, 1, 0, 0);
    maximum = cpl_image_get_max(img_slitfunc);
    minimum = cpl_image_get_min(img_slitfunc);
    // cpl_image_subtract_scalar(img_slitfunc, cpl_image_get_min(img_slitfunc));
    cpl_image_divide_scalar(img_slitfunc, maximum);

    // Set initial guesses
    // We pass the slitfunction and the ycen arrays as pointers disguised as
    // doubles in the vector. Its not very pretty, but it beats copying them
    // into the vector
    cpl_vector_set(a, 0, maximum - minimum);
    cpl_vector_set(a, 1, minimum);
    cpl_vector_set(a, 2, (double)window / 4.);
    cpl_vector_set(a, 3, peak);
    cpl_vector_set(a, 4, 0);
    cpl_vector_set(a, 5, 0);
    cpl_vector_set(a, 6, height);
    cpl_vector_set(a, 7, (double)(size_t)cpl_image_get_data(img_slitfunc));
    cpl_vector_set(a, 8, (double)(size_t)cpl_vector_get_data(ycen));

    // Fit the model
    error = cpl_fit_lvmq(x, NULL, y, NULL, a, ia, fmodel, dmodel_da, 
        CPL_FIT_LVMQ_TOLERANCE, CPL_FIT_LVMQ_COUNT, CPL_FIT_LVMQ_MAXITER,
        NULL, NULL, NULL);

    cpl_image_delete(img_slitfunc);

    // Check results
    if (error != CPL_ERROR_NONE){
        cpl_errorstate_dump(NULL, CPL_FALSE,
                        cpl_errorstate_dump_one_debug);
        cpl_error_reset();
        
        *value_b = *value_c = 0.;
        return -1;
    } else {
        *value_b = cpl_vector_get(a, 4);
        *value_c = cpl_vector_get(a, 5);
        return 0;
    }
}

/*----------------------------------------------------------------------------*/
/**
  @brief Fit the curvature of all single peak
  @param peaks    The peak positions
  @param img_rect The rectified image of the order
  @param ycen     The central line of the order (in global frame)
  @param window   The width of each peak to look at
  @param fit_second_order Whether to fit the second order curvature
  @param vec_b  [out] fitted first order coefficients of the curvature
  @param vec_c  [out] fitted second order coefficients of the curvature
  @return   0 if ok, -1 in error case
  Loops over the peaks, fitting the curvature of each using
  cr2res_slit_curv_single_peak.
 */
/*----------------------------------------------------------------------------*/
static int cr2res_slit_curv_all_peaks(
    cpl_vector * peaks,
    cpl_image  * img_rect,
    cpl_vector * ycen,
    int window,
    int fit_second_order,
    cpl_vector ** vec_b,
    cpl_vector ** vec_c)
{
    cpl_image * img_model;
    cpl_image * img_peak;
    int width, height, npeaks;
    double pos[2];
    double peak, value_b, value_c, result;
    cpl_size i, j, k;

    width = 2 * window + 1;
    height = cpl_image_get_size_y(img_rect);
    npeaks = cpl_vector_get_size(peaks);

    cpl_matrix * x;
    cpl_vector * y;
    cpl_vector * a;
    int N, D, M;

    int ia[] = {1, 1, 1, 1, 1, fit_second_order, 0, 0, 0};

    if (peaks == NULL || img_rect == NULL || ycen == NULL || 
        vec_b == NULL || vec_c == NULL) return -1;

    N = width * height;
    D = 2;
    M = 9;

    if (cpl_msg_get_level() == CPL_MSG_DEBUG){
        img_model = cpl_image_new(2 * window + 1, height, CPL_TYPE_DOUBLE);
    } else img_model = NULL;
    x = cpl_matrix_new(N, D);
    y = cpl_vector_new(N);
    a = cpl_vector_new(M);

    *vec_b = cpl_vector_new(cpl_vector_get_size(peaks));
    *vec_c = cpl_vector_new(cpl_vector_get_size(peaks));

    for (i = 0; i < npeaks; i++){
        peak = cpl_vector_get(peaks, i);
        img_peak  = cpl_image_extract(img_rect, peak - window, 1,
            peak + window, height);

        cr2res_slit_curv_single_peak(img_peak, ycen, peak, x, y, a, ia,
            &value_b, &value_c);
        cpl_vector_set(*vec_b, i, value_b);
        cpl_vector_set(*vec_c, i, value_c);

        if (cpl_msg_get_level() == CPL_MSG_DEBUG){
            for (j = 0; j < 2 * window + 1; j++){
                for (k = 0; k < height; k++){
                    pos[0] = (int)peak - window + j;
                    pos[1] = cpl_vector_get(ycen, pos[0]) - height / 2 + k;
                    fmodel(pos, cpl_vector_get_data(a), &result);
                    cpl_image_set(img_model, j+1, k+1, result);
                }
            }
            cpl_image_save(img_peak, "debug_image_curv.fits",
                CPL_TYPE_DOUBLE, NULL, CPL_IO_CREATE);
            cpl_image_save(img_model, "debug_model_curv.fits",
                CPL_TYPE_DOUBLE, NULL, CPL_IO_CREATE);
        }
        cpl_image_delete(img_peak);
    }
    cpl_matrix_delete(x);
    cpl_vector_delete(y);
    cpl_vector_delete(a);
    if (cpl_msg_get_level() == CPL_MSG_DEBUG){
        cpl_image_delete(img_model);
    }

    return 0;
}

