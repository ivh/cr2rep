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

static int cr2res_slit_curv_get_position(
        cpl_polynomial  *   trace,
        cpl_polynomial  *   wave,
        double              ref_wl,
        double          *   xpos,
        double          *   ypos) ;

static int fmodel(
    const double x[],
    const double a[],
    double *result) CPL_ATTR_NONNULL;
static int dmodel_da(
    const double x[],
    const double a[],
    double *result) CPL_ATTR_NONNULL;
static int cr2res_slit_curv_remove_peaks_at_edge(
    cpl_vector ** peaks,
    const int width,
    const int ncols) ;
static int cr2res_slit_curv_smooth_image_median(
    hdrl_image * hdrl_out,
    const cpl_image * img_in,
    const cpl_size kernel_size
);
static int cr2res_slit_curv_remove_outliers(
    cpl_vector * peaks,
    cpl_vector * vec_a,
    cpl_vector * vec_b,
    cpl_vector * vec_c,
    const int fit_second_order,
    const double divergence
);
static int cr2res_slit_curv_single_peak(
    const cpl_image  * img_peak,
    const cpl_vector * ycen,
    const double       peak,
    cpl_matrix       * x,
    cpl_vector       * y,
    cpl_vector       * a,
    const int        * ia,
    double           * value_a,
    double           * value_b,
    double           * value_c) ;
static int cr2res_slit_curv_all_peaks(
    const cpl_vector * peaks,
    const cpl_image  * img_rect,
    const cpl_vector * ycen,
    const int          window,
    const int          fit_second_order,
    cpl_vector      ** vec_a,
    cpl_vector      ** vec_b,
    cpl_vector      ** vec_c) ;

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
  the spectrum. Next for each peak, fit a Gaussian with a shifted peak (as a
  function of the vertical offset y) to the image. The shift of the peak
  determines the curvature. Finally fit the curvature coefficients to all
  peaks along the order.

  All coefficients are given in the global reference frame
  The transformation to the local frame is given by:
    a(x) += - x + yc * b(x) + yc * yc * c(x)
    b(x) += 2 * yc * c(x)

 */
/*----------------------------------------------------------------------------*/
int cr2res_slit_curv_compute_order_trace(
        const hdrl_image    *   img,
        const cpl_table     *   trace_wave,
        const int               order,
        const int               trace,
        const int               height,
        const int               window,
        const cpl_size          change_degree,
        const int               slit_degree,
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
    cpl_vector          *   vec_a;
    cpl_vector          *   vec_b;
    cpl_vector          *   vec_c;
    cpl_size                power;
    cpl_matrix          *   samppos;
    hdrl_image          *   hdrl_other;
    int                     fit_second_order;

    if (img == NULL || trace_wave == NULL || slit_poly_a == NULL ||
        slit_poly_b == NULL || slit_poly_c == NULL) return -1;
    if (slit_degree == 1)
        fit_second_order = 0;
    else if (slit_degree == 2)
        fit_second_order = 1;
    else {
        cpl_msg_error(__func__, "Only degree 1 or 2 are valid") ;
        return -1 ;
    }
    

    img_in = hdrl_image_get_image_const(img);
    const int ncols = cpl_image_get_size_x(img_in);

    // Median filter the image
    // to remove outliers, which would mess with the peak detection
    hdrl_other = hdrl_image_new(cpl_image_get_size_x(img_in),
        cpl_image_get_size_y(img_in));
    if (cr2res_slit_curv_smooth_image_median(hdrl_other, img_in, 3) != 0){
        return -1;
    }
    img_in = hdrl_image_get_image_const(hdrl_other);
    
    // Determine the peaks and remove peaks at the edges of the order
    if (cr2res_extract_sum_vert(hdrl_other, trace_wave, order, trace,
            height, &sfunc, &spec_bi, &model) != 0){
        return -1;
    }
    peaks = cr2res_etalon_find_peaks(cpl_bivector_get_x(spec_bi), 
        cpl_vector_get_mean(cpl_bivector_get_x(spec_bi)), 3);
    cr2res_slit_curv_remove_peaks_at_edge(&peaks, window, ncols);
    cpl_bivector_delete(spec_bi);
    cpl_vector_delete(sfunc);
    hdrl_image_delete(model);
    

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

    // Determine the curvature of all peaks
    if (cr2res_slit_curv_all_peaks(peaks, img_rect, ycen, window,
            fit_second_order, &vec_a, &vec_b, &vec_c) != 0){
        cpl_vector_delete(peaks);
        cpl_vector_delete(ycen);
        cpl_image_delete(img_rect);
        return -1;
    }
    cpl_image_delete(img_rect);
    cpl_vector_delete(ycen);
    hdrl_image_delete(hdrl_other);

    // Discard outliers
    // It would be better to make the fit itself more robust,
    // however this is a lot easier and works as long as there is no
    // strong variation within the order
    cr2res_slit_curv_remove_outliers(peaks, vec_a, vec_b,
        vec_c, fit_second_order, 5);

    if (cpl_msg_get_level() == CPL_MSG_DEBUG){
        cpl_vector_save(peaks, "debug_peaks.fits",
            CPL_TYPE_DOUBLE, NULL, CPL_IO_CREATE);
        cpl_vector_save(vec_a, "debug_vector_a.fits",
            CPL_TYPE_DOUBLE, NULL, CPL_IO_CREATE);
        cpl_vector_save(vec_b, "debug_vector_b.fits",
            CPL_TYPE_DOUBLE, NULL, CPL_IO_CREATE);
        cpl_vector_save(vec_c, "debug_vector_c.fits",
            CPL_TYPE_DOUBLE, NULL, CPL_IO_CREATE);
    }

    // Fit the curvature to the whole order
    samppos = cpl_matrix_wrap(1, cpl_vector_get_size(peaks), 
                    cpl_vector_get_data(peaks));

    *slit_poly_a = cpl_polynomial_new(1);
    *slit_poly_b = cpl_polynomial_new(1);
    *slit_poly_c = cpl_polynomial_new(1);

    cpl_polynomial_fit(*slit_poly_a, samppos, NULL, vec_a, NULL,
        CPL_TRUE, NULL, &change_degree);
    cpl_polynomial_fit(*slit_poly_b, samppos, NULL, vec_b, NULL,
        CPL_TRUE, NULL, &change_degree);
    cpl_polynomial_fit(*slit_poly_c, samppos, NULL, vec_c, NULL,
        CPL_TRUE, NULL, &change_degree);
    cpl_matrix_unwrap(samppos);

    // Add 1 to the linear coefficient of the a polynomial
    // So we can subtract it later in the extraction
    power = 1;
    cpl_polynomial_set_coeff(*slit_poly_a, &power,
        1 + cpl_polynomial_get_coeff(*slit_poly_a, &power));

    // Clean up memory
    cpl_vector_delete(vec_a);
    cpl_vector_delete(vec_b);
    cpl_vector_delete(vec_c);
    cpl_vector_delete(peaks);

    if (cpl_error_get_code() != CPL_ERROR_NONE){
        cpl_errorstate_dump(NULL, CPL_FALSE,
                        cpl_errorstate_dump_one);
        // Probably the fit failed
        // Maybe not enough peaks?
        return -1;
    }

    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief Smooth an image with a median filter
  @param hdrl_out   [out] an existing hdrl image of the correct size
  @param img_in     The image to smooth
  @param kernel_size The size of the median kernel in both directions
  @return   0 if ok, -1 in error case

 */
/*----------------------------------------------------------------------------*/
static int cr2res_slit_curv_smooth_image_median(
    hdrl_image * hdrl_out,
    const cpl_image * img_in,
    const cpl_size kernel_size
){
    cpl_image * img_other;
    cpl_mask * kernel;

    img_other = hdrl_image_get_image(hdrl_out);
    kernel = cpl_mask_new(kernel_size, kernel_size);
    cpl_mask_not(kernel);
    cpl_image_filter_mask(img_other, img_in, kernel,
            CPL_FILTER_MEDIAN, CPL_BORDER_FILTER);
    cpl_mask_delete(kernel);

    if (cpl_error_get_code() != CPL_ERROR_NONE){
        return -1;
    }

    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief Remove extreme outliers in the vectors
  @param peaks   [in/out] vector with peak positions
  @param vec_a   [in/out] vector of a coefficient values
  @param vec_b   [in/out] vector of b coefficient values
  @param vec_c   [in/out] vector of c coefficient values
  @param fit_second_order Whether the c coefficient was fitted or not
  @return   0 if ok, -1 in error case

  The outliers are detected using the median and the median absolute
  deviation (mad). Extreme outlier in this case means that the value diverges
  by 100 mad from the median. Something that should not happen even, for
  strong dependance on x.

  If fit_second_order is True, the c coefficients are not checked for outliers,
  but values are still removed, if outliers are present in other vectors.

  Note that it is assumed that all vectors are of the same number of elements.

 */
/*----------------------------------------------------------------------------*/
static int cr2res_slit_curv_remove_outliers(
    cpl_vector * peaks,
    cpl_vector * vec_a,
    cpl_vector * vec_b,
    cpl_vector * vec_c,
    const int fit_second_order,
    const double divergence
)
{
    cpl_vector * remove_peaks;
    cpl_image  * img_other;
    cpl_size npeaks;
    cpl_size i, j;
    double median, mad;

    // First step is to find the peaks to remove
    // We use the median and the median absolute deviation (mad)

    npeaks = cpl_vector_get_size(peaks);
    remove_peaks = cpl_vector_new(npeaks);
    for (i=0; i < npeaks; i++) cpl_vector_set(remove_peaks, i, 0);

    // Unfortunately we have to wrap the vector data in an image
    // since the relevant function is only available for images 
    img_other = cpl_image_wrap_double(1, npeaks, cpl_vector_get_data(vec_a));
    median = cpl_image_get_mad(img_other, &mad);
    cpl_image_unwrap(img_other);
    for (i = 0; i < npeaks; i++)
    {
        if (fabs(cpl_vector_get(vec_a, i) - median) > mad * divergence)
            cpl_vector_set(remove_peaks, i, 1);
    }
    img_other = cpl_image_wrap_double(1, npeaks, cpl_vector_get_data(vec_b));
    median = cpl_image_get_mad(img_other, &mad);
    cpl_image_unwrap(img_other);
    for (i = 0; i < npeaks; i++)
    {
        if (fabs(cpl_vector_get(vec_b, i) - median) > mad * divergence)
            cpl_vector_set(remove_peaks, i, 1);
    }
    if (fit_second_order){
        img_other = cpl_image_wrap_double(1, npeaks,
            cpl_vector_get_data(vec_c));
        median = cpl_image_get_mad(img_other, &mad);
        cpl_image_unwrap(img_other);
        for (i = 0; i < npeaks; i++)
        {
            if (fabs(cpl_vector_get(vec_c, i) - median) > mad * divergence)
                cpl_vector_set(remove_peaks, i, 1);
        }
    }

    // Remove the outlier peaks from the vectors
    // Overwrite the removed peaks, with data from the next entries
    j = 0;
    for (i = 0; i < npeaks; i++){
        if (cpl_vector_get(remove_peaks, i) != 0){
            j++;
        } else {
            cpl_vector_set(peaks, i-j, cpl_vector_get(peaks, i));
            cpl_vector_set(vec_a, i-j, cpl_vector_get(vec_a, i));
            cpl_vector_set(vec_b, i-j, cpl_vector_get(vec_b, i));
            cpl_vector_set(vec_c, i-j, cpl_vector_get(vec_c, i));
        }
    }
    npeaks -= j;
    cpl_vector_set_size(peaks, npeaks);
    cpl_vector_set_size(vec_a, npeaks);
    cpl_vector_set_size(vec_b, npeaks);
    cpl_vector_set_size(vec_c, npeaks);

    cpl_vector_delete(remove_peaks);

    if (cpl_error_get_code() != CPL_ERROR_NONE){
        return -1;
    }
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
cpl_polynomial * cr2res_slit_curv_build_poly(
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

/**@}*/

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
    const double height = a[0];
    const double bottom = a[1];
    const double sigma = a[2];
    const double center = a[3];
    const double tilt = a[4];
    const double shear = a[5];
    const double nrows = a[6];
    // We pass the slitfunc and ycen arrays by reference
    // but since a is a double, we have to get tricky with the casting
    const double * slitfunc = (const double*)(intptr_t)a[7];
    const double * ycen = (const double*)(intptr_t)a[8];
    const double yc = ycen[(size_t)x[0]]; 

    const double pos_x = x[0];
    const double pos_y = x[1];
    const double pos_y_rel = x[1] - yc;

    // Shift the center of the gaussian
    // The curvature should be 0, when pos_y = ycen (in the global frame)
    // i.e. we fix the offset a = - tilt * yc - shear * yc * yc
    // Which means we are using x**2 - yc**2 and NOT (x-yc)**2
    const double curvature =
        tilt * pos_y_rel + shear * (pos_y * pos_y - yc * yc);
    // Calculate the Gaussian
    const double b1 = pos_x - center - curvature;
    *result = bottom + height * exp(-(b1 * b1) / (2 * sigma * sigma));
    *result *= slitfunc[(size_t)(pos_y_rel + nrows/2)];
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
    const double height = a[0];
    const double bottom = a[1];
    const double sigma = a[2];
    const double center = a[3];
    const double tilt = a[4];
    const double shear = a[5];
    const double nrows = a[6];
    const double * ycen = (double*)(intptr_t)a[8];
    const double yc = ycen[(size_t)x[0]];

    const double pos_x = x[0];
    const double pos_y = x[1];
    const double pos_y_rel = x[1] - yc;
    const double pos_y_squared = pos_y * pos_y - yc * yc;

    const double curvature = tilt * pos_y_rel + shear * pos_y_squared;
    const double b1 = pos_x - center - curvature;
    const double b2 = sigma * sigma;
    const double factor = exp(-(b1 * b1) / (2 * sigma * sigma));

    result[0] = factor;
    result[1] = 1.0;
    result[2] = height * factor * (b1 * b1 / b2 - 1) / sigma;
    result[3] = height * factor * b1 / b2;
    result[4] = result[3] * pos_y_rel;
    result[5] = result[3] * pos_y_squared;
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
  Peaks is assumed to be sorted, for this to work!
  Note that this may change the size of the peaks vector
 */
/*----------------------------------------------------------------------------*/
static int cr2res_slit_curv_remove_peaks_at_edge(
    cpl_vector ** peaks,
    const int width,
    const int ncols)
{
    cpl_size i, j, npeaks;
    double peak;

    // Loop through the vector and shift peaks to the beginning of the vector,
    // overwriting existing values, of peaks at the edge
    j = 0;
    npeaks = cpl_vector_get_size(*peaks);
    for (i = 0; i < npeaks; i++){
        peak = cpl_vector_get(*peaks, i);
        if (peak <= width || peak >= ncols - width) continue;
        cpl_vector_set(*peaks, j, cpl_vector_get(*peaks, i));
        j++;
    }
    cpl_vector_set_size(*peaks, j);

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
  @param value_a  [out] fitted constant coeffcient of the curvature
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

  Note that the coefficients are fitted for the global reference frame of the
  image. In the local frame value_a is forced to 0.

  The transformation to the local frame is given by:
    a = a + yc * b + yc * yc * c
    b += 2 * yc * c

 */
/*----------------------------------------------------------------------------*/
static int cr2res_slit_curv_single_peak(
    const cpl_image  * img_peak,
    const cpl_vector * ycen,
    const double       peak,
    cpl_matrix       * x,
    cpl_vector       * y,
    cpl_vector       * a,
    const int        * ia,
    double           * value_a,
    double           * value_b,
    double           * value_c)
{
    cpl_size n, j, k;
    int window, badpix;
    cpl_error_code error;
    cpl_image * img_slitfunc;
    cpl_image * img_model;
    cpl_image * img_spec;
    cpl_matrix * x_extract;
    cpl_vector * y_extract;
    double minimum, maximum;
    double yc, result;
    double pos[2];
    double pix_value;

    const int width = cpl_image_get_size_x(img_peak);
    const int height = cpl_image_get_size_y(img_peak);
    window = width / 2;

    // Set x and y matrix/vector
    n = 0;
    for (j = 0; j < width; j++){
        for (k = 0; k < height; k++){
            // We want to use the absolute reference frame
            // as thats what is desired in the output
            pix_value = cpl_image_get(img_peak, j+1, k+1, &badpix);
            // Filter out bad pixels
            if (!badpix && !isnan(pix_value)){
                pos[0] = peak - window + j;
                pos[1] = cpl_vector_get(ycen, (cpl_size)pos[0]) - height / 2 + k;
                cpl_matrix_set(x, n, 0, pos[0]);
                cpl_matrix_set(x, n, 1, pos[1]);
                cpl_vector_set(y, n, pix_value);
                n++;
            }
        }
    }

    // Collapse image along y axis, to get slit illumination
    // normalized to the peak maximum
    img_slitfunc = cpl_image_collapse_median_create(img_peak, 1, 0, 0);
    maximum = cpl_image_get_max(img_slitfunc);
    // minimum = cpl_image_get_min(img_slitfunc);
    if (maximum <= 0){
        // Abort before division by 0, and weird peaks
        cpl_image_delete(img_slitfunc);
        return -1;
    }

    cpl_image_divide_scalar(img_slitfunc, maximum);

    // Collapse image along x axis, to get spectrum
    img_spec = cpl_image_collapse_median_create(img_peak, 0, 0, 0);
    maximum = cpl_image_get(img_spec, cpl_image_get_size_x(img_spec) / 2, 1,
        &badpix);
    minimum = cpl_image_get(img_spec, 1, 1, &badpix);
    cpl_image_delete(img_spec);

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
    cpl_vector_set(a, 7, (double)(intptr_t)cpl_image_get_data(img_slitfunc));
    cpl_vector_set(a, 8, (double)(intptr_t)cpl_vector_get_data_const(ycen));

    // We want to reuse the memory in every peak, but we only need a fraction of
    // it in every iteration. To adjust the size we wrap the larger memory
    // in a smaller matrix/vector, that only has the good pixels in it
    x_extract = cpl_matrix_wrap(n, 2, cpl_matrix_get_data(x));
    y_extract = cpl_vector_wrap(n, cpl_vector_get_data(y));

    // Fit the model
    error = cpl_fit_lvmq(x_extract, NULL, y_extract, NULL, a, ia, fmodel, dmodel_da, 
        CPL_FIT_LVMQ_TOLERANCE, CPL_FIT_LVMQ_COUNT, CPL_FIT_LVMQ_MAXITER,
        NULL, NULL, NULL);

    cpl_matrix_unwrap(x_extract);
    cpl_vector_unwrap(y_extract);

    if (cpl_msg_get_level() == CPL_MSG_DEBUG){
        img_model = cpl_image_new(width, height, CPL_TYPE_DOUBLE);
        for (j = 0; j < width; j++){
            for (k = 0; k < height; k++){
                pos[0] = peak - window + j;
                pos[1] = cpl_vector_get(ycen, (cpl_size)pos[0])
                    - height / 2 + k;
                fmodel(pos, cpl_vector_get_data(a), &result);
                cpl_image_set(img_model, j+1, k+1, result);
            }
        }
        cpl_image_save(img_peak, "debug_image_curv.fits",
            CPL_TYPE_DOUBLE, NULL, CPL_IO_CREATE);
        cpl_image_save(img_model, "debug_model_curv.fits",
            CPL_TYPE_DOUBLE, NULL, CPL_IO_CREATE);
        cpl_image_delete(img_model);
    }

    cpl_image_delete(img_slitfunc);

    // Check results
    if (error != CPL_ERROR_NONE){
        cpl_msg_debug(__func__, "%s", cpl_error_get_message());
        cpl_error_reset();
        
        *value_a = *value_b = *value_c = 0.;
        return -1;
    } else {
        *value_b = cpl_vector_get(a, 4);
        *value_c = cpl_vector_get(a, 5);
        // The offset a was fixed so that is 0 in the local frame
        // i.e. so that the curvature at the central line is 0
        // a = - tilt * ycen - shear * ycen * ycen
        yc = cpl_vector_get(ycen, (cpl_size)peak);
        *value_a = -(*value_b * yc + *value_c * yc * yc);
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
  @param vec_a  [out] fitted constant coefficients of the curvature
  @param vec_b  [out] fitted first order coefficients of the curvature
  @param vec_c  [out] fitted second order coefficients of the curvature
  @return   0 if ok, -1 in error case
  Loops over the peaks, fitting the curvature of each using
  cr2res_slit_curv_single_peak.
  
  All coefficients are given in the global reference frame of the detector
  The transformation to the local frame is given by:
    a(x) += yc * b(x) + yc * yc * c(x)
    b(x) += 2 * yc * c(x)
 */
/*----------------------------------------------------------------------------*/
static int cr2res_slit_curv_all_peaks(
    const cpl_vector * peaks,
    const cpl_image  * img_rect,
    const cpl_vector * ycen,
    const int          window,
    const int          fit_second_order,
    cpl_vector      ** vec_a,
    cpl_vector      ** vec_b,
    cpl_vector      ** vec_c)
{
    cpl_image * img_peak;
    double peak, value_a, value_b, value_c;
    cpl_size i;

    const int width = 2 * window + 1;
    const int height = cpl_image_get_size_y(img_rect);
    const int npeaks = cpl_vector_get_size(peaks);

    cpl_matrix * x;
    cpl_vector * y;
    cpl_vector * a;

    // Fit parameters are: 
    // The peak height of the Gaussian
    // The bottom level of the Gaussian
    // The width (sigma) of the Gaussian
    // The center position (as ofset from the peak pixel)
    // The tilt (first order curvature)
    // The shear (second order curvature), only if fit_second_order is not 0
    // Fixed parameters are:
    // Number of rows in the image
    // The slit illumination function array (as a pointer)
    // The central row (ycen) array (as a pointer)
    const int N = width * height;
    const int D = 2;
    const int M = 9;

    const int ia[] = {1, 1, 1, 1, 1, fit_second_order, 0, 0, 0};

    if (peaks == NULL || img_rect == NULL || ycen == NULL || 
        vec_a == NULL || vec_b == NULL || vec_c == NULL) return -1;

    x = cpl_matrix_new(N, D);
    y = cpl_vector_new(N);
    a = cpl_vector_new(M);

    *vec_a = cpl_vector_new(cpl_vector_get_size(peaks));
    *vec_b = cpl_vector_new(cpl_vector_get_size(peaks));
    *vec_c = cpl_vector_new(cpl_vector_get_size(peaks));

    for (i = 0; i < npeaks; i++){
        // Loop over the individual lines
        // and fit each of them seperately
        peak = cpl_vector_get(peaks, i);
        img_peak  = cpl_image_extract(img_rect, peak - window, 1,
            peak + window, height);

        cr2res_slit_curv_single_peak(img_peak, ycen, peak, x, y, a, ia,
            &value_a, &value_b, &value_c);
        cpl_vector_set(*vec_a, i, value_a);
        cpl_vector_set(*vec_b, i, value_b);
        cpl_vector_set(*vec_c, i, value_c);
        cpl_image_delete(img_peak);
    }
    cpl_matrix_delete(x);
    cpl_vector_delete(y);
    cpl_vector_delete(a);

    return 0;
}

