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
#include <time.h>

#include <cpl.h>

#include "irplib_wcs.h"

#include "cr2res_utils.h"
#include "cr2res_io.h"
#include "cr2res_pfits.h"
#include "cr2res_dfs.h"
#include "cr2res_extract.h"
#include "cr2res_trace.h"


/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_utils     Miscellaneous Utilities
 */
/*----------------------------------------------------------------------------*/
#define max(a,b) (((a)>(b))?(a):(b))
/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the current MJD-OBS
  @return   The MJD-OBS
 */
/*----------------------------------------------------------------------------*/
double cr2res_mjd_obs_now(void)
{
    double          mjd_obs ;
    int             hours, minutes, seconds, day, month, year;
    time_t          now;

    /* Get the time for MJD-OBS */
    time(&now);
    struct tm *local = localtime(&now);
    hours = local->tm_hour;          // get hours since midnight (0-23)
    minutes = local->tm_min;         // get minutes passed after the hour (0-59)
    seconds = local->tm_sec;         // get seconds passed after minute (0-59)
    day = local->tm_mday;            // get day of month (1 to 31)
    month = local->tm_mon + 1;       // get month of year (0 to 11)
    year = local->tm_year + 1900;    // get year since 1900
    /* printf("%d %d %d %d %d %d\n", hours, minutes, seconds, day,
     * month, year) ; */
    irplib_wcs_mjd_from_iso8601(&mjd_obs, year, month, day, hours,
            minutes, seconds) ;

    return mjd_obs ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Convert the order_idx into order_real
  @param    order_idx   Order index
  @param    order_zp    Order Zero Point
  @return   order_real
 */
/*----------------------------------------------------------------------------*/
int cr2res_order_idx_to_real(int order_idx, int order_zp)
{
    return order_zp + order_idx - 1;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Convert the order_real into order_idx
  @param    order_real  Order Real
  @param    order_zp    Order Zero Point
  @return   order_idx
 */
/*----------------------------------------------------------------------------*/
int cr2res_order_real_to_idx(int order_real, int order_zp)
{
    return order_real - order_zp + 1;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Convert the RA from hh mm ss to degrees 
  @param    hh      hours
  @param    mm      minutes 
  @param    ss      seconds 
  @return   RA in degrees 

  An arc-hour is 15 degrees,
  60 arc-minutes is one arc-hour and
  60 arc-seconds is one arc-minute.
  
 */
/*----------------------------------------------------------------------------*/
double cr2res_ra_hms2deg(int hh, int mm, double ss)
{
    return 15.0 * cr2res_dec_hms2deg(hh, mm, ss);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Convert the DEC from dd mm ss to degrees
  @param    dd      degrees
  @param    mm      minutes 
  @param    ss      seconds  
  @return   DEC in degrees

  60 arc-minutes is one degree and
  60 arc-seconds is one arc-minute.

 */
/*----------------------------------------------------------------------------*/
double cr2res_dec_hms2deg(int dd, int mm, double ss)
{
    return ((double)ss/60.0 + (double)mm)/60.0 + dd;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Format the setting
  @param    Setting
  @return   0 if ok, -1 in error case
    replace / by _ in the setting string
 */
/*----------------------------------------------------------------------------*/
int cr2res_format_setting(char * setting_id)
{
    int     i, len ;

    /* Check entries */
    if (setting_id == NULL) return -1 ;

    len = strlen(setting_id) ;
    for (i=0 ; i<len ; i++) if (setting_id[i] == '/') setting_id[i] = '_' ;
    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Identify Short Wavelength
  @param    Setting
  @return   1 if SW, 0 if LW, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_is_short_wavelength(char * setting_id)
{
    /* Check entries */
    if (setting_id == NULL) return -1 ;

    if (setting_id[0] == 'H' || setting_id[0] == 'J' || 
            setting_id[0] == 'K' || setting_id[0] == 'Y')
        return 1;
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Format the setting
  @param    Setting
  @return   0 if ok, -1 in error case
    replace _ by / in the setting string
 */
/*----------------------------------------------------------------------------*/
int cr2res_format_setting2(char * setting_id)
{
    int     i, len ;

    /* Check entries */
    if (setting_id == NULL) return -1 ;

    len = strlen(setting_id) ;
    for (i=0 ; i<len ; i++) if (setting_id[i] == '_') setting_id[i] = '/' ;
    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get MAD from a vector
  @param    invec   input vecotr
  @param    mad*    median absolute deviation (output)
  @return   median

  Wrap a vector into an image because the cpl does not have this for vectors
 */
/*----------------------------------------------------------------------------*/
double cr2res_vector_get_mad(cpl_vector * invec, double *mad)
{
    cpl_image * img;
    double median;
    img = cpl_image_wrap_double(1,cpl_vector_get_size(invec),
                cpl_vector_get_data(invec));
    median = cpl_image_get_mad(img, mad);
    cpl_image_unwrap(img);
    return median;
}

/*----------------------------------------------------------------------------*/
/**
  @brief
  @param    
  @return
 */
/*----------------------------------------------------------------------------*/
int * cr2res_vector_get_int(
    const cpl_vector    * ycen)
{
    int         * ycen_int;
    int         i, lenx;

    if (ycen == NULL) return NULL;

    lenx = cpl_vector_get_size(ycen);
    ycen_int = cpl_malloc(lenx*sizeof(int));
    for (i=0 ; i<lenx ; i++){
        ycen_int[i] = (int)cpl_vector_get(ycen,i);
    }
   return ycen_int;
}

/*----------------------------------------------------------------------------*/
/**
  @brief
  @param    ycen
  @return
 */
/*----------------------------------------------------------------------------*/
double * cr2res_vector_get_rest(
    const cpl_vector    * ycen)
{
    double      * ycen_rest;
    int         i, lenx ;

    if (ycen == NULL) return NULL;

    lenx = cpl_vector_get_size(ycen);
    ycen_rest = cpl_malloc(lenx*sizeof(double));
    for (i=0 ; i<lenx ; i++){
         ycen_rest[i] = fmod( cpl_vector_get(ycen,i), 1.0) ;
    }
   return ycen_rest;
}

 
/*----------------------------------------------------------------------------*/
/**
  @brief    Combine two extracted spectra by resampling one to the other's WL.
  @param    extracta        A position extracted spectrum
  @param    extractb        B position extracted spectrum
  @return   extract_combined table
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_combine_extracted(
        const cpl_table     *  extracta,
        const cpl_table     *  extractb)
{
    cpl_table                   *   extractc ;
    cpl_array                   *   col_names ;
    const char                  *   col_name ;
    char                        *   err_col ;
    char                        *   col_type ;
    char                        *   wave_col ;
    hdrl_spectrum1D_wave_scale      scale ;
    hdrl_spectrum1D             *   a_spec ;
    hdrl_spectrum1D             *   b_spec ;
    hdrl_spectrum1D             *   c_spec ;
    hdrl_spectrum1D_wavelength      spec_wav ;
    hdrl_parameter              *   params ;
    double                      *   p_flux ;
    double                      *   p_err ;
    cpl_size                        ncols, i, j, sz, sz_a, sz_b ;
    int                             trace_nb, order, increasing_values ;

    /* Check Inputs */
    if (extracta == NULL || extractb == NULL) return NULL ;

    /* Initialise */
    scale = hdrl_spectrum1D_wave_scale_linear;
    extractc = cpl_table_duplicate(extracta) ;

    /* Get the column names */
    col_names = cpl_table_get_column_names(extractc);
    ncols = cpl_table_get_ncol(extractc) ;

    /* Loop on the columns */
    for (i=0 ; i<ncols ; i++) {
        col_name = cpl_array_get_string(col_names, i);
        col_type = cr2res_dfs_SPEC_colname_parse(col_name, &order,
                &trace_nb) ;
        if (col_type != NULL && !strcmp(col_type, CR2RES_COL_SPEC_SUFFIX)) {
            /* This is a SPEC column, let's interpolate */

            /* Get the wavelength column name */
            wave_col = cr2res_dfs_WAVELENGTH_colname(order,trace_nb) ;
            err_col = cr2res_dfs_SPEC_ERR_colname(order,trace_nb) ;

            /* Get the A spectrum */
            a_spec = hdrl_spectrum1D_convert_from_table(extracta,
                    col_name, wave_col, err_col, NULL, scale);

            /* Get the B spectrum */
            b_spec = hdrl_spectrum1D_convert_from_table(extractb,
                    col_name, wave_col, err_col, NULL, scale);

            /* Check */
            if (a_spec == NULL || b_spec == NULL) {
                cpl_msg_error(__func__, "Cannot create HDRL spectra - abort") ;
                cpl_free(wave_col) ;
                cpl_free(err_col);
                if (col_type != NULL) cpl_free(col_type) ;
                cpl_array_delete(col_names) ;
                cpl_table_delete(extractc) ;
                return NULL ;
            }
            sz_a = hdrl_spectrum1D_get_size(a_spec) ;
            sz_b = hdrl_spectrum1D_get_size(b_spec) ;
            /* Check */
            if (sz_a != sz_b) {
                cpl_msg_error(__func__, 
                        "a and b spectra have diff.  sizes - abort") ;
                hdrl_spectrum1D_delete(&a_spec);
                hdrl_spectrum1D_delete(&b_spec);
                cpl_free(wave_col) ;
                cpl_free(err_col);
                if (col_type != NULL) cpl_free(col_type) ;
                cpl_array_delete(col_names) ;
                cpl_table_delete(extractc) ;
                return NULL ;
            }
   
            /* Check if the Wavelenghs are increasing */
            increasing_values = 1 ;
            for (j=1 ; j<sz_a ; j++) {
                hdrl_data_t vala1 =
                    hdrl_spectrum1D_get_wavelength_value(a_spec, j-1, NULL) ;
                hdrl_data_t vala2 =
                    hdrl_spectrum1D_get_wavelength_value(a_spec, j, NULL) ;
                hdrl_data_t valb1 =
                    hdrl_spectrum1D_get_wavelength_value(b_spec, j-1, NULL) ;
                hdrl_data_t valb2 =
                    hdrl_spectrum1D_get_wavelength_value(b_spec, j, NULL) ;
                if (vala1 >= vala2 || valb1 >=valb2) {
                    increasing_values = 0 ;
                    break ;
                }
            }
            if (increasing_values == 0) {
                cpl_msg_warning(__func__, 
                        "Column %s - Values should increase - abort", 
                        wave_col) ;
                hdrl_spectrum1D_delete(&a_spec);
                hdrl_spectrum1D_delete(&b_spec);
                cpl_free(wave_col) ;
                cpl_free(err_col);
                if (col_type != NULL) cpl_free(col_type) ;
                /* Set the column to 0 */
                p_flux = cpl_table_get_data_double(extractc, col_name) ;
                sz = cpl_table_get_nrow(extractc) ;
                for (j = 0; j < sz; j++) p_flux[j]= 0.0;
                continue ;
            }

            /* Resample B on A wavelengths */
            spec_wav = hdrl_spectrum1D_get_wavelength(a_spec);
            params = hdrl_spectrum1D_resample_interpolate_parameter_create(
                    hdrl_spectrum1D_interp_akima);
            c_spec = hdrl_spectrum1D_resample(b_spec, &spec_wav, params);
            hdrl_parameter_delete(params);    
            hdrl_spectrum1D_delete(&b_spec);

             /* Check */
            if (c_spec == NULL) {
                cpl_msg_error(__func__,"Cannot resample HDRL spectra - abort") ;
                cpl_free(wave_col) ;
                cpl_free(err_col);
                if (col_type != NULL) cpl_free(col_type) ;
                cpl_array_delete(col_names) ;
                hdrl_spectrum1D_delete(&a_spec);
                cpl_table_delete(extractc) ;
                return NULL ;
            }

            /* Sum up */
            hdrl_spectrum1D_add_spectrum(c_spec, a_spec) ;
            hdrl_spectrum1D_delete(&a_spec);

            /* Update the table with the result */
            // construct the name of the error column based on the spec column
            p_flux = cpl_table_get_data_double(extractc, col_name) ;
            if ((p_err = cpl_table_get_data_double(extractc, err_col))
                == NULL) {
                // Add the column if it doesn't exists
                cpl_error_reset();
                cpl_table_new_column(extractc, err_col, CPL_TYPE_DOUBLE);
                p_err = cpl_table_get_data_double(extractc, err_col) ;
            }
            sz = cpl_table_get_nrow(extractc) ;
            if (sz != hdrl_spectrum1D_get_size(c_spec)) {
                cpl_msg_error(__func__, "Wrong size - abort") ;
                hdrl_spectrum1D_delete(&c_spec);
                cpl_free(wave_col) ;
                cpl_free(err_col);
                if (col_type != NULL) cpl_free(col_type) ;
                cpl_array_delete(col_names) ;
                cpl_table_delete(extractc) ;
                return NULL ;
            }
            for (j = 0; j < sz; j++) {
                p_flux[j] = hdrl_spectrum1D_get_flux_value(c_spec, j, NULL).data;
                p_err[j] = hdrl_spectrum1D_get_flux_value(c_spec, j, NULL).error;
            }
            hdrl_spectrum1D_delete(&c_spec);
            cpl_free(wave_col) ;
            cpl_free(err_col);
        }
        if (col_type != NULL) cpl_free(col_type) ;
    }
    cpl_array_delete(col_names) ;

    return extractc ;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Fit a 2D Polynomial to the interorder regions
  @param    img         Image, Image with the noise and traces to fit (e.g. a 
                        observation)
  @param    trace_wave  Trace Wave Table
  @param    order_x     maximum order of the polynomial in x direction
  @param    order_y     maximum order of the polynomial in y direction
  @return   the fitted polynomial if succesful, NULL on error
 */
/*----------------------------------------------------------------------------*/
cpl_polynomial * cr2res_fit_interorder(
        cpl_image   *   img, 
        cpl_table   *   trace_wave,
        cpl_size        order_x, 
        cpl_size        order_y)
{

    if (img == NULL || trace_wave == NULL) return NULL;
    if (order_x < 0 || order_y < 0) return NULL;

    //Step 1: identify areas inbetween traces

    cpl_size order_idx;
    cpl_size trace1 = 1;
    cpl_size trace2 = 2; // if decker, set to 2, otherwise 1

    const cpl_array *lower;
    const cpl_array *upper;
    cpl_size power = 0;
    cpl_polynomial *poly_lower = NULL; //lower border of image
    cpl_polynomial *poly_upper = NULL;

    int nb_order_idx_values;
    int * order_idx_values = cr2res_trace_get_order_idx_values(trace_wave,
            &nb_order_idx_values);
    int * porder_idx_values;
    int * ptraces;

    double y_lower, y_upper;
    int c = 0; //counter for number of data points

    cpl_matrix *samppos = cpl_matrix_new(2, cpl_image_get_size_x(img) * 
            cpl_image_get_size_y(img));
    cpl_vector *fitvals = cpl_vector_new(cpl_image_get_size_x(img) * 
            cpl_image_get_size_y(img));
    cpl_vector *x = cpl_vector_new(1);

    double value;
    int pis_rejected = 0;

    // get array for order
    // find location of order and trace
    int nrows = cpl_table_get_nrow(trace_wave);
    int k;
    porder_idx_values = cpl_table_get_data_int(trace_wave, CR2RES_COL_ORDER);
    ptraces = cpl_table_get_data_int(trace_wave, CR2RES_COL_TRACENB);

    for(int m = 0; m <= nb_order_idx_values; m++) {
        // the last time is above the topmost border
        if (m != nb_order_idx_values) order_idx = order_idx_values[m];
        else {
            // find maximum order
            order_idx = order_idx_values[0];
            for(int t = 1; t < nb_order_idx_values; t++) {
                if (order_idx_values[t] > order_idx) 
                    order_idx = order_idx_values[t];
            }
            order_idx ++ ;
        }

        // lower = upper boundary of order_idx-1, or 0 if order_idx-1 = 0
        // upper = lower boundary of order_idx, or max_img of order_idx = max
        for (k=0 ; k<nrows ; k++) {
            /* If order found */
            if (porder_idx_values[k] == order_idx-1 && ptraces[k] == trace2) {
                /* Get the polynomial*/
                lower = cpl_table_get_array(trace_wave, CR2RES_COL_UPPER, k);
                poly_lower = cr2res_convert_array_to_poly(lower);
                break;
            }
        }
        for (k=0 ; k<nrows ; k++) {
            if (porder_idx_values[k] == order_idx && ptraces[k] == trace1) {
                /* Get the polynomial*/
                upper = cpl_table_get_array(trace_wave,CR2RES_COL_LOWER, k);
                poly_upper = cr2res_convert_array_to_poly(upper);
                break;
            }
        }

        // if no order found use bottom of image
        if (poly_lower == NULL) {
            poly_lower = cpl_polynomial_new(1);
            cpl_polynomial_set_coeff(poly_lower, &power, 1);
        }
        // if no order found, use top of image
        if (poly_upper == NULL) {
            poly_upper = cpl_polynomial_new(1);
            cpl_polynomial_set_coeff(poly_upper, &power, 
                    cpl_image_get_size_y(img));
        }

        // loop over image and set data points outside of traces

        for(cpl_size i = 1; i < cpl_image_get_size_x(img)-1; i++) {
            cpl_vector_set(x, 0, (double)i);
            y_lower = cpl_polynomial_eval(poly_lower, x);
            y_upper = cpl_polynomial_eval(poly_upper, x);

            for(cpl_size j = y_lower; j < y_upper; j++) {
                value = cpl_image_get(img, i, j, &pis_rejected);
                if (pis_rejected == 0){
                    cpl_matrix_set(samppos, 0, c, (double)i);
                    cpl_matrix_set(samppos, 1, c, (double)j);

                    cpl_vector_set(fitvals, c, value);
                    c++;
                }
            }
        }
        cpl_polynomial_delete(poly_lower);
        poly_lower = NULL;
        cpl_polynomial_delete(poly_upper);
        poly_upper = NULL;
    }

    // readjust size, to number of data points
    cpl_matrix_set_size(samppos, 2, c);
    cpl_vector_set_size(fitvals, c);

    //Step 2: fit 2d polynomial
    // 2d result polynomial
    cpl_polynomial *fit = cpl_polynomial_new(2); 
    //Matrix of p sample positions, with d rows and p columns
    //const cpl_matrix *samppos = mat; 
    //NULL, or d booleans, true iff the sampling is symmetric
    const cpl_boolean *sampsym = NULL; 
    //cpl_vector *fitvals = vec; //Vector of the p values to fit
    //Uncertainties of the sampled values, or NULL for all ones
    const cpl_vector *fitsigm = NULL; 
    //True iff there is a fitting degree per dimension
    const cpl_boolean dimdeg = TRUE; 
    //Pointer to 1 or d minimum fitting degree(s), or NULL
    const cpl_size * mindeg = NULL; 
    //Pointer to 1 or d maximum fitting degree(s), at least mindeg
    const cpl_size maxdeg[] = {order_x, order_y}; 

    cpl_error_code ec = cpl_polynomial_fit(fit, samppos, sampsym, fitvals, 
            fitsigm, dimdeg, mindeg, maxdeg);

    cpl_free(order_idx_values);
    cpl_matrix_delete(samppos);
    cpl_vector_delete(fitvals);
    cpl_vector_delete(x);
    if (ec == CPL_ERROR_NONE) return fit;
    else {
        cpl_free(fit);
        return NULL;
    }
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Create the polynomials needed to calculate the slit pos and 
            wavelength at any point x, y within an order
  @param    trace_wave input trace wave table with wavelength polynomials for
                       all orders
  @param    coef_slit  [out] slit pos polynomials as a function of the 
                       wavelength and the pixel row y
  @param    coef_wave  [out] wavelength polynomials as a function of the pixel
                       position x, y
  @param    size       [out] the number of orders (polynomials) in coef_slit and
                       coef_wave
  @return   0 on success, -1 on fail
 */
/*----------------------------------------------------------------------------*/
int cr2res_slit_pos(
        const cpl_table *    trace_wave, 
        cpl_polynomial  ***  coef_slit, 
        cpl_polynomial  ***  coef_wave,
        int             *    size)
{

    if (trace_wave == NULL || coef_slit == NULL || 
            coef_wave == NULL || size == NULL) return -1;

    cpl_vector      *  x;
    cpl_matrix      *  matrix_xy;
    cpl_matrix      *  matrix_wd;
    cpl_vector      *  vec_w;
    cpl_vector      *  vec_s;
    cpl_polynomial  *  wave;
    cpl_polynomial  *  line[3];
    cpl_polynomial  *  poly;
    const cpl_array *  slit;
    const cpl_array *  curv_a;
    const cpl_array *  curv_b;
    const cpl_array *  curv_c;
    cpl_polynomial * curv_poly_a;
    cpl_polynomial * curv_poly_b;
    cpl_polynomial * curv_poly_c;
    int             *  order_idx_values;
    int             *  traces;
    const cpl_size maxdeg = 2;
    int i, j, k, n, m, row;
    int order_idx, trace;
    double px, py, pw, ps;
    double pca, pcb, pcc, pc_offset;
    int nb_order_idx_values, nb_traces;
    cpl_error_code errcode;

    // pixels x, only one because thats the same for all traces
    x = cpl_vector_new(CR2RES_DETECTOR_SIZE);
    for (i = 0; i < CR2RES_DETECTOR_SIZE; i++) 
        cpl_vector_set(x, i, i+1);

    order_idx_values = cr2res_trace_get_order_idx_values(trace_wave, 
            &nb_order_idx_values);
    
    // Initialize the new polynomials
    *coef_wave = cpl_malloc(nb_order_idx_values * sizeof(cpl_polynomial*));
    *coef_slit = cpl_malloc(nb_order_idx_values * sizeof(cpl_polynomial*));
    for (i=0; i < nb_order_idx_values; i++) {
        (*coef_wave)[i] = cpl_polynomial_new(2);
        (*coef_slit)[i] = cpl_polynomial_new(2);
    }
    *size = nb_order_idx_values;

    for (i = 0; i < nb_order_idx_values; i++) {
        order_idx = order_idx_values[i];
        // For each trace of this order
        traces = cr2res_get_trace_numbers(trace_wave, order_idx, &nb_traces);

        row = -1;
        matrix_xy = cpl_matrix_new(2, CR2RES_DETECTOR_SIZE * nb_traces * 3);
        matrix_wd = cpl_matrix_new(2, CR2RES_DETECTOR_SIZE * nb_traces * 3);
        vec_w = cpl_vector_new(CR2RES_DETECTOR_SIZE * nb_traces * 3);
        vec_s = cpl_vector_new(CR2RES_DETECTOR_SIZE * nb_traces * 3);

        for (j = 0; j < nb_traces; j++) {
            trace = traces[j];
            k = cr2res_get_trace_table_index(trace_wave, order_idx, trace);
            line[0] = cr2res_get_trace_wave_poly(trace_wave, CR2RES_COL_LOWER, 
                    order_idx, trace);
            line[1] = cr2res_get_trace_wave_poly(trace_wave, CR2RES_COL_ALL, 
                    order_idx, trace);
            line[2] = cr2res_get_trace_wave_poly(trace_wave, CR2RES_COL_UPPER, 
                    order_idx, trace);

            wave = cr2res_get_trace_wave_poly(trace_wave, CR2RES_COL_WAVELENGTH,
                    order_idx, trace);
            slit = cpl_table_get_array(trace_wave, CR2RES_COL_SLIT_FRACTION, k);

            curv_a = cpl_table_get_array(trace_wave, CR2RES_COL_SLIT_CURV_A, k);
            curv_b = cpl_table_get_array(trace_wave, CR2RES_COL_SLIT_CURV_A, k);
            curv_c = cpl_table_get_array(trace_wave, CR2RES_COL_SLIT_CURV_A, k);

            curv_poly_a = cr2res_convert_array_to_poly(curv_a);
            curv_poly_b = cr2res_convert_array_to_poly(curv_b);
            curv_poly_c = cr2res_convert_array_to_poly(curv_c);

            // calculate polynomials for all traces
            for (n = 0; n < CR2RES_DETECTOR_SIZE; n++) {
                // For each of the three edges (upper, all, lower) of a trace
                for (m = 0; m < 3; m++){
                    row++;
                    px = cpl_vector_get(x, n);
                    py = cpl_polynomial_eval_1d(line[m], px, NULL);

                    // Calculate the curvature offset at this point
                    // The parameters are in the global frame so we just
                    // use them directly
                    pca = cpl_polynomial_eval_1d(curv_poly_a, px, NULL);
                    pcb = cpl_polynomial_eval_1d(curv_poly_b, px, NULL);
                    pcc = cpl_polynomial_eval_1d(curv_poly_c, px, NULL);
                    pc_offset = (pca - px) + pcb * py + pcc * py * py;

                    pw = cpl_polynomial_eval_1d(wave, px + pc_offset, NULL);
                    ps = cpl_array_get_double(slit, m, NULL);

                    cpl_matrix_set(matrix_xy, 0, row, px);
                    cpl_matrix_set(matrix_xy, 1, row, py);
                    cpl_matrix_set(matrix_wd, 0, row, pw);
                    cpl_matrix_set(matrix_wd, 1, row, py);
                    cpl_vector_set(vec_w, row, pw);
                    cpl_vector_set(vec_s, row, ps);
                }
            }

            cpl_polynomial_delete(line[0]);
            cpl_polynomial_delete(line[1]);
            cpl_polynomial_delete(line[2]);
            cpl_polynomial_delete(curv_poly_a);
            cpl_polynomial_delete(curv_poly_b);
            cpl_polynomial_delete(curv_poly_c);
            cpl_polynomial_delete(wave);
        }

        // fit 2D wavelengths
        poly = (*coef_wave)[i];
        errcode = cpl_polynomial_fit(poly, matrix_xy, NULL, vec_w, 
                NULL, FALSE, NULL, &maxdeg);
        if (errcode != CPL_ERROR_NONE){
            cpl_msg_error(__func__, "Polynomial fit failed with error: %s", cpl_error_get_message());
            cpl_error_reset();
            for (i=0; i < nb_order_idx_values; i++) {
                cpl_polynomial_delete((*coef_wave)[i]);
                cpl_polynomial_delete((*coef_slit)[i]);
            }
            cpl_free(*coef_wave);
            cpl_free(*coef_slit);
            cpl_matrix_delete(matrix_xy);
            cpl_matrix_delete(matrix_wd);
            cpl_vector_delete(vec_s);
            cpl_vector_delete(vec_w);
            cpl_free(traces);
            cpl_vector_delete(x);
            cpl_free(order_idx_values);
            *coef_wave = NULL;
            *coef_slit = NULL;
            *size = 0;
            return -1;
        }
        
        poly = (*coef_slit)[i];
        errcode = cpl_polynomial_fit(poly, matrix_wd, NULL, vec_s, 
                NULL, FALSE, NULL, &maxdeg);
        if (errcode != CPL_ERROR_NONE){
            cpl_msg_error(__func__, "Polynomial fit failed with error: %s", cpl_error_get_message());
            cpl_error_reset();
            for (i=0; i < nb_order_idx_values; i++) {
                cpl_polynomial_delete((*coef_wave)[i]);
                cpl_polynomial_delete((*coef_slit)[i]);
            }
            cpl_free(*coef_wave);
            cpl_free(*coef_slit);
            cpl_matrix_delete(matrix_xy);
            cpl_matrix_delete(matrix_wd);
            cpl_vector_delete(vec_s);
            cpl_vector_delete(vec_w);
            cpl_free(traces);
            cpl_vector_delete(x);
            cpl_free(order_idx_values);
            *coef_wave = NULL;
            *coef_slit = NULL;
            *size = 0;
            return -1;
        }
        cpl_matrix_delete(matrix_xy);
        cpl_matrix_delete(matrix_wd);
        cpl_vector_delete(vec_s);
        cpl_vector_delete(vec_w);
        cpl_free(traces);
    }

    // delete cpl pointers
    cpl_vector_delete(x);
    cpl_free(order_idx_values);

    // Just check if anything went wrong in general
    if (cpl_error_get_code() != CPL_ERROR_NONE){
        for (i=0; i < nb_order_idx_values; i++) {
                cpl_polynomial_delete((*coef_wave)[i]);
                cpl_polynomial_delete((*coef_slit)[i]);
            }
        cpl_free(*coef_wave);
        cpl_free(*coef_slit);
        *coef_wave = NULL;
        *coef_slit = NULL;
        *size = 0;
        cpl_msg_error(__func__, "ERROR: %s", cpl_error_get_message());
        cpl_error_reset();
        return -1;
    }

    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    get a image of the slitposition (and wavelength) along the slit
  @param    tw_decker_1     tracewave table for traces with decker 1
  @param    tw_decker_2     tracewave_table for traces with decker 2
  @param    slitpos         output image of slit positions
  @param    wavelength      output image of wavelength
  @return   return 0 if successful, -1 otherwise

  Uses the polynomials from cr2res_slit_pos to calculate the slit position and wavelength of each pixel 
 */
/*----------------------------------------------------------------------------*/
int cr2res_slit_pos_image(
        const cpl_table   *   trace_wave, 
        cpl_image   **  slitpos, 
        cpl_image   **  wavelength)
{
    if (trace_wave == NULL || slitpos == NULL || wavelength == NULL) return -1;

    double w, s;
    int i, k, x, y, nb_order_idx_values;
    cpl_vector * vec_xy;
    cpl_vector * vec_wd;
    cpl_polynomial ** coef_slit;
    cpl_polynomial ** coef_wave;

    if (cr2res_slit_pos(trace_wave, &coef_slit, &coef_wave, &nb_order_idx_values)){
        return -1;
    }

    // Create outut images
    *slitpos = cpl_image_new(CR2RES_DETECTOR_SIZE, CR2RES_DETECTOR_SIZE,
                                                         CPL_TYPE_DOUBLE);
    *wavelength = cpl_image_new(CR2RES_DETECTOR_SIZE, CR2RES_DETECTOR_SIZE,
                                                        CPL_TYPE_DOUBLE);

    vec_xy = cpl_vector_new(2);
    vec_wd = cpl_vector_new(2);

    for (k = 0; k < nb_order_idx_values; k++) {
        for (x=1; x <= CR2RES_DETECTOR_SIZE; x++) {
            for (y=1; y <= CR2RES_DETECTOR_SIZE; y++) {
                cpl_vector_set(vec_xy, 0, (double)x);
                cpl_vector_set(vec_xy, 1, (double)y);
                w = cpl_polynomial_eval(coef_wave[k], vec_xy);

                cpl_vector_set(vec_wd, 0, w);
                cpl_vector_set(vec_wd, 1, (double)y);
                s = cpl_polynomial_eval(coef_slit[k], vec_wd);

                if ((s >= 0) && (s <= 1)) {
                    cpl_image_set(*slitpos, x, y, s);
                    cpl_image_set(*wavelength, x, y, w);
                }
            }
        }
    }

    for (i=0; i < nb_order_idx_values; i++){
        cpl_polynomial_delete(coef_wave[i]);
        cpl_polynomial_delete(coef_slit[i]);
    }
    cpl_free(coef_slit);
    cpl_free(coef_wave);
    cpl_vector_delete(vec_wd);
    cpl_vector_delete(vec_xy);
    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Cut a bent order into a rectangle, shifting columns
  @param    img_in
  @param    ycen
  @param    height
  @return   img_out
 */
/*----------------------------------------------------------------------------*/
cpl_image * cr2res_image_cut_rectify(
        const cpl_image     * img_in,
        const cpl_vector    * ycen,
        int                   height)
{
    cpl_image       * img_out;
    cpl_image       * img_1d;
    cpl_type        imtyp;
    cpl_size        lenx, leny;
    int             * ycen_int;
    int             i, j, ymin, ymax;
    int             empty_bottom = 0;

    if (img_in == NULL || ycen == NULL || height < 1) return NULL;

    imtyp = cpl_image_get_type(img_in);
    lenx = cpl_image_get_size_x(img_in);
    leny = cpl_image_get_size_y(img_in);
    ycen_int = cr2res_vector_get_int(ycen);
    img_out = cpl_image_new(lenx, height, imtyp);

    /* Loop over columns, cut out around ycen, insert into img_out*/
    for (i=1;i<=lenx;i++){ // All image indx start at 1!

        /* treat edge cases, summing over shorter column where needed*/
        ymin = ycen_int[i-1]-(height/2);
        ymax = ycen_int[i-1]+(height/2) + height%2 ;
        if ((ymax <= 1) || (ymin > leny)) {
            // Trace is completely out of bounds, skip this column
            // set it to 0 and mark as bad pixels
            for (j = 1; j <= height; j++){
                cpl_image_set(img_out, i, j, NAN);
                cpl_image_reject(img_out, i, j);
            }
            continue;
        }
        if (ymin < 1) {
            empty_bottom = 1 - ymin; // save for later insertion
            ymin = 1;
        }
        if (ymax > leny)
            ymax = leny; // Simply stop when we reach the top.
        if (ymax <= ymin) {
            cpl_msg_error(__func__,"Unreasonable borders in column %i",i);
            cpl_free(ycen_int);
            cpl_image_delete(img_out);
            return NULL;
        }

        /* Cut out and insert */
        img_1d = cpl_image_extract(img_in, i, ymin, i, ymax);
        cpl_image_copy(img_out, img_1d, i, 1+empty_bottom);
        if (cpl_error_get_code() != CPL_ERROR_NONE) {
            cpl_msg_error(__func__,
                    "Cannot extract and copy column %d, %d %d %d, %s",
                    i, ymin, ymax, empty_bottom, cpl_error_get_where());
            cpl_free(ycen_int);
            cpl_image_delete(img_out);
            if (img_1d != NULL) cpl_image_delete(img_1d);
            cpl_error_reset();
            return NULL;
        }
        cpl_image_delete(img_1d);
    }
    cpl_free(ycen_int);
    return img_out;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Re-insert a rectangular cut-out of an order into the full frame
  @param    rect_in
  @param    ycen
  @return   img_out
 */
/*----------------------------------------------------------------------------*/
int cr2res_image_insert_rect(
        const cpl_image     * rect_in,
        const cpl_vector    * ycen,
        cpl_image           * img_out)
{
    cpl_image       * img_1d;
    cpl_size        lenx, leny, height;
    int             * ycen_int;
    int             i, j, ymin, ymax;
    int             empty_bottom;

    if (rect_in == NULL || ycen == NULL || img_out == NULL) return -1;

    lenx = cpl_image_get_size_x(img_out);
    leny = cpl_image_get_size_y(img_out);
    height = cpl_image_get_size_y(rect_in);
    if (cpl_image_get_size_x(rect_in) != lenx) {
        cpl_msg_error(__func__, "Length of rect and img need to be the same");
        return -1;
    }

    ycen_int = cr2res_vector_get_int(ycen);

    for (i=1;i<=lenx;i++){ // All image indices start at 1!
        empty_bottom = 0;
        /* treat edge cases, shorten column where needed*/
        ymin = ycen_int[i-1]-(height/2);
        ymax = ycen_int[i-1]+(height/2) + height%2 ;
        if ((ymax <= 1) || (ymin > leny)) {
            // Trace is completely out of bounds, skip this column
            // set it to 0 and mark as bad pixels 257 1754
            for (j = 1; j <= height; j++){
                cpl_image_set(img_out, i, j, NAN);
                cpl_image_reject(img_out, i, j);
            }
            continue;
        }
        if (ymin < 1) {
            empty_bottom = 1 - ymin; // save for later insertion
            ymin = 1;
        }
        if (ymax > leny)
            ymax = leny; // Simply stop when we reach the top.
        if (ymax <= ymin) {
            cpl_msg_error(__func__, "Unreasonable borders in column %i", i);
            cpl_free(ycen_int);
            return -1;
        }

        img_1d = cpl_image_extract(rect_in, i, empty_bottom+1, i, height);
        cpl_image_copy(img_out, img_1d, i, ymin);
        if (cpl_error_get_code() != CPL_ERROR_NONE) {
            cpl_msg_error(__func__, "Cannot re-insert column %d, %d %d %d, %s",
                            i, ymin, ymax, empty_bottom, cpl_error_get_where());
            cpl_free(ycen_int);
            cpl_error_reset();
            if (img_1d != NULL) cpl_image_delete(img_1d);
            return -1;
        }
        cpl_image_delete(img_1d);
    }
    cpl_free(ycen_int);
    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Evaluate a polynomial on a vector
  @param    poly
  @param    vec
  @return   Vector with evaluation result.
            Caller needs to deallocate.
 */
/*----------------------------------------------------------------------------*/
cpl_vector * cr2res_polynomial_eval_vector(
        const cpl_polynomial * poly,
        const cpl_vector     * vec)
{
    int i;
    cpl_size nx;
    cpl_vector * outvec;

    if (poly == NULL || vec == NULL) return NULL;

    nx = cpl_vector_get_size(vec);
    outvec = cpl_vector_new(nx);
    for (i=0; i<nx; i++){
        cpl_vector_set(outvec, i,
            cpl_polynomial_eval_1d(poly, cpl_vector_get(vec,i), NULL));
    }
    return outvec;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Find the regions with over-average values in a vector
  @param    invector    The vector to be analyzed
  @param    smooth      The size of the boxcar smoothing kernel
  @return   Vector derived as (invector-smoothed_vector - thresh),
            meaning that positive values are at least thresh larger than
            the smoothed vector.
            The returned vector needs to be deallocated by the caller.
 */
/*----------------------------------------------------------------------------*/
cpl_vector * cr2res_threshold_spec(
        const cpl_vector * invector,
        int smooth,
        double thresh)
{
    cpl_vector * smoothed;

    if (invector == NULL || smooth < 0) return NULL;

    smoothed = cpl_vector_filter_median_create(invector, (smooth/2)+1);
    cpl_vector_subtract(smoothed, invector);
    cpl_vector_add_scalar(smoothed, thresh);
    cpl_vector_multiply_scalar(smoothed, -1.0);
    return smoothed;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Find out the base name of a file (i.e. without prefix path)
  @param    filename    Full path name to scan.
  @return   Pointer to char within the input string.
 */
/*----------------------------------------------------------------------------*/
char * cr2res_get_base_name(const char *filename)
{
    char *p ;
    if (filename == NULL) return NULL;

    p = strrchr (filename, '/');
    return p ? p + 1 : (char *) filename;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Find out the root part of a basename (name without extension).
  @param    filename    File name to scan.
  @return   Pointer to statically allocated string.
 */
/*----------------------------------------------------------------------------*/
char * cr2res_get_root_name(const char * filename)
{
    static char path[4096+1];
    char * lastdot ;
    if (filename == NULL) return NULL;

    if (strlen(filename)>4096) return NULL ;
    memset(path, 4096, 0);
    strcpy(path, filename);
    lastdot = strrchr(path, '.');
    if (lastdot == NULL) return path ;
    if ((!strcmp(lastdot, ".fits")) || (!strcmp(lastdot, ".FITS")) ||
        (!strcmp(lastdot, ".dat")) || (!strcmp(lastdot, ".DAT")) ||
        (!strcmp(lastdot, ".paf")) || (!strcmp(lastdot, ".PAF")) ||
        (!strcmp(lastdot, ".txt")) || (!strcmp(lastdot, ".TXT")) ||
        (!strcmp(lastdot, ".ascii")) || (!strcmp(lastdot, ".ASCII")))
    {
        lastdot[0] = (char)0;
    }
    return path ;
}

/*----------------------------------------------------------------------------*/
/**
   @brief   Extract the filename for the first frame of the given tag
   @param   in      A non-empty frameset
   @param   tag     The tag of the requested file
   @return  Pointer to the file
 */
/*----------------------------------------------------------------------------*/
const char * cr2res_extract_filename(
        const cpl_frameset  *   in,
        const char          *   tag)
{
    const cpl_frame *   cur_frame ;

    /* Get the frame  */
    if ((cur_frame = cpl_frameset_find_const(in, tag)) == NULL) return NULL ;
    return cpl_frame_get_filename(cur_frame) ;
}

/*----------------------------------------------------------------------------*/
/**
   @brief   Extract the frames with the given tag from a frameset
   @param   in      A non-empty frameset
   @param   tag     The tag of the requested frames
   @return  The newly created frameset or NULL on error

   The returned frameset must be de allocated with cpl_frameset_delete()
 */
/*----------------------------------------------------------------------------*/
cpl_frameset * cr2res_extract_frameset(
        const cpl_frameset  *   in,
        const char          *   tag)
{
    cpl_frameset    *   out ;
    const cpl_frame *   cur_frame ;
    cpl_frame       *   loc_frame ;
    int                 nbframes;
    int                 i ;

    /* Test entries */
    if (in == NULL) return NULL ;
    if (tag == NULL) return NULL ;

    /* Initialise */
    nbframes = cpl_frameset_get_size(in) ;

    /* Count the frames with the tag */
    if ((cpl_frameset_count_tags(in, tag)) == 0) return NULL ;

    /* Create the output frameset */
    out = cpl_frameset_new() ;

    /* Loop on the requested frames and store them in out */
    for (i=0 ; i<nbframes ; i++) {
        cur_frame = cpl_frameset_get_position_const(in, i) ;
        if (!strcmp(cpl_frame_get_tag(cur_frame), tag)) {
            loc_frame = cpl_frame_duplicate(cur_frame) ;
            cpl_frameset_insert(out, loc_frame) ;
        }
    }
    return out ;
}

/*----------------------------------------------------------------------------*/
/**
   @brief   Extract the frames whose tag is within a list from a frameset
   @param   in      A non-empty frameset
   @param   tags    The list of allowed tags of the requested frames
   @param   ntags   The number of tags
   @return  The newly created frameset or NULL on error

   The returned frameset must be de allocated with cpl_frameset_delete()
 */
/*----------------------------------------------------------------------------*/
cpl_frameset * cr2res_extract_frameset_several_tags(
        const cpl_frameset  *   in,
        const char          **  tags,
        int                     ntags)
{
    cpl_frameset    *   out ;
    const cpl_frame *   cur_frame ;
    cpl_frame       *   loc_frame ;
    int                 nbframes;
    int                 match ;
    int                 i, j ;

    /* Test entries */
    if (in == NULL) return NULL ;
    if (tags == NULL) return NULL ;
    if (ntags < 1) return NULL ;

    /* Initialise */
    nbframes = cpl_frameset_get_size(in) ;

    /* Create the output frameset */
    out = cpl_frameset_new() ;

    /* Loop on the requested frames and store them in out */
    for (i=0 ; i<nbframes ; i++) {
        cur_frame = cpl_frameset_get_position_const(in, i) ;
        match = 0 ;
        for (j=0 ; j<ntags ; j++) {
            if (!strcmp(cpl_frame_get_tag(cur_frame), tags[j])) 
                match = 1 ;
        }
        if (match) {
            loc_frame = cpl_frame_duplicate(cur_frame) ;
            cpl_frameset_insert(out, loc_frame) ;
        }
    }

    /* If no frame is valid, return NULL rather than an empty list */
    if (cpl_frameset_get_size(out) == 0) {
        cpl_frameset_delete(out) ;
        return NULL ;
    }
    return out ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the decker position string for display
  @param    dpos	The decker position
  @return  	the newly allocated string
 */
/*----------------------------------------------------------------------------*/
char * cr2res_decker_print_position(cr2res_decker dpos)
{
    char    *   out ;

    /* Initialise */
    out = NULL ;

    if (dpos == CR2RES_DECKER_INVALID) {
        out = cpl_strdup("INVALID") ;
    } else if (dpos == CR2RES_DECKER_NONE) {
        out = cpl_strdup("NONE") ;
    } else if (dpos == CR2RES_DECKER_1_3) {
        out = cpl_strdup("1_3") ;
    } else if (dpos == CR2RES_DECKER_2_4) {
        out = cpl_strdup("2_4") ;
    } else {
        out = cpl_strdup("Unknown Decker Code") ;
    }
    return out ;
}

/*----------------------------------------------------------------------------*/
/**
   @brief   Convert an array to polynomial
   @param  	arr		An array
   @return  The newly created polynomial or NULL
   The returned object must be de allocated with cpl_polynomial_delete()
 */
/*----------------------------------------------------------------------------*/
cpl_polynomial * cr2res_convert_array_to_poly(const cpl_array * arr)
{
    cpl_polynomial  *   out ;
    double              val ;
    cpl_size            i ;

    /* Test entries */
    if (arr == NULL) return NULL ;

    /* Create Output poly */
	out = cpl_polynomial_new(1) ;

    /* Fill it  */
	for (i=0 ; i<cpl_array_get_size(arr) ; i++) {
        val = cpl_array_get(arr, i, NULL) ;
        if (isnan(val)) {
            cpl_polynomial_delete(out) ;
            return NULL ;
        } 
		cpl_polynomial_set_coeff(out, &i, cpl_array_get(arr, i, NULL)) ;
    }

    return out ;
}

/*----------------------------------------------------------------------------*/
/**
   @brief   Convert a  polynomial to array
   @param  	poly    A polynomial
   @param   size    The requested array size
   @return  The newly created array or NULL
   The returned object must be de allocated with cpl_array_delete()
 */
/*----------------------------------------------------------------------------*/
cpl_array * cr2res_convert_poly_to_array(
        const cpl_polynomial    *   poly,
        int                         size)
{
    cpl_array   *   out ;
    cpl_size        degree, i ;

    /* Test entries */
    if (poly == NULL || size < 1) return NULL ;

    /* Initialise */
    degree = cpl_polynomial_get_degree(poly) ;
                        
    /* Check */
    if (size < degree+1) {
        cpl_msg_error(__func__,
                "The requested array size is too small for the polynomial") ;
        return NULL ;
    }

    /* Create Output array */
	out = cpl_array_new(size, CPL_TYPE_DOUBLE) ;
    cpl_array_fill_window(out, 0, size, 0.0) ;

    /* Fill it  */
    for (i=0 ; i<=degree ; i++) {
        cpl_array_set(out, i, cpl_polynomial_get_coeff(poly, &i)) ;
    }
    return out ;
}

// Remove an element from a vector, the data is modified and resized
int cr2res_vector_erase_element(cpl_vector * vector, cpl_size pos)
{
    cpl_size i, n;

    if (vector == NULL) return -1;

    n = cpl_vector_get_size(vector);
    if (pos >= n || pos < 0) return -1;

    // we shift all elements past pos one step to the left
    for (i = pos; i < n - 1; i++){
        cpl_vector_set(vector, i, cpl_vector_get(vector, i + 1));
    }
    // and then remove the last element in the vector
    cpl_vector_set_size(vector, n-1);
    return 0;
}

// applies the absolute to each element
int cr2res_vector_abs(cpl_vector * vector){
    cpl_size i;

    if (vector == NULL) return 0;

    for (i = 0; i < cpl_vector_get_size(vector); i++){
        cpl_vector_set(vector, i, 
            fabs(cpl_vector_get(vector, i)));
    }
    return 0;
}

/* This function is copied from HDRLDEMO -> should not be changed */
/* It could be added in HDRL */
/*----------------------------------------------------------------------------*/
/**
  @brief   compute photon count error in [ADU]
  @param   ima_data in [ADU]
  @param   gain detector's gain in [e- / ADU]
  @param   ron  detector's read out noise in [ADU]
  @param   ima_errs output error image in [ADU]
  @return  cpl_error_code
  @note ima_errs need to be deallocated
        ima_data must contain the photon counts with no offsets
        this usually means the image must be overscan and bias corrected
        Then the shot noise can be calculated from the poissonian distribution
        as sqrt(electron-counts). To this (transformed back into ADUs) the
        readout noise is added in quadrature.
  @doc
  error is computed with standard formula

  \f$ err_{ADU} = \sqrt{ \frac{ counts }{ gain } + ron^{ 2 } } \f$

  If an image value is negative the associated error is set to RON
 */
/*----------------------------------------------------------------------------*/
cpl_error_code cr2res_detector_shotnoise_model(
        const cpl_image *   ima_data,
        const double        gain,
        const double        ron,
        cpl_image       **  ima_errs)
{
    cpl_ensure_code(ima_data, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(ima_errs, CPL_ERROR_NULL_INPUT);
    cpl_ensure_code(gain > 0., CPL_ERROR_ILLEGAL_INPUT);
    cpl_ensure_code(ron > -1e-5, CPL_ERROR_ILLEGAL_INPUT);

    *ima_errs = cpl_image_duplicate(ima_data);
    /* set negative values (= zero measurable electrons) to read out noise */
    cpl_image_threshold(*ima_errs, 0., DBL_MAX, ron, ron);
    cpl_image_divide_scalar(*ima_errs, gain);
    cpl_image_add_scalar(*ima_errs, ron * ron);
    cpl_image_power(*ima_errs, 0.5);

    return cpl_error_get_code();
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Plot the spectrum with the catalog
  @param    extracted_spec  The extracted spectrum
  @param    catalog         The lines catalog
  @param    title           The title for the plot
  @param    wmin            Min Wavelength to display or -1
  @param    wmax            Max Wavelength to display or -1
  @return   0 if ok
 */
/*----------------------------------------------------------------------------*/
int cr2res_plot_wavecal_result(
        const cpl_bivector      *   extracted_spec,
        const cpl_bivector      *   catalog,
        const char              *   title,
        double                      wmin,
        double                      wmax)
{
    cpl_bivector    **  bivectors ;
    double          *   p0x ;
    double          *   p0y ;
    double          *   p1x ;
    double          *   p1y ;
    double          *   ptmp_biv0x ;
    double          *   ptmp_biv0y ;
    double          *   ptmp_biv1x ;
    double          *   ptmp_biv1y ;
    cpl_bivector    *   tmp_biv0 ;
    cpl_bivector    *   tmp_biv1 ;
    double              rate ;
    int                 n0, n1, i ;

    /* Check entries */
    if (extracted_spec==NULL||catalog==NULL||title==NULL)
        return -1 ;

    /* Initialise */
    n0 = n1 = 0 ;

    /* Create bivectors */
    bivectors = cpl_malloc(2*sizeof(cpl_bivector*)) ;
    bivectors[0] = cpl_bivector_duplicate(extracted_spec) ;
    bivectors[1] = cpl_bivector_duplicate(catalog) ;

    /* Sort the bivectors - Should not be necessary */
    cpl_bivector_sort(bivectors[0], bivectors[0], CPL_SORT_ASCENDING,
            CPL_SORT_BY_X);
    p0x = cpl_bivector_get_x_data(bivectors[0]) ;
    p0y = cpl_bivector_get_y_data(bivectors[0]) ;
    cpl_bivector_sort(bivectors[1], bivectors[1], CPL_SORT_ASCENDING,
            CPL_SORT_BY_X);
    p1x = cpl_bivector_get_x_data(bivectors[1]) ;
    p1y = cpl_bivector_get_y_data(bivectors[1]) ;

    /* Shrink bivectors to wmin-wmax only if requested */
    if (wmin > 0.0 && wmax > 0.0 && wmax > wmin) {

        /* Count the bins falling in the wished interval  */
        for (i=0 ; i<cpl_bivector_get_size(bivectors[0]) ; i++)
            if (p0x[i] >= wmin && p0x[i] <= wmax) n0++ ;
        for (i=0 ; i<cpl_bivector_get_size(bivectors[1]) ; i++)
            if (p1x[i] >= wmin && p1x[i] <= wmax) n1++ ;

        if (n0 > 0 && n1 > 0) {
            /* Allocate the smaller bivectors */
            tmp_biv0 = cpl_bivector_new(n0) ;
            ptmp_biv0x = cpl_bivector_get_x_data(tmp_biv0) ;
            ptmp_biv0y = cpl_bivector_get_y_data(tmp_biv0) ;
            n0 = 0 ;
            for (i=0 ; i<cpl_bivector_get_size(bivectors[0]) ; i++)
                if (p0x[i] >= wmin && p0x[i] <= wmax) {
                    ptmp_biv0x[n0] = p0x[i] ;
                    ptmp_biv0y[n0] = p0y[i] ;
                    n0++ ;
                }
            tmp_biv1 = cpl_bivector_new(n1) ;
            ptmp_biv1x = cpl_bivector_get_x_data(tmp_biv1) ;
            ptmp_biv1y = cpl_bivector_get_y_data(tmp_biv1) ;
            n1 = 0 ;
            for (i=0 ; i<cpl_bivector_get_size(bivectors[1]) ; i++)
                if (p1x[i] >= wmin && p1x[i] <= wmax) {
                    ptmp_biv1x[n1] = p1x[i] ;
                    ptmp_biv1y[n1] = p1y[i] ;
                    n1++ ;
                }
            cpl_bivector_delete(bivectors[0]) ;
            cpl_bivector_delete(bivectors[1]) ;

            bivectors[0] = tmp_biv0 ;
            bivectors[1] = tmp_biv1 ;
        }
    }

    /* Adjust the signal */
    rate = -1* fabs(10*cpl_vector_get_mean(cpl_bivector_get_y(bivectors[0])) /
        cpl_vector_get_mean(cpl_bivector_get_y(bivectors[1]))) ;
    cpl_vector_multiply_scalar(cpl_bivector_get_y(bivectors[1]), rate) ;

    /* Plot */
    if (cpl_bivector_get_size(bivectors[0]) > 0 &&
                cpl_bivector_get_size(bivectors[1]) > 0) {

        char ** options = cpl_malloc(2*sizeof(char*)) ;
        options[0] = cpl_sprintf("t '1-Extracted %s' w lines", title) ;
        options[1] = "t '2-Catalog' w impulses" ;
        cpl_plot_bivectors("set grid;set xlabel 'Wavelength (nm)';",
            (const char **)options,
            "", (const cpl_bivector **)bivectors, 2);
        cpl_free(options[0]) ;
        cpl_free(options) ;
    }

    /* Free */
    cpl_bivector_delete(bivectors[0]) ;
    cpl_bivector_delete(bivectors[1]) ;
    cpl_free(bivectors) ;
    return 0 ;
}

/*
   opt_filter_1d performs optimal filtering of a 1D double array. The mandatory parameters
   are the data array Yarg, the output array Result, the filtering parameter Lam1 and the
   array of Options.
   Options is a 3-element integer array indicating the presence of optional parameters Xarg
   (x argument of the data array), Weights (weights of the data points) and Lam2 (filtering
   parameter for the 2nd derivatives). Xarg, if present, must be sorted in ascending or descending
   order. Absent parameters can by replaced with NULL at the calling.   
*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Apply the optimal filter to a 1D vector
  @param    Yarg            Input vector of type double
  @param    Lam1            Regularization parameter for the 1st derivatives
  @param    Result          Result of the filtering (same size and type as Yarg)
  @param    n               Size of all vectors
  @param    Options         Flags for optional parameters (see below)
  @param    Xarg            X coordinates of Yarg for non-equispaced arrays
  @param    Weights         Weights assigned to data points (same size as Yarg)
  @param    Lam2            Regularization parameter for 2nd derivatives
  @return   CPL_ERROR_NONE  if OK, error code otherwise

  The optimal filter is a smoothing algorithm minimizing the absolut values of the
  1st and (optionally) 2nd derivatives.
 */
/*----------------------------------------------------------------------------*/
int cr2res_util_optimal_filter_1d(
    double     * Yarg,
    double       Lam1,
    double     * Result,
    int          n,
    int          Options[],
    double     * Xarg,
    double     * Weights,
    double       Lam2)
{
   int i, nzero, flag_x, flag_w, flag_lam2;
   double *dx, *ddx2, *aij, *bj;
   double dddd;

   flag_x   =(Options[0])?1:0;  // Flags for optional parameters
   flag_w   =(Options[1])?1:0;
   flag_lam2=(Options[2])?1:0;

   if(flag_x)                          // Xarg is present
   {
     dx=(double *)cpl_malloc((n-1)*sizeof(double));
     nzero=0;
     for(i=0; i<n-1; i++)
     {
       dx[i]=Xarg[i+1]-Xarg[i];
       if(dx[i]>0.)      nzero++;       // Check that Xarg is sorted one way
       else if(dx[i]<0.) nzero--;       // or the other
     }

     if(abs(nzero)!=n-1)                // Xarg is present but not sorted
     {
        cpl_free(dx);
        return 8;
     }

     if(flag_lam2)
     {
        ddx2=(double *)cpl_malloc((n-2)*sizeof(double));
        for(i=0; i<n-2; i++)
        {
          ddx2[i]=(Xarg[i+2]-Xarg[i])*0.5;
          ddx2[i]*=ddx2[i];
        }

        aij=(double *)cpl_malloc(n*5*sizeof(double));
        bj =(double *)cpl_malloc(n*  sizeof(double));
//
// 2nd lower subdiagonal
        aij[0]=aij[1]=0.;
        for(i=0; i<n-2; i++) aij[i+2]=Lam2/(dx[i]*dx[i+1]*ddx2[i]);
//
// Lower subdiagonal
        aij[n]=0.;
        for(i=0; i<n-1; i++) aij[n+i+1] =-Lam1/(dx[i]*dx[i]);
        for(i=0; i<n-2; i++)
        {
          aij[n+i+1]-= Lam2/(dx[i  ]*dx[i  ]*ddx2[i])
                      +Lam2/(dx[i  ]*dx[i+1]*ddx2[i]);
          aij[n+i+2]-= Lam2/(dx[i+1]*dx[i+1]*ddx2[i])
                      +Lam2/(dx[i  ]*dx[i+1]*ddx2[i]);
        }

//
// Main diagonal
        if(flag_w) for(i=0; i<n; i++) aij[2*n+i]=Weights[i];
        else       for(i=0; i<n; i++) aij[2*n+i]=1.;

        for(i=0; i<n-1; i++) aij[2*n+i  ]+=Lam1/(dx[i]*dx[i]);

        for(i=0; i<n-1; i++) aij[2*n+i+1]+=Lam1/(dx[i]*dx[i]);

        for(i=0; i<n-2; i++)
        {
          aij[2*n+i  ]+=Lam2/(dx[i  ]*dx[i  ]*ddx2[i]);
          dddd=1./dx[i]+1./dx[i+1]; dddd*=dddd;
          aij[2*n+i+1]+=Lam2*dddd/ddx2[i];
          aij[2*n+i+2]+=Lam2/(dx[i+1]*dx[i+1]*ddx2[i]);
        }

//
// Upper subdiagonal
        for(i=0; i<n-1; i++) aij[3*n+i]=-Lam1/(dx[i]*dx[i]);
        aij[4*n-1]=0.;

        for(i=0; i<n-2; i++)
        {
         aij[3*n+i  ]-= Lam2/(dx[i  ]*dx[i  ]*ddx2[i])
                       +Lam2/(dx[i  ]*dx[i+1]*ddx2[i]);
         aij[3*n+i+1]-= Lam2/(dx[i+1]*dx[i+1]*ddx2[i])
                       +Lam2/(dx[i  ]*dx[i+1]*ddx2[i]);
        }

//
// 2nd upper subdiagonal
        for(i=0; i<n-2; i++) aij[4*n+i]=Lam2/(dx[i+1]*dx[i]*ddx2[i]);
        aij[5*n-2]=aij[5*n-1]=0.;

//
// RHS
        if(flag_w) for(i=0; i<n; i++) bj[i]=Yarg[i]*Weights[i];
        else       for(i=0; i<n; i++) bj[i]=Yarg[i];

        i=cr2res_extract_slitdec_bandsol(aij, bj, n, 5, Lam1); // Solve the band diagonal SLE
        for(i=0; i<n; i++) Result[i]=bj[i]; // Copy results

        cpl_free(bj);
        cpl_free(aij);
        cpl_free(ddx2);
        cpl_free(dx);
        return CPL_ERROR_NONE;
     }
     else                // No 2nd derivative filtering
     {
        aij=(double *)cpl_malloc(n*3*sizeof(double));
        bj =(double *)cpl_malloc(n*  sizeof(double));
//
// Lower subdiagonal
        aij[0]=0.;
        for(i=0; i<n-1; i++) aij[i+1]=-Lam1/(dx[i]*dx[i]);

//
// Main diagonal
        if(flag_w) for(i=0; i<n; i++) aij[n+i]=Weights[i];
        else       for(i=0; i<n; i++) aij[n+i]=1.;

        for(i=0; i<n-1; i++) aij[n+i  ]+=Lam1/(dx[i]*dx[i]);
        for(i=0; i<n-1; i++) aij[n+i+1]+=Lam1/(dx[i]*dx[i]);

//
// Upper subdiagonal
        for(i=0; i<n-1; i++) aij[2*n+i]=-Lam1/(dx[i]*dx[i]);
        aij[3*n-1]=0.;

//
// RHS
        if(flag_w) for(i=0; i<n; i++) bj[i]=Yarg[i]*Weights[i];
        else       for(i=0; i<n; i++) bj[i]=Yarg[i];

        i=cr2res_extract_slitdec_bandsol(aij, bj, n, 3, Lam1); // Solve the band diagonal SLE
        for(i=0; i<n; i++) Result[i]=bj[i]; // Copy results

        cpl_free(bj);
        cpl_free(aij);
        cpl_free(dx);
        return CPL_ERROR_NONE;
     }
   }
   else          // No Xarg is set
   {
     if(flag_lam2)
     {
        aij=(double *)cpl_malloc(n*5*sizeof(double));
        bj =(double *)cpl_malloc(n*  sizeof(double));
//
// 2nd lower subdiagonal
        aij[0]=aij[1]=0.;
        for(i=0; i<n-2; i++) aij[i+2]=Lam2;

//
// Lower subdiagonal
        aij[n]=0.;
        for(i=0; i<n-1; i++) aij[n+i+1]=-Lam1;
        for(i=0; i<n-2; i++) aij[n+i+1]-=2*Lam2;
        for(i=0; i<n-2; i++) aij[n+i+2]-=2*Lam2;

//
// Main diagonal
        if(flag_w) for(i=0; i<n; i++) aij[2*n+i]=Weights[i];
        else       for(i=0; i<n; i++) aij[2*n+i]=1.;
        for(i=0; i<n-1; i++) aij[2*n+i  ]+=Lam1;
        for(i=0; i<n-1; i++) aij[2*n+i+1]+=Lam1;
        for(i=0; i<n-2; i++) aij[2*n+i  ]+=Lam2;
        for(i=0; i<n-2; i++) aij[2*n+i+1]+=Lam2*4;
        for(i=0; i<n-2; i++) aij[2*n+i+2]+=Lam2;

//
// Upper subdiagonal
        for(i=0; i<n-1; i++) aij[3*n+i]=-Lam1;
        aij[4*n-1]=0.;
        for(i=0; i<n-2; i++) aij[3*n+i  ]-=Lam2*2;
        for(i=0; i<n-2; i++) aij[3*n+i+1]-=Lam2*2;

//
// 2nd lower subdiagonal
        for(i=0; i<n-2; i++) aij[4*n+i]=Lam2;
        aij[5*n-2]=aij[5*n-1]=0.;

//
// RHS
        if(flag_w) for(i=0; i<n; i++) bj[i]=Yarg[i]*Weights[i];
        else       for(i=0; i<n; i++) bj[i]=Yarg[i];

        i=cr2res_extract_slitdec_bandsol(aij, bj, n, 5, Lam1); // Solve the band diagonal SLE
        for(i=0; i<n; i++) Result[i]=bj[i]; // Copy results

        cpl_free(bj);
        cpl_free(aij);
        return CPL_ERROR_NONE;
     }
     else                // No 2nd derivative filtering
     {
        aij=(double *)cpl_malloc(n*3*sizeof(double));
        bj =(double *)cpl_malloc(n*  sizeof(double));
//
// Lower subdiagonal
        aij[0]=0.;
        for(i=0; i<n-1; i++) aij[i+1]=-Lam1;

//
// Main diagonal
        if(flag_w) for(i=0; i<n; i++) aij[n+i]=Weights[i];
        else       for(i=0; i<n; i++) aij[n+i]=1.;

        for(i=0; i<n-1; i++) aij[n+i  ]+=Lam1;
        for(i=0; i<n-1; i++) aij[n+i+1]+=Lam1;

//
// Upper subdiagonal
        for(i=0; i<n-1; i++) aij[2*n+i]=-Lam1;
        aij[3*n-1]=0.;

//
// RHS
        if(flag_w) for(i=0; i<n; i++) bj[i]=Yarg[i]*Weights[i];
        else       for(i=0; i<n; i++) bj[i]=Yarg[i];

        i=cr2res_extract_slitdec_bandsol(aij, bj, n, 3, Lam1); // Solve the band diagonal SLE
        for(i=0; i<n; i++) Result[i]=bj[i]; // Copy filtered vector to Result

        cpl_free(bj);
        cpl_free(aij);
        return CPL_ERROR_NONE;
     }
   }
}

#define aij_index(x, y) ((y) * n) + (x)
#define w_index(x, y) ((y) * nx) + (x)

/*----------------------------------------------------------------------------*/
/**
  @brief    Apply the optimal filter in the 2D case
  @param    img             Input image to filter
  @param    weight          Weights for each point, should be between 0 and 1
  @param    lam_x           Regularization factor in x direction
  @param    lam_y           Regularization factor in y direction
  @return   filtered image if ok, NULL otherwise

  The optimal filter applies a restricition on the 1st derivatives in the x and
  y directions.
 */
/*----------------------------------------------------------------------------*/
cpl_image * cr2res_util_optimal_filter_2d(
    const cpl_image * img,
    const cpl_image * weight,
    double lam_x, 
    double lam_y)
{
    cpl_size nx, ny, i, j, k;
    int n, ndiag ;
    double * aij;
    const double * weight_data;
    cpl_image * rhs_image;
    double * rhs;

    nx = cpl_image_get_size_x(img);
    ny = cpl_image_get_size_y(img);
    n = nx * ny;
    ndiag = 2 * nx + 1;

    lam_x = fabs(lam_x);
    lam_y = fabs(lam_y);

    weight_data = cpl_image_get_data_double_const(weight);
    // aij=dblarr(n, ndiag)
    aij = cpl_malloc(n * ndiag * sizeof(double));
    // aij[0, nx] = weight[0, 0] + lam_x + lam_y
    aij[aij_index(0, nx)] = weight_data[w_index(0, 0)] + lam_x + lam_y;
    // aij[1:nx-1, nx] = weight[1:nx-1, 0] +2 * lam_x + lam_y
    for (i = 1; i < nx; i++)
    {
        aij[aij_index(i, nx)] = weight_data[w_index(i, 0)] + 2 * lam_x + lam_y;
    }
    // aij[nx:n-nx-1,nx]=reform(weight[*,1:ny-2],nx*(ny-2L))+2.d0*(lam_x+lam_y)
    k = nx;
    for (i = 0; i < nx; i++)
    {
        for (j = 1; j < ny-2; j++){
            aij[aij_index(k, nx)] = weight_data[w_index(i, j)] 
                                        + 2 * lam_x + 2 * lam_y;
            k++;
        }
    }
    // aij[n-nx:n-2,nx]= weight[0:nx-2,ny-1] + 2 * lam_x + lam_y
    j = 0;
    for (i = n - nx; i < n - 1; i++)
    {
        aij[aij_index(i, nx)] = weight_data[w_index(j, ny-1)] 
                                    + 2 * lam_x + lam_y;
        j++;
    }
    // aij[n-1,nx] = weight[nx-1,ny-1]+lam_x+lam_y
    aij[aij_index(n-1, nx)] = weight_data[w_index(nx-1, ny-1)] + lam_x + lam_y;
    for (i = 1; i < n; i++)
    {
        // aij[1:n-1, nx-1]=-lam_x
        aij[aij_index(i, nx-1)] = -1 * lam_x;
        // aij[0:n-2, nx+1]=-lam_x
        aij[aij_index(i-1, nx+1)] = -1 * lam_x;
    }

    for (i = 0; i < ny-1; i++)
    {
        // ind = lindgen(ny-1) * nx + nx + nx * n
        j = i * nx + nx + nx * n;
        // aij[ind-1] = aij[ind-1] - lam_x
        aij[j-1] = aij[j - 1] - lam_x;
        // aij[ind  ] = aij[ind] - lam_x
        aij[j] = aij[j] - lam_x;
    }
    for (i = 0; i < ny - 1; i++){
        // ind = lindgen(ny-1) * nx + nx
        j = i * nx + nx;
        // aij[ind-1,nx+1] = 0
        aij[aij_index(j - 1, nx + 1)] = 0;
        // aij[ind, nx-1] = 0
        aij[aij_index(j, nx - 1)] = 0;
    }
    // aij[nx:n-1, 0] = -lam_y
    for (i = nx; i < n; i++)
    {
        aij[aij_index(i, 0)] = -lam_y;
    }
    // aij[0:n-nx-1, ndiag-1] = -lam_y
    for (i = 0; i < n - nx; i++)
    {
        aij[aij_index(i, ndiag - 1)] = -lam_y;
    }

    rhs_image = cpl_image_multiply_create(img, weight);
    rhs = cpl_image_get_data_double(rhs_image);

    cr2res_extract_slitdec_bandsol(aij, rhs, n, ndiag, max(lam_x, lam_y));
    cpl_free(aij);

    return rhs_image;
}

#undef aij_index
#undef w_index

/*

def _scale(x, y):
    # Normalize x and y to avoid huge numbers
    # Mean 0, Variation 1
    offset_x, offset_y = np.mean(x), np.mean(y)
    norm_x, norm_y = np.std(x), np.std(y)
    if norm_x == 0:
        norm_x = 1
    if norm_y == 0:
        norm_y = 1
    x = (x - offset_x) / norm_x
    y = (y - offset_y) / norm_y
    return x, y, (norm_x, norm_y), (offset_x, offset_y)


def _unscale(x, y, norm, offset):
    x = x * norm[0] + offset[0]
    y = y * norm[1] + offset[1]
    return x, y

def polyfit2d(
    x, y, z, degree=1, max_degree=None, scale=True, plot=False, plot_title=None
):
    """A simple 2D plynomial fit to data x, y, z
    The polynomial can be evaluated with numpy.polynomial.polynomial.polyval2d

    Parameters
    ----------
    x : array[n]
        x coordinates
    y : array[n]
        y coordinates
    z : array[n]
        data values
    degree : int, optional
        degree of the polynomial fit (default: 1)
    max_degree : {int, None}, optional
        if given the maximum combined degree of the coefficients is limited to this value
    scale : bool, optional
        Wether to scale the input arrays x and y to mean 0 and variance 1, to avoid numerical overflows.
        Especially useful at higher degrees. (default: True)
    plot : bool, optional
        wether to plot the fitted surface and data (slow) (default: False)

    Returns
    -------
    coeff : array[degree+1, degree+1]
        the polynomial coefficients in numpy 2d format, i.e. coeff[i, j] for x**i * y**j
    """
    # Flatten input
    x = np.asarray(x).ravel()
    y = np.asarray(y).ravel()
    z = np.asarray(z).ravel()

    # Removed masked values
    mask = ~(np.ma.getmask(z) | np.ma.getmask(x) | np.ma.getmask(y))
    x, y, z = x[mask].ravel(), y[mask].ravel(), z[mask].ravel()

    if scale:
        x, y, norm, offset = _scale(x, y)

    # Create combinations of degree of x and y
    # usually: [(0, 0), (1, 0), (0, 1), (1, 1), (2, 0), ....]
    if np.isscalar(degree):
        degree = (int(degree), int(degree))
    assert len(degree) == 2, "Only 2D polynomials can be fitted"
    degree = [int(degree[0]), int(degree[1])]
    # idx = [[i, j] for i, j in product(range(degree[0] + 1), range(degree[1] + 1))]
    coeff = np.zeros((degree[0] + 1, degree[1] + 1))
    idx = _get_coeff_idx(coeff)

    # Calculate elements 1, x, y, x*y, x**2, y**2, ...
    A = polyvander2d(x, y, degree)

    # We only want the combinations with maximum order COMBINED power
    if max_degree is not None:
        mask = idx[:, 0] + idx[:, 1] <= int(max_degree)
        idx = idx[mask]
        A = A[:, mask]

    # Do least squares fit
    C, *_ = lstsq(A, z)

    # Reorder coefficients into numpy compatible 2d array
    for k, (i, j) in enumerate(idx):
        coeff[i, j] = C[k]

    # # Backup copy of coeff
    if scale:
        coeff = polyscale2d(coeff, *norm, copy=False)
        coeff = polyshift2d(coeff, *offset, copy=False)

    if plot:  # pragma: no cover
        if scale:
            x, y = _unscale(x, y, norm, offset)
        plot2d(x, y, z, coeff, title=plot_title)

    return coeff
*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Perform a 2D polynomial fit f(x, y) = z
  @param    x coordinates
  @param    y coordinates
  @param    z values to fit
  @param    degrees to fit with shape (ndegrees, 2) with one entry per 
            xy combination, x degrees in the 1st column, and y degrees in the 2nd,
            Note that you need to include 0, 0 as well for the constant offset
  @return   Fitted polynomail, NULL for error

    Perform a 2D polynomial fit, where the fit degrees are explicitly specified.

 */
/*----------------------------------------------------------------------------*/
cpl_polynomial * cr2res_polyfit_2d_loc(
    const cpl_vector * x, 
    const cpl_vector * y, 
    const cpl_vector * z,
    const cpl_matrix * degree)
{
    cpl_size npoints, ndegrees;
    double offset_x, offset_y;
    double norm_x, norm_y;
    double coef;
    cpl_size i, j;

    cpl_polynomial * poly;
    cpl_vector * xhat, * yhat, *zhat;
    cpl_matrix * mh, * mz;
    cpl_matrix * mcoef;
    cpl_size power[2];

    if (x == NULL || y == NULL || z == NULL || degree == NULL) return NULL;

    npoints = cpl_vector_get_size(x);
    if (cpl_vector_get_size(y) != npoints || cpl_vector_get_size(z) != npoints){
        // All vectors need to be the same size
        cpl_msg_error(__func__, "Vectors in polyfit2d are not of the same size");
        return NULL;
    }
    ndegrees = cpl_matrix_get_nrow(degree);
    if (ndegrees <= 0){
        cpl_msg_error(__func__, "No fit degrees passed to polyfit2d");
        return NULL;
    }

    // Normalize x and y to avoid huge numbers
    // Mean 0, Variation 1
    // offset_x = cpl_vector_get_mean(x);
    // offset_y = cpl_vector_get_mean(y);
    // norm_x = cpl_vector_get_stdev(x);
    // norm_y = cpl_vector_get_stdev(y);
    // We normalize to 0 and 1
    // this only works well when the regular scale is 0 to N
    // which we do have for a detector
    // TODO: figure out how to correct for an offset
    norm_x = cpl_vector_get_max(x);
    norm_y = cpl_vector_get_max(y);
    if (norm_x == 0) norm_x = 1;
    if (norm_y == 0) norm_y = 1;

    xhat = cpl_vector_duplicate(x);
    yhat = cpl_vector_duplicate(y);

    // cpl_vector_subtract_scalar(xhat, offset_x);
    cpl_vector_divide_scalar(xhat, norm_x);

    // cpl_vector_subtract_scalar(yhat, offset_y);
    cpl_vector_divide_scalar(yhat, norm_y);

    // Create the Vandermode Matrix A
    // Calculate elements 1, x, y, x*y, x**2, y**2, ...
    const double * xdata, *ydata;
    double xdegree, ydegree;
    mh = cpl_matrix_new(npoints, ndegrees);
    
    xdata = cpl_vector_get_data_const(xhat);
    ydata = cpl_vector_get_data_const(yhat);


    for (i = 0; i < ndegrees; i++)
    {
        xdegree = cpl_matrix_get(degree, i, 0);
        ydegree = cpl_matrix_get(degree, i, 1);

        for (j = 0; j < npoints; j++)
        {
            cpl_matrix_set(mh, j, i, pow(xdata[j], xdegree) * pow(ydata[j], ydegree));
        }
    }

    // Wrap the z vector in a matrix
    // Copy it to avoid the const classifier
    zhat = cpl_vector_duplicate(z);
    mz = cpl_matrix_wrap(npoints, 1, cpl_vector_get_data(zhat));

    /* Solve XA=B by a least-square solution (aka pseudo-inverse). */
    mcoef = cpl_matrix_solve_svd(mh, mz);

    cpl_matrix_unwrap(mz);
    cpl_matrix_delete(mh);

    // reset the scale
    for (i = 0; i < ndegrees; i++){
        xdegree = cpl_matrix_get(degree, i, 0);
        ydegree = cpl_matrix_get(degree, i, 1);

        coef = cpl_matrix_get(mcoef, i, 0);
        coef /= pow(norm_x, xdegree) * pow(norm_y, ydegree);
        cpl_matrix_set(mcoef, i, 0, coef);
    }


    // turn coefficients into a polynomial
    poly = cpl_polynomial_new(2);
    for (i = 0; i < ndegrees; i++)
    {
        power[0] = cpl_matrix_get(degree, i, 0);
        power[1] = cpl_matrix_get(degree, i, 1);
        coef = cpl_matrix_get(mcoef, i, 0);
        cpl_polynomial_set_coeff(poly, power, coef);
    }

    cpl_vector_delete(xhat);
    cpl_vector_delete(yhat);
    cpl_vector_delete(zhat);
    cpl_matrix_delete(mcoef);
    

    return poly;
}

cpl_polynomial * cr2res_polyfit_2d(
    const cpl_vector * x, 
    const cpl_vector * y,
    const cpl_vector * z,
    const cpl_size degree[2]
){

    cpl_matrix * deg;
    cpl_polynomial * result;
    cpl_size ndegrees = degree[0] * degree[1] + 1;
    cpl_size max_degree = max(degree[0], degree[1]);
    cpl_size i, j, n;

    deg = cpl_matrix_new(ndegrees, 2);

    n = 0;
    for (i = 0; i < degree[0]; i++)
    {
        for (j = 0; j < degree[1]; j++)
        {
            if (i + j > max_degree) continue;
            cpl_matrix_set(deg, n, 0, i);
            cpl_matrix_set(deg, n, 1, j);
            n++;
        }
    }
    
    cpl_matrix_set_size(deg, n, 2);

    result = cr2res_polyfit_2d_loc(x, y, z, deg);
    cpl_matrix_delete(deg);
    return result;
}

#define SECONDS_TO_DAYS 1.157407e-5
/*----------------------------------------------------------------------------*/
/**
  @brief    Calculate the middle of the exposures in the frameset
  @param    frameset - the input frames
  @return   The center of the exposures in MJD

  The calculation is a weighted arithmetic mean of the midpoints in mjd, 
  weigted by the exposure times.
 */
/*----------------------------------------------------------------------------*/
double cr2res_utils_get_center_mjd(const cpl_frameset * frameset)
{
    cpl_size nframes = 0;
    cpl_size i;
    const cpl_frame * frame;
    cpl_propertylist * pl;
    double exptime, mjd;
    double total, total_weight;

    const char * fname;
    const char exptime_name[] = "EXPTIME";
    const char mjd_mid_name[] = "MJD-OBS";

    if (frameset == NULL) return -1;


    total = 0;
    total_weight = 0;

    // Iterate over all frames in the framelist
    nframes = cpl_frameset_get_size(frameset);
    for ( i = 0; i < nframes; i++)
    {
        // load the header from the frame
        frame = cpl_frameset_get_position_const(frameset, i);
        fname = cpl_frame_get_filename(frame);
        pl = cpl_propertylist_load(fname, 0);
        if (pl == NULL){
            cpl_msg_error(__func__, "Could not read file %s", fname);
            cpl_error_reset();
            return -1;
        }

        // load the two properties from the header we care about
        // abort the whole function if we can't find them
        mjd = cpl_propertylist_get_double(pl, mjd_mid_name);
        if (mjd == 0){
            cpl_msg_error(__func__, "Property %s not found in file %s", mjd_mid_name, fname);
            cpl_propertylist_delete(pl);
            cpl_error_reset();
            return -1;
        }
        exptime = cpl_propertylist_get_double(pl, exptime_name);
        if (exptime <= 0){
            cpl_msg_error(__func__, "Property %s not found in file %s", exptime_name, fname);
            cpl_propertylist_delete(pl);
            cpl_error_reset();
            return -1;
        }
        cpl_propertylist_delete(pl);

        // MJD is the start of the observation, but we want the middle
        mjd += exptime / 2.0 * SECONDS_TO_DAYS;

        // weighted mean maths
        total += mjd * exptime;
        total_weight += exptime;
    }

    // normalize by total weights
    total /= total_weight;
    return total;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the pipeline copyright and license
  @return   The copyright and license string

  The function returns a pointer to the statically allocated license string.
  This string should not be modified using the returned pointer.
 */
/*----------------------------------------------------------------------------*/
const char * cr2res_get_license(void)
{
    const char  *   cr2res_license =
        "This file is part of the CR2RES Instrument Pipeline\n"
        "Copyright (C) 2002,2003 European Southern Observatory\n"
        "\n"
        "This program is free software; you can redistribute it and/or modify\n"
        "it under the terms of the GNU General Public License as published by\n"
        "the Free Software Foundation; either version 2 of the License, or\n"
        "(at your option) any later version.\n"
        "\n"
        "This program is distributed in the hope that it will be useful,\n"
        "but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
        "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
        "GNU General Public License for more details.\n"
        "\n"
        "You should have received a copy of the GNU General Public License\n"
        "along with this program; if not, write to the Free Software\n"
        "Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, \n"
        "MA  02111-1307  USA" ;
    return cr2res_license ;
}

/**@}*/
