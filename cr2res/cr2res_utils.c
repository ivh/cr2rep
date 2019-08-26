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

#include "cr2res_utils.h"
#include "cr2res_pfits.h"
#include "cr2res_dfs.h"
#include "cr2res_extract.h"
#include "cr2res_trace.h"

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_utils     Miscellaneous Utilities
 */
/*----------------------------------------------------------------------------*/

/**@{*/


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
    int         i, lenx, val;

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
  @brief    Fit a 2D Polynomial to the interorder noise
  @param    img         Image, Image with the noise and traces to fit (e.g. a 
                        observation)
  @param    trace_wave  Trace Wave Table
  @param    order_x     maximum order of the polynomial in x direction
  @param    order_y     maximum order of the polynomial in y direction
  @return   the fitted polynomial if succesful, NULL on error
 */
/*----------------------------------------------------------------------------*/
cpl_polynomial * cr2res_fit_noise(
        cpl_image   *   img, 
        cpl_table   *   trace_wave,
        cpl_size        order_x, 
        cpl_size        order_y)
{

    if (img == NULL || trace_wave == NULL) return NULL;
    if (order_x < 0 || order_y < 0) return NULL;

    //Step 1: identify areas inbetween traces

    cpl_size order;
    cpl_size trace1 = 1;
    cpl_size trace2 = 2; // if decker, set to 2, otherwise 1

    const cpl_array *lower;
    const cpl_array *upper;
    cpl_size power = 0;
    cpl_polynomial *poly_lower = NULL; //lower border of image
    cpl_polynomial *poly_upper = NULL;

    int nb_orders;
    int * orders = cr2res_trace_get_order_numbers(trace_wave, &nb_orders);
    int * porders;
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
    porders = cpl_table_get_data_int(trace_wave, CR2RES_COL_ORDER);
    ptraces = cpl_table_get_data_int(trace_wave, CR2RES_COL_TRACENB);

    for(int m = 0; m <= nb_orders; m++) {
        // the last time is above the topmost border
        if (m != nb_orders) order = orders[m];
        else {
            // find maximum order
            order = orders[0];
            for(int t = 1; t < nb_orders; t++) {
                if (orders[t] > order) order = orders[t];
            }
            order = order + 1;
        }

        // lower = upper boundary of order-1, or 0 if order-1 = 0
        // upper = lower boundary of order, or max_img of order = max
        for (k=0 ; k<nrows ; k++) {
            /* If order found */
            if (porders[k] == order-1 && ptraces[k] == trace2) {
                /* Get the polynomial*/
                lower = cpl_table_get_array(trace_wave, CR2RES_COL_UPPER, k);
                poly_lower = cr2res_convert_array_to_poly(lower);
                break;
            }
        }
        for (k=0 ; k<nrows ; k++) {
            if (porders[k] == order && ptraces[k] == trace1) {
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

    cpl_free(orders);
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
  @brief    Get a picture of the slit position (and wavelength?) depend on x, y
  @param
  @param
  @param
  @param
  @return   0 on success, -1 on fail
 */
/*----------------------------------------------------------------------------*/
int cr2res_slit_pos(
        const cpl_table *    trace_wave, 
        cpl_polynomial  ***  coef_slit, 
        cpl_polynomial  ***  coef_wave)
{

    if (trace_wave == NULL || coef_slit == NULL || 
            coef_wave == NULL) return -1;
    if (*coef_slit == NULL || *coef_wave == NULL) return -1;

    // load data

    cpl_vector      *  x;
    cpl_matrix      *  matrix_xy;
    cpl_matrix      *  matrix_wd;
    cpl_vector      *  vec_w;
    cpl_vector      *  vec_s;
    cpl_polynomial  *  wave;
    cpl_polynomial  *  line[3];
    const cpl_array *  slit;
    int             *  orders;
    int             *  traces;
    const cpl_size maxdeg = 2;
    int i, j, k, row;
    int order, trace;
    double px, py, pw, ps;
    int nb_orders, nb_traces;
    cpl_error_code errcode;


    // pixels x, only one because thats the same for all traces
    x = cpl_vector_new(CR2RES_DETECTOR_SIZE);
    for (i = 0; i < CR2RES_DETECTOR_SIZE; i++) 
        cpl_vector_set(x, i, i+1);

    orders = cr2res_trace_get_order_numbers(trace_wave, &nb_orders);
    for (i = 0; i < nb_orders; i++) {
        order = orders[i];
        // For each trace of this order
        traces = cr2res_get_trace_numbers(trace_wave, order, &nb_traces);

        row = -1;
        matrix_xy = cpl_matrix_new(2, CR2RES_DETECTOR_SIZE * nb_traces * 3);
        matrix_wd = cpl_matrix_new(2, CR2RES_DETECTOR_SIZE * nb_traces * 3);
        vec_w = cpl_vector_new(CR2RES_DETECTOR_SIZE * nb_traces * 3);
        vec_s = cpl_vector_new(CR2RES_DETECTOR_SIZE * nb_traces * 3);

        for (j = 0; j < nb_traces; j++) {
            trace = traces[j];
            k = cr2res_get_trace_table_index(trace_wave, order, trace);
            line[0] = cr2res_get_trace_wave_poly(trace_wave, CR2RES_COL_LOWER, order, trace);
            line[1] = cr2res_get_trace_wave_poly(trace_wave, CR2RES_COL_ALL, order, trace);
            line[2] = cr2res_get_trace_wave_poly(trace_wave, CR2RES_COL_UPPER, order, trace);

            wave = cr2res_get_trace_wave_poly(trace_wave, CR2RES_COL_WAVELENGTH, order, trace);
            slit = cpl_table_get_array(trace_wave, CR2RES_COL_SLIT_FRACTION, k);

            // calculate polynomials for all traces
            for (i = 0; i < CR2RES_DETECTOR_SIZE; i++) {
                // For each of the three edges (upper, all, lower) of a trace
                for (k = 0; k < 3; k++){
                    row++;
                    px = cpl_vector_get(x, i);
                    py = cpl_polynomial_eval_1d(line[k], px, NULL);
                    pw = cpl_polynomial_eval_1d(wave, px, NULL);
                    ps = cpl_array_get_double(slit, k, NULL);

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
            cpl_polynomial_delete(wave);
        }

        // fit 2D wavelengths
        errcode = cpl_polynomial_fit((*coef_wave)[j], matrix_xy, NULL, vec_w, 
                NULL, FALSE, NULL, &maxdeg);
        if (errcode != CPL_ERROR_NONE){
            // TODO: What to do in case of error?
            cpl_error_reset();
        }
        errcode = cpl_polynomial_fit((*coef_slit)[j], matrix_wd, NULL, vec_s, 
                NULL, FALSE, NULL, &maxdeg);
        if (errcode != CPL_ERROR_NONE){
            cpl_error_reset();
        }
        cpl_matrix_delete(matrix_xy);
        cpl_matrix_delete(matrix_wd);
        cpl_vector_delete(vec_s);
        cpl_vector_delete(vec_w);
        cpl_free(traces);
    }

    // delete cpl pointers
    cpl_vector_delete(x);
    cpl_free(orders);
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
    if (trace_wave == NULL | slitpos == NULL | 
            wavelength == NULL) return -1;
    if (*slitpos == NULL | *wavelength == NULL) return -1;

    double w, s;
    int i, k, x, y, nb_orders;
    cpl_vector * vec_xy;
    cpl_vector * vec_wd;
    cpl_polynomial ** coef_slit;
    cpl_polynomial ** coef_wave;
    int *orders;

    orders = cr2res_trace_get_order_numbers(trace_wave, &nb_orders);
    cpl_free(orders);

    coef_wave = cpl_malloc(nb_orders * sizeof(cpl_polynomial*));
    coef_slit = cpl_malloc(nb_orders * sizeof(cpl_polynomial*));
    for (i=0; i < nb_orders; i++) {
        coef_wave[i] = cpl_polynomial_new(2);
        coef_slit[i] = cpl_polynomial_new(2);
    }

    if (-1 == cr2res_slit_pos(trace_wave, &coef_slit, &coef_wave)){
        for (i=0; i < nb_orders; i++){
            cpl_polynomial_delete(coef_wave[i]);
            cpl_polynomial_delete(coef_slit[i]);
        }
        cpl_free(coef_wave);
        cpl_free(coef_slit);
        return -1;
    }

    vec_xy = cpl_vector_new(2);
    vec_wd = cpl_vector_new(2);

    for (k = 0; k < nb_orders; k++)
    {
        for (x=1; x <= CR2RES_DETECTOR_SIZE; x++)
        {
            for (y=1; y <= CR2RES_DETECTOR_SIZE; y++)
            {
                cpl_vector_set(vec_xy, 0, (double)x);
                cpl_vector_set(vec_xy, 1, (double)y);
                w = cpl_polynomial_eval(coef_wave[k], vec_xy);

                cpl_vector_set(vec_wd, 0, w);
                cpl_vector_set(vec_wd, 1, (double)y);
                s = cpl_polynomial_eval(coef_slit[k], vec_wd);

                if ((s > 0) && (s < 10))
                {
                    cpl_image_set(*slitpos, x, y, s);
                    cpl_image_set(*wavelength, x, y, w);
                }
            }
        }
    }

    for (i=0; i < nb_orders; i++){
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
    int             i, ymin, ymax;
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
        img_1d = cpl_image_extract(img_in,i,ymin, i, ymax);
        cpl_image_copy(img_out, img_1d, i, 1+empty_bottom);
        if (cpl_error_get_code() != CPL_ERROR_NONE) {
            cpl_msg_error(__func__,
                    "Cannot extract and copy column %d, %d %d %d, %s",
                    i, ymin, ymax, empty_bottom, cpl_error_get_where());
            cpl_free(ycen_int);
            cpl_image_delete(img_out);
            if (img_1d != NULL) cpl_image_delete(img_1d);
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
        cpl_image           * img_out  )
{
    cpl_image       * img_1d;
    cpl_size        lenx, leny, height;
    int             * ycen_int;
    int             i, ymin, ymax;
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

    for (i=1;i<=lenx;i++){ // All image indx start at 1!
        empty_bottom = 0;
        /* treat edge cases, shorten column where needed*/
        ymin = ycen_int[i-1]-(height/2);
        ymax = ycen_int[i-1]+(height/2) + height%2 ;
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
