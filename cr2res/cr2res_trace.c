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
#include "cr2res_dfs.h"
#include "cr2res_trace.h"
#include "cr2res_pfits.h"
#include "cr2res_io.h"
#include "cr2res_wave.h"
#include "cr2res_utils.h"
#include "cr2res_cluster.h"

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static cpl_array * cr2res_trace_get_slit_frac(
        cpl_table       *   traces,
        cpl_size            idx,
        cr2res_decker       decker_pos) ;
static cpl_mask * cr2res_trace_split_traces(
        cpl_mask    *   mask,
        cpl_table   *   trace_table) ;
static cpl_mask * cr2res_trace_signal_detect(
        const cpl_image *   image,
        int                 trace_sep,
        double              smoothfactor,
        double              thresh) ;
static cpl_table * cr2res_trace_fit_traces(
        cpl_table   *   clustertable,
        int             degree) ;
static cpl_array * cr2res_trace_fit_trace(
        cpl_table   *   table,
        int             degree) ;
static cpl_image * cr2res_trace_convert_cluster_to_labels(
        cpl_table   *   cluster,
        int             nx,
        int             ny) ;
static cpl_table * cr2res_trace_convert_labels_to_cluster(cpl_image * labels) ;
static cpl_mask * cr2res_trace_clean_blobs(
        cpl_mask    *   mask,
        int             min_cluster) ;
static int cr2res_trace_extract_edges(
        cpl_table   *   pixels_table,
        cpl_table   **  edge_lower_table,
        cpl_table   **  edge_upper_table) ;

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_trace       Trace Functions
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief  Main function for running all parts of the trace algorithm
  @param    ima             input image
  @param    smoothfactor    Used for detection
  @param    opening         Used for cleaning the mask
  @param    degree          Fitted polynomial degree
  @param    min_cluster     A trace must be bigger - discarded otherwise
  @return The newly allocated trace table or NULL in error case

  A detection is applied to create a mask. This one is labelised.
  The function converts the label image in the proper cluster table in
  trace to call the traces fitting function.
  The cluster table contains the label image information in the form of
  a table. One column per pixel. The columns are xs (pixel x position),
  ys (pixel y position) and cluster (label number).
  The returned table contains 1 line per trace. Each line has 3 polynomials
  (All, Upper and Lower).
    For example with degree 1 :
                 All|               Upper|               Lower|
  24.3593, 0.0161583|  34.6822, 0.0164165|  14.0261, 0.0159084|
  225.479, 0.0167469|  236.604, 0.0168986|  214.342, 0.0166058|
   436.94, 0.0173438|   448.436, 0.017493|   425.423, 0.017203|
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_trace(
        cpl_image       *   ima,
        double              smoothfactor,
        int                 opening,
        int                 degree,
        int                 min_cluster)
{
    cpl_mask        *   mask ;
    cpl_mask        *   mask_clean ;
    cpl_image       *   labels ;
    cpl_apertures   *   aperts ;
    cpl_table       *   clustertable ;
    cpl_table       *   trace_table ;
    cpl_size            nlabels ;

    /* Check Entries */
    if (ima == NULL) return NULL ;

    /* Initialise */

    /* TODO This needs to come from a static calibration, each band */
    /* Alternative: remove altogether and have single kernel size value */
    int                     trace_sep=80;

    /* TODO Should be set to noise level*/
    /* Best to calcule from input image, not relying on external value */
    /* E.g. mask with thresh=0 once, measure noise in non-signal regions */
    /* and do a second pass with this thresh.  */
    double                  thresh=0.0;

    /* Apply detection */
    cpl_msg_info(__func__, "Detect the signal") ;
    if ((mask = cr2res_trace_signal_detect(ima, trace_sep, smoothfactor,
                    thresh)) == NULL) {
        cpl_msg_error(__func__, "Detection failed") ;
        return NULL ;
    }

    if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
        cpl_mask_save(mask, "debug_mask_before_cleaning.fits", NULL,
                    CPL_IO_CREATE);
        }
    /* Clean the traces in the image */
    cpl_msg_info(__func__, "Traces cleaning") ;
    if ((mask_clean = cr2res_trace_clean(mask, opening, min_cluster)) == NULL) {
        cpl_msg_error(__func__, "Cannot clean the traces") ;
        cpl_mask_delete(mask) ;
        return NULL ;
    }
    cpl_mask_delete(mask) ;

    /* Labelization */
    cpl_msg_info(__func__, "Labelise the traces") ;
    if ((labels=cpl_image_labelise_mask_create(mask_clean, &nlabels)) == NULL) {
        cpl_msg_error(__func__, "Cannot labelise") ;
        cpl_mask_delete(mask_clean);
        return NULL ;
    }

    /* Analyse and dump traces */
    if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
        aperts = cpl_apertures_new_from_image(ima, labels);
        cpl_apertures_dump(aperts, stdout) ;
        cpl_apertures_delete(aperts) ;
    }

    /* Create cluster table needed for fitting */
    clustertable = cr2res_trace_convert_labels_to_cluster(labels) ;

    /* Debug Saving */
    if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
		cpl_image_save(labels, "debug_labels.fits",
				CPL_TYPE_INT, NULL, CPL_IO_CREATE);
        cpl_table_save(clustertable, NULL, NULL, "debug_cluster_table.fits",
                CPL_IO_CREATE);
    }
    cpl_image_delete(labels) ;

    /* Alternative with cr2res_cluster_detect() */
    /*
    clustertable = cr2res_cluster_detect(mask, min_cluster) ;
    labels = cr2res_trace_convert_cluster_to_labels(clustertable,
            cpl_image_get_size_x(ima), cpl_image_get_size_y(ima)) ;
    */

    /* Fit the traces */
    cpl_msg_info(__func__, "Fit the trace edges") ;
    if ((trace_table = cr2res_trace_fit_traces(clustertable, degree)) == NULL) {
        cpl_msg_error(__func__, "Failed to Fit") ;
        cpl_table_delete(clustertable);
        cpl_mask_delete(mask_clean);
        return NULL ;
    }
    cpl_table_delete(clustertable);
    cpl_mask_delete(mask_clean);

    /* Debug Saving */
    if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
        cpl_table_save(trace_table, NULL, NULL, "debug_trace_table.fits",
                CPL_IO_CREATE);
    }
    return trace_table ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Clean small blobs
  @param    mask        input mask with small blobs
  @param    opening     Flag to apply opening filtering to the traces
  @param    min_cluster Remove all clusters smaller than this
  @return A newly allocated mask or NULL in error case
 */
/*----------------------------------------------------------------------------*/
cpl_mask * cr2res_trace_clean(
        cpl_mask    *   mask,
        int             opening,
        int             min_cluster)
{
    cpl_mask    *   mask_kernel ;
    cpl_mask    *   new_mask ;
    cpl_mask    *   diff_mask ;
    cpl_mask    *   clean_mask ;

    /* Check entries */
    if (mask == NULL) return NULL ;

    /* Apply a opening to join horizontally the close clusters */
    if (opening) {
        cpl_msg_info(__func__, "Apply Opening to cleanup the traces") ;
        mask_kernel = cpl_mask_new(5, 1) ;
        cpl_mask_not(mask_kernel);
        new_mask = cpl_mask_duplicate(mask) ;
        cpl_mask_filter(new_mask, mask, mask_kernel, CPL_FILTER_OPENING,
                CPL_BORDER_COPY) ;
        cpl_mask_delete(mask_kernel) ;

        /* Compute the difference */
        diff_mask = cpl_mask_duplicate(mask) ;
        cpl_mask_xor(diff_mask, new_mask) ;
        if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
            cpl_mask_save(diff_mask, "debug_diff_mask.fits", NULL,
                    CPL_IO_CREATE);
        }
        cpl_mask_delete(diff_mask) ;
    } else {
        new_mask = cpl_mask_duplicate(mask) ;
    }

    /* Clean the small blobs */
    cpl_msg_info(__func__, "Remove the small blobs (<= %d pixels)",
            min_cluster) ;
    if ((clean_mask = cr2res_trace_clean_blobs(new_mask, min_cluster)) == NULL){
        cpl_msg_error(__func__, "Cannot clean the blobs") ;
        cpl_mask_delete(new_mask) ;
        return NULL ;
    }
    cpl_mask_delete(new_mask) ;

    /* Debug Saving */
    if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
        cpl_mask_save(clean_mask, "debug_mask.fits", NULL, CPL_IO_CREATE);
    }
    return clean_mask ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Make an image out of the trace solution
  @param trace  The trace table
  @param nx     X size of the produced image
  @param ny     Y size of the produced image
  @return   A newly allocated image or NULL in error case
  The returned INT image is of size nx x ny, is filled with -1. The
  polynomials of the different trace edges are used to fill the traces with
  the value of the order.
 */
/*----------------------------------------------------------------------------*/
cpl_image * cr2res_trace_gen_image(
        cpl_table   *   trace,
        int             nx,
        int             ny)
{
    cpl_image       *   out ;
    int             *   pout ;
    const cpl_array *   coeffs_upper ;
    const cpl_array *   coeffs_lower ;
    cpl_polynomial  *   poly_upper ;
    cpl_polynomial  *   poly_lower ;
    int                 order, i ;
    cpl_size            j, k, y_pos_lower, y_pos_upper ;

    /* Check entries */
    if (trace == NULL) return NULL ;
    if (nx < 1 || ny < 1) return NULL ;

    /* Create the empty image */
    out = cpl_image_new(nx, ny, CPL_TYPE_INT) ;
    cpl_image_add_scalar(out, -1.0) ;
    pout = cpl_image_get_data_int(out) ;

    /* Loop on the traces */
    for (i=0 ; i<cpl_table_get_nrow(trace) ; i++) {
        /* Get the Order - Use fix value when no order column */
        if (cpl_table_has_column(trace, CR2RES_COL_ORDER))
            order = cpl_table_get(trace, CR2RES_COL_ORDER, i, NULL) ;
        else
            order = 100 ;

        /* Get the Upper polynomial*/
        coeffs_upper = cpl_table_get_array(trace, CR2RES_COL_UPPER, i) ;
        poly_upper = cr2res_convert_array_to_poly(coeffs_upper) ;

        /* Get the Lower polynomial*/
        coeffs_lower = cpl_table_get_array(trace, CR2RES_COL_LOWER, i) ;
        poly_lower = cr2res_convert_array_to_poly(coeffs_lower) ;

        /* Draw It  */
        for (j=0 ; j<nx ; j++) {
            y_pos_lower = (cpl_size)cpl_polynomial_eval_1d(poly_lower,
                    (double)j+1, NULL) ;
            y_pos_upper = (cpl_size)cpl_polynomial_eval_1d(poly_upper,
                    (double)j+1, NULL) ;
            for (k = y_pos_lower-1 ; k < y_pos_upper ; k++)
                if (k < ny && k >=0)
                    pout[j+k*nx] = order ;

        }
        cpl_polynomial_delete(poly_upper) ;
        cpl_polynomial_delete(poly_lower) ;
    }
    return out;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Count and return the order numbers in a trace table
  @param    trace       trace table
  @param    nb_orders   [output] number of orders
  @return   newly allocated int array

  The int array will need to be freed by the caller. Its size iѕ
  nb_orders. It contains the list of orders found in the trace table.
 */
/*----------------------------------------------------------------------------*/
int * cr2res_trace_get_order_numbers(
        cpl_table   *   trace,
        int         *   nb_orders)
{
    int     *   porders ;
    int         nrows, count_orders, new_order ;
    int     *   tmp_orders_list ;
    int     *   out ;
    int         i, j ;

    /* Check entries */
    if (trace == NULL || nb_orders == NULL) return NULL ;

    /* Initialise */
    nrows = cpl_table_get_nrow(trace) ;
    porders = cpl_table_get_data_int(trace, CR2RES_COL_ORDER);

    /* Allocate orders list */
    tmp_orders_list = cpl_malloc(nrows * sizeof(int)) ;

    /* Count the different orders */
    count_orders = 0 ;
    for (i=0 ; i<nrows ; i++) {
        /* Is the current order a new one ? */
        new_order = 1 ;
        for (j=0 ; j<count_orders ; j++)
            if (tmp_orders_list[j] == porders[i])
                new_order = 0 ;

        /* Current order not yet stored */
        if (new_order) {
            tmp_orders_list[count_orders] = porders[i] ;
            count_orders ++ ;
        }
    }

    /* Allocate and fill output array */
    out = cpl_malloc(count_orders * sizeof(int)) ;
    for (i=0 ; i<count_orders ; i++) out[i] = tmp_orders_list[i] ;

    /* Free and return */
    cpl_free(tmp_orders_list) ;
    *nb_orders = count_orders ;
    return out ;
}

/*----------------------------------------------------------------------------*/
/**
   @brief   Get the number of traces for a specified order
   @param   tab         A TRACE_WAVE table
   @param   order       the order number
   @return  the number or traces in this order or -1 in error case
 */
/*----------------------------------------------------------------------------*/
cpl_size cr2res_get_nb_traces(
        const cpl_table     *   trace_wave,
        int                     order)
{
    cpl_size        nrows, i, count ;

    /* Check Entries */
    if (trace_wave == NULL) return -1 ;

    /* Initialise */
    count = 0 ;
    nrows = cpl_table_get_nrow(trace_wave) ;

    /* Loop on the table rows */
    for (i=0 ; i<nrows ; i++)
        if (cpl_table_get(trace_wave, CR2RES_COL_ORDER,i,NULL)==order)
            count ++ ;
    return count ;
}

/*----------------------------------------------------------------------------*/
/**
   @brief   Get the number of traces for a specified order that have a WL
   @param   tab         A TRACE_WAVE table
   @param   order       the order number
   @return  the number or traces in this order or -1 in error case
   Only Count the traces that have a valid Wavelength solution
 */
/*----------------------------------------------------------------------------*/
cpl_size cr2res_get_nb_traces_with_wavelength(
        const cpl_table     *   trace_wave,
        int                     order)
{
    const cpl_array *   wave_array ;
    cpl_size            nrows, i, count ;

    /* Check Entries */
    if (trace_wave == NULL) return -1 ;

    /* Initialise */
    count = 0 ;
    nrows = cpl_table_get_nrow(trace_wave) ;

    /* Loop on the table rows */
    for (i=0 ; i<nrows ; i++) {
        wave_array = cpl_table_get_array(trace_wave, CR2RES_COL_WAVELENGTH, i) ;
        if (cpl_table_get(trace_wave, CR2RES_COL_ORDER, i, NULL) == order &&
                wave_array != NULL)
            count ++ ;
    }
    return count ;
}

/* TODO */
/*----------------------------------------------------------------------------*/
/**
   @brief   
   @param  
   @return
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_trace_merge(
        const cpl_table     *   trace_wave1,
        const cpl_table     *   trace_wave2)
{
    return NULL ;
}

/*----------------------------------------------------------------------------*/
/**
   @brief   Get the index in a TRACE_WAVE table
   @param   tab         A TRACE_WAVE table
   @param   order       the order number
   @param   trace_nb    the trace number
   @return  the row index or -1 in error case
 */
/*----------------------------------------------------------------------------*/
cpl_size cr2res_get_trace_table_index(
        const cpl_table     *   trace_wave,
        int                     order,
        int                     trace_nb)
{
    cpl_size        nrows, i ;

    /* Check Entries */
    if (trace_wave == NULL) return -1 ;

    /* Initialise */
    nrows = cpl_table_get_nrow(trace_wave) ;

    /* Loop on the table rows */
    for (i=0 ; i<nrows ; i++) {
        if (cpl_table_get(trace_wave, CR2RES_COL_ORDER,i,NULL)==order &&
                cpl_table_get(trace_wave, CR2RES_COL_TRACENB,i,NULL)==trace_nb)
            return i ;
    }
    return -1 ;
}

/*----------------------------------------------------------------------------*/
/**
   @brief   Get the Wavelength polynomial from a TRACE_WAVE table
   @param   tab     A TRACE_WAVE table
   @param   poly_column CR2RES_COL_WAVELENGTH, CR2RES_COL_UPPER,
                        CR2RES_COL_LOWER or CR2RES_COL_ALL
   @return  The newly created polynomial or NULL in error case
   The returned object must be de allocated with cpl_polynomial_delete()
 */
/*----------------------------------------------------------------------------*/
cpl_polynomial * cr2res_get_trace_wave_poly(
        const cpl_table     *   trace_wave,
        const char          *   poly_column,
        int                     order,
        int                     trace_nb)
{
    const cpl_array     *   wave_arr ;
    cpl_polynomial      *   wave_poly ;
    cpl_size                index ;

    /* Check Entries */
    if (trace_wave == NULL) return NULL ;
    if (strcmp(poly_column, CR2RES_COL_WAVELENGTH) &&
            strcmp(poly_column, CR2RES_COL_UPPER) &&
            strcmp(poly_column, CR2RES_COL_LOWER) &&
            strcmp(poly_column, CR2RES_COL_ALL) &&
            strcmp(poly_column, CR2RES_COL_SLIT_CURV_A) &&
            strcmp(poly_column, CR2RES_COL_SLIT_CURV_B) &&
            strcmp(poly_column, CR2RES_COL_SLIT_CURV_C)) return NULL ;

    /* Get Table index from order and trace */
    index = cr2res_get_trace_table_index(trace_wave, order, trace_nb) ;

    /* Read the Table */
    wave_arr = cpl_table_get_array(trace_wave, poly_column, index) ;
    if (wave_arr == NULL) return NULL ;
        
    /* Convert to Polynomial */
    wave_poly = cr2res_convert_array_to_poly(wave_arr) ;

    return wave_poly ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Retrieves the middle (All) polynomial from trace table and evaluates
  @param    trace   TRACE table
  @param    order_nb   Wished order
  @param    trace_nb   Wished trace
  @param    size    Output vector size
  @return
  The returned vector contains the poly evaluation reszult on vector from 1 to
  size. It needs to be destryed by the caller.
 */
/*----------------------------------------------------------------------------*/
cpl_vector * cr2res_trace_get_ycen(
            cpl_table   *   trace,
            cpl_size        order_nb,
            cpl_size        trace_nb,
            int             size)
{
    cpl_vector      *    out ;
    cpl_polynomial  *    poly;
    int                  nrows, i, found  ;
    int             *    porders ;
    int             *    ptraces ;
    const cpl_array *    coeffs;

    /* Initialize */
    /* poly is not necessarly initialized in the for loop !!! */
    /* If not initialized here, it might be used uninitialised  */
    /* by cpl_vector_fill_polynomial() and cause a seg fault */
    poly = NULL ;

    if (trace == NULL || size < 1) return NULL ;

    nrows = cpl_table_get_nrow(trace) ;

    // prevent reading garbage if invalid data in table
    cpl_table_fill_invalid_int(trace, CR2RES_COL_ORDER, -1);
    cpl_table_fill_invalid_int(trace, CR2RES_COL_TRACENB, -1);
    porders = cpl_table_get_data_int(trace, CR2RES_COL_ORDER);
    ptraces = cpl_table_get_data_int(trace, CR2RES_COL_TRACENB);

    /* Loop on the orders */
    for (i=0 ; i<nrows ; i++) {
        /* If order found */
        if (porders[i] == order_nb && ptraces[i] == trace_nb) {
            /* Get the polynomial*/
            coeffs = cpl_table_get_array(trace, CR2RES_COL_ALL, i) ;
            if (coeffs == NULL) {
                cpl_msg_warning(__func__,
                    "Row %d should have our array, but error %d", 
                    i, cpl_error_get_code());
                cpl_error_reset();
                continue;
            }
            poly = cr2res_convert_array_to_poly(coeffs) ;
            break; // The first one is enough.
        }
    }

    // if no order found
    if (poly == NULL){
        cpl_msg_warning(__func__,
                "Cannot find trace %"CPL_SIZE_FORMAT" of order %"CPL_SIZE_FORMAT,
                trace_nb, order_nb);
        return NULL;
    }

    out = cpl_vector_new(size) ;
    if (cpl_vector_fill_polynomial(out, poly, 1, 1) != CPL_ERROR_NONE) {
        cpl_vector_delete(out);
        return NULL;
    };

    cpl_polynomial_delete(poly);

    return out ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the average height (pix) of an order, from trace polys.
  @param    trace   TRACE table
  @param    order_nb   Wished order
  @param    trace_nb   Wished trace
  @return   height in pixels
 */
/*----------------------------------------------------------------------------*/
int cr2res_trace_get_height(
            cpl_table   *   trace,
            cpl_size        order_nb,
            cpl_size        trace_nb)
{
    int                 height;
    cpl_polynomial  *   poly_upper ;
    cpl_polynomial  *   poly_lower ;

    /* Check entries */
    if (trace == NULL) return -1 ;

    /* Get the trace edges */
    poly_upper = cr2res_get_trace_wave_poly(trace, CR2RES_COL_UPPER, order_nb, 
            trace_nb); 
    poly_lower = cr2res_get_trace_wave_poly(trace, CR2RES_COL_LOWER, order_nb, 
            trace_nb); 
    if (poly_upper == NULL || poly_lower == NULL) return -1;

    /* Compute the height */
    height = cr2res_trace_compute_height(poly_upper, poly_lower, 
            CR2RES_DETECTOR_SIZE);

    cpl_polynomial_delete(poly_upper);
    cpl_polynomial_delete(poly_lower);
    return height;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the positions between 2 trace polynomials
  @param    poly1   First trace
  @param    poly2   Second trace
  @param    size    Output vector size
  @return
  The returned vector contains the pixel positions of the middle of the
  2 traces.
  The nth vector value is trace1(n) + trace2(n) / 2
  n=1 for the first value
 */
/*----------------------------------------------------------------------------*/
cpl_vector * cr2res_trace_compute_middle(
        cpl_polynomial  *   trace1,
        cpl_polynomial  *   trace2,
        int                 vector_size)
{
    cpl_vector  *   out ;
    cpl_polynomial * tmp;

    if (trace1 == NULL || trace2 == NULL || vector_size < 1){
        return NULL;
    }

    out = cpl_vector_new(vector_size) ;
    tmp =  cpl_polynomial_new(1);
    cpl_polynomial_add(tmp, trace1, trace2);
    cpl_polynomial_multiply_scalar(tmp, tmp, 0.5);

    cpl_vector_fill_polynomial(out, tmp, 1, 1);
    cpl_polynomial_delete(tmp);

    return out ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes extraction height between 2 trace polynomials
  @param trace1         First trace
  @param trace2         Second trace
  @param vector_size    detector x size
  @return   The average height between 2 polynomials or -1 in error case

  The returned int is the rounded-up mean difference between the two
  input polynomials, evaluated on a vector from 1 to vector_size.
 */
/*----------------------------------------------------------------------------*/
int cr2res_trace_compute_height(
        cpl_polynomial  *   trace1,
        cpl_polynomial  *   trace2,
        int                 vector_size)
{
    cpl_polynomial  *   diff_poly;
    cpl_vector      *   diff_vec ;
    int                 height ;

    /* Check Entries */
    if (trace1 == NULL || trace2 == NULL || vector_size < 1) return -1 ;

    diff_poly =  cpl_polynomial_new(1);
    diff_vec =  cpl_vector_new(vector_size);
    cpl_polynomial_subtract(diff_poly, trace1, trace2);
    cpl_vector_fill_polynomial(diff_vec, diff_poly, 1, 1);
    height = (int)ceil(fabs( cpl_vector_get_mean(diff_vec) ));
    cpl_polynomial_delete(diff_poly) ;

    if (cpl_vector_get_stdev(diff_vec) > 10){ // TODO: make this not hardcoded?
        cpl_msg_warning(__func__, "Stdev of extraction height is large: %.1f",
                    cpl_vector_get_stdev(diff_vec));
    }
    cpl_vector_delete(diff_vec) ;

    cpl_msg_debug(__func__, "Computed height is %d pix.", height);
    return height;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Compute the y position of the trace
  @param traces     The traces table
  @param idx        The index of the trace row
  @return   The y position of the center of the trace
 */
/*----------------------------------------------------------------------------*/
double cr2res_trace_get_trace_ypos(
        cpl_table   *   traces,
        int             idx)
{
    const cpl_array *   coeffs ;
    cpl_polynomial  *   poly ;
    double              ypos ;
    cpl_size            i ;

    /* Check Entries */
    if (traces == NULL) return -1.0 ;
    if (idx < 0) return -1.0 ;
    if (idx >= cpl_table_get_nrow(traces)) return -1.0 ;

    /* Get the trace polynomial*/
    coeffs = cpl_table_get_array(traces, CR2RES_COL_ALL, idx) ;
    poly = cr2res_convert_array_to_poly(coeffs) ;

    /* Evaluate the central pixel */
    ypos = cpl_polynomial_eval_1d(poly, 1024, NULL) ;
    cpl_polynomial_delete(poly) ;

    return ypos ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief   Add extra columns to the plain trace 
            table
  @param    traces      The plain traces table
  @param    infile      File used for WL infos and Decker position
  @param    det_nr      Detector
  @return   0 if ok
  
  infile is the input file from which we use the header to get the
  Wavelength and the Decker position.

  The added columns are:
    CR2RES_COL_ORDER
    CR2RES_COL_TRACENB
    CR2RES_COL_WAVELENGTH
    CR2RES_COL_SLIT_CURV_A
    CR2RES_COL_SLIT_CURV_B
    CR2RES_COL_SLIT_CURV_C
    CR2RES_COL_SLIT_FRACTION
  The values are set as defaults
  The resulting table is a proper TRACE_WAVE table
 */
/*----------------------------------------------------------------------------*/
int cr2res_trace_add_extra_columns(
        cpl_table           *   traces,
        const char          *   infile,
        int                     det_nr)
{
    cpl_propertylist    *   ext_plist ;
    cpl_propertylist    *   main_plist ;
    int                 *   orders ;
    cpl_array           *   wl_array ;
    cpl_array           *   array_id ;
    cpl_array           *   array_zero ;
    cpl_array           *   array_neg ;
    cpl_array           *   slit_frac ;
    double                  y_pos ;
    cr2res_decker           decker_pos ;
    int                     order, trace_nb, nb_orders, trace_id, i, j ;

    if (traces == NULL) return -1;

    /* Load the Plist */
    main_plist = cpl_propertylist_load(infile, 0) ;
    ext_plist = cpl_propertylist_load(infile,
            cr2res_io_get_ext_idx(infile, det_nr, 1)) ;
    if (ext_plist == NULL || main_plist == 0) {
        cpl_error_reset();
        return -1;
    }
 
    /* Get the decker position */
    decker_pos = cr2res_pfits_get_decker_position(main_plist) ;
    cpl_propertylist_delete(main_plist) ;

    /* Add The Order column using the header */
    cpl_table_new_column(traces, CR2RES_COL_ORDER, CPL_TYPE_INT) ;

    /* Loop on the traces */
    for (i=0 ; i<cpl_table_get_nrow(traces) ; i++) {
        /* Get the current trace Y position */
        y_pos = cr2res_trace_get_trace_ypos(traces, i) ;

        /* Compute the trace order from the header */
        order = cr2res_pfits_get_order(ext_plist, y_pos) ;

        /* Store the Order in the table */
        cpl_table_set(traces, CR2RES_COL_ORDER, i, order);
    }
    cpl_propertylist_delete(ext_plist) ;

    /* Add The TraceNb column */
    cpl_table_new_column(traces, CR2RES_COL_TRACENB, CPL_TYPE_INT);

    orders = cr2res_trace_get_order_numbers(traces, &nb_orders) ;
    for (i=0 ; i<nb_orders ; i++) {
        /* Initialise */
        trace_nb = 1 ;
        /* Loop on the traces */
        for (j=0 ; j<cpl_table_get_nrow(traces) ; j++) {
            if (cpl_table_get(traces, CR2RES_COL_ORDER, j, NULL) == orders[i]) {
                cpl_table_set(traces, CR2RES_COL_TRACENB, j, trace_nb);
                trace_nb ++ ;
            }
        }
    }
    cpl_free(orders) ;

    /* Add The Wavelength column using the header */
    cpl_table_new_column_array(traces, CR2RES_COL_WAVELENGTH,CPL_TYPE_DOUBLE,2);

    /* Loop on the traces */
    for (i=0 ; i<cpl_table_get_nrow(traces) ; i++) {
        /* Get the Order number */
        order = cpl_table_get(traces, CR2RES_COL_ORDER, i, NULL) ;
        trace_id = cpl_table_get(traces, CR2RES_COL_TRACENB, i, NULL) ;

        /* Get the Wavelength estimates from the header */
        if ((wl_array = cr2res_wave_get_estimate(infile, det_nr,
                        order)) == NULL) {
            cpl_msg_warning(__func__,
                    "No Wavelength estimate for Detector %d / order %d",
                    det_nr, order) ;
            cpl_error_reset() ;
            cpl_table_set_array(traces, CR2RES_COL_WAVELENGTH, i, NULL);
        } else {
            /* Store the Wavelength in the table */
            cpl_table_set_array(traces, CR2RES_COL_WAVELENGTH, i, wl_array);
            cpl_array_delete(wl_array) ;
        }
    }

    /* Add the Slit Curvature Coefficients                                  */
    /* A(X) = X, B(X) = 0, C(X) = 0, X:1->2048                              */
    /* This means that for X=1024, the Slit curvature polynomial is:        */
    /*           x = A(1024) + B(1024) * y + C(1024) *y^2                   */
    /*      i.e. x = 1024                                                   */
    /*      where x, y are the detector coordinates, (1,1) as Lower left    */
    cpl_table_new_column_array(traces,CR2RES_COL_SLIT_CURV_A,CPL_TYPE_DOUBLE,3);
    cpl_table_new_column_array(traces,CR2RES_COL_SLIT_CURV_B,CPL_TYPE_DOUBLE,3);
    cpl_table_new_column_array(traces,CR2RES_COL_SLIT_CURV_C,CPL_TYPE_DOUBLE,3);
    array_id = cpl_array_new(3, CPL_TYPE_DOUBLE) ;
    cpl_array_set(array_id, 0, 0.0) ;
    cpl_array_set(array_id, 1, 1.0) ;
    cpl_array_set(array_id, 2, 0.0) ;
    array_zero = cpl_array_new(3, CPL_TYPE_DOUBLE) ;
    cpl_array_set(array_zero, 0, 0.0) ;
    cpl_array_set(array_zero, 1, 0.0) ;
    cpl_array_set(array_zero, 2, 0.0) ;

    /* Loop on the traces */
    for (i=0 ; i<cpl_table_get_nrow(traces) ; i++) {
        cpl_table_set_array(traces, CR2RES_COL_SLIT_CURV_A, i, array_id);
        cpl_table_set_array(traces, CR2RES_COL_SLIT_CURV_B, i, array_zero);
        cpl_table_set_array(traces, CR2RES_COL_SLIT_CURV_C, i, array_zero);
    }
    cpl_array_delete(array_id) ;
    cpl_array_delete(array_zero) ;

    /* Add the Slit Fraction position */
    /* 0 is the Bottom of the slit, 1 is the Top. */
    /* The 3 numbers correspond to the Lower, All, Upper trace positions */
    cpl_table_new_column_array(traces,CR2RES_COL_SLIT_FRACTION, 
            CPL_TYPE_DOUBLE,3);
    array_neg = cpl_array_new(3, CPL_TYPE_DOUBLE) ;
    cpl_array_set(array_neg, 0, -1.0) ;
    cpl_array_set(array_neg, 1, -1.0) ;
    cpl_array_set(array_neg, 2, -1.0) ;
    for (i=0 ; i<cpl_table_get_nrow(traces) ; i++) {
        /* Get the slit fraction array */
        slit_frac = cr2res_trace_get_slit_frac(traces, i, decker_pos) ;
        if (slit_frac == NULL) {
            cpl_table_set_array(traces, CR2RES_COL_SLIT_FRACTION, i, array_neg);
        } else { 
            cpl_table_set_array(traces, CR2RES_COL_SLIT_FRACTION, i, slit_frac);
            cpl_array_delete(slit_frac) ;
        }
    }
    cpl_array_delete(array_neg) ;
    return 0 ;
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the slit fraction position for a given trace
  @param traces     The trace_wave file
  @param idx        The index of the trace 
  @param decker_pos The decker position
  @return   the slit_frac array or NULL

  The returned array contains 3 numbers between 0 and 1 identifing the low,
  middle and high positions of the trace along the slit.
  0 is for the slit bottom, 1 for the top.
 */
/*----------------------------------------------------------------------------*/
static cpl_array * cr2res_trace_get_slit_frac(
        cpl_table       *   traces,
        cpl_size            idx,
        cr2res_decker       decker_pos)
{
    cpl_array       *   slit_frac ;
    int                 order ;
    cpl_size            nb_traces, other_idx, i ;
    double              center_pos_curr, center_pos_other ;

    /* Get the current trace order */
    order = cpl_table_get(traces, CR2RES_COL_ORDER, idx, NULL) ;

    /* Get the number of traces of this order */
    nb_traces = cr2res_get_nb_traces(traces, order) ;

    /* The number of traces and the decker positions need to be consistent */
    if (nb_traces == 1 && decker_pos == CR2RES_DECKER_NONE) {
        /* TRACE_OPEN */
        slit_frac = cpl_array_new(3, CPL_TYPE_DOUBLE) ;
        cpl_array_set(slit_frac, 0, 0.0) ;
        cpl_array_set(slit_frac, 1, 0.5) ;
        cpl_array_set(slit_frac, 2, 1.0) ;
    } else if (nb_traces == 2 && 
            (decker_pos==CR2RES_DECKER_1_3 || decker_pos==CR2RES_DECKER_2_4)) {
        /* DECKER - Identify if the trace is the lower or upper one */
        center_pos_curr = cr2res_trace_get_trace_ypos(traces, idx) ;

        /* Search the other trace */
        other_idx = -1 ;
        for (i=0 ; i<cpl_table_get_nrow(traces) ; i++) {
            if (cpl_table_get(traces, CR2RES_COL_ORDER, i, NULL)==order &&
                    i != idx) other_idx = i ; 
        }
        if (other_idx < 0) return NULL ;

        /* Get the position of the other trace */
        center_pos_other = cr2res_trace_get_trace_ypos(traces, other_idx) ;
        
        /* Write the result */
        slit_frac = cpl_array_new(3, CPL_TYPE_DOUBLE) ;
        if (center_pos_curr < center_pos_other) {
            /* Lower Trace */
            if (decker_pos==CR2RES_DECKER_1_3) {
                /* DECKER_1_3 */
                cpl_array_set(slit_frac, 0, 0.250) ;
                cpl_array_set(slit_frac, 1, 0.375) ;
                cpl_array_set(slit_frac, 2, 0.500) ;
            } else {
                /* DECKER_2_4 */
                cpl_array_set(slit_frac, 0, 0.000) ;
                cpl_array_set(slit_frac, 1, 0.125) ;
                cpl_array_set(slit_frac, 2, 0.250) ;
            }
        } else {
            /* Upper Trace */
            if (decker_pos==CR2RES_DECKER_1_3) {
                /* DECKER_1_3 */
                cpl_array_set(slit_frac, 0, 0.750) ;
                cpl_array_set(slit_frac, 1, 0.875) ;
                cpl_array_set(slit_frac, 2, 1.000) ;
            } else {
                /* DECKER_2_4 */
                cpl_array_set(slit_frac, 0, 0.500) ;
                cpl_array_set(slit_frac, 1, 0.625) ;
                cpl_array_set(slit_frac, 2, 0.750) ;
            }
        }
    } else {
        /* Inconsistent */
        slit_frac = NULL ;
    }
    return slit_frac ;
}

/* TODO */
/*----------------------------------------------------------------------------*/
/**
  @brief
  @param
  @return
 */
/*----------------------------------------------------------------------------*/
static cpl_mask * cr2res_trace_split_traces(
        cpl_mask    *   mask,
        cpl_table   *   trace_table)
{
    cpl_mask   		*   out ;
    cpl_binary      *   pout ;
    cpl_polynomial  *   poly_all ;
    cpl_polynomial  *   poly_upper ;
    cpl_polynomial  *   poly_lower ;
    int                 i, j, upper_y, lower_y, center_y, start_y, stop_y ;
    cpl_size            nx, ny, k, l;

    /* Check entries */
    if (mask == NULL || trace_table == NULL) return NULL ;

    /* Create out label image */
    out = cpl_mask_duplicate(mask) ;
    nx = cpl_mask_get_size_x(out) ;
    ny = cpl_mask_get_size_y(out) ;
    pout = cpl_mask_get_data(out) ;

    /* Loop on the traces */
    for (i=0 ; i<cpl_table_get_nrow(trace_table) ; i++) {
        poly_all = cr2res_convert_array_to_poly(cpl_table_get_array(
                    trace_table, CR2RES_COL_ALL, i)) ;
        poly_upper = cr2res_convert_array_to_poly(cpl_table_get_array(
                    trace_table, CR2RES_COL_UPPER, i)) ;
        poly_lower = cr2res_convert_array_to_poly(cpl_table_get_array(
                    trace_table, CR2RES_COL_LOWER, i)) ;

        /* For each x */
        for (k=0 ; k<nx ; k++) {
            lower_y = cpl_polynomial_eval_1d(poly_lower, k+1, NULL) ;
            upper_y = cpl_polynomial_eval_1d(poly_upper, k+1, NULL) ;
            center_y = cpl_polynomial_eval_1d(poly_all, k+1, NULL) ;

            start_y = (int)(center_y - (upper_y-lower_y)/3.0) ;
            stop_y = (int)(center_y + (upper_y-lower_y)/3.0) ;

            if (start_y < 1) start_y = 1 ;
            if (stop_y > ny) stop_y = ny ;

            for (l=start_y-1 ; l<stop_y ; l++) pout[k+l*nx] = CPL_BINARY_0 ;
        }
        cpl_polynomial_delete(poly_all) ;
        cpl_polynomial_delete(poly_upper) ;
        cpl_polynomial_delete(poly_lower) ;
    }
    return out ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Detect the Traces signal
  @param image          The input image with the traces
  @param trace_sep       The approximate number of pixels between 2 traces
  @param smoothfactor   Mult. factor for the low pass filter kernel size
  @param thresh         The threshold used for detection
  @return   A newly allocated mask or NULL in error case.

  The returned mask identifies the pixels belonging to a trace
  The input image is smoothed, subtracted to the result, and a simple
  thresholding is applied. The smoothing kernel size is
  trace_sep*smoothfactor x 1
 */
/*----------------------------------------------------------------------------*/
static cpl_mask * cr2res_trace_signal_detect(
        const cpl_image *   image,
        int                 trace_sep,
        double              smoothfactor,
        double              thresh)
{
    cpl_image       *   smimage ;
    int                 trace_sep_loc ;
    cpl_matrix      *   kernel ;
    cpl_mask        *   mask ;

    if (image == NULL) return NULL;
    if (trace_sep < 0 || smoothfactor < 0) return NULL;

    /* Prepare the kernel used for median filtering */
    trace_sep_loc = (int) (trace_sep*smoothfactor);
    if (trace_sep_loc % 2 == 0) trace_sep_loc +=1;
    cpl_msg_debug(__func__, "Traces separation: %d", trace_sep_loc);
    kernel = cpl_matrix_new(trace_sep_loc, 1);

    /* kernel should have normalized values */
    cpl_matrix_add_scalar(kernel, 1.0/((double)trace_sep_loc)) ;

    /* Median filtering */
    smimage = cpl_image_duplicate(image);
    if (cpl_image_filter(smimage, image, kernel, CPL_FILTER_LINEAR,
                CPL_BORDER_FILTER) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot filter the image") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        cpl_matrix_delete(kernel);
        cpl_image_delete(smimage) ;
        return NULL ;
    }
    cpl_matrix_delete(kernel);

    /* The pixels we want are the ones with values below -thresh */
    cpl_image_subtract(smimage, image);

    if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
        cpl_image_save(smimage, "debug_smimage.fits", CPL_TYPE_FLOAT, NULL,
                CPL_IO_CREATE);
    }

    mask=cpl_mask_new(cpl_image_get_size_x(image),cpl_image_get_size_y(image));
    cpl_mask_not(mask) ;
    cpl_mask_threshold_image(mask,smimage,-1*thresh,DBL_MAX,CPL_BINARY_0);
    cpl_image_delete(smimage) ;

    return mask ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief Fit a polynomial on the different traces (center and edges)
  @param clustertable   The table holding the traces pixels with their labels
  @param degree         Fitting polynomial degree
  @return   A newly allocated table or NULL in error case

  The function loops over the traces labels. For each of them, it
  identifies the upper and lower edges and fits a polynomial to them. It
  also fits a polynomial using all the pixels of the trace.

  The returned table contains 1 line per trace. Each line has 3 polynomials
  (All, Upper and Lower).
 */
/*----------------------------------------------------------------------------*/
static cpl_table * cr2res_trace_fit_traces(
        cpl_table   *   clustertable,
        int             degree)
{
    cpl_array   *    fitparams_all;
    cpl_array   *    fitparams_upper;
    cpl_array   *    fitparams_lower;
    cpl_table   *    traces_table;
    cpl_table   *    trace_table;
    cpl_table   *    edge_upper_table;
    cpl_table   *    edge_lower_table;
    cpl_size         nclusters_cur;
    int              i, nclusters;

    /* Check entries */
    if (clustertable == NULL) return NULL ;

    /* Create the output table */
    nclusters = cpl_table_get_column_max(clustertable, CR2RES_COL_CLUSTERS);
    traces_table = cpl_table_new(nclusters);
    cpl_table_new_column_array(traces_table, CR2RES_COL_ALL, CPL_TYPE_DOUBLE,
            degree+1) ;
    cpl_table_new_column_array(traces_table, CR2RES_COL_UPPER, CPL_TYPE_DOUBLE,
            degree+1) ;
    cpl_table_new_column_array(traces_table, CR2RES_COL_LOWER, CPL_TYPE_DOUBLE,
            degree+1) ;

    /* Loop on the clusters */
    for (i=1 ; i<=nclusters ; i++) {
        /* Select the pixels of the current cluster */
        nclusters_cur = cpl_table_and_selected_int(clustertable,
                CR2RES_COL_CLUSTERS, CPL_EQUAL_TO, i);
        cpl_msg_debug(__func__, "Cluster %d has %"CPL_SIZE_FORMAT" pixels",
                i, nclusters_cur);

        /* Extract the table with the current trace pixels */
        trace_table = cpl_table_extract_selected(clustertable);

        /* Fit the current trace */
        fitparams_all = cr2res_trace_fit_trace(trace_table, degree);
        cpl_table_set_array(traces_table, "All", i-1, fitparams_all);
        cpl_array_delete(fitparams_all);

        /* Extract the edges of the current trace pixels */
        cr2res_trace_extract_edges(trace_table, &edge_lower_table,
                &edge_upper_table) ;

        /* Fit the upper edge of the current trace */
        fitparams_upper = cr2res_trace_fit_trace(edge_upper_table, degree);
        cpl_table_delete(edge_upper_table);
        cpl_table_set_array(traces_table, CR2RES_COL_UPPER, i-1,
                fitparams_upper);
        cpl_array_delete(fitparams_upper);

        /* Fit the lower edge of the current trace */
        fitparams_lower = cr2res_trace_fit_trace(edge_lower_table, degree);
        cpl_table_delete(edge_lower_table);
        cpl_table_set_array(traces_table, CR2RES_COL_LOWER, i-1,
                fitparams_lower);
        cpl_array_delete(fitparams_lower);

        cpl_table_delete(trace_table);

        /* Reset the selection */
        cpl_table_select_all(clustertable);
    }
    return traces_table;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Simple fitting
  @param table  Table containing all the pixels to fit
  @param degree Fitting polynomial degree
  @return   A newly allocated array or NULL in error case

  The pixels in the input table (columns are xs, ys, clusters) are all
  used for the fitting of a polynomial of degree degree. The clusters
  column is IGNORED.
  If the x range of pixels does not exceed 1500 pixels, a linear fit is
  applied.
  The polynomial coefficients are returned in an array.
 */
/*----------------------------------------------------------------------------*/
static cpl_array * cr2res_trace_fit_trace(
        cpl_table   *   table,
        int             degree)
{
    cpl_matrix      *   x ;
    cpl_vector      *   y ;
    cpl_polynomial  *   poly1 ;
    cpl_array       *   result ;
    int             *   xs;
    int             *   ys;
    int                 x_min, x_max ;
    cpl_size            i, degree_local, n ;

    /* Check Entries */
    if (table == NULL || degree < 0) return NULL ;

    /* Initialise */
    n = cpl_table_get_nrow(table) ;
    degree_local = (cpl_size)degree ;
    x_min = x_max = -1 ;

    /* Create Objects */
    x = cpl_matrix_new(1, n) ;
    y = cpl_vector_new(n) ;

    xs = cpl_table_get_data_int(table, CR2RES_COL_XS);
    ys = cpl_table_get_data_int(table, CR2RES_COL_YS);
    for (i=0 ; i<n ; i++) {
        /* Compute x_min and x_max */
        if (x_min < 0 || xs[i] < x_min) x_min = xs[i] ;
        if (x_max < 0 || xs[i] > x_max) x_max = xs[i] ;

        /* Fill the objects used for fitting */
        cpl_matrix_set(x, 0, i, xs[i]) ;
        cpl_vector_set(y, i, (double)ys[i]) ;
    }

    /* If the xs range is too small, reduce the degree */
    /* This case corresponds to the traces that appear on the image corner */
    if (x_max - x_min < 1500) degree_local = 1 ;

    /* Apply the fit */
    poly1 = cpl_polynomial_new(1) ;
    if (cpl_polynomial_fit(poly1, x, NULL, y, NULL, CPL_FALSE, NULL,
                &degree_local) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot fit the data") ;
        cpl_matrix_delete(x);
        cpl_vector_delete(y);
        cpl_polynomial_delete(poly1);
        return NULL;
    }
    cpl_matrix_delete(x);
    cpl_vector_delete(y);

    /* Store the result */
    result = cr2res_convert_poly_to_array(poly1, degree+1) ;
    cpl_polynomial_delete(poly1);
    return result;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Convert cluster table to labels image
  @param cluster    A cluster table (xs, ys, cluster)
  @param nx         X size of the returned label image
  @param ny         Y size of the returned label image
  @return   A newly allocated INT image or NULL in error case

  A new label image is created from the cluster table. Each entry in the
  cluster table is used to set a pixel in the label image.
 */
/*----------------------------------------------------------------------------*/
static cpl_image * cr2res_trace_convert_cluster_to_labels(
        cpl_table   *   cluster,
        int             nx,
        int             ny)
{
    cpl_image       *   labels ;
    const int       *   pxs ;
    const int       *   pys ;
    const int       *   pclusters ;
    int                 i ;

    if (cluster == NULL || nx < 1 || ny < 1) return NULL;
    /* Create labels image  */
    labels = cpl_image_new(nx, ny, CPL_TYPE_INT);
    pxs = cpl_table_get_data_int_const(cluster, CR2RES_COL_XS) ;
    pys = cpl_table_get_data_int_const(cluster, CR2RES_COL_YS) ;
    pclusters = cpl_table_get_data_int_const(cluster, CR2RES_COL_CLUSTERS) ;
    for (i=0 ; i<cpl_table_get_nrow(cluster) ; i++)
        cpl_image_set(labels, pxs[i], pys[i], pclusters[i]);
    return labels ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Convert Labels image to the cluster table
  @param labels     The label image
  @return   A newly allocated cluster table or NULL in error case

  The cluster table contains the label image information in the form of
  a table. One column per pixel. The columns are xs (pixel x position),
  ys (pixel y position) and cluster (label number).
 */
/*----------------------------------------------------------------------------*/
static cpl_table * cr2res_trace_convert_labels_to_cluster(cpl_image * labels)
{
    cpl_image   *   non_zero_image ;
    const int   *   plabels ;
    cpl_table   *   clustertable ;
    int             nb_table_entries, nx, ny, i, j ;


    if (labels == NULL) return NULL;
    /* Count the number of pixels that are not 0 */
    nb_table_entries = 0 ;
    plabels = cpl_image_get_data_int_const(labels) ;
    nx = cpl_image_get_size_x(labels) ;
    ny = cpl_image_get_size_y(labels) ;
    for (i=0 ; i<nx*ny ; i++)
        if (plabels[i] > 0.5) nb_table_entries ++ ;
    cpl_msg_debug(__func__, "Number of table entries: %d", nb_table_entries) ;

    /* Create the output table */
    clustertable = cpl_table_new(nb_table_entries);
    cpl_table_new_column(clustertable, CR2RES_COL_XS, CPL_TYPE_INT) ;
    cpl_table_new_column(clustertable, CR2RES_COL_YS, CPL_TYPE_INT) ;
    cpl_table_new_column(clustertable, CR2RES_COL_CLUSTERS, CPL_TYPE_INT) ;
    nb_table_entries = 0 ;
    for (j=0 ; j<ny ; j++) {
        for (i=0 ; i<nx ; i++) {
            if (plabels[i+j*nx] > 0.5) {
                cpl_table_set_int(clustertable, CR2RES_COL_XS, nb_table_entries,
                        i+1) ;
                cpl_table_set_int(clustertable, CR2RES_COL_YS, nb_table_entries,
                        j+1) ;
                cpl_table_set_int(clustertable, CR2RES_COL_CLUSTERS,
                        nb_table_entries, plabels[i+j*nx]) ;
                nb_table_entries ++ ;
            }
        }
    }
    return clustertable ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Cleans small size group of pixels from a mask
  @param mask           Input mask
  @param min_cluster    Size of clusters under which they need to be removed
  @return   A newly allocated mask or NULL in error case
 */
/*----------------------------------------------------------------------------*/
static cpl_mask * cr2res_trace_clean_blobs(
        cpl_mask    *   mask,
        int             min_cluster)
{
    cpl_mask    *   new_mask ;
    cpl_image   *   labels ;
    cpl_size        nlabels ;
    cpl_binary  *   pnew_mask ;
    const int   *   plabels ;
    int             i, curr_label, pix_count, npix ;

    /* Check entries */
    if (mask == NULL) return NULL ;

    /* Labelise */
    if ((labels = cpl_image_labelise_mask_create(mask, &nlabels)) == NULL) {
        cpl_msg_error(__func__, "Failed to labelise") ;
        return NULL ;
    }
    if (nlabels > 10000) {
        cpl_msg_error(__func__, "Too many labels resulting from mask") ;
        cpl_image_delete(labels) ;
        return NULL ;
    }
    plabels = cpl_image_get_data_int_const(labels) ;

    /* Number of pixels */
    npix = cpl_mask_get_size_x(mask) * cpl_mask_get_size_y(mask) ;

    /* Create Output mask */
    new_mask=cpl_mask_new(cpl_mask_get_size_x(mask),cpl_mask_get_size_y(mask));
    pnew_mask = cpl_mask_get_data(new_mask) ;

    /* Loop on the labels */
    for (curr_label=1 ; curr_label<=nlabels ; curr_label++) {
        pix_count = 0 ;

        /* Count ocurences of current label */
        for (i=0 ; i<npix ; i++) {
            if (plabels[i] == curr_label) pix_count++ ;
            if (pix_count >= min_cluster) break;
        }

        /* Blob big enough ? */
        if (pix_count >= min_cluster) {
            /* Fill the new mask */
            for (i=0 ; i<npix ; i++)
                if (plabels[i] == curr_label) pnew_mask[i] = CPL_BINARY_1 ;
        }
    }
    cpl_image_delete(labels) ;

    return new_mask ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief  Extracts pixels on the upper and lower edges of a trace
  @param pixels_table       Input cluster table with a single trace
  @param edge_lower_table   [output] Lower edge pixels table
  @param edge_upper_table   [output] Upper edge pixels table
  @return   0 if ok, -1 otherwise

  For each found x in the input cluster table, the min and max y are
  used for adding an entry in the 2 output cluster tables.
 */
/*----------------------------------------------------------------------------*/
static int cr2res_trace_extract_edges(
        cpl_table   *   pixels_table,
        cpl_table   **  edge_lower_table,
        cpl_table   **  edge_upper_table)
{
    cpl_table   *   upper_sel ;
    cpl_table   *   lower_sel ;
    const int   *   pxs ;
    const int   *   pys ;
    int         *   min_y ;
    int         *   max_y ;
    int             max_x, i ;

    /* Check entries */
    if (pixels_table == NULL || edge_lower_table == NULL ||
            edge_upper_table == NULL) return -1 ;

    /* Initialise */
    pxs = cpl_table_get_data_int_const(pixels_table, CR2RES_COL_XS) ;
    pys = cpl_table_get_data_int_const(pixels_table, CR2RES_COL_YS) ;

    /* Get the maximum x position */
    max_x = (int)cpl_table_get_column_max(pixels_table, CR2RES_COL_XS) ;

    /* Allocate  */
    min_y = cpl_malloc(max_x * sizeof(int)) ;
    max_y = cpl_malloc(max_x * sizeof(int)) ;
    for (i=0 ; i<max_x ; i++) min_y[i] = max_y[i] = -1 ;

    /* Loop over the pixels table to compute the edges positions */
    for (i=0 ; i<cpl_table_get_nrow(pixels_table) ; i++) {
        if (pys[i] < min_y[pxs[i]-1] || min_y[pxs[i]-1] < 0)
            min_y[pxs[i]-1] = pys[i] ;
        if (pys[i] > max_y[pxs[i]-1] || max_y[pxs[i]-1] < 0)
            max_y[pxs[i]-1] = pys[i] ;
    }

    /* Upper Edge extraction */
    upper_sel = cpl_table_duplicate(pixels_table) ;
    cpl_table_unselect_all(upper_sel) ;
    pxs = cpl_table_get_data_int_const(upper_sel, CR2RES_COL_XS) ;
    pys = cpl_table_get_data_int_const(upper_sel, CR2RES_COL_YS) ;
    for (i=0 ; i<cpl_table_get_nrow(upper_sel) ; i++)
        if (max_y[pxs[i]-1] == pys[i])
            cpl_table_select_row(upper_sel, i) ;
    cpl_free(max_y) ;
    *edge_upper_table = cpl_table_extract_selected(upper_sel) ;
    cpl_table_delete(upper_sel) ;

    /* Lower Edge extraction */
    lower_sel = cpl_table_duplicate(pixels_table) ;
    cpl_table_unselect_all(lower_sel) ;
    pxs = cpl_table_get_data_int_const(lower_sel, CR2RES_COL_XS) ;
    pys = cpl_table_get_data_int_const(lower_sel, CR2RES_COL_YS) ;
    for (i=0 ; i<cpl_table_get_nrow(lower_sel) ; i++)
        if (min_y[pxs[i]-1] == pys[i])
            cpl_table_select_row(lower_sel, i) ;
    cpl_free(min_y) ;
    *edge_lower_table = cpl_table_extract_selected(lower_sel) ;
    cpl_table_delete(lower_sel) ;

    return CPL_ERROR_NONE ;
}
