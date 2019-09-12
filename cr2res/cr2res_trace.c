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

#define any(arr, f) ({  \
    int isError = FALSE;\
    for (cpl_size i = 0; i < cpl_array_get_size(arr); i++) \
    {\
        isError = f(cpl_array_get(arr, i, NULL));\
        if (isError) break;\
    }\
    isError;\
    })

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_trace_apply_shift(
        cpl_table       *   traces,
        double              shift) ;
static double cr2res_trace_compute_shift(
        const cpl_table *   traces1,
        const cpl_table *   traces2) ;
static int cr2res_trace_new_trace(
        const cpl_array     *   slit_fraction_in,
        const cpl_array     **  trace_in,
        const int               nb_traces,
        const cpl_array     *   slit_fraction_wished,
        cpl_array           **  trace_all_new,
        cpl_array           **  trace_upper_new,
        cpl_array           **  trace_lower_new) ;
static int cr2res_trace_check_slit_fraction(const cpl_array * slit_fraction) ;
static cpl_array * cr2res_trace_get_slit_frac(
        const cpl_table *   traces,
        cpl_size            idx,
        cr2res_decker       decker_pos) ;
static cpl_mask * cr2res_trace_signal_detect(
        const cpl_image *   image,
        int                 smooth_x,
        int                 smooth_y,
        double              thresh) ;
static cpl_table * cr2res_trace_fit_traces(
        cpl_table   *   clustertable,
        int             degree) ;
static cpl_array * cr2res_trace_fit_trace(
        cpl_table   *   table,
        int             degree) ;
static cpl_table * cr2res_trace_convert_labels_to_cluster(cpl_image * labels) ;
static cpl_mask * cr2res_trace_clean_blobs(
        cpl_mask    *   mask,
        int             min_cluster) ;
static int cr2res_trace_extract_edges(
        cpl_table   *   pixels_table,
        cpl_table   **  edge_lower_table,
        cpl_table   **  edge_upper_table) ;
static int cr2res_trace_get_subtrace(
        cpl_table   *   trace_wave, 
        double          slit_pos, 
        double          height, 
        int             order,
        cpl_array   **  bottom, 
        cpl_array   **  center, 
        cpl_array   **  top,
        cpl_array   **  fraction, 
        cpl_array   **  wave) ;

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
  @param    smooth_x        Low pass filter kernel size in x
  @param    smooth_y        Low pass filter kernel size in y
  @param    threshold       The threshold used for detection
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
        int                 smooth_x,
        int                 smooth_y,
        double              threshold,
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

    /* Apply detection */
    cpl_msg_info(__func__, "Detect the signal") ;
    if ((mask = cr2res_trace_signal_detect(ima, smooth_x, smooth_y,
                    threshold)) == NULL) {
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
        if (cpl_image_get_max(labels) != 0){
            aperts = cpl_apertures_new_from_image(ima, labels);
            cpl_apertures_dump(aperts, stdout) ;
            cpl_apertures_delete(aperts) ;
        } else {
            cpl_msg_debug(__func__, "No labels found, can not create aperture");
        }
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
    if ((opening != 0) & (opening != 1)) return NULL;

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
        const cpl_table   *   trace,
        int         *   nb_orders)
{
    const int     *   porders ;
    int         nrows, count_orders, new_order ;
    int     *   tmp_orders_list ;
    int     *   out ;
    int         i, j ;

    /* Check entries */
    if (trace == NULL || nb_orders == NULL) return NULL ;

    /* Initialise */
    nrows = cpl_table_get_nrow(trace) ;
    porders = cpl_table_get_data_int_const(trace, CR2RES_COL_ORDER);

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
        if (cpl_table_get(trace_wave, CR2RES_COL_ORDER, i, NULL)==order){
            count ++ ;
        }
    return count ;
}

/*----------------------------------------------------------------------------*/
/**
   @brief   Get the trace numbers for a specified order
   @param   tab         A TRACE_WAVE table
   @param   order       the order number
   @param   nb_traces   [out] number of traces
   @return  the trace numbers, or NULL in error cases

  return value must be freed with cpl_free

 */
/*----------------------------------------------------------------------------*/
int * cr2res_get_trace_numbers(
        const cpl_table * trace_wave,
        int               order,
        int *             nb_traces)
{
    cpl_size nrows, i, k;
    int number_traces;
    int * traces;
    /* Check Entries */
    if (trace_wave == NULL) return NULL ;

    number_traces = cr2res_get_nb_traces(trace_wave, order);
    if (number_traces == -1) return NULL ;


    /* Initialise */
    k = 0;
    nrows = cpl_table_get_nrow(trace_wave) ;
    traces = cpl_malloc(number_traces * sizeof(int));

    /* Loop on the table rows */
    for (i=0 ; i<nrows ; i++)
        if (cpl_table_get(trace_wave, CR2RES_COL_ORDER, i, NULL)==order)
        {
            traces[k] = cpl_table_get(trace_wave, CR2RES_COL_TRACENB, i, NULL);
            k++;
        }

    if (nb_traces != NULL) *nb_traces = number_traces;

    if (k != number_traces){
        cpl_msg_warning(__func__, "Not all expected traces found");
    }

    return traces;
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

/*----------------------------------------------------------------------------*/
/**
   @brief   Merge 2 trace_wave tables 
   @param  trace_wave1  Table to merge
   @param  trace_wave2  Table to merge
   @return  a newly allocated trace_wave table or NULL in error case
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_trace_merge(
        const cpl_table     *   trace_wave1,
        const cpl_table     *   trace_wave2)
{
    cpl_table       *   merged;
    int             *   orders ;
    int                 norders, my_order, j, cur_trace_id ;
    cpl_size            i, nrows ;
    
    /* Check Entries */
    if (trace_wave1 == NULL || trace_wave2 == NULL) return NULL ;
    nrows = cpl_table_get_nrow(trace_wave1) ;

    /* Create the merged table */
    merged = cpl_table_duplicate(trace_wave1) ;

    /* Merge the tables */
    if (cpl_table_insert(merged, trace_wave2, nrows) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Failed merging the tables") ;
        cpl_table_delete(merged) ;
        return NULL ;
    }

    /* Recompute the merged table size */
    nrows = cpl_table_get_nrow(merged) ;

    /* Get the orders information */
	orders = cr2res_trace_get_order_numbers(merged, &norders) ;

    /* Update the trace_id Information */
    /* Loop on the orders */
    for (j=0 ; j<norders ; j++) {
        /* Initialise the trace ID */
        cur_trace_id = 1 ;

        /* Loop on the table */
        for (i=0 ; i<nrows ; i++) {
            /* Each time the curent order is encountered, update the trace id */
            if (cpl_table_get(merged, CR2RES_COL_ORDER, i, NULL)==orders[j]) {
                cpl_table_set(merged, CR2RES_COL_TRACENB, i, cur_trace_id);
                cur_trace_id++ ;
            }
        }
    }
    cpl_free(orders) ;
    return merged ;
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
   @brief   Get a polynomial from a TRACE_WAVE table
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
    if (index == -1) return NULL;

    /* Read the Table */
    wave_arr = cpl_table_get_array(trace_wave, poly_column, index) ;
    if (wave_arr == NULL) return NULL ;
        
    /* Convert to Polynomial */
    wave_poly = cr2res_convert_array_to_poly(wave_arr) ;

    return wave_poly ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief   Get the Wavelength vector from a TRACE_WAVE table
  @param   trace_wave  A TRACE_WAVE table
  @param   order       Wished order
  @param   trace_nb    Wished trace number
  @param   size    Output vector size
  @return  The newly created vector or NULL in error case
  The returned object must be de allocated with cpl_vector_delete()
 */
/*----------------------------------------------------------------------------*/
cpl_vector * cr2res_trace_get_wl(
        const cpl_table *   trace_wave,
        int                 order,
        int                 trace_nb,
        int                 size)
{
    cpl_vector      *    out ;
    cpl_polynomial  *    wl_poly;

    /* Check Inputs */
    if (trace_wave == NULL || size < 1) return NULL ;

    /* Get the Wavelength polynomial */
    wl_poly = cr2res_get_trace_wave_poly(trace_wave, CR2RES_COL_WAVELENGTH, 
            order, trace_nb) ;
    if (wl_poly == NULL) return NULL ;

    /* Allocate output */
    out = cpl_vector_new(size) ;
    if (cpl_vector_fill_polynomial(out, wl_poly, 1, 1) != CPL_ERROR_NONE) {
        cpl_vector_delete(out);
        cpl_polynomial_delete(wl_poly);
        return NULL;
    }
    cpl_polynomial_delete(wl_poly);
    return out ;
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
        const cpl_table *   trace,
        cpl_size            order_nb,
        cpl_size            trace_nb,
        int                 size)
{
    cpl_vector      *    out ;
    cpl_polynomial  *    pos_poly;

    if (trace == NULL || size < 1) return NULL ;

    /* Get the Positions polynomial */
    pos_poly = cr2res_get_trace_wave_poly(trace, CR2RES_COL_ALL, 
            order_nb, trace_nb) ;
    if (pos_poly == NULL) return NULL ;

    /* Allocate output */
    out = cpl_vector_new(size) ;
    if (cpl_vector_fill_polynomial(out, pos_poly, 1, 1) != CPL_ERROR_NONE) {
        cpl_vector_delete(out);
        cpl_polynomial_delete(pos_poly);
        return NULL;
    }
    cpl_polynomial_delete(pos_poly);
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
        const cpl_table *   trace,
        cpl_size            order_nb,
        cpl_size            trace_nb)
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
        const cpl_table *   traces,
        int                 idx)
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
    ypos = cpl_polynomial_eval_1d(poly, CR2RES_DETECTOR_SIZE / 2, NULL) ;
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
    CR2RES_COL_WAVELENGTH_ERROR
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
        cpl_propertylist_delete(main_plist);
        cpl_propertylist_delete(ext_plist);
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

    /* Add The Wavelength(_Error) column using the header */
    cpl_table_new_column_array(traces, CR2RES_COL_WAVELENGTH,
            CPL_TYPE_DOUBLE, 2) ;
    cpl_table_new_column_array(traces, CR2RES_COL_WAVELENGTH_ERROR, 
            CPL_TYPE_DOUBLE, 2) ;

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
            cpl_table_set_array(traces, CR2RES_COL_WAVELENGTH_ERROR, i, NULL);
        } else {
            /* Store the Wavelength in the table */
            cpl_table_set_array(traces, CR2RES_COL_WAVELENGTH, i, wl_array);
            cpl_array_delete(wl_array) ;

            wl_array = cpl_array_new(2, CPL_TYPE_DOUBLE) ;
            cpl_array_set(wl_array, 0, CR2RES_WAVELENGTH_ERROR_DEFAULT) ;
            cpl_array_set(wl_array, 1, CR2RES_WAVELENGTH_ERROR_DEFAULT) ;
            cpl_table_set_array(traces, CR2RES_COL_WAVELENGTH_ERROR,i,wl_array);
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
            cpl_msg_warning(__func__, "Cannot assign slit fraction.");
        } else { 
            cpl_table_set_array(traces, CR2RES_COL_SLIT_FRACTION, i, slit_frac);
            cpl_array_delete(slit_frac) ;
        }
    }
    cpl_array_delete(array_neg) ;
    return 0 ;
}

/* TODO : NOT WORKING - needs review / Fix */
/*----------------------------------------------------------------------------*/
/**
  @brief    Recompute the traces at a newly specified slit fraction
  @param    traces              The input traces
  @param    new_slit_fraction   The newly wishded slit fraction
  @return   The newly computed trace 
  @see  cr2res_trace_new_trace()
  
  For each order in the input tracewave, the function will produce a new row 
  in the output trace table corresponding to the passed slit fraction, based 
  on all available traces (rows).

  - CR2RES_COL_SLIT_FRACTION is filled with the input fraction
  - CR2RES_COL_ORDER and CR2RES_COL_TRACENB are copied from input trace file
  - CR2RES_COL_UPPER CR2RES_COL_LOWER CR2RES_COL_ALL are computed by
        cr2res_trace_new_trace()
  - CR2RES_COL_WAVELENGTH and CR2RES_COL_WAVELENGTH_ERROR are computed
    like this:
    TODO
  - CR2RES_COL_SLIT_CURV_A CR2RES_COL_SLIT_CURV_B CR2RES_COL_SLIT_CURV_C
    are computed like this:
    TODO

    COMMENTS / TODO
        Currently only calculates the new position of the trace
        Wavelength and slit curvature remain unchanged
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_trace_new_slit_fraction(
        const cpl_table     *   traces,
        const cpl_array     *   new_slit_fraction)
{
    cpl_table       *   out ;
    cpl_array *   slit_frac_old_trace ;
    const cpl_array *   slit_frac_old ;
    const cpl_array *   trace_all_old ;
    const cpl_array *   trace_upper_old ;
    const cpl_array *   trace_lower_old ;
    const cpl_array **  trace_old;
    cpl_array       *   trace_all_new ;
    cpl_array       *   trace_upper_new ;
    cpl_array       *   trace_lower_new ;
    cpl_array       *   wave ;
    const cpl_array *   const_wave ;
    cpl_array       *   wave_err ;
    cpl_array       *   slit_curv_a ;
    cpl_array       *   slit_curv_b ;
    cpl_array       *   slit_curv_c ;
    const cpl_array *   const_slit_curv_a;
    const cpl_array *   const_slit_curv_b;
    const cpl_array *   const_slit_curv_c;
    cpl_polynomial  *   poly_a ;
    cpl_polynomial  *   poly_b ;
    cpl_polynomial  *   poly_c ;
    cpl_polynomial  *   poly_tmp;
    cpl_size            i, j, k, m, nrows ;
    int nb_orders, nb_traces;
    int * orders, * trace_numbers;
    double a, b, c;
    double pix_lower, pix_upper;
    double sf_lower, sf_upper, sf_all, sf_new;
    double pix_shift;
    int isError = FALSE;

    /* Check entries */
    if (traces == NULL || new_slit_fraction == NULL) return NULL ;
    if (cr2res_dfs_check_traces_table(traces) != 1) {
        cpl_msg_error(__func__, "The traces table is incomplete") ;
        return NULL ;
    }
    if (cpl_array_get_size(new_slit_fraction) != 3){
        cpl_msg_error(__func__, 
            "New slitfraction must have 3 components (Upper, Middle, Lower)");
        return NULL;
    }

    /* Initialise */
    nrows = cpl_table_get_nrow(traces) ;
    orders = cr2res_trace_get_order_numbers(traces, &nb_orders);
    out = cpl_table_new(nb_orders);
    cpl_table_copy_structure(out, traces);
    cpl_table_select_all(out);

    /* Loop on the orders */
    for (i=0 ; i<nb_orders; i++) {
        trace_numbers = cr2res_get_trace_numbers(traces, orders[i], &nb_traces);
        trace_old = cpl_malloc(nb_traces * 3 * sizeof(cpl_array*));
        for (j = 0; j < nb_traces * 3; j++) trace_old[j] = NULL;
        slit_frac_old_trace = cpl_array_new(nb_traces * 3, CPL_TYPE_DOUBLE);        

        for(j = m = 0; j < nb_traces; j++ + m++) {
            k=cr2res_get_trace_table_index(traces, orders[i], trace_numbers[j]);
            /* Check if the input trace slit_fraction is available */
            slit_frac_old = cpl_table_get_array(traces,
                    CR2RES_COL_SLIT_FRACTION, k) ;
            trace_all_old = cpl_table_get_array(traces, CR2RES_COL_ALL, k) ;
            trace_upper_old = cpl_table_get_array(traces, CR2RES_COL_UPPER, k) ;
            trace_lower_old = cpl_table_get_array(traces, CR2RES_COL_LOWER, k) ;

            /* Unselect rows with wrong slit_fraction or without trace */
            /* to be erased below */
            if (cpl_table_is_selected(out, i) == 0 ||
                    cr2res_trace_check_slit_fraction(slit_frac_old) != 1 ||
                    trace_upper_old == NULL || trace_lower_old == NULL ||
                    trace_all_old == NULL) {
                cpl_table_unselect_row(out, i) ;
                m--;
                continue ;
            }
            trace_old[m + 0] = trace_lower_old;
            trace_old[m + 1] = trace_all_old;
            trace_old[m + 2] = trace_upper_old;
            cpl_array_set(slit_frac_old_trace, m + 0, cpl_array_get(slit_frac_old, 0, NULL));
            cpl_array_set(slit_frac_old_trace, m + 1, cpl_array_get(slit_frac_old, 1, NULL));
            cpl_array_set(slit_frac_old_trace, m + 2, cpl_array_get(slit_frac_old, 2, NULL));
        }

        if (m <= 0){
            cpl_msg_error(__func__, 
                    "No valid traces found for order %i", orders[i]) ;
            cpl_array_delete(slit_frac_old_trace);
            cpl_free(trace_numbers);
            cpl_free(trace_old);
            continue;
        }

        trace_old = cpl_realloc(trace_old, m * 3 * sizeof(cpl_array*));
        cpl_array_set_size(slit_frac_old_trace, m * 3);

        /* Fill out fields */
        cpl_table_set_array(out, CR2RES_COL_SLIT_FRACTION, i,new_slit_fraction);
        cpl_table_set_int(out, CR2RES_COL_ORDER, i, orders[i]);
        cpl_table_set_int(out, CR2RES_COL_TRACENB, i, 1);

        /* Compute the new trace */
        if (cr2res_trace_new_trace(slit_frac_old_trace, trace_old, nb_traces, 
                    new_slit_fraction, &trace_all_new, &trace_upper_new, 
                    &trace_lower_new) == -1) {
            cpl_msg_error(__func__, 
                    "Cannot compute the new trace for order %i", orders[i]) ;
            cpl_free(trace_old);
            cpl_free(trace_numbers);
            cpl_array_delete(slit_frac_old_trace);
            continue;
        }
        cpl_array_delete(slit_frac_old_trace);

        /* Set new trace */
        cpl_table_set_array(out, CR2RES_COL_ALL, i, trace_all_new);
        cpl_table_set_array(out, CR2RES_COL_UPPER, i, trace_upper_new);
        cpl_table_set_array(out, CR2RES_COL_LOWER, i, trace_lower_new);
        cpl_array_delete(trace_all_new) ;
        cpl_array_delete(trace_upper_new) ;
        cpl_array_delete(trace_lower_new) ;

        // Read data arrays from existing table
        // TODO this only works if trace_numbers 0 is a valid trace
        k = cr2res_get_trace_table_index(traces, orders[i], trace_numbers[0]);
        const_wave = cpl_table_get_array(traces, CR2RES_COL_WAVELENGTH, k) ;
        if (any(const_wave, isnan)){
            cpl_msg_error(__func__, 
                    "Wavelength polynomial is not set for order %i", orders[i]) ;
            cpl_free(trace_old);
            cpl_free(trace_numbers);
            continue;
        }

        wave_err = cpl_array_duplicate(
                cpl_table_get_array(traces, CR2RES_COL_WAVELENGTH_ERROR, k)) ; 

        const_slit_curv_a = cpl_table_get_array(traces, CR2RES_COL_SLIT_CURV_A,
                k) ; 
        const_slit_curv_b = cpl_table_get_array(traces, CR2RES_COL_SLIT_CURV_B,
                k) ; 
        const_slit_curv_c = cpl_table_get_array(traces, CR2RES_COL_SLIT_CURV_C,
                k) ; 

        /* Compute the new wavelength */
        // Calculate the horizontal pixel shift of the new trace, using
        // the slit curvature first determine the vertical shift
        slit_frac_old = cpl_table_get_array(traces, CR2RES_COL_SLIT_FRACTION,k);
        trace_lower_old = trace_old[0];
        trace_all_old = trace_old[1];
        trace_upper_old = trace_old[2];

        poly_tmp = cr2res_convert_array_to_poly(trace_lower_old);
        pix_lower = cpl_polynomial_eval_1d(poly_tmp, 
                (double)(CR2RES_DETECTOR_SIZE/2.0), NULL);
        cpl_polynomial_delete(poly_tmp);

        poly_tmp = cr2res_convert_array_to_poly(trace_upper_old);
        pix_upper = cpl_polynomial_eval_1d(poly_tmp, 
                (double)(CR2RES_DETECTOR_SIZE/2.0), NULL);
        cpl_polynomial_delete(poly_tmp);

        sf_lower = cpl_array_get_double(slit_frac_old, 0, NULL);
        sf_all = cpl_array_get_double(slit_frac_old, 1, NULL);
        sf_upper = cpl_array_get_double(slit_frac_old, 2, NULL);
        sf_new = cpl_array_get_double(new_slit_fraction, 1, NULL);

        // then use the slit curvature to translate that into a horizontal shift
        poly_a = cr2res_convert_array_to_poly(const_slit_curv_a);
        poly_b = cr2res_convert_array_to_poly(const_slit_curv_b);
        poly_c = cr2res_convert_array_to_poly(const_slit_curv_c);

        a = cpl_polynomial_eval_1d(poly_a,(CR2RES_DETECTOR_SIZE/2.0) + 1, NULL);
        b = cpl_polynomial_eval_1d(poly_b,(CR2RES_DETECTOR_SIZE/2.0) + 1, NULL);
        c = cpl_polynomial_eval_1d(poly_c,(CR2RES_DETECTOR_SIZE/2.0) + 1, NULL);

        // vertical pixel shift
        pix_shift = (pix_upper - pix_lower) / 
            (sf_upper - sf_lower) * (sf_all - sf_new);
        // horizontal pixel shift
        pix_shift = (a - CR2RES_DETECTOR_SIZE/2. - 1) + 
            b * pix_shift + c * pix_shift * pix_shift;

        poly_tmp = cr2res_convert_array_to_poly(const_wave);
        cpl_polynomial_shift_1d(poly_tmp, 0, pix_shift);
        wave = cr2res_convert_poly_to_array(poly_tmp, 
                cpl_array_get_size(const_wave));
        cpl_polynomial_delete(poly_tmp);

        /* Compute new Curvature */
        cpl_polynomial_shift_1d(poly_a, 0, pix_shift);
        cpl_polynomial_shift_1d(poly_b, 0, pix_shift);
        cpl_polynomial_shift_1d(poly_c, 0, pix_shift);

        slit_curv_a = cr2res_convert_poly_to_array(poly_a, 
                cpl_array_get_size(const_slit_curv_a));
        slit_curv_b = cr2res_convert_poly_to_array(poly_b, 
                cpl_array_get_size(const_slit_curv_b));
        slit_curv_c = cr2res_convert_poly_to_array(poly_c, 
                cpl_array_get_size(const_slit_curv_c));

        cpl_polynomial_delete(poly_a);
        cpl_polynomial_delete(poly_b);
        cpl_polynomial_delete(poly_c);


        /* Set new Wavelength  */
        cpl_table_set_array(out, CR2RES_COL_WAVELENGTH, i, wave);
        cpl_table_set_array(out, CR2RES_COL_WAVELENGTH_ERROR, i, wave_err);
        cpl_array_delete(wave) ;
        cpl_array_delete(wave_err) ;

        /* Set the new slit curvature */
        cpl_table_set_array(out, CR2RES_COL_SLIT_CURV_A, i, slit_curv_a);
        cpl_table_set_array(out, CR2RES_COL_SLIT_CURV_B, i, slit_curv_b);
        cpl_table_set_array(out, CR2RES_COL_SLIT_CURV_C, i, slit_curv_c);
        cpl_array_delete(slit_curv_a) ;
        cpl_array_delete(slit_curv_b) ;
        cpl_array_delete(slit_curv_c) ;

        cpl_free(trace_numbers);
        cpl_free(trace_old);
    }
    cpl_free(orders);

    /* Clean the selected invalid slit_fractions entries */
    cpl_table_not_selected(out) ;
    cpl_table_erase_selected(out) ;

    /* Check if there is at lease one remaining trace */
    if (cpl_table_get_nrow(out) == 0) {
        cpl_table_delete(out) ;
        return NULL ;
    }
    cpl_table_not_selected(out) ;
    return out ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief Adjust the traces spatial positions using a new FLAT frameset
  @param trace_wave     The trace_wave to adjust
  @param flat_raw       The FLAT frameset
  @param det_nr         The detector number
  @return   the new adjusted trace_wave or NULL
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_trace_adjust(
        const cpl_table     *   trace_wave,
        const cpl_frameset  *   flat_raw,
        int                     det_nr)
{
    hdrl_imagelist      *   imlist ;
    hdrl_image          *   collapsed ;
    cpl_image           *   contrib ;
    cpl_table           *   new_traces ;
    const char          *   first_file ;
    cpl_table           *   corrected_traces ;
    cpl_table           *   traces ;
    int                     trace_opening, trace_degree,
                            trace_min_cluster, trace_smooth_x, trace_smooth_y ;
    double                  trace_threshold, traces_shift ;
 
    /* Check Entries */
    if (trace_wave == NULL || flat_raw == NULL || 
            det_nr < 1 || det_nr > CR2RES_NB_DETECTORS) 
        return NULL ;

    /* Initialise */
    trace_smooth_x = 200 ;
    trace_smooth_y = 11 ;
    trace_threshold = 5.0 ;
    trace_opening = 1 ;
    trace_degree = 5 ;
    trace_min_cluster = 10000 ;

    /* Load the image list */
    imlist = cr2res_io_load_image_list_from_set(flat_raw, det_nr) ;
    if (imlist == NULL) {
        cpl_msg_error(__func__, "Failed to load the RAW flat frames") ;
        return NULL ;
    }

    /* Collapse the inputs */
    cpl_msg_info(__func__, "Collapse the input images") ;
    if (hdrl_imagelist_collapse_mean(imlist, &collapsed, &contrib) !=
            CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Failed to Collapse") ;
        hdrl_imagelist_delete(imlist) ;
        cpl_msg_indent_less() ;
        return NULL ;
    }
    hdrl_imagelist_delete(imlist) ;
    cpl_image_delete(contrib) ;

    /* Compute the trace */
    cpl_msg_info(__func__, "Compute the traces") ;
    if ((new_traces = cr2res_trace(hdrl_image_get_image(collapsed),
                    trace_smooth_x, trace_smooth_y, trace_threshold,
                    trace_opening, trace_degree, trace_min_cluster)) == NULL) {
        cpl_msg_error(__func__, "Failed compute the traces") ;
        hdrl_image_delete(collapsed) ;
        cpl_msg_indent_less() ;
        return NULL ;
    }
    hdrl_image_delete(collapsed) ;

    /* Add The remaining Columns to the new trace table */
    first_file = cpl_frame_get_filename(
            cpl_frameset_get_position_const(flat_raw, 0)) ;
    cr2res_trace_add_extra_columns(new_traces, first_file, det_nr) ;

    /* Compute the shift */
    cpl_msg_info(__func__, "Compute the Shift between 2 traces tables") ;
    traces_shift = cr2res_trace_compute_shift(trace_wave, new_traces) ;
    cpl_table_delete(new_traces) ;

    /* Apply the shift */
    corrected_traces = cpl_table_duplicate(trace_wave) ;
    cpl_msg_error(__func__, "Apply correction shift of %g pixels",
            traces_shift) ;
    cr2res_trace_apply_shift(corrected_traces, traces_shift) ;

    /* Free and return */
    return corrected_traces ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get a standard slit fraction information
  @param    slit_frac   Slit Fraction array
  @param    up_or_down  [out]   0 (if CR2RES_DECKER_NONE), 1 for up, 2 for down
  @return   CR2RES_DECKER_NONE, _1_3 or _2_4 or CR2RES_DECKER_INVALID
 */
/*----------------------------------------------------------------------------*/
cr2res_decker cr2res_trace_slit_fraction_info(
        const cpl_array *   slit_frac,
        int             *   up_or_down)
{
    double          up, down, mid ;
    cr2res_decker   decker_pos ;

    /* Initialise */
    if (up_or_down != NULL) *up_or_down=-1 ;

    /* Check entries */
    if (slit_frac == NULL) return CR2RES_DECKER_INVALID ;

    /* Read the positions */
    down = cpl_array_get(slit_frac, 0, NULL) ;
    mid = cpl_array_get(slit_frac, 1, NULL) ;
    up = cpl_array_get(slit_frac, 2, NULL) ;
    if (cpl_error_get_code()) return CR2RES_DECKER_INVALID ;

    if (fabs((down-0.0)<1e-3) && fabs((up-1.0)<1e-3)) {
        decker_pos = CR2RES_DECKER_NONE ;
        if (up_or_down != NULL) *up_or_down=0 ;
    } else if (fabs((down-0.750)<1e-3) && fabs((up-1.000)<1e-3)) {
        decker_pos = CR2RES_DECKER_1_3 ;
        if (up_or_down != NULL) *up_or_down=1 ;
    } else if (fabs((down-0.250)<1e-3) && fabs((up-0.500)<1e-3)) {
        decker_pos = CR2RES_DECKER_1_3 ;
        if (up_or_down != NULL) *up_or_down=2 ;
    } else if (fabs((down-0.500)<1e-3) && fabs((up-0.750)<1e-3)) {
        decker_pos = CR2RES_DECKER_2_4 ;
        if (up_or_down != NULL) *up_or_down=1 ;
    } else if (fabs((down-0.000)<1e-3) && fabs((up-0.250)<1e-3)) {
        decker_pos = CR2RES_DECKER_2_4 ;
        if (up_or_down != NULL) *up_or_down=2 ;
    } else {
        decker_pos = CR2RES_DECKER_INVALID ;
        if (up_or_down != NULL) *up_or_down=-1 ;
    }
    return decker_pos ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get a standard slit fraction
  @param    decker_position     CR2RES_DECKER_NONE, _1_3 or _2_4
  @param    up_or_down          0 (if CR2RES_DECKER_NONE), 1 for up, 2 for down
  @return   A newly allocated mask or NULL in error case.

 */
/*----------------------------------------------------------------------------*/
cpl_array * cr2res_trace_slit_fraction_create(
        cr2res_decker   decker_position,
        int             up_or_down)
{
    cpl_array   *   slit_frac ;

    /* Check entries */
    if (decker_position != CR2RES_DECKER_NONE &&
            decker_position != CR2RES_DECKER_1_3 &&
            decker_position != CR2RES_DECKER_2_4) 
        return NULL ;
    if (up_or_down != 0 && up_or_down != 1 && up_or_down != 2)
        return NULL ;

    /* Allocate */
    slit_frac = cpl_array_new(3, CPL_TYPE_DOUBLE) ;

    if (decker_position == CR2RES_DECKER_NONE) {
        /* TRACE_OPEN */
        cpl_array_set(slit_frac, 0, 0.0) ;
        cpl_array_set(slit_frac, 1, 0.5) ;
        cpl_array_set(slit_frac, 2, 1.0) ;
    } else if (decker_position == CR2RES_DECKER_1_3) {
        /* DECKER 1_3 */
        if (up_or_down == 1) {
            /* Upper trace */
            cpl_array_set(slit_frac, 0, 0.750) ;
            cpl_array_set(slit_frac, 1, 0.875) ;
            cpl_array_set(slit_frac, 2, 1.000) ;
        } else if (up_or_down == 2) {
            /* Lower trace */
            cpl_array_set(slit_frac, 0, 0.250) ;
            cpl_array_set(slit_frac, 1, 0.375) ;
            cpl_array_set(slit_frac, 2, 0.500) ;
        } else {
            cpl_array_delete(slit_frac) ;
            slit_frac = NULL ;
        }
    } else if (decker_position == CR2RES_DECKER_2_4) {
        /* DECKER 2_4 */
        if (up_or_down == 1) {
            /* Upper trace */
            cpl_array_set(slit_frac, 0, 0.500) ;
            cpl_array_set(slit_frac, 1, 0.625) ;
            cpl_array_set(slit_frac, 2, 0.750) ;
        } else if (up_or_down == 2) {
            /* Lower trace */
            cpl_array_set(slit_frac, 0, 0.000) ;
            cpl_array_set(slit_frac, 1, 0.125) ;
            cpl_array_set(slit_frac, 2, 0.250) ;
        } else {
            cpl_array_delete(slit_frac) ;
            slit_frac = NULL ;
        }
    }
    return slit_frac ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Splits full slit traces into several sub traces
  @param    trace_wave       trace wave table
  @param    order            order to split (-1 for all)
  @param    nb_subtraces     number of subtraces to create
  @return   the new trace_wave table or NULL in error case

  All input traces that have a full slit_fraction are splitted.

  The returned trace wave table only contains the newly created traces.
  
  Sub traces are created by evenly splitting the existing trace, so that no 
  subtraces overlap.
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_trace_split(
        cpl_table   *   trace_wave, 
        int             order, 
        int             nb_subtraces) 
{
    cpl_table       *   sub_trace_wave ;
    cpl_size            i, table_index;
    cpl_array       *bottom, *center, *top, *fraction, *wave, *tmp, *tmp2;
    const cpl_array *wave_err, *slit_a, *slit_b, *slit_c;
    double              height = 1. / (nb_subtraces * 2.);
    double              pos;
    int                 res;

    if (trace_wave == NULL || nb_subtraces <= 0) return NULL ;

    /* TMP */
    return cpl_table_duplicate(trace_wave) ;
    /* end TMP */

    // Expand the trace table to fit the new data
    sub_trace_wave = cpl_table_new(nb_subtraces);
    cpl_table_copy_structure(sub_trace_wave, trace_wave);

    // Prepare data to fill the other columns of the table
    i = cr2res_get_trace_table_index(trace_wave, order, 1);
    wave_err = cpl_table_get_array(trace_wave, CR2RES_COL_WAVELENGTH_ERROR, i);
    slit_a = cpl_table_get_array(trace_wave, CR2RES_COL_SLIT_CURV_A, i);
    slit_b = cpl_table_get_array(trace_wave, CR2RES_COL_SLIT_CURV_B, i);
    slit_c = cpl_table_get_array(trace_wave, CR2RES_COL_SLIT_CURV_C, i);

    for (i = 0; i < nb_subtraces; i++){
        // center of first subtrace is one height above bottom
        pos = height + 2 * height * i; 
        res = cr2res_trace_get_subtrace(trace_wave, pos, height, order, &bottom,
                &center, &top, &fraction, &wave);
        if (res == -1) break; // if something went wrong stop here

        // attach new trace information to the bottom of the table
        table_index = i;
        cpl_table_set_int(sub_trace_wave, CR2RES_COL_ORDER, table_index,order);
        // Traces start counting at 1
        cpl_table_set_int(sub_trace_wave, CR2RES_COL_TRACENB, table_index, 
                i + 1);

        // Assign interpolated data
        cpl_table_set_array(sub_trace_wave, CR2RES_COL_ALL, table_index, 
                center);
        cpl_table_set_array(sub_trace_wave, CR2RES_COL_LOWER, table_index, 
                bottom);
        cpl_table_set_array(sub_trace_wave, CR2RES_COL_UPPER, table_index, 
                top);
        cpl_table_set_array(sub_trace_wave, CR2RES_COL_SLIT_FRACTION, 
                table_index, fraction);
        cpl_table_set_array(sub_trace_wave, CR2RES_COL_WAVELENGTH, 
                table_index, wave);

        // keep wavelength error from original order
        cpl_table_set_array(sub_trace_wave, CR2RES_COL_WAVELENGTH_ERROR, 
                table_index, wave_err);

        // set slit polynomials to default values for now
        cpl_table_set_array(sub_trace_wave, CR2RES_COL_SLIT_CURV_A, 
                table_index, slit_a);
        cpl_table_set_array(sub_trace_wave, CR2RES_COL_SLIT_CURV_B, 
                table_index, slit_b);
        cpl_table_set_array(sub_trace_wave, CR2RES_COL_SLIT_CURV_C, 
                table_index, slit_c);

        // Clean up data
        cpl_array_delete(bottom);
        cpl_array_delete(center);
        cpl_array_delete(top);
        cpl_array_delete(fraction);
    }
    cpl_array_unwrap(tmp);
    cpl_array_unwrap(tmp2);

    if (res == -1) {
        cpl_table_delete(sub_trace_wave);
        return NULL ;
    }
    return sub_trace_wave ;
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief  Compares 2 traces files, and returns the estimated shift between both
  @param    traces1     A trace wave file
  @param    traces2     A trace wave file
  @return   The shift between both. Return 0.0 in case of problem

  The function tries to detect the shift caused by reproducibility issue
  of the instrument. Both set of traces should have the same spacing,
  they should be consistently shifted. 
  If the function does not recognise the same pattern, it should return
  0.0.
 */
/*----------------------------------------------------------------------------*/
static double cr2res_trace_compute_shift(
        const cpl_table *   traces1,
        const cpl_table *   traces2)
{
    /* TODO */
    int             order, trace_id ;
    cpl_size        i ;

    for (i=0 ; i<cpl_table_get_nrow(traces1) ; i++) {
        /* Only the Open slit */

        /* Get the order */
        order = cpl_table_get(traces1, CR2RES_COL_ORDER, i, NULL) ;
        if (order >= 0) {

        }
    }
    return 0.0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief  Apply a shift to the traces only (Upper, Lower and All)
  @param    traces  The trace table to shift
  @param    shift   The shift that needs to be applied
  @return   0 if ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
static int cr2res_trace_apply_shift(
        cpl_table       *   traces,
        double              shift)
{
    /* TODO */
    return -1 ;
}

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

  The function assumes that the slit is open if there is 1 trace in the
  order, and that the decker is used if there are 2 traces.
 */
/*----------------------------------------------------------------------------*/
static cpl_array * cr2res_trace_get_slit_frac(
        const cpl_table *   traces,
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
        slit_frac = cr2res_trace_slit_fraction_create(decker_pos, 0) ;
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
        if (center_pos_curr < center_pos_other) {
            /* Lower Trace */
            slit_frac = cr2res_trace_slit_fraction_create(decker_pos, 2) ;
        } else {
            /* Upper Trace */
            slit_frac = cr2res_trace_slit_fraction_create(decker_pos, 1) ;
        }
    } else {
        /* Inconsistent */
        slit_frac = NULL ;
    }
    return slit_frac ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Detect the Traces signal
  @param image          The input image with the traces
  @param smooth_x       Low pass filter kernel size in x
  @param smooth_y       Low pass filter kernel size in y
  @param thresh         The threshold used for detection
  @return   A newly allocated mask or NULL in error case.

  The returned mask identifies the pixels belonging to a trace
  The input image is smoothed, subtracted to the result, and a simple
  thresholding is applied. 
 */
/*----------------------------------------------------------------------------*/
static cpl_mask * cr2res_trace_signal_detect(
        const cpl_image *   image,
        int                 smooth_x,
        int                 smooth_y,
        double              thresh)
{
    cpl_image       *   smx_image ;
    cpl_image       *   smxy_image ;
    int                 kernel_x, kernel_y ;
    cpl_mask        *   kernel ;
    cpl_mask        *   mask ;

    /* Check Entries */
    if (image == NULL) return NULL;
    if (smooth_x < 0 || smooth_y < 0) return NULL;

    /* Prepare kernel */
    kernel_x = smooth_x ;
    if (kernel_x % 2 == 0) kernel_x++ ;
    kernel_y = smooth_y ;
    if (kernel_y % 2 == 0) kernel_y++ ;

    // Start with smooth in X
    kernel = cpl_mask_new(kernel_x, 1);
    cpl_mask_not(kernel);
    
    smx_image = cpl_image_duplicate(image);
    if (cpl_image_filter_mask(smx_image, image, kernel, CPL_FILTER_AVERAGE_FAST,
                CPL_BORDER_FILTER) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot filter the image") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        cpl_mask_delete(kernel);
        cpl_image_delete(smx_image) ;
        return NULL ;
    }
    cpl_mask_delete(kernel);

    // Now smooth Y
    kernel = cpl_mask_new(1, kernel_y);
    cpl_mask_not(kernel);
    
    smxy_image = cpl_image_duplicate(smx_image);
    if (cpl_image_filter_mask(smxy_image, smx_image, kernel, CPL_FILTER_AVERAGE_FAST,
                CPL_BORDER_FILTER) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot filter the image") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        cpl_mask_delete(kernel);
        cpl_image_delete(smx_image) ;
        cpl_image_delete(smxy_image) ;
        return NULL ;
    }
    cpl_mask_delete(kernel);
    
    /* Subtract x-smoothed image from xy-smoothed */
    cpl_image_subtract(smxy_image, smx_image);

    if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
        cpl_image_save(smxy_image, "debug_smimage.fits", CPL_TYPE_DOUBLE, NULL,
                CPL_IO_CREATE);
        cpl_msg_debug(__func__, "Smooth X: %d, Y: %d, Threshold: %.1f",
                kernel_x, kernel_y, thresh);
    }

    /* Wanted pixels are where input image exceeds sm_image by thresh */
    mask = cpl_mask_threshold_image_create(smxy_image,-DBL_MAX, thresh);
    cpl_image_delete(smx_image) ;
    cpl_image_delete(smxy_image) ;

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
    if (degree < 0) return NULL;

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
    if (min_cluster < 0) return NULL;

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

/*----------------------------------------------------------------------------*/
/**
  @brief    Check if the passed slit fraction is valid
  @param    slit_frac   Input slit fraction
  @return   1 if valid, 0 if not, -1 in error case
  The slit fraction values should be in [0,1] and increasing
 */
/*----------------------------------------------------------------------------*/
static int cr2res_trace_check_slit_fraction(const cpl_array * slit_fraction)
{
    double      low, mid, up ;

    /* Check entries */
    if (slit_fraction == NULL) return -1; 

    /* Get the values */
    low = cpl_array_get(slit_fraction, 0, NULL) ;
    mid = cpl_array_get(slit_fraction, 1, NULL) ;
    up = cpl_array_get(slit_fraction, 2, NULL) ;

    if (cpl_error_get_code() != CPL_ERROR_NONE) return -1 ;

    if (up<0.0 || low<0.0 || mid<0.0 || up>1.0 || low>1.0 || mid>1.0)
        return 0 ;
    if (low>mid || mid > up) 
        return 0 ;
    return 1 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief   Compute new trace polynomials by using slit_fraction specs
  @param slit_fraction_in     Input slit fraction
  @param trace_in             Array os size [3 * nb_traces]
                              Traces are ordered in blocks of (lower, all, 
                              upper), then the next trace.
  @param nb_traces            number of traces in the order
  @param slit_fraction_wished
  @param trace_all_new
  @param trace_upper_new
  @param trace_lower_new
  @return 0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
static int cr2res_trace_new_trace(
        const cpl_array     *   slit_fraction_in,
        const cpl_array     **  trace_in,
        const int               nb_traces,
        const cpl_array     *   slit_fraction_wished,
        cpl_array           **  trace_all_new,
        cpl_array           **  trace_upper_new,
        cpl_array           **  trace_lower_new)
{
    cpl_polynomial  **  poly_in ;
    cpl_polynomial  *   poly_slit;
    cpl_polynomial  *   poly_upper_out ;
    cpl_polynomial  *   poly_lower_out ;
    cpl_polynomial  *   poly_all_out ;
    const cpl_size 		power = 0 ;
    cpl_array       *   trace_all_out ;
    cpl_array       *   trace_upper_out ;
    cpl_array       *   trace_lower_out ;
    cpl_matrix      *   sf_matrix;
    cpl_error_code      error;
    cpl_size            fitdeg;
    const double    *   sf_in;
    cpl_vector      *   pix_in;
    double              sf_out_l, sf_out_m, sf_out_u, pix_out_l, pix_out_u, 
                        pix_out_m, sf_wished, new_coeff, b, c, divisor,
                        x2, x3, y2, y3 ;

    /* Check entries */
    if (slit_fraction_in==NULL || trace_in==NULL || slit_fraction_wished==NULL
            || trace_all_new==NULL || trace_upper_new==NULL
            || trace_lower_new==NULL || nb_traces <= 0) return -1; 

    // Both in and output slitfractions need to have 3 components 
    // (3 per trace for input)
    if (cpl_array_get_size(slit_fraction_in) != 3 * nb_traces || 
                    cpl_array_get_size(slit_fraction_wished) != 3) return -1;



    /* Get wished slit positions in arcseconds */
    sf_out_l = cpl_array_get(slit_fraction_wished, 0, NULL) ;
    sf_out_m = cpl_array_get(slit_fraction_wished, 1, NULL) ;
    sf_out_u = cpl_array_get(slit_fraction_wished, 2, NULL) ;

    // Allocate memory
    poly_in = cpl_malloc(3 * nb_traces * sizeof(cpl_polynomial*));
    for (cpl_size i = 0; i < nb_traces * 3; i++) poly_in[i] = NULL;
    pix_in = cpl_vector_new(3 * nb_traces);

    /* Get input slit positions in pixels from middle of the traces */
    for(cpl_size i = 0; i < nb_traces * 3; i++)
    {
        poly_in[i] = cr2res_convert_array_to_poly(trace_in[i]);
        cpl_vector_set(pix_in, i, cpl_polynomial_eval_1d(poly_in[i], 
            (double)(CR2RES_DETECTOR_SIZE/2.0), NULL));
    }
    
    // Fit a 2nd order polynomial to all points
    poly_slit = cpl_polynomial_new(1);
    /* Get input slit positions in arcseconds */
    sf_in = cpl_array_get_data_double_const(slit_fraction_in);
    sf_matrix = cpl_matrix_wrap(1, nb_traces * 3, (double*)sf_in);
    fitdeg = 2;
    error = cpl_polynomial_fit(poly_slit, sf_matrix, NULL, pix_in, NULL, 
                                CPL_FALSE, NULL, &fitdeg);
    cpl_matrix_unwrap(sf_matrix);

    // in case something went wrong
    if (error != CPL_ERROR_NONE){
        cpl_polynomial_delete(poly_slit);
        for(cpl_size i = 0; i < nb_traces*3; i++)
        {
            cpl_polynomial_delete(poly_in[i]);
        }
        cpl_free(poly_in);
        cpl_vector_delete(pix_in);
        cpl_error_reset();
        return -1;
    }

    /* Compute the Pixel positions on the wished slit fraction */
    pix_out_l = cpl_polynomial_eval_1d(poly_slit, sf_out_l, NULL);
    pix_out_m = cpl_polynomial_eval_1d(poly_slit, sf_out_m, NULL);
    pix_out_u = cpl_polynomial_eval_1d(poly_slit, sf_out_u, NULL);

    /* Debug Message */
    cpl_msg_debug(__func__, 
    "Slit fraction: [%g-%g-%g]->[%g-%g-%g] - Pixel positions: [%g-%g-%g]->[%g-%g-%g]",
    sf_in[0], sf_in[1], sf_in[2], sf_out_l, sf_out_m, sf_out_u, 
    cpl_vector_get(pix_in, 0), cpl_vector_get(pix_in, 1),
    cpl_vector_get(pix_in, 2), pix_out_l, pix_out_m, pix_out_u) ;

    /* Correct the polynomials */
    poly_upper_out = cpl_polynomial_duplicate(poly_in[2]) ;
    poly_lower_out = cpl_polynomial_duplicate(poly_in[0]) ;
    poly_all_out = cpl_polynomial_duplicate(poly_in[1]) ;

    new_coeff = cpl_polynomial_get_coeff(poly_upper_out, &power) + pix_out_u 
                - cpl_vector_get(pix_in, 2) ;
    cpl_polynomial_set_coeff(poly_upper_out, &power, new_coeff) ;
    new_coeff = cpl_polynomial_get_coeff(poly_lower_out, &power) + pix_out_l 
                - cpl_vector_get(pix_in, 0) ;
    cpl_polynomial_set_coeff(poly_lower_out, &power, new_coeff) ;
    new_coeff = cpl_polynomial_get_coeff(poly_all_out, &power) + pix_out_m 
                - cpl_vector_get(pix_in, 1) ;
    cpl_polynomial_set_coeff(poly_all_out, &power, new_coeff) ;

    for(cpl_size i = 0; i < nb_traces*3; i++)
    {
        cpl_polynomial_delete(poly_in[i]);
    }
    cpl_free(poly_in);
    cpl_vector_delete(pix_in);
    cpl_polynomial_delete(poly_slit);
    

    /* Convert to arrays */
    trace_all_out = cr2res_convert_poly_to_array(poly_all_out, 
            cpl_array_get_size(trace_in[1]));
    trace_upper_out = cr2res_convert_poly_to_array(poly_upper_out, 
            cpl_array_get_size(trace_in[2]));
    trace_lower_out = cr2res_convert_poly_to_array(poly_lower_out, 
            cpl_array_get_size(trace_in[0]));
    cpl_polynomial_delete(poly_all_out) ;
    cpl_polynomial_delete(poly_upper_out) ;
    cpl_polynomial_delete(poly_lower_out) ;
   
    /* Check */
    if (cpl_error_get_code() != CPL_ERROR_NONE) {
        cpl_array_delete(trace_all_out) ;
        cpl_array_delete(trace_upper_out) ;
        cpl_array_delete(trace_lower_out) ;
        return -1 ;
    }

    /* Return  */
    *trace_all_new = trace_all_out ;
    *trace_upper_new = trace_upper_out ;
    *trace_lower_new = trace_lower_out ;
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Create a trace at given slit position with given height, based 
            on table data
  @param    trace_wave     trace wave table
  @param    slit_pos       fractional slit position (0-1)
  @param    height         fractional height (0-1)
  @param    order          order to split
  @param    trace          trace to split
  @param    bottom         [out] subtrace bottom polynomial coefficients
  @param    center         [out] subtrace center polynomial coefficients
  @param    top            [out] subtrace top polynomial coefficients
  @param    fraction       [out] subtrace polynomial fractions
  @return   0 on success, -1 on failure

  fits a 2nd order polynomial to the three coefficients of LOWER, ALL, 
  UPPER order trace polynomials and evaluates them at slit_pos - height, 
  slit_pos, and slit_pos + height

 */
/*----------------------------------------------------------------------------*/
static int cr2res_trace_get_subtrace(
        cpl_table   *   trace_wave, 
        double          slit_pos, 
        double          height, 
        int             order,
        cpl_array   **  bottom, 
        cpl_array   **  center, 
        cpl_array   **  top,
        cpl_array   **  fraction, 
        cpl_array   **  wave)
{
   
    // check input values
    if (slit_pos < 0 | slit_pos > 1) return -1;
    // get trace numbers
    int nb_traces;
    int * traces;
    traces = cr2res_get_trace_numbers(trace_wave, order, &nb_traces);
    if (traces == NULL) return -1;
    
    const cpl_array * bounds[3], *old_fraction, *old_wave;
    double res;
    cpl_size power, i, j, k, ndegree;
    
    cpl_matrix * samppos;
    cpl_vector *fitvals;
    cpl_polynomial * poly;

    cpl_size mindeg = 2;
    cpl_size maxdeg = 5;
    cpl_error_code error;

    
    // interpolate order tracing
    j = cr2res_get_trace_table_index(trace_wave, order, traces[0]);
    ndegree = cpl_table_get_column_dimension(trace_wave, CR2RES_COL_ALL, j);

    *bottom = cpl_array_new(ndegree, CPL_TYPE_DOUBLE);
    *top = cpl_array_new(ndegree, CPL_TYPE_DOUBLE);
    *center = cpl_array_new(ndegree, CPL_TYPE_DOUBLE);
    *fraction = cpl_array_new(ndegree, CPL_TYPE_DOUBLE);

    samppos = cpl_matrix_new(3 * nb_traces, 1);
    fitvals = cpl_vector_new(3 * nb_traces);
    poly = cpl_polynomial_new(1);


    for (i = 0; i < ndegree; i++){
        for (k = 0; k < nb_traces; k++){
            j = cr2res_get_trace_table_index(trace_wave, order, traces[k]);

            bounds[0] = cpl_table_get_array(trace_wave, CR2RES_COL_LOWER, j);
            bounds[1] = cpl_table_get_array(trace_wave, CR2RES_COL_ALL, j);
            bounds[2] = cpl_table_get_array(trace_wave, CR2RES_COL_UPPER, j);
            old_fraction = cpl_table_get_array(trace_wave, 
                    CR2RES_COL_SLIT_FRACTION, j);

            for (j = 0; j < 3; j++){
                res = cpl_array_get_double(old_fraction, 0, NULL);
                if (res == -1) res = 0.5 * j; // bottom: 0, center: 0.5, top: 1
                cpl_matrix_set(samppos, k * 3 + j, 0, res);
                cpl_vector_set(fitvals, k * 3 + j, 
                        cpl_array_get_double(bounds[j], i, NULL));
            }
        }

        error = cpl_polynomial_fit(poly, samppos, NULL, fitvals, NULL, 1, 
                &mindeg, &maxdeg);
        
        res = cpl_polynomial_eval_1d(poly, slit_pos, NULL);
        cpl_array_set_double(*center, i, res);

        res = cpl_polynomial_eval_1d(poly, slit_pos + height, NULL);
        cpl_array_set_double(*top, i, res);

        res = cpl_polynomial_eval_1d(poly, slit_pos - height, NULL);
        cpl_array_set_double(*bottom, i, res);
        
    }
    cpl_matrix_delete(samppos);
    cpl_vector_delete(fitvals);
    cpl_polynomial_delete(poly);

    cpl_array_set_double(*fraction, 0, slit_pos-height);
    cpl_array_set_double(*fraction, 1, slit_pos);
    cpl_array_set_double(*fraction, 2, slit_pos+height);

    // interpolate wavelength solution
    j = cr2res_get_trace_table_index(trace_wave, order, traces[0]);
    ndegree = cpl_table_get_column_dimension(trace_wave, CR2RES_COL_WAVELENGTH,
            j);
    *wave = cpl_array_new(ndegree, CPL_TYPE_DOUBLE);

    samppos = cpl_matrix_new(nb_traces, 1);
    fitvals = cpl_vector_new(nb_traces);
    poly = cpl_polynomial_new(1);

    if (nb_traces == 1){
        old_wave = cpl_table_get_array(trace_wave, CR2RES_COL_WAVELENGTH, j);
        cpl_array_copy_data_double(*wave, 
                cpl_array_get_data_double_const(old_wave));
    } else {
        for (i = 0; i < ndegree; i++){
            for (k = 0; k < nb_traces; k++){
                j = cr2res_get_trace_table_index(trace_wave, order, traces[k]);
                old_fraction = cpl_table_get_array(trace_wave, 
                        CR2RES_COL_SLIT_FRACTION, j);
                old_wave = cpl_table_get_array(trace_wave, 
                        CR2RES_COL_WAVELENGTH, j);

                res = cpl_array_get_double(old_fraction, 0, NULL);
                if (res == -1) res = 0.5 * j; // bottom: 0, center: 0.5, top: 1

                cpl_matrix_set(samppos, k, 0, res);
                cpl_vector_set(fitvals, k, cpl_array_get_double(old_wave, i, 
                            NULL));
            }
            error = cpl_polynomial_fit(poly, samppos, NULL, fitvals, NULL, 1, 
                    &mindeg, &maxdeg);
            res = cpl_polynomial_eval_1d(poly, slit_pos, NULL);
            cpl_array_set_double(*wave, i, res);
        }
    }
    cpl_matrix_delete(samppos);
    cpl_vector_delete(fitvals);
    cpl_polynomial_delete(poly);

    if (error != CPL_ERROR_NONE){
        cpl_array_delete(*bottom);
        cpl_array_delete(*center);
        cpl_array_delete(*top);
        cpl_array_delete(*fraction);
        cpl_array_delete(*wave);
        return -1;
    }
    return 0;
}

