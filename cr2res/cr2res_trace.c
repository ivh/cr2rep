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

#include <cpl.h>
#include "cr2res_trace.h"
#include "cr2res_utils.h"
#include "cr2res_cluster.h"

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static cpl_mask * cr2res_trace_signal_detect(
        const cpl_image *   image,
        int                 ordersep,
        double              smoothfactor,
        double              thresh) ;
static cpl_table * cr2res_trace_orders_fit(
        cpl_table   *   clustertable,
        int             degree) ;
static cpl_array * cr2res_trace_order_fit(
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
  @param    ima     input image
  @param    decker  slit layout
  @return  
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_trace_cpl(
        cpl_image       *   ima, 
        cr2res_decker       decker, 
        double              smoothfactor,
        int                 opening,
        int                 degree,
        int                 min_cluster)
{
    cpl_mask        *   mask ;
    cpl_image       *   labels ;
    cpl_apertures   *   aperts ;
    cpl_table       *   clustertable ;
    cpl_table       *   trace_table ;

    /* Check Entries */
    if (ima == NULL) return NULL ;

    /* Initialise */

    /* Detect the orders in the image */
    cpl_msg_info(__func__, "Orders detection") ;
    if ((mask = cr2res_trace_detect(ima, smoothfactor, opening,
                    min_cluster)) == NULL) {
        cpl_msg_error(__func__, "Cannot detect the orders") ;
        return NULL ;
    }
    
    /* Labelization */
    cpl_msg_info(__func__, "Labelise the orders") ;
    if ((labels = cr2res_trace_labelize(mask)) == NULL) {
        cpl_msg_error(__func__, "Cannot labelise") ;
        cpl_mask_delete(mask);
        return NULL ;
    }
    cpl_mask_delete(mask);
   
    /* Analyse and dump orders */
    if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
        aperts = cpl_apertures_new_from_image(ima, labels);
        cpl_apertures_dump(aperts, stdout) ;
        cpl_apertures_delete(aperts) ;
    }

    /* Fit Labels Edges */
    cpl_msg_info(__func__, "Fit the order edges") ;
    if ((trace_table = cr2res_trace_fit(labels, degree)) == NULL) {
        cpl_msg_error(__func__, "Cannot fit the order edges") ;
        cpl_image_delete(labels) ;
        return NULL ;
    }
    cpl_image_delete(labels) ;

    return trace_table ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief  Main function for running all parts of the trace algorithm 
  @param    ima     input image
  @param    decker  slit layout
  @param    npolys  [out] the number of trace polynomials determined
  @return   Set of trace polynomials that describe the orders
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_trace_nocpl(
        cpl_image       *   ima, 
        cr2res_decker       decker, 
        double              smoothfactor,
        int                 opening,
        int                 degree,
        int                 min_cluster)
{
    cpl_mask        *   mask ;
    cpl_table       *   clustertable ;
    cpl_image       *   labels ;
    cpl_table       *   trace_table ;

    /* Check Entries */
    if (ima == NULL) return NULL ;

    /* Initialise */

    /* Detect the orders in the image */
    mask = cr2res_trace_detect(ima, smoothfactor, opening, min_cluster) ;
    
    /* Detect the clusters */
    clustertable = cr2res_cluster_detect(mask, min_cluster) ;

    /* Create the labels image */
    labels = cr2res_trace_convert_cluster_to_labels(clustertable, 
            cpl_image_get_size_x(ima), cpl_image_get_size_y(ima)) ;
    cpl_table_delete(clustertable);

    /* Fit Labels Edges */
    if ((trace_table = cr2res_trace_fit(labels, degree)) == NULL) {
        cpl_msg_error(__func__, "Fit the order edges") ;
        cpl_image_delete(labels) ;
        return NULL ;
    }
    cpl_image_delete(labels) ;

    return trace_table ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief  Determine which pixels belong to an order
  @param    ima             input image
  @param    smoothfactor    
  @param    opening         Flag to apply opening to the orders
  @return   mask for pixels in orders
 */
/*----------------------------------------------------------------------------*/
cpl_mask * cr2res_trace_detect(
        cpl_image   *   ima,
        double          smoothfactor,
        int             opening,
        int             min_cluster)
{
    cpl_mask    *   mask ;
    cpl_mask    *   mask_kernel ;
    cpl_mask    *   new_mask ;
    cpl_mask    *   diff_mask ;
    cpl_mask    *   clean_mask ;

    /* Check entries */
    if (ima == NULL) return NULL ;
 
    /* TODO This needs to come from a static calibration, each band */
    int                     ordersep=80;
    /* TODO Set to read-noise later, also input-para */
    double                  thresh=10;

    /* Apply detection */
    cpl_msg_info(__func__, "Detect the signal") ;
    if ((mask = cr2res_trace_signal_detect(ima, ordersep, smoothfactor,
                    thresh)) == NULL) {
        cpl_msg_error(__func__, "Detection failed") ;
        return NULL ;
    }

    /* Apply a opening to join horizontally the close clusters */
    if (opening) {
        cpl_msg_info(__func__, "Apply Opening to cleanup the orders") ;
        mask_kernel = cpl_mask_new(5, 1) ;
        cpl_mask_not(mask_kernel);
        new_mask = cpl_mask_duplicate(mask) ;
        cpl_mask_filter(new_mask, mask, mask_kernel, CPL_FILTER_OPENING,
                CPL_BORDER_NOP) ;
        cpl_mask_delete(mask_kernel) ;

        /* Compute the difference */
        diff_mask = cpl_mask_duplicate(mask) ;
        cpl_mask_xor(diff_mask, new_mask) ;
        if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
            cpl_mask_save(diff_mask, "debug_diff_mask.fits", NULL, 
                    CPL_IO_CREATE);
        }
        cpl_mask_delete(diff_mask) ;

        cpl_mask_delete(mask) ;
        mask = new_mask ;
    }

    /* Clean the small blobs */
    cpl_msg_info(__func__, "Remove the small blobs (<= %d pixels)",
            min_cluster) ;
    if ((clean_mask = cr2res_trace_clean_blobs(mask, min_cluster)) == NULL) {
        cpl_msg_error(__func__, "Cannot clean the blobs") ;
        cpl_mask_delete(mask) ;
        return NULL ;
    }
    cpl_mask_delete(mask) ;

    /* Debug Saving */
    if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
        cpl_mask_save(clean_mask, "debug_mask.fits", NULL, CPL_IO_CREATE);
    }

    return clean_mask ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Find out which pixels belong to the same order, label orders
  @param    bin_ima     input binary image
  @return   image with labelled pixels.
 */
/*----------------------------------------------------------------------------*/
cpl_image * cr2res_trace_labelize(cpl_mask * mask)
{
    cpl_image   *   labels ;
    cpl_size        nlabels ;

    /* Check entries */
    if (mask == NULL) return NULL ;

    /* Call CPL function to labelise */
    if ((labels = cpl_image_labelise_mask_create(mask, &nlabels)) == NULL) {
        cpl_msg_error(__func__, "Failed to labelise") ;
        return NULL ;
    }
    
    /* Debug Saving */
    if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
        /* Analyse and dump orders */
        cpl_image_save(labels, "debug_labels.fits", CPL_TYPE_INT, NULL,
                CPL_IO_CREATE);
    }
    return labels ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Fit polynomials to pixel coordinates in each order
  @param    
  @param   
  @return 
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_trace_fit(
        cpl_image   *   labels, 
        int             degree)
{
    cpl_table       *   clustertable ;
    cpl_table       *   trace_table ;

    /* Create cluster table needed for fitting */
    clustertable = cr2res_trace_convert_labels_to_cluster(labels) ;

    /* Debug Saving */
    if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
        cpl_table_save(clustertable, NULL, NULL, "debug_cluster_table.fits", 
                CPL_IO_CREATE);
    }
    
    /* Fit the orders */
    if ((trace_table = cr2res_trace_orders_fit(clustertable, degree)) == NULL) {
        cpl_msg_error(__func__, "Failed to Fit") ;
        cpl_table_delete(clustertable);
        return NULL ;
    }
    cpl_table_delete(clustertable);

    /* Debug Saving */
    if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
        cpl_table_save(trace_table, NULL, NULL, "debug_trace_table.fits", 
                CPL_IO_CREATE);
    }
 
    /* Free and return */
    return trace_table ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Compare two traces
  @param    trace1  first trace
  @param    trace2  second trace
  @return   sqrt(sum(distances along x-axis ^2))
 */
/*----------------------------------------------------------------------------*/
cpl_vector * cr2res_trace_compare(
        cpl_table   *   trace1, 
        cpl_table   *   trace2)
{
    return NULL;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Combine two decker traces into an open trace
  @param    td_13   1-3-decker trace
  @param    td_24   2-4-decker trace
  @return   open trace
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_trace_combine(
        cpl_table   *   td_13, 
        cpl_table   *   td_24)
{
    return NULL;
}

/*----------------------------------------------------------------------------*/
/**
  @brief Make an image out of the trace solution
  @param trace_open (or NULL)
  @param trace_decker_1_3 (or NULL)
  @param trace_decker_2_4 (or NULL)
  @return image for trace visualization
 */
/*----------------------------------------------------------------------------*/
cpl_image * cr2res_trace_gen_image(
        cpl_table   *   trace, 
        int             nx,
        int             ny)
{
    cpl_image       *   out ;
    const cpl_array *   coeffs ;
    cpl_polynomial  *   poly ;
    int             *   pout ;
    int                 order, i ;
    cpl_size            j, y_pos ;

    /* Check entries */
    if (trace == NULL) return NULL ;
    if (nx < 1 || ny < 1) return NULL ;

    /* Create the empty image */
    out = cpl_image_new(nx, ny, CPL_TYPE_INT) ;
    pout = cpl_image_get_data_int(out) ;

    /* Loop on the traces */
    for (i=0 ; i<cpl_table_get_nrow(trace) ; i++) {
        /* Get the Order */
        order = cpl_table_get(trace, "Order", i, NULL) ;
        
        /* Get the Upper polynomial*/
        coeffs = cpl_table_get_array(trace, "Upper", i) ;
        poly = cpl_polynomial_new(1) ;
        for (j=0 ; j<cpl_array_get_size(coeffs) ; j++) 
            cpl_polynomial_set_coeff(poly, &j, cpl_array_get(coeffs, j, NULL)) ;

        /* Draw It  */
        for (j=0 ; j<nx ; j++) {
            y_pos = (cpl_size)cpl_polynomial_eval_1d(poly, (double)j+1, NULL) ;
            pout[j+(y_pos-1)*nx] = order ;
        }
        cpl_polynomial_delete(poly) ;

        /* Get the Lower polynomial*/
        coeffs = cpl_table_get_array(trace, "Lower", i) ;
        poly = cpl_polynomial_new(1) ;
        for (j=0 ; j<cpl_array_get_size(coeffs) ; j++) 
            cpl_polynomial_set_coeff(poly, &j, cpl_array_get(coeffs, j, NULL)) ;

        /* Draw It  */
        for (j=0 ; j<nx ; j++) {
            y_pos = (cpl_size)cpl_polynomial_eval_1d(poly, (double)j+1, NULL) ;
            pout[j+(y_pos-1)*nx] = order ;
        }
        cpl_polynomial_delete(poly) ;
    }
    return out;
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Detect the Orders signal
  @param
  @return
 */
/*----------------------------------------------------------------------------*/
static cpl_mask * cr2res_trace_signal_detect(
        const cpl_image *   image,
        int                 ordersep,
        double              smoothfactor,
        double              thresh)
{
    cpl_image       *   smimage ;
    int                 ordersep_loc ;
    cpl_matrix      *   kernel ;
    cpl_mask        *   mask ;

    /* Prepare the kernel used for median filtering */
    ordersep_loc = (int) (ordersep*smoothfactor);
    if (ordersep_loc % 2 == 0) ordersep_loc +=1;
    cpl_msg_debug(__func__, "Order separation: %d", ordersep_loc);
    kernel = cpl_matrix_new(ordersep_loc, 1);
    
    /* kernel should have normalized values */
    cpl_matrix_add_scalar(kernel, 1.0/((double)ordersep_loc)) ;

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

    /* the pixels we want are the ones with values below -thresh */
    cpl_image_subtract(smimage, image);
    mask=cpl_mask_new(cpl_image_get_size_x(image),cpl_image_get_size_y(image));
    cpl_mask_threshold_image(mask,smimage,-1*thresh,DBL_MAX,CPL_BINARY_0);
    cpl_image_delete(smimage) ;

    return mask ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Go through the orders and initiate the polynomial fit for each
  @param
  @return
 */
/*----------------------------------------------------------------------------*/
static cpl_table * cr2res_trace_orders_fit(
        cpl_table   *   clustertable,
        int             degree)
{
    cpl_array   *    fitparams_all;
    cpl_array   *    fitparams_upper;
    cpl_array   *    fitparams_lower;
    cpl_table   *    trace_table;
    cpl_table   *    edge_upper_table;
    cpl_table   *    edge_lower_table;
    cpl_table   *    order_table;
    cpl_size         nclusters_cur;
    int              i, nclusters;

    /* Check entries */
    if (clustertable == NULL) return NULL ;

    /* Create the output table */
    nclusters = cpl_table_get_column_max(clustertable, "clusters");
    trace_table = cpl_table_new(nclusters);
    cpl_table_new_column_array(trace_table, "All", CPL_TYPE_DOUBLE,degree+1) ;
    cpl_table_new_column_array(trace_table, "Upper", CPL_TYPE_DOUBLE,degree+1) ;
    cpl_table_new_column_array(trace_table, "Lower", CPL_TYPE_DOUBLE,degree+1) ;
    cpl_table_new_column(trace_table, "Order", CPL_TYPE_INT) ;

    /* Loop on the clusters */
    for (i=1 ; i<=nclusters ; i++) {
        cpl_table_set(trace_table, "Order", i-1, i);

        /* Select the pixels of the current cluster */
        nclusters_cur = cpl_table_and_selected_int(clustertable, "clusters",
                CPL_EQUAL_TO, i);
        cpl_msg_debug(__func__, "Cluster %d has %"CPL_SIZE_FORMAT" pixels",
                i, nclusters_cur);

        /* Extract the table with the current order pixels */
        order_table = cpl_table_extract_selected(clustertable);

        /* Fit the current order */
        fitparams_all = cr2res_trace_order_fit(order_table, degree);
        cpl_table_set_array(trace_table, "All", i-1, fitparams_all);
        cpl_array_delete(fitparams_all);

        /* Extract the edges of the current order pixels */
        cr2res_trace_extract_edges(order_table, &edge_lower_table,
                &edge_upper_table) ;

        /* Fit the upper edge of the current order */
        fitparams_upper = cr2res_trace_order_fit(edge_upper_table, degree);
        cpl_table_delete(edge_upper_table);
        cpl_table_set_array(trace_table, "Upper", i-1, fitparams_upper);
        cpl_array_delete(fitparams_upper);

        /* Fit the lower edge of the current order */
        fitparams_lower = cr2res_trace_order_fit(edge_lower_table, degree);
        cpl_table_delete(edge_lower_table);
        cpl_table_set_array(trace_table, "Lower", i-1, fitparams_lower);
        cpl_array_delete(fitparams_lower);

        cpl_table_delete(order_table);

        /* Reset the selection */
        cpl_table_select_all(clustertable);
    }

    return trace_table;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Fit a single order
  @param
  @return
 */
/*----------------------------------------------------------------------------*/
static cpl_array * cr2res_trace_order_fit(
        cpl_table   *   table,
        int             degree)
{
    cpl_matrix      *   x ;
    cpl_vector      *   y ;
    cpl_polynomial  *   poly1 ;
    cpl_array       *   result ;
    int             *   xs;
    int             *   ys;
    cpl_size            i, degree_local, n ;

    /* Check Entries */
    if (table == NULL) return NULL ;

    /* Initialise */
    n = cpl_table_get_nrow(table) ;
    degree_local = (cpl_size)degree ;

    /* Create Objects */
    x = cpl_matrix_new(1, n) ;
    y = cpl_vector_new(n) ;
    poly1 = cpl_polynomial_new(1) ;
    result = cpl_array_new(degree_local+1, CPL_TYPE_DOUBLE) ;

    xs = cpl_table_get_data_int(table,"xs");
    ys = cpl_table_get_data_int(table,"ys");
    for (i=0 ; i<n ; i++) {
        cpl_matrix_set(x, 0, i, xs[i]) ;
        cpl_vector_set(y, i, (double)ys[i]) ;
    }

    /* Apply the fit */
    if (cpl_polynomial_fit(poly1, x, NULL, y, NULL, CPL_FALSE, NULL, 
                &degree_local) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot fit the data") ;
        cpl_matrix_delete(x);
        cpl_vector_delete(y);
        cpl_polynomial_delete(poly1);
        cpl_array_delete(result);
        return NULL;
    }
    cpl_matrix_delete(x);
    cpl_vector_delete(y);

    /* Store the result */
    for (i=0 ; i<=degree_local ; i++) {
        cpl_array_set(result, i, cpl_polynomial_get_coeff(poly1, &i)) ;
    }
    cpl_polynomial_delete(poly1);
    return result;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Convert cluster table to labels image
  @param
  @return
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

    /* Create labels image  */
    labels = cpl_image_new(nx, ny, CPL_TYPE_INT);
    pxs = cpl_table_get_data_int_const(cluster, "xs") ;
    pys = cpl_table_get_data_int_const(cluster, "ys") ;
    pclusters = cpl_table_get_data_int_const(cluster, "clusters") ;
    for (i=0 ; i<cpl_table_get_nrow(cluster) ; i++)
        cpl_image_set(labels, pxs[i], pys[i], pclusters[i]);

    return labels ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Convert Labels image to the cluster table
  @param
  @return
 */
/*----------------------------------------------------------------------------*/
static cpl_table * cr2res_trace_convert_labels_to_cluster(cpl_image * labels) 
{
    cpl_image   *   non_zero_image ;
    const int   *   plabels ;
    cpl_table   *   clustertable ;
    int             nb_table_entries, nx, ny, i, j ;

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
    cpl_table_new_column(clustertable, "xs", CPL_TYPE_INT) ;
    cpl_table_new_column(clustertable, "ys", CPL_TYPE_INT) ;
    cpl_table_new_column(clustertable, "clusters", CPL_TYPE_INT) ;
    nb_table_entries = 0 ;
    for (j=0 ; j<ny ; j++) {
        for (i=0 ; i<nx ; i++) {
            if (plabels[i+j*nx] > 0.5) {
                cpl_table_set_int(clustertable, "xs", nb_table_entries, i+1) ;
                cpl_table_set_int(clustertable, "ys", nb_table_entries, j+1) ;
                cpl_table_set_int(clustertable, "clusters", nb_table_entries,
                        plabels[i+j*nx]) ;
                nb_table_entries ++ ;
            }
        }
    }
    return clustertable ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief  
  @param
  @return
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
        for (i=0 ; i<npix ; i++) 
            if (plabels[i] == curr_label) pix_count++ ;

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
  @brief  
  @param
  @return
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
    pxs = cpl_table_get_data_int_const(pixels_table, "xs") ;
    pys = cpl_table_get_data_int_const(pixels_table, "ys") ;
   
    /* Get the maximum x position */
    max_x = (int)cpl_table_get_column_max(pixels_table, "xs") ;

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
    pxs = cpl_table_get_data_int_const(upper_sel, "xs") ;
    pys = cpl_table_get_data_int_const(upper_sel, "ys") ;
    for (i=0 ; i<cpl_table_get_nrow(upper_sel) ; i++) 
        if (max_y[pxs[i]-1] == pys[i])
            cpl_table_select_row(upper_sel, i) ;
    cpl_free(max_y) ;
    *edge_upper_table = cpl_table_extract_selected(upper_sel) ;
    cpl_table_delete(upper_sel) ;

    /* Lower Edge extraction */
    lower_sel = cpl_table_duplicate(pixels_table) ;
    cpl_table_unselect_all(lower_sel) ;
    pxs = cpl_table_get_data_int_const(lower_sel, "xs") ;
    pys = cpl_table_get_data_int_const(lower_sel, "ys") ;
    for (i=0 ; i<cpl_table_get_nrow(lower_sel) ; i++) 
        if (min_y[pxs[i]-1] == pys[i])
            cpl_table_select_row(lower_sel, i) ;
    cpl_free(min_y) ;
    *edge_lower_table = cpl_table_extract_selected(lower_sel) ;
    cpl_table_delete(lower_sel) ;

    return CPL_ERROR_NONE ;
}

