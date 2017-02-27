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

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_trace		Trace Functions
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief  Main function for running all parts of the trace algorithm 
  @param    ima     input image
  @param    decker  slit layout
  @param    npolys  [out] the number of trace polynomials determined
  @return   Set of trace polynomials that describe the orders
 */
/*----------------------------------------------------------------------------*/
cpl_polynomial ** cr2res_trace(
        cpl_image       *   ima, 
        cr2res_decker       decker, 
        double              smoothfactor,
        int                 closing,
        int                 cpl_lab,
        int             *   npolys)
{
    cpl_mask        *   mask ;
    cpl_mask        *   mask_kernel ;
    cpl_mask        *   joined_mask ;
    cpl_image       *   labels ;
	int				*	plabels ;
    cpl_apertures   *   aperts ;
    cpl_table       *   clustertable ;
	int				*	pxs ;
	int				*	pys ;
	int				*	pclusters ;
    cpl_table       *   fittable ;
    int                 i, j, nlabels, nx, ny, counts, qc_nclusters ;

    /* Check Entries */
    if (ima == NULL) return NULL ;

    /* Initialise */

    /* TODO This needs to come from a static calibration, each band */
    int                     ordersep=80;
    /* TODO Set to read-noise later, also input-para */
    double                  thresh=10;

    /* Detect the orders in the image */
    mask = cr2res_signal_detect(image, ordersep, smoothfactor, thresh) ;

    /* Apply a closing to join horizontally the close clusters */
    if (closing) {
        mask_kernel = cpl_mask_new(3, 1) ;
        cpl_mask_not(mask_kernel);
        joined_mask = cpl_mask_duplicate(mask) ;
        cpl_mask_filter(joined_mask, mask, mask_kernel, CPL_FILTER_CLOSING,
                CPL_BORDER_COPY) ;
        cpl_mask_delete(mask_kernel) ;
        cpl_mask_delete(mask) ;
        mask = joined_mask ;
    }
    /* cpl_mask_save(mask, "mask.fits", NULL, CPL_IO_CREATE); */

    /* Labelization */
    if (cpl_lab) {
        /* Use CPL method */
        labels = cpl_image_labelise_mask_create(mask, &nlabels);
        aperts = cpl_apertures_new_from_image(image, labels);
        cpl_apertures_dump(aperts, stdout) ;
        cpl_apertures_delete(aperts) ;

        /* Create cluster table */
        clustertable = cpl_table_new(cpl_mask_count(mask));
        cpl_table_new_column(clustertable, "xs", CPL_TYPE_INT) ;
        cpl_table_new_column(clustertable, "ys", CPL_TYPE_INT) ;
        cpl_table_new_column(clustertable, "clusters", CPL_TYPE_INT) ;
        plabels = cpl_image_get_data_int_const(labels) ;
        count = 0 ;
        nx = cpl_image_get_size_x(image) ;
        ny = cpl_image_get_size_y(image) ;
        for (j=0 ; j<ny ; j++) {
            for (i=0 ; i<nx ; i++) {
                if (plabels[i+j*nx] > 0) {
                    cpl_table_set_int(clustertable, "xs", count, i+1) ;
                    cpl_table_set_int(clustertable, "ys", count, j+1) ;
                    cpl_table_set_int(clustertable, "clusters", count,
                            plabels[i+j*nx]) ;
                    count ++ ;
                }
            }
        }
    } else {
        /* Detect the clusters */
        clustertable = cr2re_cluster_detect(mask, mincluster) ;

        /* Create labels image  */
        labels = cpl_image_new(cpl_image_get_size_x(image),
                cpl_image_get_size_y(image), CPL_TYPE_INT);
        pxs = cpl_table_get_data_int_const(clustertable, "xs") ;
        pys = cpl_table_get_data_int_const(clustertable, "ys") ;
        pclusters = cpl_table_get_data_int_const(clustertable, "clusters") ;
        for (i=0 ; i<cpl_table_get_nrow(clustertable) ; i++)
            cpl_image_set(labels, pxs[i], pys[i], pclusters[i]);
    }
    cpl_mask_delete(mask);

    /* Save Labels image */
    cpl_image_save(labels, "labels.fits", CPL_TYPE_INT, NULL, CPL_IO_CREATE);
    cpl_image_delete(labels) ;

    /* Nb of clusters */
    qc_nclusters = cpl_table_get_column_max(clustertable, "clusters");
    cpl_msg_info(__func__, "Number of clusters: %d", qc_nclusters);
    cpl_table_save(clustertable, NULL, NULL, "clustertable.fits",CPL_IO_CREATE);

    /* Fit ys */
    fittable = cr2res_orders_fit(clustertable);
    cpl_table_delete(clustertable);
    cpl_table_save(fittable, NULL, NULL, "fittable.fits", CPL_IO_CREATE);



    return NULL;
}

/*----------------------------------------------------------------------------*/
/**
  @brief  Determine which pixels belong to an order
  @param    ima                 input image
  @param    min_cluster_size    smaller clusters are disregarded
  @return   binary image with 1s for pixels in orders
 */
/*----------------------------------------------------------------------------*/
cpl_image * cr2res_trace_detect(
        cpl_image   *   ima, 
        int             min_cluster_size)
{
    return NULL;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Find out which pixels belong to the same order, label orders
  @param    bin_ima     input binary image
  @return   image with labelled pixels.
 */
/*----------------------------------------------------------------------------*/
cpl_image * cr2res_trace_labelize(cpl_image * bin_ima)
{
    return NULL;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Fit polynomials to pixel coordinates in each order
  @param    ima     input image with labels
  @param    degree  polynomial degree
  @return   fit result polynomials
 */
/*----------------------------------------------------------------------------*/
cpl_polynomial ** cr2res_trace_fit(
        cpl_image   *   ima, 
        int             degree)
{
    return NULL;
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
        cpl_table   *   trace_open,
        cpl_table   *   trace_decker_1_3,
        cpl_table   *   trace_decker_2_4)
{
    return NULL;
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Detect the Orders signal
  @param
  @return
 */
/*----------------------------------------------------------------------------*/
static cpl_mask * cr2res_signal_detect(
        const cpl_image     *   image,
        int                     ordersep,
        double                  smoothfactor,
        double                  thresh)
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

    cpl_image_save(smimage, "smimage.fits", CPL_TYPE_DOUBLE, NULL,
            CPL_IO_CREATE);

    /* save in smimage since image is static */
    /* tis means the pixels we want are the ones with values below -thresh */
    cpl_image_subtract(smimage,image);

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
static cpl_table * cr2res_orders_fit(
    cpl_table         *    clustertable)
{
    int                    i, n;
    cpl_array         *    fitparams;
    cpl_table         *    fittable;
    cpl_table         *    seltable;
    cpl_size               nclusters_cur;

    n = cpl_table_get_column_max(clustertable,"clusters");
    fittable = cpl_table_new(n);
    cpl_table_new_column_array(fittable,"fitparams",CPL_TYPE_DOUBLE,4);

    for (i=1;i <= n ;i++){
        nclusters_cur = cpl_table_and_selected_int(clustertable,"clusters",
                CPL_EQUAL_TO,i);
        cpl_msg_debug(__func__, "Cluster %d has %"CPL_SIZE_FORMAT" pixels",
                i, nclusters_cur);
        seltable = cpl_table_extract_selected(clustertable);
        fitparams = cr2res_order_fit(seltable,nclusters_cur);
        cpl_table_set_array(fittable,"fitparams",i-1,fitparams);
        cpl_array_delete(fitparams);
        cpl_table_delete(seltable);
        cpl_table_select_all(clustertable);
    }

    return fittable;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Fit a single order
  @param
  @return
 */
/*----------------------------------------------------------------------------*/
static cpl_array * cr2res_order_fit(
        cpl_table             *    table,
        cpl_size                   n)
{
    cpl_size                   i;
    int                   *    xs;
    int                   *    ys;
    const cpl_size             deg=3;
    cpl_matrix            *    x = cpl_matrix_new(1,n);
    cpl_vector            *    y = cpl_vector_new(n);
    cpl_polynomial        *    poly1 = cpl_polynomial_new(1);
    cpl_array             *    result = cpl_array_new(deg+1,CPL_TYPE_DOUBLE);

    xs = cpl_table_get_data_int(table,"xs");
    ys = cpl_table_get_data_int(table,"ys");
    for (i=0;i<n;i++){
        //cpl_msg_debug(__func__, "i,xs: %d %d %d", i,xs[i]);
        cpl_matrix_set(x,0,i,xs[i]);
    }
    for (i=0;i<n;i++){
        cpl_vector_set(y,i,(double)ys[i]);
    }

    if ( cpl_polynomial_fit(poly1, x, NULL, y, NULL,
            CPL_FALSE, NULL, &deg) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot fit the data") ;
        cpl_error_set(__func__, CPL_ERROR_CONTINUE) ;
        cpl_polynomial_delete(poly1);
        cpl_matrix_delete(x);
        cpl_vector_delete(y);
        cpl_array_delete(result);
        return NULL;
    }

    for (i=0;i<=deg;i++){
        cpl_array_set(result,i,cpl_polynomial_get_coeff(poly1,&i));
    }
    cpl_polynomial_delete(poly1);
    cpl_matrix_delete(x);
    cpl_vector_delete(y);
    return result;
}



