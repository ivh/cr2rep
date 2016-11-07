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

#include "cr2res_utils.h"

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_utils     Miscellaneous Utilities
 */
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/**
  @brief    Detect the Orders signal
  @param
  @return
 */
/*----------------------------------------------------------------------------*/
cpl_mask * cr2res_signal_detect(
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
  @brief    Fit a single order
  @param
  @return
 */
/*----------------------------------------------------------------------------*/
cpl_array * cr2res_order_fit(
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


/*----------------------------------------------------------------------------*/
/**
  @brief    Go through the orders and initiate the polynomial fit for each
  @param
  @return
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_orders_fit(
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
        nclusters_cur = cpl_table_and_selected_int(clustertable,"clusters",CPL_EQUAL_TO,i);
        cpl_msg_debug(__func__, "Cluster %d has %d pixels", i, nclusters_cur);
        seltable = cpl_table_extract_selected(clustertable);
        fitparams = cr2res_order_fit(seltable,nclusters_cur);
        cpl_table_set_array(fittable,"fitparams",i-1,fitparams);
        cpl_array_delete(fitparams);
        cpl_table_delete(seltable);
        cpl_table_select_all(clustertable);
    }


    return fittable;
}

/**@{*/

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
        "This file is part of the CR2RE Instrument Pipeline\n"
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
