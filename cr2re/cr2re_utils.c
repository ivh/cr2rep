/*
 * This file is part of the CR2RE Pipeline
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

#include "cr2re_utils.h"

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2re_utils     Miscellaneous Utilities
 */
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/**
  @brief    Detect the Orders signal
  @param
  @return
 */
/*----------------------------------------------------------------------------*/
cpl_mask * cr2re_signal_detect(
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
    cpl_image_save(smimage, "smimage.fits", CPL_TYPE_DOUBLE, NULL, CPL_IO_CREATE);

    /* save in smimage since image is static */
    /* tis means the pixels we want are the ones with values below -thresh */
    cpl_image_subtract(smimage,image);

    mask = cpl_mask_threshold_image_create(image,-thresh,DBL_MAX);
    cpl_mask_save(mask, "mask.fits", NULL, CPL_IO_CREATE);
    cpl_image_delete(smimage) ;


    return mask ;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Fit polynomial to the detected orders
  @param
  @return
 */
/*----------------------------------------------------------------------------*/


cpl_error_code cr2re_order_fit(int *xs, int *ys){
/* fit a single order */


}

cpl_error_code cr2re_orders_fit(cpl_table *table, cpl_table *fittable){

    int i, n;
    cpl_size nclusters_cur;

    n=cpl_table_get_column_max(table,"clusters");
    for (i=1;i <= n ;i++){
        nclusters_cur = cpl_table_and_selected_int(table,"clusters",CPL_EQUAL_TO,i);
        cpl_msg_debug(__func__, "Cluster %d has %d pixels", i, nclusters_cur);

    }

    fittable = cpl_table_new(55);

    return 0;
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
const char * cr2re_get_license(void)
{
    const char  *   cr2re_license = 
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
    return cr2re_license ;
}

/**@}*/

