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
#include <cpl.h>

#include "cr2res_utils.h"

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_utils     Miscellaneous Utilities
 */
/*----------------------------------------------------------------------------*/

/**@{*/

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
    cpl_ensure_code(ron > 0., CPL_ERROR_ILLEGAL_INPUT);

    *ima_errs = cpl_image_duplicate(ima_data);
    /* set negative values (= zero measurable electrons) to read out noise */
    cpl_image_threshold(*ima_errs, 0., DBL_MAX, ron, ron);

    /* err_ADU = sqrt(counts/gain + ron * ron)*/

    cpl_image_divide_scalar(*ima_errs, gain);
    cpl_image_add_scalar(*ima_errs, ron * ron);
    cpl_image_power(*ima_errs, 0.5);

    return cpl_error_get_code();
}

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
    double      *   pout ;
    int             i ;

    /* Allocate the output vector */
    out = cpl_vector_new(vector_size) ;
    pout = cpl_vector_get_data(out) ;

    /* Loop on the vector */
    for (i=0 ; i<vector_size ; i++) {
        pout[i] = (cpl_polynomial_eval_1d(trace1, (double)(i+1), NULL) +
            cpl_polynomial_eval_1d(trace1, (double)(i+1), NULL)) / 2.0 ;
    }
    return out ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Turn a single table column into a polynomial
  @param    trace table
  @param    column name
  @return   cpl_polynomial

  Read a table column as doubles and assign them as coefficients to a
  cpl_polynomial.
  The index of the column corresponds to the degree of the coefficient.

  The allocated cpl_polynomial will need to be destroyed by the caller.
 */
/*----------------------------------------------------------------------------*/
cpl_polynomial * cr2res_trace_column_to_polynomial(
            cpl_table   *   trace,
            char        *   col_name) {
    cpl_size i;
    double coeff;
    int flag;
    int polyorder = cpl_table_get_nrow(trace);
    cpl_polynomial * poly = cpl_polynomial_new(1);

    for (i=0;i<polyorder;i++){
        coeff = cpl_table_get_double(trace, col_name, i, &flag) ;
        /* TODO: check return flag?*/
        cpl_polynomial_set_coeff(poly, &i, coeff) ;
    }
    return poly;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Select the columns that belong to the same order as polynomials
  @param    open trace cpl_table
  @param    order number
  @return   array of two polynomials

  The polynomials will need to be destroyed by the caller.

 */
/*----------------------------------------------------------------------------*/
cpl_polynomial ** cr2es_trace_open_get_polynomials(
            cpl_table   *   trace,
            cpl_size        order_nb ) {

    cpl_polynomial ** polys ;
    char * col_name;

    /* Allocate the returned pointer */
    polys = cpl_malloc(2 * sizeof(cpl_polynomial*)) ;

    col_name = cpl_sprintf("%02d_Upper",(int)order_nb);
    polys[0] = cr2res_trace_column_to_polynomial(trace, col_name);
    cpl_free(col_name);

    col_name = cpl_sprintf("%02d_Lower",(int)order_nb);
    polys[1] = cr2res_trace_column_to_polynomial(trace, col_name);
    cpl_free(col_name);

    return polys;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Count and return the order numbers in a trace table
  @param    trace cpl_table
  @param    number of orders (output)
  @return   int array of order numbers

The int array will need to be freed by the caller.

 */
/*----------------------------------------------------------------------------*/
int * cr2res_trace_get_order_numbers(
        cpl_table * trace, int * nb_orders) {

    cpl_array * col_names;
    cpl_size ncols;
    const char * col_name;
    char * numstr;
    numstr = cpl_malloc(2*sizeof(char)) ;
    int i,j;
    int * order_numbers;
    int * order_indices;
    order_indices = cpl_calloc(128, sizeof(int));

    col_names = cpl_table_get_column_names(trace);
    ncols = cpl_array_get_size(col_names);
    for (i=0;i<ncols;i++){
        col_name = cpl_array_get_string(col_names, i);
        memcpy(numstr,col_name,2*sizeof(char));
        j = atoi(numstr);
        order_indices[j] = 1;
    }
    *nb_orders=0;
    for (i=0;i<128;i++){
        *nb_orders += order_indices[i];
    }
    order_numbers = cpl_malloc(*nb_orders * sizeof(int));
    j=0;
    for (i=0;i<128;i++){
        if (order_indices[i]==1){
            order_numbers[j]=i;
            j++;
        }
    }

    cpl_free(numstr);
    cpl_free(order_indices);
    cpl_array_delete(col_names);
    return order_numbers;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Select the columns that belong to the same order as polynomials
  @param    decker trace cpl_table
  @param    order number
  @return   array of four polynomials
 */
/*----------------------------------------------------------------------------*/
cpl_polynomial ** cr2es_trace_decker_get_polynomials(
            cpl_table * trace, cpl_size order_nb ) {


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
