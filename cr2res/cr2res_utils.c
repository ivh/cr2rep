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
  @param    poly1   First trace
  @param    poly2   Second trace
  @return
  The returned int is the rounded-up mean difference between the two
  input polynomials, evaluated on a vector from 1 to vector_size.

 */
/*----------------------------------------------------------------------------*/
int cr2res_trace_compute_height(
        cpl_polynomial  *   trace1,
        cpl_polynomial  *   trace2,
        int                 vector_size)
{
    int height;
    cpl_polynomial * diff_poly;
    cpl_vector  *   diff_vec ;

    diff_poly =  cpl_polynomial_new(1);
    diff_vec =  cpl_vector_new(vector_size);
    cpl_polynomial_subtract(diff_poly, trace2, trace1);
    cpl_vector_fill_polynomial(diff_vec, diff_poly, 1, 1);
    height = (int)ceil(fabs( cpl_vector_get_mean(diff_vec) ));

    if (cpl_vector_get_stdev(diff_vec) > 5){ // TODO: make this not hardcoded?
        cpl_msg_warning(__func__, "Stdev of extraction height is large.");
    }

    return height;
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
