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
#include "cr2res_dfs.h"

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_utils     Miscellaneous Utilities
 */
/*----------------------------------------------------------------------------*/

/**@{*/


/*----------------------------------------------------------------------------*/
/**
  @brief
  @param    ycen
  @return
 */
/*----------------------------------------------------------------------------*/
int * cr2res_vector_get_int(
    const cpl_vector    * ycen)
{
    int         * ycen_int;
    int         i, lenx;

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

    lenx = cpl_vector_get_size(ycen);
    ycen_rest = cpl_malloc(lenx*sizeof(double));
    for (i=0 ; i<lenx ; i++){
         ycen_rest[i] = fmod( cpl_vector_get(ycen,i), 1.0) ;
    }

   return ycen_rest;

}
/*----------------------------------------------------------------------------*/
/**
  @brief    Cut a bent order into a rectangle, shifting columns
  @param    img_in
  @param    ycen
  @param    height
  @param    yecen_int
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
            cpl_msg_error(__func__,"Cannot extract and copy column %d, %d %d, %s",
                            i, ymin, ymax, cpl_error_get_where());
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
   @brief   Get the TRACE_WAVE table orders list
   @param   tab         A TRACE_WAVE table
   @param   nb_orders   The output array size
   @return  the array of orders or NULL in error case
    Needs to be deallocated with cpl_free()
 */
/*----------------------------------------------------------------------------*/
int * cr2res_get_trace_table_orders(
        const cpl_table     *   trace_wave,
        int                 *   nb_orders)
{
    int         *   orders_big ;
    int         *   orders ;
    int             nb_orders_loc, nb_orders_max, new_order, cur_order ;
    cpl_size        nrows, i, j ;


    /* Check Entries */
    if (trace_wave == NULL) return NULL ;

    /* Initialise */
    nrows = cpl_table_get_nrow(trace_wave) ;
    nb_orders_max = 100 ;

    /* Allocate the orders_big array */
    orders_big = cpl_malloc(nb_orders_max * sizeof(int)) ;

    /* Loop on the table rows */
    nb_orders_loc = 0 ;
    for (i=0 ; i<nrows ; i++) {
        cur_order = cpl_table_get(trace_wave, CR2RES_COL_ORDER, i, NULL) ;

        /* Is this order already there ? */
        new_order = 1 ;
        for (j=0 ; j<nb_orders_loc ; j++)
            if (cur_order == orders_big[j])
                new_order = 0 ;

        /* Store the new order */
        if (new_order) {
            orders_big[nb_orders_loc] = cur_order ;
            nb_orders_loc ++ ;
        }
    }

    /* Resize */
    orders = cpl_malloc(nb_orders_loc * sizeof(int)) ;
    for (i=0 ; i<nb_orders_loc ; i++) orders[i] = orders_big[i] ;
    cpl_free(orders_big) ;

    /* Return */
    *nb_orders = nb_orders_loc ;
    return orders ;
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
            strcmp(poly_column, CR2RES_COL_ALL)) return NULL ;

    /* Get Table index from order and trace */
    index = cr2res_get_trace_table_index(trace_wave, order, trace_nb) ;

    /* Read the Table */
    wave_arr = cpl_table_get_array(trace_wave, poly_column, index) ;

    /* Convert to Polynomial */
	wave_poly = cr2res_convert_array_to_poly(wave_arr) ;

    return wave_poly ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Compute the polynomial from boundaries
  @param    wmin    First pixel wavelength
  @param    wmax    Last pixel wavelength
  @return   the polynomial or NULL in error case

  The returned polynomial must be deallocated with cpl_polynomial_delete()
 */
/*----------------------------------------------------------------------------*/
cpl_polynomial * cr2res_wlestimate_compute(
        double          wmin,
        double          wmax)
{
    cpl_polynomial      *   poly ;
    cpl_size                power ;
    double                  a, b ;
    int                     nbpix = CR2RES_DETECTOR_SIZE ;

    /* Test entries */
    if (wmin <= 0.0) return NULL ;
    if (wmax <= 0.0) return NULL ;
    if (wmax <= wmin) return NULL ;

    /* Compute polynomial coeffs */
    a = (wmax - wmin) / (nbpix-1) ;
    b = wmin - a ;

    /* Create polynomial */
    poly = cpl_polynomial_new(1) ;
    power = 0 ;
    cpl_polynomial_set_coeff(poly, &power, b) ;
    power = 1 ;
    cpl_polynomial_set_coeff(poly, &power, a) ;

    return poly ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Convert the order to the keyword index
  @param    order   Order (-49 to 50)
  @return   the order index or a negative value in error case
            (00 to 99)
 */
/*----------------------------------------------------------------------------*/
int cr2res_convert_order_to_idx(int order)
{
    /* Check entries */
    if (order < -49 || order > 50) return -1 ;

    /* Conversion order <-> keyword Index */
    if (order < 0)  return order + 100 ;
    else            return order ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Convert the keyword index to the order
  @param    order_idx   the order index (00 to 99)
  @return   Order (-50 to 50)
 */
/*----------------------------------------------------------------------------*/
int cr2res_convert_idx_to_order(int order_idx)
{
    /* Check entries */
    if (order_idx < 0 || order_idx > 99) return -1 ;

    /* Conversion order <-> keyword Index */
    if (order_idx > 50) return order_idx - 100 ;
    else                return order_idx ;
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
    cpl_size            i ;

    /* Test entries */
    if (arr == NULL) return NULL ;

    /* Create Output poly */
	out = cpl_polynomial_new(1) ;

    /* Fill it  */
	for (i=0 ; i<cpl_array_get_size(arr) ; i++)
		cpl_polynomial_set_coeff(out, &i, cpl_array_get(arr, i, NULL)) ;

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
    if (poly == NULL) return NULL ;

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
