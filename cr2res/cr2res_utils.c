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

/* TODO 1 : Move/Rename to cr2res_io.c */
/*----------------------------------------------------------------------------*/
/**
  @brief    Get the DITS from a frame set
  @param    set     Input frame set
  @return   the DITS or NULL in error case
 */
/*----------------------------------------------------------------------------*/
cpl_vector * cr2res_read_dits(const cpl_frameset * in)
{
    cpl_vector          *   dits ;
    cpl_propertylist    *   plist ;
    cpl_size                i ;

    /* Check entries */
    if (in == NULL) return NULL ;

    /* Allocate the vector */
    dits = cpl_vector_new(cpl_frameset_get_size(in)) ;

    /* Loop on the frames */
    for (i=0 ; i< cpl_vector_get_size(dits) ; i++) {
        plist = cpl_propertylist_load(cpl_frame_get_filename(
                    cpl_frameset_get_position_const(in, i)), 0) ;
        cpl_vector_set(dits, i, cr2res_pfits_get_dit(plist)) ;
        cpl_propertylist_delete(plist) ;
    }

    return dits ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the decker positions from a frame set
  @param    set     Input frame set
  @return   the DECKER positions or NULL in error case
 */
/*----------------------------------------------------------------------------*/
cr2res_decker * cr2res_decker_read_positions(const cpl_frameset * in)
{
    cr2res_decker  		*   out ;
    cpl_propertylist    *   plist ;
    const char          *   fname ;
    cpl_size                nframes, i ;

    /* Check entries */
    if (in == NULL) return NULL ;

    /* Initialise */
    nframes = cpl_frameset_get_size(in) ;

    /* Allocate the vector */
    out = cpl_malloc(nframes * sizeof(cr2res_decker)) ;

    /* Loop on the frames */
    for (i=0 ; i< nframes ; i++) {
        plist = cpl_propertylist_load(cpl_frame_get_filename(
                    cpl_frameset_get_position_const(in, i)), 0) ;
        out[i] = cr2res_pfits_get_decker_position(plist) ;
        cpl_propertylist_delete(plist) ;
    }
    return out ;
}

/*----------------------------------------------------------------------------*/
/**
   @brief   Extract the frames with the given tag and Decker position
   @param   in      A non-empty frameset
   @param   tag     The tag of the requested frames
   @param   decker  CR2RES_DECKER_NONE,CR2RES_DECKER_1_3 or CR2RES_DECKER_2_4
   @return  The newly created frameset or NULL on error

   The returned frameset must be de allocated with cpl_frameset_delete()
 */
/*----------------------------------------------------------------------------*/
cpl_frameset * cr2res_extract_decker_frameset(
        const cpl_frameset  *   in,
        const char          *   tag,
        cr2res_decker           decker)
{
    cpl_frameset        *   out ;
    const cpl_frame     *   cur_frame ;
    cpl_frame           *   loc_frame ;
    int                     nbframes;
    cpl_propertylist    *   plist ;
    int                     i ;

    /* Test entries */
    if (in == NULL) return NULL ;
    if (tag == NULL) return NULL ;
    if (decker != CR2RES_DECKER_NONE && decker != CR2RES_DECKER_1_3 && 
            decker != CR2RES_DECKER_2_4)  
        return NULL ;

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
            /* Get the propertylist */
            plist = cpl_propertylist_load(cpl_frame_get_filename(cur_frame), 0);
            if (cr2res_pfits_get_decker_position(plist) == decker) {
                loc_frame = cpl_frame_duplicate(cur_frame) ;
                cpl_frameset_insert(out, loc_frame) ;
            }
            cpl_propertylist_delete(plist) ;
        }
    }

    /* No matching frame found */
    if (cpl_frameset_get_size(out) == 0){
        cpl_frameset_delete(out) ;
        return NULL ;
    }
    return out ;
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

/* End TODO 1 */

/* TODO 2 : Move/Rename to cr2res_wave.c */


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

/* End TODO 2 */


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
            /* cpl_polynomial_delete(out) ; */
            /* return NULL ; */
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

/* TODO 3 : Remove those unused functions ? */

/*----------------------------------------------------------------------------*/
/**
  @brief    Extract both traces from 2 images for demodulation
  @param im1                Image at angle 1
  @param im2                Image at angle 2
  @param trace_wave         Trace Wave Table
  @param order              Order to extract
  @param trace_a_angle_1    OUTPUT data of Trace A from observation at angle 1
  @param trace_b_angle_1    OUTPUT data of Trace B from observation at angle 1
  @param trace_a_angle_2    OUTPUT data of Trace A from observation at angle 2
  @param trace_b_angle_2    OUTPUT data of Trace B from observation at angle 2
  @return   0 if successful, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
int cr2res_pol_demod_extract(
        hdrl_image      *   im1, 
        hdrl_image      *   im2,
        cpl_table       *   trace_wave, 
        int                 order, 
        cpl_bivector    **  trace_a_angle_1,
        cpl_bivector    **  trace_b_angle_1, 
        cpl_bivector    **  trace_a_angle_2,
        cpl_bivector    **  trace_b_angle_2)
{
    const int Trace_A = 1;
    const int Trace_B = 2;
    int extr_res = 0;

    // Define extraction parameters
    const int height = 10;
    const int swath = 64;
    const int oversample = 10;
    const double smooth_slit = 1.0;

    // We dont care about those, but need to delete them
    cpl_vector *slit_func;
    hdrl_image *model;

    if ((cr2res_get_trace_table_index(trace_wave, order, Trace_A) == -1) \
        || (cr2res_get_trace_table_index(trace_wave, order, Trace_B) == -1))
    {
        return -1;
    } else {
        // extract values for the current order
        extr_res = cr2res_extract_slitdec_vert(im1, trace_wave, order, 
                Trace_A, height, swath, oversample, smooth_slit, &slit_func, 
                trace_a_angle_1, &model);
        cpl_vector_delete(slit_func);
        hdrl_image_delete(model);
        extr_res = cr2res_extract_slitdec_vert(im1, trace_wave, order, Trace_B,
                height, swath, oversample, smooth_slit, &slit_func, 
                trace_b_angle_1, &model);
        cpl_vector_delete(slit_func);
        hdrl_image_delete(model);
        extr_res = cr2res_extract_slitdec_vert(im2, trace_wave, order, Trace_A,
                height, swath, oversample, smooth_slit, &slit_func, 
                trace_a_angle_2, &model);
        cpl_vector_delete(slit_func);
        hdrl_image_delete(model);
        extr_res = cr2res_extract_slitdec_vert(im2, trace_wave, order, Trace_B,
                height, swath, oversample, smooth_slit, &slit_func, 
                trace_b_angle_2, &model);
        cpl_vector_delete(slit_func);
        hdrl_image_delete(model);
    }
    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Demodulate one order
  @param    trace_a_angle_1     data of Trace A from observation at angle 1
  @param    trace_b_angle_1     data of Trace B from observation at angle 1
  @param    trace_a_angle_2     data of Trace A from observation at angle 2
  @param    trace_b_angle_2     data of Trace B from observation at angle 2
  @param    sx                  OUTPUT Stokes X, e.g. Stokes V
  @param    si                  OUTPUT Stokes I, i.e. Intensity
  @param    sn                  OUTPUT Stokes N, i.e. the null signal
  @return   0 if successful, -1 otherwise

       R^(1/2) - 1           Trace_B(angl1)   Trace_A(angl2)
  x =  ----------- where R = -------------- * --------------
       R^(1/2) + 1           Trace_A(angl1)   Trace_B(angl2)
    
  and x is V/I or Q/I or U/I. The Stokes I is as before:
    
  I = 1/4 * [Trace_A(angl1) + Trace_B(angl1) + Trace_A(angl2) + Trace_B(angl2)]
    
      RN^(1/4) - 1
  N = -----------
      RN^(1/4) + 1
    
             Fiber_B(angl1)   Fiber_A(angl2)   Fiber_A(angl3)   Fiber_B(angl4)
  where RN = -------------- * -------------- * -------------- * --------------
             Fiber_A(angl1)   Fiber_B(angl2)   Fiber_B(angl3)   Fiber_A(angl4)
 */
/*----------------------------------------------------------------------------*/
int cr2res_pol_demod_order(
        cpl_bivector    *   trace_a_angle_1,
        cpl_bivector    *   trace_b_angle_1, 
        cpl_bivector    *   trace_a_angle_2,
        cpl_bivector    *   trace_b_angle_2, 
        cpl_vector      **  sx, 
        cpl_vector      **  si,
        cpl_vector      **  sn)
{
    cpl_vector  *   R, 
                *   RN ;

    /* Check Inputs */

    /* Initialise */

    /* Calculate R */
    R = cpl_vector_duplicate(cpl_bivector_get_x(trace_b_angle_1));
    cpl_vector_multiply(R, cpl_bivector_get_x(trace_a_angle_2));
    cpl_vector_divide(R, cpl_bivector_get_x(trace_a_angle_1));
    cpl_vector_divide(R, cpl_bivector_get_x(trace_b_angle_2));
    cpl_vector_power(R, 0.5);
     
    /* Calculate Stokes X */
    *sx = cpl_vector_duplicate(R);
    cpl_vector_subtract_scalar(*sx, 1.);
    cpl_vector_add_scalar(R, 1.);
    cpl_vector_divide(*sx, R);

    /* Calculate Stokes I */
    *si = cpl_vector_duplicate(cpl_bivector_get_x(trace_a_angle_1));
    cpl_vector_add(*si, cpl_bivector_get_x(trace_b_angle_1));
    cpl_vector_add(*si, cpl_bivector_get_x(trace_a_angle_2));
    cpl_vector_add(*si, cpl_bivector_get_x(trace_b_angle_2));
    cpl_vector_multiply_scalar(*si, 0.25);

    /* Calculate RN */
    RN = cpl_vector_duplicate(cpl_bivector_get_x(trace_b_angle_1));
    cpl_vector_multiply(RN, cpl_bivector_get_x(trace_a_angle_2));
    cpl_vector_divide(  RN, cpl_bivector_get_x(trace_a_angle_1));
    cpl_vector_divide(  RN, cpl_bivector_get_x(trace_b_angle_2));
    cpl_vector_power(   RN, 0.5);

    /* Calculate Stokes N */
    *sn = cpl_vector_duplicate(RN);
    cpl_vector_subtract_scalar(*sn, 1.);
    cpl_vector_add_scalar(RN, 1.);
    cpl_vector_divide(*sn, RN);

    /* Delete temp values */
    cpl_vector_delete(R);
    cpl_vector_delete(RN);
    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Demodulate frames
  @param    sum1        Image 1, Trace_A and Trace_B at angle 1
  @param    sum2        Image 2, Trace_A and Trace_B at angle 2
  @param    trace_wave  Trace Wave Table
  @return   Table with columns "POL_X_1" for each order and stokes parameter, 
            NULL on error
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_demod(
        hdrl_image  *   sum1, 
        hdrl_image  *   sum2, 
        cpl_table   *   trace_wave)
{
    // Check input
    if (sum1 == NULL || sum2 == NULL || trace_wave == NULL) return NULL;
    if (hdrl_image_get_size_x(sum1) != hdrl_image_get_size_x(sum2)) return NULL;
    if (hdrl_image_get_size_y(sum1) != hdrl_image_get_size_y(sum2)) return NULL;

    // Step 0: Define variables
    // vectors for each combination of fiber and angle
    cpl_bivector *trace_b_angle_1;
    cpl_bivector *trace_a_angle_1;
    cpl_bivector *trace_b_angle_2;
    cpl_bivector *trace_a_angle_2;

    // Get number of orders
    cpl_size order;
    int nb_orders;
    int * orders = cr2res_trace_get_order_numbers(trace_wave, &nb_orders);

    // temporary Stokes parameters
    cpl_vector * sx, *si, *sn;
    int res = 0;

    // output stokes table
    // number of elements = number of orders * number of chips * 
    //      number of stokes parameters
    cpl_size size = hdrl_image_get_size_x(sum1);
    cpl_table * stokes = cpl_table_new(size);
    char column_name[10]; //for constructing the names of columns

    // step 2: extraction
    for(int i = 0; i < nb_orders; i++) {
        order = orders[i];

        // Step 2: Extract individual traces
        res = cr2res_pol_demod_extract(sum1, sum2, trace_wave, order, 
                &trace_a_angle_1,
            &trace_b_angle_1, &trace_a_angle_2, &trace_b_angle_2);
        if (res == -1) {
            sx = cpl_vector_new(size);
            cpl_vector_fill(sx, 0.);
            si = cpl_vector_new(size);
            cpl_vector_fill(si, 0.);
            sn = cpl_vector_new(size);
            cpl_vector_fill(sn, 0.);
        } else {
            // Step 3: Calculate Stokes spectra
            cr2res_pol_demod_order(trace_a_angle_1, trace_b_angle_1, 
                    trace_a_angle_2, trace_b_angle_2, &sx, &si, &sn);
            cpl_bivector_delete(trace_a_angle_1);
            cpl_bivector_delete(trace_b_angle_1);
            cpl_bivector_delete(trace_a_angle_2);
            cpl_bivector_delete(trace_b_angle_2);
        }
        sprintf(column_name, "POL_X_%d", (int)order);
        cpl_table_wrap_double(stokes, cpl_vector_get_data(sx), column_name);
        cpl_vector_unwrap(sx);

        sprintf(column_name, "POL_I_%d", (int)order);
        cpl_table_wrap_double(stokes, cpl_vector_get_data(si), column_name);
        cpl_vector_unwrap(si);

        sprintf(column_name, "POL_N_%d", (int)order);
        cpl_table_wrap_double(stokes, cpl_vector_get_data(sn), column_name);
        cpl_vector_unwrap(sn);
    }

    // delete all temporary values
    cpl_free(orders);

    return stokes;
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
        cpl_table       *   tw_decker1, 
        cpl_table       *   tw_decker2, 
        cpl_polynomial  **  coef_slit, 
        cpl_polynomial  **  coef_wave)
{

    if (tw_decker1 == NULL || tw_decker2 == NULL || coef_slit == NULL || 
            coef_wave == NULL) return -1;

    // load data
    int i, j, k;
    int order;
    double _x, _y, _w;
    cpl_array *tmp;
    // orders and traces in decker 1
    int nb_orders_1;
    int *unique_orders_1 = cr2res_trace_get_order_numbers(tw_decker1, 
            &nb_orders_1);

    // orders and traces in decker 2
    int nb_orders_2;
    int *unique_orders_2 = cr2res_trace_get_order_numbers(tw_decker1, 
            &nb_orders_2);

    if (nb_orders_1 != nb_orders_2){
        cpl_free(unique_orders_1);
        cpl_free(unique_orders_2);
        return -1;
    }

    volatile cpl_error_code errcode;

    // fixed slit positions of each of the four trace per order
    double slit[] = {8.75, 6.25, 3.75, 1.25};

    // pixels x, only one because thats the same for all traces
    cpl_vector *x = cpl_vector_new(CR2RES_DETECTOR_SIZE);
    for (i = 0; i < CR2RES_DETECTOR_SIZE; i++) 
        cpl_vector_set(x, i, (double)i+1);

    for (i=0; i < nb_orders_1; i++) {
        coef_wave[i] = cpl_polynomial_new(2);
        coef_slit[i] = cpl_polynomial_new(2);
    }

    cpl_matrix *matrix_xy = cpl_matrix_new(2, CR2RES_DETECTOR_SIZE * 4);
    cpl_matrix *matrix_wd = cpl_matrix_new(2, CR2RES_DETECTOR_SIZE * 4);

    cpl_vector *vec_w = cpl_vector_new(CR2RES_DETECTOR_SIZE*4);
    cpl_vector *vec_s = cpl_vector_new(CR2RES_DETECTOR_SIZE*4);

    const cpl_size maxdeg = 2;

    cpl_polynomial *wave;
    cpl_polynomial *line;

    for (j = 0; j < nb_orders_1; j++) {
        // iterate over orders of trace wave 1, has to be the same number 
        // of orders as trace wave 2
        order = unique_orders_1[j];
        // TODO better criterion
        if ((order == 0) || (order == 8)) continue;

        for (k = 1; k <= 4; k++) {
            // load wavelenght polynomials, one per trace, note the
            // numbering depends on decker
            if (k == 1) {
                line = cr2res_get_trace_wave_poly(tw_decker1, 
                        CR2RES_COL_ALL, order, 1);
                wave = cr2res_get_trace_wave_poly(tw_decker1, 
                        CR2RES_COL_WAVELENGTH, order, 1);
            } else if (k == 2) {
                line = cr2res_get_trace_wave_poly(tw_decker2, 
                        CR2RES_COL_ALL, order, 1);
                wave = cr2res_get_trace_wave_poly(tw_decker2, 
                        CR2RES_COL_WAVELENGTH, order, 1);
            } else if (k == 3) {
                line = cr2res_get_trace_wave_poly(tw_decker1, 
                        CR2RES_COL_ALL, order, 2);
                wave = cr2res_get_trace_wave_poly(tw_decker1, 
                        CR2RES_COL_WAVELENGTH, order, 2);
            } else if (k == 4) {
                line = cr2res_get_trace_wave_poly(tw_decker2, 
                        CR2RES_COL_ALL, order, 2);
                wave = cr2res_get_trace_wave_poly(tw_decker2, 
                        CR2RES_COL_WAVELENGTH, order, 2);
            }
            // calculate polynomials for all traces
            for (i = 0; i < CR2RES_DETECTOR_SIZE; i++) {
                // Center line
                _x = cpl_vector_get(x, i);
                _y = cpl_polynomial_eval_1d(line, _x, NULL);
                _w = cpl_polynomial_eval_1d(wave, _x, NULL);

                cpl_matrix_set(matrix_xy, 0, CR2RES_DETECTOR_SIZE*(k-1) + i,_x);
                cpl_matrix_set(matrix_xy, 1, CR2RES_DETECTOR_SIZE*(k-1) + i,_y);

                cpl_matrix_set(matrix_wd, 0, CR2RES_DETECTOR_SIZE*(k-1) + i,_w);
                cpl_matrix_set(matrix_wd, 1, CR2RES_DETECTOR_SIZE*(k-1) + i,_y);

                cpl_vector_set(vec_w, CR2RES_DETECTOR_SIZE*(k-1) + i, _w);
                cpl_vector_set(vec_s, CR2RES_DETECTOR_SIZE*(k-1) + i,slit[k-1]);
            }
            cpl_polynomial_delete(line);
            cpl_polynomial_delete(wave);
        }

        // fit 2D wavelengths
        errcode = cpl_polynomial_fit(coef_wave[j], matrix_xy, NULL, vec_w, 
                NULL, FALSE, NULL, &maxdeg);
        errcode = cpl_polynomial_fit(coef_slit[j], matrix_wd, NULL, vec_s, 
                NULL, FALSE, NULL, &maxdeg);
    }

    // delete cpl pointers
    cpl_vector_delete(x);
    cpl_matrix_delete(matrix_xy);
    cpl_matrix_delete(matrix_wd);
    cpl_vector_delete(vec_s);
    cpl_vector_delete(vec_w);

    cpl_free(unique_orders_1);
    cpl_free(unique_orders_2);
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
        cpl_table   *   tw_decker1, 
        cpl_table   *   tw_decker2, 
        cpl_image   *   slitpos, 
        cpl_image   *   wavelength)
{
    if (tw_decker1 == NULL | tw_decker2 == NULL | slitpos == NULL | 
            wavelength == NULL) return -1;
    double w, s;
    int nb_orders;
    int *orders = cr2res_trace_get_order_numbers(tw_decker1, &nb_orders);
    cpl_free(orders);

    cpl_polynomial *coef_slit[nb_orders];
    cpl_polynomial *coef_wave[nb_orders];

    if (-1 == cr2res_slit_pos(tw_decker1, tw_decker2, coef_slit, coef_wave)){
        return -1;
    }

    //slitpos = cpl_image_new(CR2RES_DETECTOR_SIZE, CR2RES_DETECTOR_SIZE, 
    //              CPL_TYPE_DOUBLE);
    //wavelength = cpl_image_new(CR2RES_DETECTOR_SIZE, CR2RES_DETECTOR_SIZE, 
    //              CPL_TYPE_DOUBLE);

    cpl_vector *vec_xy = cpl_vector_new(2);
    cpl_vector *vec_wd = cpl_vector_new(2);

    for (int k = 0; k < nb_orders; k++)
    {
        for (int x=1; x <= CR2RES_DETECTOR_SIZE; x++)
        {
            for (int y=1; y <= CR2RES_DETECTOR_SIZE; y++)
            {
                cpl_vector_set(vec_xy, 0, (double)x);
                cpl_vector_set(vec_xy, 1, (double)y);
                w = cpl_polynomial_eval(coef_wave[k], vec_xy);

                cpl_vector_set(vec_wd, 0, w);
                cpl_vector_set(vec_wd, 1, (double)y);
                s = cpl_polynomial_eval(coef_slit[k], vec_wd);

                if ((s > 0) && (s < 10))
                {
                    cpl_image_set(slitpos, x, y, s);
                    cpl_image_set(wavelength, x, y, w);
                }
            }
        }
    }

    for (int i=0; i < nb_orders; i++){
        cpl_polynomial_delete(coef_wave[i]);
        cpl_polynomial_delete(coef_slit[i]);
    }
    cpl_vector_delete(vec_wd);
    cpl_vector_delete(vec_xy);
    return 0;
}

/* End TODO 3  */

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
