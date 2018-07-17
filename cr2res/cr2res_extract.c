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
#include <cpl.h>
#include "cr2res_dfs.h"
#include "cr2res_trace.h"
#include "cr2res_extract.h"

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

typedef unsigned char byte;
#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))
#define signum(a) (((a)>0)?1:((a)<0)?-1:0)

typedef struct
        {
          int x; int y; /* Coordinates of target pixel x,y  */
          double w;     /* Contribution weight <= 1/osample */
        } xi_ref;

typedef struct
        {
          int x; int iy;/* Contributing subpixel  x,iy      */
          double w;     /* Contribution weight <= 1/osample */
        } zeta_ref;

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_extract_slit_func_vert(
         int         ncols,
         int         nrows,
         int         osample,
         double  *   im,
         int     *   mask,
         double  *   ycen,
         double  *   sL,
         double  *   sP,
         double  *   model,
         double      lambda_sP,
         double      lambda_sL,
         double      sP_stop,
         int         maxiter) ;
static int cr2res_extract_slit_func_curved(
         int         ncols,
         int         nrows,
         int         osample,
         double  *   im,
         double  *   pix_unc,
         int     *   mask,
         double  *   ycen,
         int     *   ycen_offset,
         int         y_lower_lim,
         xi_ref  *   xi,
         zeta_ref *  zeta,
         int     *   m_zeta,
         double  *   sL,
         double  *   sP,
         double  *   model,
         double  *   unc,
         double      lambda_sP,
         double      lambda_sL,
         double      sP_stop,
         int         maxiter);
static int cr2res_extract_xi_zeta_tensors(
                    int ncols,
                    int nrows,
                    double * ycen,
                    int * ycen_offset,
                    int y_lower_lim,


                    int osample,
                    double * PSF_curve);
static int cr2res_extract_slitdec_bandsol(double *, double *, int, int) ;
static int cr2res_extract_slitdec_adjust_swath(int sw, int nx, cpl_vector *bins_begin, cpl_vector *bins_end);

/*----------------------------------------------------------------------------*/
/**
  @defgroup cr2res_slitdec  Extraction routines (Slit Decomposition,...)
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Simple extraction function
  @param    img_in      full detector image
  @param    trace_tab   The traces table
  @param    order       The order to extract
  @param    trace_id    The Trace to extract
  @param    height      number of pix above and below mid-line or -1
  @param    slit_func   the returned slit function, normalized to sum=1
  @param    spec        the returned spectrum, sum of rows
  @param    model       the reconstructed image
  @return   0 if ok, -1 otherwise

  This func takes a single image (containing many orders), a trace table,
  the order and trace numbers to extract (one each), and the extraction
  height.
  It uses the trace to shift the image columns by integer values,
  thereby straightening the order. The output spectra and slit
  function are then the collapsed image in x and y, respectively.
  Also returned is the "model", i.e. an image reconstruction from
  the two extracted vectors.

 */
/*----------------------------------------------------------------------------*/
int cr2res_extract_sum_vert(
        cpl_image   *   img_in,
        cpl_table   *   trace_tab,
        int             order,
        int             trace_id,
        int             height,
        cpl_vector  **  slit_func,
        cpl_vector  **  spec,
        hdrl_image  **  model)
{
    int             *   ycen_int;
    cpl_vector      *   ycen ;
    cpl_image       *   img_tmp;
    cpl_image       *   img_1d;
    cpl_vector      *   spc;
    cpl_vector      *   slitfu;
    cpl_size            lenx, leny;
    int                 i, j;
    int                 ymin, ymax;
    int                 empty_bottom = 0;
    cpl_type            imtyp;

    /* Check Entries */
    if (img_in == NULL || trace_tab == NULL) return -1 ;

    /* use the same type as input for temp images below */
    imtyp = cpl_image_get_type(img_in);
    lenx = cpl_image_get_size_x(img_in);
    leny = cpl_image_get_size_y(img_in);

    /* Compute height if not given */
    if (height <= 0) {
        height = cr2res_trace_get_height(trace_tab, order, trace_id);
        if (height <= 0) {
            cpl_msg_error(__func__, "Cannot compute height");
            return -1;
        }
    }
    /* Get ycen */
    if ((ycen = cr2res_trace_get_ycen(trace_tab, order,
                    trace_id, lenx)) == NULL) {
        cpl_msg_error(__func__, "Cannot get ycen");
        return -1 ;
    }

    img_tmp = cr2res_image_cut_rectify(img_in, ycen, height);
    if (img_tmp == NULL) {
        cpl_msg_error(__func__, "Cannot rectify order");
        cpl_vector_delete(ycen);
        return -1;
    }

    if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
        cpl_image_save(img_tmp, "debug_rectorder.fits", imtyp,
                NULL, CPL_IO_CREATE);
    }
    img_1d = cpl_image_collapse_create(img_tmp, 0); // sum of rows
    spc = cpl_vector_new_from_image_row(img_1d, 1);
    cpl_image_delete(img_1d);

    img_1d = cpl_image_collapse_create(img_tmp, 1); // sum of columns
    slitfu = cpl_vector_new_from_image_column(img_1d, 1);
    cpl_vector_divide_scalar(slitfu, cpl_vector_get_sum(slitfu));
    cpl_image_delete(img_1d);
    cpl_image_delete(img_tmp);

    // reconstruct the 2d image with the "model"
    img_tmp = cpl_image_new(lenx, leny, imtyp);
    ycen_int = cr2res_vector_get_int(ycen);
    for (i=1;i<=lenx;i++){
        for (j=1;j<=height;j++){
            cpl_image_set(img_tmp, i, ycen_int[i-1]-(height/2)+j,
                cpl_vector_get(spc,i-1)*cpl_vector_get(slitfu,j-1) );
        }
    }
    cpl_vector_delete(ycen);
    cpl_free(ycen_int);

    if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
        cpl_image_save(img_tmp, "debug_model.fits", imtyp,
                NULL, CPL_IO_CREATE);
    }


    *slit_func = slitfu;
    *spec = spc;
    *model = hdrl_image_create(img_tmp, NULL);
    cpl_image_delete(img_tmp);

    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Create the extract 1D table to be saved
  @param    spectrum    The extracted spectra of the different orders
  @param    trace_table The trace table
  @return   the extract_1D table or NULL
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_extract_EXTRACT1D_create(
        cpl_vector      **  spectrum,
        cpl_table       *   trace_table)
{
    cpl_table       *   out ;
    char            *   col_name ;
    const double    *   pspec ;
    int                 nrows, all_null, i, order, trace_id, nb_traces ;

    /* Check entries */
    if (spectrum == NULL || trace_table == NULL) return NULL ;

    /* Initialise */
    nb_traces = cpl_table_get_nrow(trace_table) ;

    /* Check if all vectorѕ are not null */
    all_null = 1 ;
    for (i=0 ; i<nb_traces ; i++)
        if (spectrum[i] != NULL) {
            nrows = cpl_vector_get_size(spectrum[i]) ;
            all_null = 0 ;
        }
    if (all_null == 1) return NULL ;

    /* Check the sizes */
    for (i=0 ; i<nb_traces ; i++)
        if (spectrum[i] != NULL && cpl_vector_get_size(spectrum[i]) != nrows)
            return NULL ;

    /* Create the table */
    out = cpl_table_new(nrows);
    for (i=0 ; i<nb_traces ; i++) {
        order = cpl_table_get(trace_table, CR2RES_COL_ORDER, i, NULL) ;
        trace_id = cpl_table_get(trace_table, CR2RES_COL_TRACENB, i, NULL) ;
        col_name = cr2res_dfs_SPEC_colname(order, trace_id) ;
        cpl_table_new_column(out, col_name, CPL_TYPE_DOUBLE);
        cpl_free(col_name) ;
    }

    /* Fill the table */
    for (i=0 ; i<nb_traces ; i++) {
        if (spectrum[i] != NULL) {
            order = cpl_table_get(trace_table, CR2RES_COL_ORDER, i, NULL) ;
            trace_id = cpl_table_get(trace_table, CR2RES_COL_TRACENB, i, NULL) ;
            pspec = cpl_vector_get_data_const(spectrum[i]) ;
            col_name = cr2res_dfs_SPEC_colname(order, trace_id) ;
            cpl_table_copy_data_double(out, col_name, pspec) ;
            cpl_free(col_name) ;
        }
    }
    return out ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Create the slit functions table to be saved
  @param    slit_func   The slit functions of the different orders
  @param    trace_table The trace table
  @return   the slit_func table or NULL
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_extract_SLITFUNC_create(
        cpl_vector      **  slit_func,
        cpl_table       *   trace_table)
{
    cpl_table       *   out ;
    char            *   col_name ;
    const double    *   pslit ;
    int                 nrows, all_null, i, order, trace_id, nb_traces ;

    /* Check entries */
    if (slit_func == NULL || trace_table == NULL) return NULL ;

    /* Initialise */
    nb_traces = cpl_table_get_nrow(trace_table) ;

    /* Check the all vectorѕ are not null */
    all_null = 1 ;
    for (i=0 ; i<nb_traces ; i++)
        if (slit_func[i] != NULL) {
            nrows = cpl_vector_get_size(slit_func[i]) ;
            all_null = 0 ;
        }
    if (all_null == 1) return NULL ;

    /* Check the sizes */
    for (i=0 ; i<nb_traces ; i++)
        if (slit_func[i] != NULL && cpl_vector_get_size(slit_func[i]) != nrows)
            return NULL ;

    /* Create the table */
    out = cpl_table_new(nrows);
    for (i=0 ; i<nb_traces ; i++) {
        order = cpl_table_get(trace_table, CR2RES_COL_ORDER, i, NULL) ;
        trace_id = cpl_table_get(trace_table, CR2RES_COL_TRACENB, i, NULL) ;
        col_name = cr2res_dfs_SLIT_FUNC_colname(order, trace_id) ;
        cpl_table_new_column(out, col_name, CPL_TYPE_DOUBLE);
        cpl_free(col_name) ;
    }

    /* Fill the table */
    for (i=0 ; i<nb_traces ; i++) {
        if (slit_func[i] != NULL) {
            order = cpl_table_get(trace_table, CR2RES_COL_ORDER, i, NULL) ;
            trace_id = cpl_table_get(trace_table, CR2RES_COL_TRACENB, i, NULL) ;
            pslit = cpl_vector_get_data_const(slit_func[i]) ;
            col_name = cr2res_dfs_SLIT_FUNC_colname(order, trace_id) ;
            cpl_table_copy_data_double(out, col_name, pslit) ;
            cpl_free(col_name) ;
        }
    }
    return out ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get a Spectrum from the EXTRACT_1D table
  @param    tab         the EXTRACT_1D table
  @param    order       the order
  @param    trace_nb    the wished trace
  @return   the spectrum or NULL
 */
/*----------------------------------------------------------------------------*/
cpl_vector * cr2res_extract_EXTRACT1D_get_spectrum(
        cpl_table   *   tab,
        int             order,
        int             trace_nb)
{
    cpl_vector  *   out ;
    char        *   col_name ;
    double      *   pcol ;
    double      *   pout ;
    int             i, tab_size ;

    /* Check entries */
    if (tab == NULL) return NULL ;

    /* Col name */
    col_name = cr2res_dfs_SPEC_colname(order, trace_nb) ;

    /* Get the column */
    if ((pcol = cpl_table_get_data_double(tab, col_name)) == NULL) {
        cpl_msg_error(__func__, "Cannot find the extracted spectrum") ;
        cpl_free(col_name) ;
        return NULL ;
    }
    cpl_free(col_name) ;

    /* Create the output vector */
    tab_size = cpl_table_get_nrow(tab) ;
    out = cpl_vector_new(tab_size) ;
    pout = cpl_vector_get_data(out) ;
    for (i=-0 ; i<tab_size ; i++)
        pout[i] = pcol[i] ;

    return out ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief   Adjust the swath width to match the length of detector
  @param    sw Swath width to start from
  @param    nx number of pixel columns to match

  Returns the new swath width. CURRENTLY UNIMPLEMENTED, except it ensures
  an even number. TODO

 */
/*----------------------------------------------------------------------------*/
static int cr2res_extract_slitdec_adjust_swath(int sw, int nx, cpl_vector* bins_begin, cpl_vector *bins_end)
{
    if (sw == 0  || nx == 0) return -1;
    if (bins_begin == NULL || bins_end == NULL) return -1;

    int nbin, i = 0;
    double step = 0;

    // Calculate number of bins
    nbin =  (int) round((double)nx / (double)sw);
    if (nbin < 1) nbin = 1;
    
    step = (double)nx / (double)nbin / 2.;
    double bins[nbin * 2 + 1];
    cpl_vector_set_size(bins_begin, 2*nbin-1);
    cpl_vector_set_size(bins_end, 2*nbin-1);

    // boundaries of bins
    for(i = 0; i < nbin*2 + 1; i++)
    {
        bins[i] = i * step;
        if (i < nbin*2 + 1 - 2){
            cpl_vector_set(bins_begin, i, floor(bins[i]));
        }
        if (i >= 2){
            cpl_vector_set(bins_end, i-2, ceil(bins[i]));
        }
    }

    sw = (int) ceil(step * 2. + 1);
    return sw;
}

/*----------------------------------------------------------------------------*/
/**
  @brief   Main vertical slit decomposition wrapper with swath loop
  @param    img_in      full detector image
  @param    trace_tab   The traces table
  @param    order       The order to extract
  @param    trace_id    The Trace to extract
  @param    height      number of pix above and below mid-line or -1
  @param    swath       width per swath
  @param    oversample  factor for oversampling
  @param    smooth_slit
  @param    slit_func   the returned slit function
  @param    spec        the returned spectrum
  @param    model       the returned model
  @return   0 if ok, -1 otherwise

  This func takes a single image (contining many orders), and a *single*
  order definition in the form of central y-corrds., plus the height.
  Swath widht and oversampling are passed through.

  The task of this function then is to

    -cut out the relevant pixels of the order
    -shift im in y integers, so that nrows becomes minimal,
        adapt ycen accordingly
    -loop over swaths, in half-steps
    -derive a good starting guess for the spectrum, by median-filter
        average along slit, beware of cosmics
    -run cr2res_extract_slit_func_vert()
    -merge overlapping swath results by linear weights from swath-width to edge.
    -return re-assembled model image, slit-fu, spectrum, new mask.
    -calculate the errors and return them. This is done by comparing the
        variance of (im-model) to the poisson-statistics of the spectrum.

 */
/*----------------------------------------------------------------------------*/
int cr2res_extract_slitdec_vert(
        cpl_image   *   img_in,
        cpl_table   *   trace_tab,
        int             order,
        int             trace_id,
        int             height,
        int             swath,
        int             oversample,
        double          smooth_slit,
        cpl_vector  **  slit_func,
        cpl_vector  **  spec,
        hdrl_image  **  model)
{
    cpl_polynomial  **  traces ;
    int             *   ycen_int;
    double          *   ycen_rest;
    double          *   ycen_sw;
    double          *   img_sw_data;
    double          *   spec_sw_data;
    double          *   slitfu_sw_data;
    double          *   model_sw;
    int             *   mask_sw;
    cpl_image       *   img_sw;
    cpl_image       *   img_rect;
    cpl_image       *   model_rect;
    cpl_vector      *   ycen ;
    cpl_image       *   img_tmp;
    cpl_image       *   img_out;
    cpl_vector      *   spec_sw;
    cpl_vector      *   slitfu_sw;
    cpl_vector      *   spc;
    cpl_vector      *   slitfu;
    cpl_vector      *   weights_sw;
    cpl_vector      *   tmp_vec;
    cpl_vector      *   bins_begin;
    cpl_vector      *   bins_end;
    cpl_size            lenx, leny;
    cpl_type            imtyp;
    double              pixval, img_median;
    int                 i, j, nswaths, row, col, x, y, ny_os,
                        sw_start, sw_end, badpix;



    /* Check Entries */
    if (img_in == NULL || trace_tab == NULL) return -1 ;

    /* Compute height if not given */
    if (height <= 0) {
        height = cr2res_trace_get_height(trace_tab, order, trace_id);
        if (height <= 0) {
            cpl_msg_error(__func__, "Cannot compute height");
            return -1;
        }
    }
    /* Initialise */
    imtyp = cpl_image_get_type(img_in);
    lenx = cpl_image_get_size_x(img_in);
    leny = cpl_image_get_size_y(img_in);

    /* Get ycen */
    if ((ycen = cr2res_trace_get_ycen(trace_tab, order,
                    trace_id, lenx)) == NULL) {
        cpl_msg_error(__func__, "Cannot get ycen");
        return -1 ;
    }

    bins_begin = cpl_vector_new(10);
    bins_end =   cpl_vector_new(lenx / swath);

    /* Number of rows after oversampling */
    ny_os = oversample*(height+1) +1;
    if ((swath = cr2res_extract_slitdec_adjust_swath(swath, lenx, bins_begin, bins_end)) == -1){
        cpl_msg_error(__func__, "Cannot calculate swath size");
        cpl_vector_delete(ycen);
        cpl_vector_delete(bins_begin);
        cpl_vector_delete(bins_end);
        return -1;
    }
    nswaths = cpl_vector_get_size(bins_begin);

    // Get cut-out rectified order
    img_rect = cr2res_image_cut_rectify(img_in, ycen, height);
    if (img_rect == NULL){
        cpl_msg_error(__func__, "Cannot rectify order");
        cpl_vector_delete(ycen);
        cpl_vector_delete(bins_begin);
        cpl_vector_delete(bins_end);
        return -1;
    }
    if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
        cpl_image_save(img_rect, "debug_rectorder.fits", imtyp,
                NULL, CPL_IO_CREATE);
    }
    ycen_rest = cr2res_vector_get_rest(ycen);

    // Work vectors
    slitfu_sw = cpl_vector_new(ny_os);
    slitfu_sw_data = cpl_vector_get_data(slitfu_sw);



    // Local versions of return data
    slitfu = cpl_vector_new(ny_os);
    spc = cpl_vector_new(lenx);
    img_out = cpl_image_new(lenx, leny, CPL_TYPE_DOUBLE);
    model_rect = cpl_image_new(lenx, height, CPL_TYPE_DOUBLE);


    for (i=0;i<nswaths-1;i++){ // TODO: Treat last swath!
        sw_start = cpl_vector_get(bins_begin, i);
        sw_end = cpl_vector_get(bins_end, i);
        swath = sw_end - sw_start;

        /* Allocate */
        mask_sw = cpl_malloc(height*swath*sizeof(int));
        model_sw = cpl_malloc(height*swath*sizeof(double));
        ycen_sw = cpl_malloc(swath*sizeof(double));
        img_sw = cpl_image_new(swath, height, CPL_TYPE_DOUBLE);

        weights_sw = cpl_vector_new(swath);

        /* Pre-calculate the weights for overlapping swaths*/
        for (i=0; i < swath/2; i++) {
            cpl_vector_set(weights_sw, i, i + 1);
            cpl_vector_set(weights_sw, swath/2 - i - 1, i+1);
        }
        cpl_vector_divide_scalar(weights_sw, i + 1); // normalize such that max(w)=1


        // Copy swath image into seperate image
        for(col=1; col<=swath; col++){      // col is x-index in swath
            x = sw_start + col;          // coords in large image
            for(y=1;y<=height;y++){
                pixval = cpl_image_get(img_rect, x, y, &badpix);
                cpl_image_set(img_sw, col, y, pixval);
                if(cpl_error_get_code() != CPL_ERROR_NONE)
                    cpl_msg_error(__func__, "%d %d %s", x, y, cpl_error_get_where());
                j = (y-1)*swath + (col-1) ; // raw index for mask, start with 0!
                if (badpix == 0) mask_sw[j] = 1;
                else mask_sw[j] = 0;
            }
        }

        // img_median = cpl_image_get_median(img_sw);
        // for (j=0;j<ny_os;j++) cpl_vector_set(slitfu_sw,j,img_median);
        // cpl_image_turn(img_sw, 1);
        img_sw_data = cpl_image_get_data_double(img_sw);
        img_tmp = cpl_image_collapse_median_create(img_sw, 0, 0, 0);
        spec_sw = cpl_vector_new_from_image_row(img_tmp,1);
        cpl_image_delete(img_tmp);
        spec_sw_data = cpl_vector_get_data(spec_sw);
        for (j=sw_start;j<sw_end;j++) ycen_sw[j-sw_start] = ycen_rest[j];


        /* Finally ready to call the slit-decomp */
        cr2res_extract_slit_func_vert(swath, height, oversample, img_sw_data,
                mask_sw, ycen_sw, slitfu_sw_data, spec_sw_data, model_sw,
                0.0, smooth_slit, 1.0e-5, 20);

        for(col=1; col<=swath; col++){   // col is x-index in cut-out
            x = sw_start + col;          // coords in large image
            for(y=1;y<=height;y++){
                j = (y-1)*swath + (col-1) ; // raw index for mask, start with 0!
                cpl_image_set(model_rect,x,y, model_sw[j]);
                if(cpl_error_get_code() != CPL_ERROR_NONE)
                    cpl_msg_error(__func__, "%d %d %s", x, y, cpl_error_get_where());
            }
        }

        // add up slit-functions, divide by nswaths below to get average
        if (i==0) cpl_vector_copy(slitfu,slitfu_sw);
        else cpl_vector_add(slitfu, slitfu_sw);

        if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
            cpl_vector_save(spec_sw, "debug_spc.fits", CPL_TYPE_DOUBLE, NULL,
                    CPL_IO_CREATE);
            tmp_vec = cpl_vector_wrap(swath, ycen_sw);
            cpl_vector_save(tmp_vec, "debug_ycen.fits", CPL_TYPE_DOUBLE, NULL,
                    CPL_IO_CREATE);
            cpl_vector_unwrap(tmp_vec);
            cpl_vector_save(slitfu_sw, "debug_slitfu.fits", CPL_TYPE_DOUBLE,
                    NULL, CPL_IO_CREATE);
            img_tmp = cpl_image_wrap_double(swath, height, model_sw);
            cpl_image_save(img_tmp, "debug_model_sw.fits", CPL_TYPE_FLOAT,
                    NULL, CPL_IO_CREATE);
            cpl_image_unwrap(img_tmp);
            cpl_image_save(img_sw, "debug_img_sw.fits", CPL_TYPE_FLOAT, NULL,
                    CPL_IO_CREATE);
        }
        /* Multiply by weights and add to output array */
        cpl_vector_multiply(spec_sw, weights_sw);
        if (i==0){
            for (j=0;j<swath/2;j++) {
                cpl_vector_set(spec_sw, j,
                    cpl_vector_get(spec_sw,j)/cpl_vector_get(weights_sw,j));
            }
        }
        if (i==nswaths-1) {
            for (j=swath/2; j<swath;j++) {
                cpl_vector_set(spec_sw, j,
                    cpl_vector_get(spec_sw,j)/cpl_vector_get(weights_sw,j));
            }
        }
        for (j=sw_start;j<sw_end;j++) {
            cpl_vector_set(spc, j,
                cpl_vector_get(spec_sw, j - sw_start) + cpl_vector_get(spc, j) );
        }        cpl_vector_delete(spec_sw);

        cpl_free(mask_sw);
        cpl_free(model_sw);
        cpl_free(ycen_sw);
        cpl_image_delete(img_sw);

    } // End loop over swaths
    cpl_vector_delete(slitfu_sw);
    cpl_vector_delete(weights_sw);
    cpl_vector_delete(bins_begin);
    cpl_vector_delete(bins_end);

    // insert model_rect into large frame
    cr2res_image_insert_rect(model_rect, ycen, img_out);

    // divide by nswaths to make the slitfu into the average over all swaths.
    cpl_vector_divide_scalar(slitfu, nswaths);

    // TODO: Update BPM in img_out
    // TODO: Calculate error and return it.

    // TODO: Deallocate return arrays in case of error, return -1
    cpl_image_delete(img_rect);
    cpl_image_delete(model_rect);
    cpl_vector_delete(ycen);
    cpl_free(ycen_rest);

    *slit_func = slitfu;
    *spec = spc;
    *model = hdrl_image_create(img_out, NULL);
    cpl_image_delete(img_out);

    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Extract optimally (slit-decomposition) with polynomial slit
  @param    img_in      full detector image
  @param    trace_tab   The traces table
  @param    order       The order to extract
  @param    trace_id    The Trace to extract
  @param    height      number of pix above and below mid-line or -1
  @param    swath       width per swath
  @param    oversample  factor for oversampling
  @param    smooth_slit
  @param    slit_func   the returned slit function
  @param    spec        the returned spectrum
  @param    model       the returned model
  @return   0 if ok, -1 otherwise

  This func takes a single image (contining many orders), and a *single*
  order definition in the form of central y-corrds., plus the height.
  Swath widht and oversampling are passed through.

  The task of this function then is to

    -cut out the relevant pixels of the order
    -shift im in y integers, so that nrows becomes minimal,
        adapt ycen accordingly
    -loop over swaths, in half-steps
    -derive a good starting guess for the spectrum, by median-filter
        average along slit, beware of cosmics
    -run cr2res_extract_slit_func_vert()
    -merge overlapping swath results by linear weights from swath-width to edge.
    -return re-assembled model image, slit-fu, spectrum, new mask.
    -calculate the errors and return them. This is done by comparing the
        variance of (im-model) to the poisson-statistics of the spectrum.

 */
/*----------------------------------------------------------------------------*/
int cr2res_extract_slitdec_curved(
        hdrl_image  *   img_hdrl,
        cpl_table   *   trace_tab,
        int             order,
        int             trace_id,
        int             height,
        int             swath,
        int             oversample,
        double          smooth_slit,
        cpl_vector  **  slit_func,
        cpl_vector  **  spec,
        hdrl_image  **  model)
{
    cpl_polynomial  **  traces ;
    int             *   ycen_int;
    double          *   ycen_rest;
    double          *   ycen_sw;
    double          *   img_sw_data;
    double          *   spec_sw_data;
    double          *   slitfu_sw_data;
    double          *   model_sw;
    int             *   mask_sw;
    const cpl_image       *   img_in;
    const cpl_image       *   err_in;
    cpl_image       *   img_sw;
    cpl_image       *   err_sw;
    cpl_image       *   img_rect;
    cpl_image       *   model_rect;
    cpl_vector      *   ycen ;
    cpl_image       *   img_tmp;
    cpl_image       *   img_out;
    cpl_vector      *   spec_sw;
    cpl_vector      *   slitfu_sw;
    cpl_vector      *   spc;
    cpl_vector      *   slitfu;
    cpl_vector      *   weights_sw;
    cpl_vector      *   tmp_vec;
    cpl_vector      *   bins_begin;
    cpl_vector      *   bins_end;
    cpl_size            lenx, leny;
    cpl_type            imtyp;
    double              pixval, img_median;
    int                 i, j, nswaths, halfswath, row, col, x, y, ny_os,
                        sw_start, sw_end, badpix;

    img_in = hdrl_image_get_image_const(img_hdrl);
    err_in = hdrl_image_get_error_const(img_hdrl);

    /* Check Entries */
    if (img_in == NULL || trace_tab == NULL) return -1 ;

    /* Compute height if not given */
    if (height <= 0) {
        height = cr2res_trace_get_height(trace_tab, order, trace_id);
        if (height <= 0) {
            cpl_msg_error(__func__, "Cannot compute height");
            return -1;
        }
    }
    /* Initialise */
    imtyp = cpl_image_get_type(img_in);
    lenx = cpl_image_get_size_x(img_in);
    leny = cpl_image_get_size_y(img_in);
    /* Get ycen */
    if ((ycen = cr2res_trace_get_ycen(trace_tab, order,
                    trace_id, lenx)) == NULL) {
        cpl_msg_error(__func__, "Cannot get ycen");
        return -1 ;
    }

    bins_begin = cpl_vector_new((int)lenx / swath - 2);
    bins_end =   cpl_vector_new((int)lenx / swath - 2);

    /* Number of rows after oversampling */
    ny_os = oversample*(height+1) +1;
    swath = cr2res_extract_slitdec_adjust_swath(swath, lenx, bins_begin, bins_end);
    halfswath = swath/2;
    nswaths = (lenx / swath) *2; // *2 because we step in half swaths!
    if (lenx%swath >= halfswath) nswaths +=1;

    // Get cut-out rectified order
    img_rect = cr2res_image_cut_rectify(img_in, ycen, height);
    if (img_rect == NULL){
        cpl_msg_error(__func__, "Cannot rectify order");
        cpl_vector_delete(ycen);
        return -1;
    }
    if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
        cpl_image_save(img_rect, "debug_rectorder.fits", imtyp,
                NULL, CPL_IO_CREATE);
    }
    ycen_rest = cr2res_vector_get_rest(ycen);

    /* Allocate */
    mask_sw = cpl_malloc(height*swath*sizeof(int));
    model_sw = cpl_malloc(height*swath*sizeof(double));
    img_sw = cpl_image_new(swath, height, CPL_TYPE_DOUBLE);
    err_sw = cpl_image_new(swath, height, CPL_TYPE_DOUBLE);
    ycen_sw = cpl_malloc(swath*sizeof(double));

    // Local versions of return data
    slitfu = cpl_vector_new(ny_os);
    spc = cpl_vector_new(lenx);
    img_out = cpl_image_new(lenx, leny, CPL_TYPE_DOUBLE);
    model_rect = cpl_image_new(lenx, height, CPL_TYPE_DOUBLE);

    // Work vectors
    slitfu_sw = cpl_vector_new(ny_os);
    slitfu_sw_data = cpl_vector_get_data(slitfu_sw);
    weights_sw = cpl_vector_new(swath);

    /* Pre-calculate the weights for overlapping swaths*/
    for (i=0;i<halfswath;i++) {
         cpl_vector_set(weights_sw,i,i+1);
         cpl_vector_set(weights_sw,swath-i-1,i+1);
    }
    cpl_vector_divide_scalar(weights_sw,i+1); // normalize such that max(w)=1


    for (i=0;i<nswaths-1;i++){ // TODO: Treat last swath!
        sw_start = i*halfswath;
        sw_end = sw_start + swath;
        for(col=1; col<=swath; col++){      // col is x-index in swath
            x = i*halfswath + col;          // coords in large image
            for(y=1;y<=height;y++){
                pixval = cpl_image_get(img_rect, x, y, &badpix);
                if(cpl_error_get_code() != CPL_ERROR_NONE)
                    cpl_msg_error(__func__, "%d %d %s", x, y, cpl_error_get_where());
                cpl_image_set(img_sw, col, y, pixval);
                j = (y-1)*swath + (col-1) ; // raw index for mask, start with 0!
                if (badpix ==0) mask_sw[j] = 1;
                else mask_sw[j] = 0;
            }
        }

        img_median = cpl_image_get_median(img_sw);
        for (j=0;j<ny_os;j++) cpl_vector_set(slitfu_sw,j,img_median);
        img_sw_data = cpl_image_get_data_double(img_sw);
        img_tmp = cpl_image_collapse_median_create(img_sw, 0, 0, 0);
        spec_sw = cpl_vector_new_from_image_row(img_tmp,1);
        cpl_image_delete(img_tmp);
        spec_sw_data = cpl_vector_get_data(spec_sw);
        for (j=sw_start;j<sw_end;j++) ycen_sw[j-sw_start] = ycen_rest[j];

        /* Finally ready to call the slit-decomp */
        //cr2res_extract_slit_func_curved(swath, height, oversample, img_sw_data,
        //        mask_sw, ycen_sw, slitfu_sw_data, spec_sw_data, model_sw,
        //        0.0, smooth_slit, 1.0e-5, 20);

        for(col=1; col<=swath; col++){      // col is x-index in cut-out
            x = i*halfswath + col;          // coords in large image
            for(y=1;y<=height;y++){
                j = (y-1)*swath + (col-1) ; // raw index for mask, start with 0!
                cpl_image_set(model_rect,x,y, model_sw[j]);
            }
        }

        // add up slit-functions, divide by nswaths below to get average
        if (i==0) cpl_vector_copy(slitfu,slitfu_sw);
        else cpl_vector_add(slitfu,slitfu_sw);

        if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
            cpl_vector_save(spec_sw, "debug_spc.fits", CPL_TYPE_DOUBLE, NULL,
                    CPL_IO_CREATE);
            tmp_vec = cpl_vector_wrap(swath, ycen_sw);
            cpl_vector_save(tmp_vec, "debug_ycen.fits", CPL_TYPE_DOUBLE, NULL,
                    CPL_IO_CREATE);
            cpl_vector_unwrap(tmp_vec);
            cpl_vector_save(slitfu_sw, "debug_slitfu.fits", CPL_TYPE_DOUBLE,
                    NULL, CPL_IO_CREATE);
            img_tmp = cpl_image_wrap_double(swath, height, model_sw);
            cpl_image_save(img_tmp, "debug_model_sw.fits", CPL_TYPE_FLOAT,
                    NULL, CPL_IO_CREATE);
            cpl_image_unwrap(img_tmp);
            cpl_image_save(img_sw, "debug_img_sw.fits", CPL_TYPE_FLOAT, NULL,
                    CPL_IO_CREATE);
        }
        /* Multiply by weights and add to output array */
        cpl_vector_multiply(spec_sw, weights_sw);
        if (i==0){
            for (j=0;j<halfswath;j++) {
                cpl_vector_set(spec_sw, j,
                    cpl_vector_get(spec_sw,j)/cpl_vector_get(weights_sw,j));
            }
        }
        if (i==nswaths-1) {
            for (j=halfswath;j<swath;j++) {
                cpl_vector_set(spec_sw, j,
                    cpl_vector_get(spec_sw,j)/cpl_vector_get(weights_sw,j));
            }
        }
        for (j=sw_start;j<sw_end;j++) {
            cpl_vector_set(spc, j,
                cpl_vector_get(spec_sw,j-sw_start) + cpl_vector_get(spc, j) );
        }        cpl_vector_delete(spec_sw);
    } // End loop over swaths
    cpl_vector_delete(slitfu_sw);
    cpl_vector_delete(weights_sw);

    // insert model_rect into large frame
    cr2res_image_insert_rect(model_rect, ycen, img_out);

    // divide by nswaths to make the slitfu into the average over all swaths.
    cpl_vector_divide_scalar(slitfu,nswaths);

    // TODO: Update BPM in img_out
    // TODO: Calculate error and return it.

    // TODO: Deallocate return arrays in case of error, return -1
    cpl_image_delete(img_rect);
    cpl_image_delete(model_rect);
    cpl_image_delete(img_sw);
    cpl_free(mask_sw) ;
    cpl_free(model_sw) ;
    cpl_vector_delete(ycen);
    cpl_free(ycen_rest);
    cpl_free(ycen_sw);

    *slit_func = slitfu;
    *spec = spc;
    *model = hdrl_image_create(img_out, NULL);
    cpl_image_delete(img_out);

    return 0;
}



/*----------------------------------------------------------------------------*/
/**
  @brief    Slit-decomposition of a single swath, assuming vertical slit
  @param    ncols       Swath width in pixels
  @param    nrows       Extraction slit height in pixels
  @param    osample     Subpixel ovsersampling factor
  @param    im          Image to be decomposed
  @param    mask        int mask of same dimension as image
  @param    ycen        Order centre line offset from pixel row boundary
  @param    sL          Slit function resulting from decomposition, start
                        guess is input, gets overwriteten with result
  @param    sP          Spectrum resulting from decomposition
  @param    model       the model reconstruction of im
  @param    lambda_sP   Smoothing parameter for the spectrum, could be zero
  @param    lambda_sL   Smoothing parameter for the slit function, usually >0
  @param    sP_stop     Fraction of spectyrum change, stop condition
  @param    maxiter     Max number of iterations
  @return
 */
/*----------------------------------------------------------------------------*/
static int cr2res_extract_slit_func_vert(
        int         ncols,
        int         nrows,
        int         osample,
        double  *   im,
        int     *   mask,
        double  *   ycen,
        double  *   sL,
        double  *   sP,
        double  *   model,
        double      lambda_sP,
        double      lambda_sL,
        double      sP_stop,
        int         maxiter)
{
    int x, y, iy, jy, iy1, iy2, ny, nd, i, j;
    double step, d1, d2, sum, norm, dev, lambda, diag_tot, sP_change, sP_max;
    int info, iter, isum;
    /* Initialise */
    nd=2*osample+1;
    ny=osample*(nrows+1)+1; /* The size of the sf array */
    step=1.e0/osample;
    double * E = cpl_malloc(ncols*sizeof(double)); // double E[ncols];
    double * sP_old = cpl_malloc(ncols*sizeof(double)); // double sP_old[ncols];
    double * Aij = cpl_malloc(ny*ny*sizeof(double)); // double Aij[ny*ny];
    double * bj = cpl_malloc(ny*sizeof(double)); // double bj[ny];
    double * Adiag = cpl_malloc(ncols*3*sizeof(double)); // double Adiag[ncols*3];
    double * omega = cpl_malloc(ny*nrows*ncols*sizeof(double)); // double omega[ny][nrows][ncols];
    // index as: [iy+(y*ny)+(x*ny*nrows)]

    /*
      Construct the omega tensor. Normally it has the dimensionality of
      ny*nrows*ncols.
      The tensor is mostly empty and can be easily compressed to ny*nx, but
      this will complicate matrix operations at later stages. I will keep
      it as it is for now.
      Note, that omega is used in in the equations for sL, sP and for the model
      but it does not involve the data, only the geometry. Thus it can be
      pre-computed once.
      */
    for(x=0; x<ncols; x++) {
        iy2=(1.e0-ycen[x])*osample;
        /*
           The initial offset should be reconsidered.
           It looks fine but needs theory.
         */
        iy1=iy2-osample;
        if(iy2==0) d1=step;
        else if(iy1==0) d1=0.e0;
        else d1=fmod(ycen[x], step);
        d2=step-d1;
        for(y=0; y<nrows; y++)
        {
            iy1+=osample;
            iy2+=osample;
            for(iy=0; iy<ny; iy++)
            {
                if(iy<iy1)                omega[iy+(y*ny)+(x*ny*nrows)] = 0.;
                else if(iy==iy1)          omega[iy+(y*ny)+(x*ny*nrows)] = d1;
                else if(iy>iy1 && iy<iy2) omega[iy+(y*ny)+(x*ny*nrows)] = step;
                else if(iy==iy2)          omega[iy+(y*ny)+(x*ny*nrows)] = d2;
                else                      omega[iy+(y*ny)+(x*ny*nrows)] = 0.;
            }
        }
    }

    /* Loop through sL , sP reconstruction until convergence is reached */
    iter=0;
    do
    {
        /* Compute slit function sL */

        /* Fill in SLE arrays */
        diag_tot=0.e0;
        for(iy=0; iy<ny; iy++)
        {
            bj[iy]=0.e0;
            for(jy=max(iy-osample,0); jy<=min(iy+osample,ny-1); jy++)
            {
            /* printf("iy=%d jy=%d %d\n", iy, jy, iy+ny*(jy-iy+osample)); */
                Aij[iy+ny*(jy-iy+osample)]=0.e0;
                for(x=0; x<ncols; x++)
                {
                    sum=0.e0;
                    for(y=0; y<nrows; y++)
                        sum+=omega[iy+(y*ny)+(x*ny*nrows)] * omega[jy+(y*ny)+(x*ny*nrows)]*mask[y*ncols+x];
                    Aij[iy+ny*(jy-iy+osample)]+=sum*sP[x]*sP[x];
                }
            }
            for(x=0; x<ncols; x++)
            {
                sum=0.e0;
                for(y=0; y<nrows; y++)
                    sum+=omega[iy+(y*ny)+(x*ny*nrows)]*
                        mask[y*ncols+x]*im[y*ncols+x];
                bj[iy]+=sum*sP[x];
            }
            diag_tot+=Aij[iy+ny*osample];
        }
        // printf("SUM : %e\n", sum);

        /* Scale regularization parameters */
        lambda=lambda_sL*diag_tot/ny;

        /* Add regularization parts for the slit function */
        Aij[ny*osample]    +=lambda;           /* Main diagonal  */
        Aij[ny*(osample+1)]-=lambda;           /* Upper diagonal */
        for(iy=1; iy<ny-1; iy++)
        {
            Aij[iy+ny*(osample-1)]-=lambda;      /* Lower diagonal */
            Aij[iy+ny*osample    ]+=lambda*2.e0; /* Main diagonal  */
            Aij[iy+ny*(osample+1)]-=lambda;      /* Upper diagonal */
        }
        Aij[ny-1+ny*(osample-1)]-=lambda;      /* Lower diagonal */
        Aij[ny-1+ny*osample]    +=lambda;      /* Main diagonal  */

        /* Solve the system of equations */
        info=cr2res_extract_slitdec_bandsol(Aij, bj, ny, nd);
        if(info) printf("info(sL)=%d\n", info);

        /* Normalize the slit function */
        norm=0.e0;
        for(iy=0; iy<ny; iy++)
        {
            sL[iy]=bj[iy];
            norm+=sL[iy];
        }
        norm/=osample;
        for(iy=0; iy<ny; iy++) sL[iy]/=norm;

        /* Compute spectrum sP */
        for(x=0; x<ncols; x++)
        {
            Adiag[x+ncols]=0.e0;
            E[x]=0.e0;
            for(y=0; y<nrows; y++)
            {
                sum=0.e0;
                for(iy=0; iy<ny; iy++)
                {
                    sum+=omega[iy+(y*ny)+(x*ny*nrows)]*sL[iy];
                }

                Adiag[x+ncols]+=sum*sum*mask[y*ncols+x];
                E[x]+=sum*im[y*ncols+x]*mask[y*ncols+x];
            }
        }
        if(lambda_sP>0.e0)
        {
            norm=0.e0;
            for(x=0; x<ncols; x++)
            {
                sP_old[x]=sP[x];
                norm+=sP[x];
            }
            norm/=ncols;
            lambda=lambda_sP*norm;
            Adiag[0        ] = 0.e0;
            Adiag[0+ncols  ]+= lambda;
            Adiag[0+ncols*2] =-lambda;
            for(x=1; x<ncols-1; x++)
            {
                Adiag[x]=-lambda;
                Adiag[x+ncols  ]+= 2.e0*lambda;
                Adiag[x+ncols*2] =-lambda;
            }
            Adiag[ncols-1        ] =-lambda;
            Adiag[ncols*2-1+ncols]+= lambda;
            Adiag[ncols*3-1+ncols] = 0.e0;

            info=cr2res_extract_slitdec_bandsol(Adiag, E, ncols, 3);
            for(x=0; x<ncols; x++) sP[x]=E[x];
        }
        else
        {
            for(x=0; x<ncols; x++)
            {
                sP_old[x]=sP[x];
                sP[x]=E[x]/Adiag[x+ncols];
            }
        }

        /* Compute the model */
        for(y=0; y<nrows; y++)
        {
            for(x=0; x<ncols; x++)
            {
                sum=0.e0;
                for(iy=0; iy<ny; iy++)
                    sum+=omega[iy+(y*ny)+(x*ny*nrows)]*sL[iy];
                model[y*ncols+x]=sum*sP[x];
            }
        }

        /* Compare model and data */
        sum=0.e0;
        isum=0;
        for(y=0; y<nrows; y++)
        {
            for(x=0;x<ncols; x++)
            {
                sum+=mask[y*ncols+x]*(model[y*ncols+x]-im[y*ncols+x]) *
                                (model[y*ncols+x]-im[y*ncols+x]);
                isum+=mask[y*ncols+x];
            }
        }
        dev=sqrt(sum/isum);

        /* Adjust the mask marking outlyers */
        for(y=0; y<nrows; y++)
        {
            for(x=0;x<ncols; x++)
            {
                if(fabs(model[y*ncols+x]-im[y*ncols+x])>6.*dev)
                    mask[y*ncols+x]=0;
                else mask[y*ncols+x]=1;
            }
        }

        /* Compute the change in the spectrum */
        sP_change=0.e0;
        sP_max=1.e0;
        for(x=0; x<ncols; x++)
        {
            if(sP[x]>sP_max) sP_max=sP[x];
            if(fabs(sP[x]-sP_old[x])>sP_change) sP_change=fabs(sP[x]-sP_old[x]);
        }
        /* Check the convergence */
    } while(iter++ < maxiter && sP_change > sP_stop*sP_max);

    cpl_free(E);
    cpl_free(sP_old);
    cpl_free(omega);
    cpl_free(Aij);
    cpl_free(bj);
    cpl_free(Adiag);

    return 0;
}




/*----------------------------------------------------------------------------*/
/**
  @brief    Helper function for cr2res_extract_slit_func_curved
  @param
  @return
 */
/*----------------------------------------------------------------------------*/

static int cr2res_extract_xi_zeta_tensors(
                    int ncols,      /* Swath width in pixels                                 */
                    int nrows,                     /* Extraction slit height in pixels                      */
                    double * ycen,            /* Order centre line offset from pixel row boundary      */
                    int * ycen_offset,        /* Order image column shift                              */
                    int y_lower_lim,               /* Number of detector pixels below the pixel containing  */
                                                   /* the central line yc.                                  */
                                                   /* ycen_offset[0]=0, absolute ycen=y_lower_lim+ycen_offset+ycen */
                    int osample,                   /* Subpixel ovsersampling factor                         */
                    double * PSF_curve    /* Parabolic fit to the slit image curvature.            */
                                                   /* For column d_x = PSF_curve[ncols *  3  + 0] +                */
                                                   /*                  PSF_curve[ncols *  3  + 1] *d_y +           */
                                                   /*                  PSF_curve[ncols *  3  + 2] *d_y^2,          */
                                                   /* where d_y is the offset from the central line ycen.   */
                                                   /* Thus central subpixel of omega[x][y'][delta_x][iy']   */
                                                   /* does not stick out of column x.                       */
)
{
  int x, xx, y, yy, ix, ix1, ix2, iy, iy1, iy2, ny;
  double step, delta, dy, w, d1, d2;
  xi_ref * xi;
  zeta_ref * zeta;
  int * m_zeta;

  ny=osample*(nrows+1)+1; /* The size of the sf array */

  /* Allocate working matrices */
  xi = cpl_malloc(ncols*ny*4*sizeof(xi_ref)); //xi_ref xi[ncols][ny][4],
                /* Convolution tensor telling the coordinates of detector*/
                /* pixels on which {x, iy} element falls and the         */
                /* corresponding projections.                            */
  zeta = cpl_malloc(ncols*nrows*3*(osample+1)*sizeof(zeta_ref)); // zeta_ref zeta[ncols][nrows][3*(osample+1)],

                /* Convolution tensor telling the coordinates*/
                /* of subpixels {x, iy} contributing to detector pixel   */
                /* {x, y}.                                               */
  m_zeta = cpl_malloc(ncols*nrows*sizeof(int)); //int m_zeta[ncols][nrows])
                    /* The actual number of controbuting elements in zeta    */
  step=1.e0/osample;

  for(x=0; x<ncols*ny*4; x++) /* Clean xi   */
  {
      xi[x].x=xi[x].y=xi[x].w=0;
  }

  for(x=0; x<ncols*nrows*3*(osample+1); x++) /* Clean zeta */
  {
    zeta[x].x =0;
    zeta[x].iy=0;
    zeta[x].w =0.;
  }
  for(x=0; x<ncols*nrows; x++){
    m_zeta[x] = 0;
  }
//    printf("%g %g %g; %g %g %g; %g %g %g\n",PSF_curve[313 *  3  + 0],PSF_curve[313 *  3  + 1],PSF_curve[313 *  3  + 2]
//                                           ,PSF_curve[314 *  3  + 0],PSF_curve[314 *  3  + 1],PSF_curve[314 *  3  + 2]
//                                           ,PSF_curve[315 *  3  + 0],PSF_curve[315 *  3  + 1],PSF_curve[315 *  3  + 2]);
/*
   Construct the xi and zeta tensors. They contain pixel references and contribution.
   values going from a given subpixel to other pixels (xi) and coming from other subpixels
   to a given detector pixel (zeta).
   Note, that xi and zeta are used in the equations for sL, sP and for the model but they
   do not involve the data, only the geometry. Thus it can be pre-computed once.
*/

  for(x=0; x<ncols; x++)
  {
   /*
     I promised to reconsider the initial offset. Here it is. For the original layout
     (no column shifts and discontinuities in ycen) there is pixel y that contains the
     central line yc. There are two options here (by construction of ycen that can be 0
     but cannot be 1): (1) yc is inside pixel y and (2) yc falls at the boundary between
     pixels y and y-1. yc cannot be at the foundary of pixels y+1 and y because we would
     select y+1 to be pixel y in that case.

     Next we need to define starting and ending indices iy for sL subpixels that contribute
     to pixel y. I call them iy1 and iy2. For both cases we assume osample+1 subpixels covering
     pixel y (wierd). So for case 1 iy1 will be (y-1)*osample and iy2 == y*osample. Special
     treatment of the boundary subpixels will compensate for introducing extra subpixel in
     case 1. In case 2 things are more logical: iy1=(yc-y)*osample+(y-1)*osample;
     iy2=(y+1-yc)*osample)+(y-1)*osample. ycen is yc-y making things simpler. Note also that
     the same pattern repeates for all rows: we only need to initialize iy1 and iy2 and keep
     incrementing them by osample.
   */

    iy2=osample-floor(ycen[x]/step)-1;
    iy1=iy2-osample;

   /*
     Handling partial subpixels cut by detector pixel rows is again tricky. Here we have three
     cases (mostly because of the decision to assume that we always have osample+1 subpixels
     per one detector pixel). Here d1 is the fraction of the subpixel iy1 inside detector pixel y.
     d2 is then the fraction of subpixel iy2 inside detector pixel y. By definition d1+d2==step.
     Case 1: ycen falls on the top boundary of each detector pixel (ycen == 1). Here we conclude
             that the first subpixel is fully contained inside pixel y and d1 is set to step.
     Case 2: ycen falls on the bottom boundary of each detector pixel (ycen == 0). Here we conclude
             that the first subpixel is totally outside of pixel y and d1 is set to 0.
     Case 3: ycen falls inside of each pixel (0>ycen>1). In this case d1 is set to the fraction of
             the first step contained inside of each pixel.
     And BTW, this also means that central line coinsides with the upper boundary of subpixel iy2
     when the y loop reaches pixel y_lower_lim. In other words:

     dy=(iy-(y_lower_lim+ycen[x])*osample)*step-0.5*step

   */

    if(iy2==0)      d1=step;                /* Case 1 */
    else if(iy1==0) d1=0.e0;                /* Case 2 */
    else            d1=fmod(ycen[x], step); /* Case 3: This is very clever */
    d2=step-d1;

   /*
     The final hurdle for 2D slit decomposition is to construct two 3D reference tensors. We proceed
     similar to 1D case except that now each iy subpixel can be shifted left or right following
     the curvature of the slit image on the detector. We assume for now that each subpixel is
     exactly 1 detector pixel wide. This may not be exactly true if the curvature changes accross
     the focal plane but will deal with it when the necessity will become apparent. For now we
     just assume that a shift delta the weight w assigned to subpixel iy is divided between
     ix1=int(delta) and ix2=int(delta)+signum(delta) as (1-|delta-ix1|)*w and |delta-ix1|*w.

     The curvature is given by a quadratic polynomial evaluated from an approximation for column
     x: delta = PSF_curve[x *  3  + 0] + PSF_curve[x *  3  + 1] * (y-yc[x]) + PSF_curve[x *  3  + 2] * (y-yc[x])^2.
     It looks easy except that y and yc are set in the global detector coordinate system rather than
     in the shifted and cropped swath passed to slit_func_2d. One possible solution I will try here
     is to modify PSF_curve before the call such as:
     delta = PSF_curve'[x][0] + PSF_curve'[x][1] * (y'-ycen[x]) + PSF_curve'[x][2] * (y'-ycen[x])^2
     where y' = y - floor(yc).
    */

//    dy=-ceil((ycen[x]+y_lower_lim)*osample)*step-step*0.5;
//    printf("dy=%g,", dy);
    dy=-(y_lower_lim*osample+floor(ycen[x]/step)+0.5)*step; /* Define initial distance from ycen       */
                                                           /* It is given by the center of the first  */
                                                           /* subpixel falling into pixel y_lower_lim */
//    printf("x=%d, dy=%g, step=%g, iy1=%d, iy2=%d, ycen=%g\n", x, dy, step, iy1, iy2, ycen[x]);

/*
   Now we go detector pixels x and y incrementing subpixels looking for their controibutions
   to the current and adjacent pixels. Note that the curvature/tilt of the projected slit
   image could be so large that subpixel iy may no contribute to column x at all. On the
   other hand, subpixels around ycen by definition must contribute to pixel x,y.
   3rd index in xi refers corners of pixel xx,y: 0:LL, 1:LR, 2:UL, 3:UR.
*/

    for(y=0; y<nrows; y++)
    {
      iy1+=osample;  // Bottom subpixel falling in row y
      iy2+=osample;  // Top subpixel falling in row y
      dy-=step;
      for(iy=iy1; iy<=iy2; iy++)
      {
        if(iy==iy1)      w=d1;
        else if(iy==iy2) w=d2;
        else             w=step;
        dy+=step;
        delta = (PSF_curve[x *  3  + 1] + PSF_curve[x *  3  + 2] * dy) * dy;
        ix1=delta;
        ix2=ix1+signum(delta);

/* Three cases: bottom boundary of row y, intermediate subpixels and top boundary */

        if(iy==iy1)   /* Subpixel iy is entering detector row y        */
        {
          if(ix1<ix2) /* Subpixel iy shifts to the right from column x */
          {
            if(x+ix1>=0 && x+ix2<ncols)
            {
              xx=x+ix1;
              yy=y+ycen_offset[x]-ycen_offset[xx];
              xi[x * (ncols * nrows) + iy * ncols + 1].x=xx;
              xi[x * (ncols * nrows) + iy * ncols + 1].y=yy;
              xi[x * (ncols * nrows) + iy * ncols + 1].w=w-fabs(delta-ix1)*w;
              if(xx>=0 && xx<ncols && yy>=0 && yy<nrows && xi[x * (ncols * nrows) + iy * ncols + 1].w>0)
              {
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].x =x;
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].iy=iy;
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].w=xi[x * (ncols * nrows) + iy * ncols + 1].w;
                m_zeta[xx *  ncols  + yy]++;
              }
              xx=x+ix2;
              yy=y+ycen_offset[x]-ycen_offset[xx];
              xi[x * (ncols * nrows) + iy * ncols + 0].x=xx;
              xi[x * (ncols * nrows) + iy * ncols + 0].y=yy;
              xi[x * (ncols * nrows) + iy * ncols + 0].w=fabs(delta-ix1)*w;
              if(xx>=0 && xx<ncols && yy>=0 && yy<nrows && xi[x * (ncols * nrows) + iy * ncols + 0].w>0)
              {
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].x =x;
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].iy=iy;
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].w=xi[x * (ncols * nrows) + iy * ncols + 0].w;
                m_zeta[xx *  ncols  + yy]++;
              }
            }
          }
          else if(ix1>ix2) /* Subpixel iy shifts to the left from column x */
          {
            if(x+ix2>=0 && x+ix1<ncols)
            {
              xx=x+ix2;
              yy=y+ycen_offset[x]-ycen_offset[xx];
              xi[x * (ncols * nrows) + iy * ncols + 1].x=xx;
              xi[x * (ncols * nrows) + iy * ncols + 1].y=yy;
              xi[x * (ncols * nrows) + iy * ncols + 1].w=fabs(delta-ix1)*w;
              if(xx>=0 && xx<ncols && yy>=0 && yy<nrows && xi[x * (ncols * nrows) + iy * ncols + 1].w>0)
              {
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].x =x;
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].iy=iy;
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].w=xi[x * (ncols * nrows) + iy * ncols + 1].w;
                m_zeta[xx *  ncols  + yy]++;
              }
              xx=x+ix1;
              yy=y+ycen_offset[x]-ycen_offset[xx];
              xi[x * (ncols * nrows) + iy * ncols + 0].x=xx;
              xi[x * (ncols * nrows) + iy * ncols + 0].y=yy;
              xi[x * (ncols * nrows) + iy * ncols + 0].w=w-fabs(delta-ix1)*w;
              if(xx>=0 && xx<ncols && yy>=0 && yy<nrows && xi[x * (ncols * nrows) + iy * ncols + 0].w>0)
              {
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].x =x;
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].iy=iy;
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].w=xi[x * (ncols * nrows) + iy * ncols + 0].w;
                m_zeta[xx *  ncols  + yy]++;
              }
            }
          }
          else             /* Subpixel iy stays inside column x */
          {
            xx=x+ix1;
            yy=y+ycen_offset[x]-ycen_offset[xx];
            xi[x * (ncols * nrows) + iy * ncols + 0].x=xx;
            xi[x * (ncols * nrows) + iy * ncols + 0].y=yy;
            xi[x * (ncols * nrows) + iy * ncols + 0].w=w;
            if(xx>=0 && xx<ncols && yy>=0 && yy<nrows && w>0)
            {
              zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].x =x;
              zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].iy=iy;
              zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].w=w;
              m_zeta[xx *  ncols  + yy]++;
            }
          }
        }
        else if(iy==iy2)   /* Subpixel iy is leaving detector row y    */
        {
          if(ix1<ix2) /* Subpixel iy shifts to the right from column x */
          {
            if(x+ix1>=0 && x+ix2<ncols)
            {
              xx=x+ix1;
              yy=y+ycen_offset[x]-ycen_offset[xx];
              xi[x * (ncols * nrows) + iy * ncols + 3].x=xx;
              xi[x * (ncols * nrows) + iy * ncols + 3].y=yy;
              xi[x * (ncols * nrows) + iy * ncols + 3].w=w-fabs(delta-ix1)*w;
              if(xx>=0 && xx<ncols && yy>=0 && yy<nrows && xi[x * (ncols * nrows) + iy * ncols + 3].w>0)
              {
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].x =x;
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].iy=iy;
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].w=xi[x * (ncols * nrows) + iy * ncols + 3].w;
                m_zeta[xx *  ncols  + yy]++;
              }
              xx=x+ix2;
              yy=y+ycen_offset[x]-ycen_offset[xx];
              xi[x * (ncols * nrows) + iy * ncols + 2].x=xx;
              xi[x * (ncols * nrows) + iy * ncols + 2].y=yy;
              xi[x * (ncols * nrows) + iy * ncols + 2].w=fabs(delta-ix1)*w;
              if(xx>=0 && xx<ncols && yy>=0 && yy<nrows && xi[x * (ncols * nrows) + iy * ncols + 2].w>0)
              {
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].x =x;
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].iy=iy;
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].w=xi[x * (ncols * nrows) + iy * ncols + 2].w;
                m_zeta[xx *  ncols  + yy]++;
              }
            }
          }
          else if(ix1>ix2) /* Subpixel iy shifts to the left from column x */
          {
            if(x+ix2>=0 && x+ix1<ncols)
            {
              xx=x+ix2;
              yy=y+ycen_offset[x]-ycen_offset[xx];
              xi[x * (ncols * nrows) + iy * ncols + 3].x=xx;
              xi[x * (ncols * nrows) + iy * ncols + 3].y=yy;
              xi[x * (ncols * nrows) + iy * ncols + 3].w=fabs(delta-ix1)*w;
              if(xx>=0 && xx<ncols && yy>=0 && yy<nrows && xi[x * (ncols * nrows) + iy * ncols + 3].w>0)
              {
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].x =x;
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].iy=iy;
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].w=xi[x * (ncols * nrows) + iy * ncols + 3].w;
                m_zeta[xx *  ncols  + yy]++;
              }
              xx=x+ix1;
              yy=y+ycen_offset[x]-ycen_offset[xx];
              xi[x * (ncols * nrows) + iy * ncols + 2].x=xx;
              xi[x * (ncols * nrows) + iy * ncols + 2].y=yy;
              xi[x * (ncols * nrows) + iy * ncols + 2].w=w-fabs(delta-ix1)*w;
              if(xx>=0 && xx<ncols && yy>=0 && yy<nrows && xi[x * (ncols * nrows) + iy * ncols + 2].w>0)
              {
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].x =x;
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].iy=iy;
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].w=xi[x * (ncols * nrows) + iy * ncols + 2].w;
                m_zeta[xx *  ncols  + yy]++;
              }
            }
          }
          else             /* Subpixel iy stays inside column x        */
          {
            xx=x+ix1;
            yy=y+ycen_offset[x]-ycen_offset[xx];
            xi[x * (ncols * nrows) + iy * ncols + 2].x=xx;
            xi[x * (ncols * nrows) + iy * ncols + 2].y=yy;
            xi[x * (ncols * nrows) + iy * ncols + 2].w=w;
            if(xx>=0 && xx<ncols && yy>=0 && yy<nrows && w>0)
            {
              zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].x =x;
              zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].iy=iy;
              zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].w=w;
              m_zeta[xx *  ncols  + yy]++;
            }
          }
        }
        else               /* Subpixel iy is fully inside detector row y */
        {
          if(ix1<ix2) /* Subpixel iy shifts to the right from column x   */
          {
            if(x+ix1>=0 && x+ix2<ncols)
            {
              xx=x+ix1;
              yy=y+ycen_offset[x]-ycen_offset[xx];
              xi[x * (ncols * nrows) + iy * ncols + 1].x=xx;
              xi[x * (ncols * nrows) + iy * ncols + 1].y=yy;
              xi[x * (ncols * nrows) + iy * ncols + 1].w=w-fabs(delta-ix1)*w;
              if(xx>=0 && xx<ncols && yy>=0 && yy<nrows && xi[x * (ncols * nrows) + iy * ncols + 1].w>0)
              {
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].x =x;
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].iy=iy;
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].w=xi[x * (ncols * nrows) + iy * ncols + 1].w;
                m_zeta[xx *  ncols  + yy]++;
              }
              xx=x+ix2;
              yy=y+ycen_offset[x]-ycen_offset[xx];
              xi[x * (ncols * nrows) + iy * ncols + 0].x=xx;
              xi[x * (ncols * nrows) + iy * ncols + 0].y=yy;
              xi[x * (ncols * nrows) + iy * ncols + 0].w=fabs(delta-ix1)*w;
              if(xx>=0 && xx<ncols && yy>=0 && yy<nrows && xi[x * (ncols * nrows) + iy * ncols + 0].w>0)
              {
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].x =x;
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].iy=iy;
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].w=xi[x * (ncols * nrows) + iy * ncols + 0].w;
                m_zeta[xx *  ncols  + yy]++;
              }
            }
          }
          else if(ix1>ix2) /* Subpixel iy shifts to the left from column x */
          {
            if(x+ix2>=0 && x+ix1<ncols)
            {
              xx=x+ix2;
              yy=y+ycen_offset[x]-ycen_offset[xx];
              xi[x * (ncols * nrows) + iy * ncols + 1].x=xx;
              xi[x * (ncols * nrows) + iy * ncols + 1].y=yy;
              xi[x * (ncols * nrows) + iy * ncols + 1].w=fabs(delta-ix1)*w;
              if(xx>=0 && xx<ncols && yy>=0 && yy<nrows && xi[x * (ncols * nrows) + iy * ncols + 1].w>0)
              {
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].x =x;
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].iy=iy;
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].w=xi[x * (ncols * nrows) + iy * ncols + 1].w;
                m_zeta[xx *  ncols  + yy]++;
              }
              xx=x+ix1;
              yy=y+ycen_offset[x]-ycen_offset[xx];
              xi[x * (ncols * nrows) + iy * ncols + 0].x=xx;
              xi[x * (ncols * nrows) + iy * ncols + 0].y=yy;
              xi[x * (ncols * nrows) + iy * ncols + 0].w=w-fabs(delta-ix1)*w;
              if(xx>=0 && xx<ncols && yy>=0 && yy<nrows && xi[x * (ncols * nrows) + iy * ncols + 0].w>0)
              {
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].x =x;
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].iy=iy;
                zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].w=xi[x * (ncols * nrows) + iy * ncols + 0].w;
                m_zeta[xx *  ncols  + yy]++;
              }
            }
          }
          else             /* Subpixel iy stays inside column x */
          {
            xx=x+ix2;
            yy=y+ycen_offset[x]-ycen_offset[xx];
            xi[x * (ncols * nrows) + iy * ncols + 0].x=xx;
            xi[x * (ncols * nrows) + iy * ncols + 0].y=yy;
            xi[x * (ncols * nrows) + iy * ncols + 0].w=w;
            if(xx>=0 && xx<ncols && yy>=0 && yy<nrows && w>0)
            {
              zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].x =x;
              zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].iy=iy;
              zeta[xx * nrows *  ncols  + yy * ncols + m_zeta[xx *  ncols  + yy]].w=w;
              m_zeta[xx *  ncols  + yy]++;
            }
          }
        }
      }
    }
  }
  return 0;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Slit decomposition of single swath with slit tilt & curvature
  @param
  @return
 */
/*----------------------------------------------------------------------------*/
static int cr2res_extract_slit_func_curved(
         int         ncols,
         int         nrows,
         int         osample,
         double  *   im,
         double  *   pix_unc,
         int     *   mask,
         double  *   ycen,
         int     *   ycen_offset,
         int         y_lower_lim,
         xi_ref  *   xi,
         zeta_ref *  zeta,
         int     *   m_zeta,
         double  *   sL,
         double  *   sP,
         double  *   model,
         double  *   unc,
         double      lambda_sP,
         double      lambda_sL,
         double      sP_stop,
         int         maxiter)

/* ORIGINAL DECLARATION
static int cr2res_extract_slit_func_curved(int ncols,      // Swath width in pixels
                     int nrows,                     // Extraction slit height in pixels
                     int nx,                        // Range of columns affected by PSF tilt: nx=2*delta_x+1
                     int ny,                        // Size of the slit function array: ny=osample(nrows+1)+1
                     double im[nrows][ncols],       // Image to be decomposed
                     double pix_unc[nrows][ncols],  // Individual pixel uncertainties. Set to zero if unknown
                     byte mask[nrows][ncols],       // Initial and final mask for the swath
                     double ycen[ncols],            // Order centre line offset from pixel row boundary
                     int ycen_offset[ncols],        // Order image column shift
                     int y_lower_lim,               // Number of detector pixels below the pixel containing
                                                    // the central line yc.
                     int osample,                   // Subpixel ovsersampling factor
                     double PSF_curve[ncols][3],    // Parabolic fit to the slit image curvature.
                                                    // For column d_x = PSF_curve[ncols *  3  + 0] +
                                                    //                  PSF_curve[ncols *  3  + 1] *d_y +
                                                    //                  PSF_curve[ncols *  3  + 2] *d_y^2,
                                                    // where d_y is the offset from the central line ycen.
                                                    // Thus central subpixel of omega[x][y'][delta_x][iy']
                                                    // does not stick out of column x.
                     double lambda_sP,              // Smoothing parameter for the spectrum, coiuld be zero
                     double lambda_sL,              // Smoothing parameter for the slit function, usually >0
                     double sP[ncols],              // Spectrum resulting from decomposition
                     double sL[ny],                 // Slit function resulting from decomposition
                     double model[nrows][ncols],    // Model constructed from sp and sf
                     double unc[ncols],             // Spectrum uncertainties based on data - model
                     xi_ref xi[ncols][ny][4],       // Convolution tensor telling the coordinates of detector
                                                    // pixels on which {x, iy} element falls and the
                                                    // corresponding projections.
                     zeta_ref zeta[ncols][nrows][3*(osample+1)],// Convolution tensor telling the coordinates
                                                    // of subpixels {x, iy} contributing to detector pixel
                                                    // {x, y}.
                     int m_zeta[ncols][nrows],      // The actual number of controbuting elements in zeta
                     double sP_old[ncols],          // Work array to control the convergence
                     double l_Aij[],                // Various LAPACK arrays (ny*ny)
                     double l_bj[ny],               // ny
                     double p_Aij[],                // Various LAPACK arrays (ncols*ncols)
                     double p_bj[ncols])            // ncols (RHS)
*/
{
  int x, xx, xxx, y, yy, iy, jy, n, m, nx, ny;
  double step, delta_x, sum, norm, dev, lambda, diag_tot, ww, www, sP_change, sP_max;
  int info, iter, isum;

  delta_x=nx/2;           /* Maximum horizontal shift in detector pixels due to slit image curvature         */
  ny=osample*(nrows+1)+1; /* The size of the sL array. Extra osample is because ycen can be between 0 and 1. */
  step=1.e0/osample;

  double * l_Aij = cpl_malloc(ny * ny * sizeof(double));
  double * l_bj = cpl_malloc(ny * sizeof(double));
  double * p_Aij = cpl_malloc(ncols * ncols * sizeof(double));
  double * p_bj = cpl_malloc(ncols * sizeof(double));
  double * sP_old = cpl_malloc(ncols*sizeof(double)); // double sP_old[ncols];

/* Loop through sL , sP reconstruction until convergence is reached */
  iter=0;
  do
  {

/*
  ====================================================================
  Compute slit function sL
*/

/* Prepare the RHS and the matrix */
    for(iy=0; iy<ny; iy++)
    {
      l_bj[iy]=0.e0;                                       /* Clean RHS                */
      for(jy=0; jy<=4*osample; jy++) l_Aij[jy+ny*jy]=0.e0; /* Clean matrix row         */
    }

/* Fill in SLE arrays for slit function */
   	diag_tot=0.e0;

    for(iy=0; iy<ny; iy++)
    {
      for(x=0; x<ncols; x++)
      {
        n=0;
        ww=xi[x * (ncols * nrows) + iy * ncols + n].w;
        if(ww>0)
        {
          xx=xi[x * (ncols * nrows) + iy * ncols + n].x;
          yy=xi[x * (ncols * nrows) + iy * ncols + n].y;
          if(m_zeta[xx *  ncols  + yy]>0 && xx>=0 && xx<ncols && yy>=0 && yy<nrows)
          {
            for(m=0; m<m_zeta[xx *  ncols  + yy]; m++)
            {
              xxx=zeta[xx * (ncols * nrows) + yy * ncols + m].x;
              jy =zeta[xx * (ncols * nrows) + yy * ncols + m].iy;
              www=zeta[xx * (ncols * nrows) + yy * ncols + m].w;
              l_Aij[iy+ny*(jy-iy+2*osample)]+=sP[xxx]*sP[x]*www*ww*mask[yy * ncols + xx];
            }
            l_bj[iy]+=im[yy * ncols + xx]*mask[yy * ncols + xx]*sP[x]*ww;
          }
        }
        n=1;
        ww=xi[x * (ncols * nrows) + iy * ncols + n].w;
        if(ww>0)
        {
          xx=xi[x * (ncols * nrows) + iy * ncols + n].x;
          yy=xi[x * (ncols * nrows) + iy * ncols + n].y;
          if(m_zeta[xx *  ncols  + yy]>0 && xx>=0 && xx<ncols && yy>=0 && yy<nrows)
          {
            for(m=0; m<m_zeta[xx *  ncols  + yy]; m++)
            {
              xxx=zeta[xx * (ncols * nrows) + yy * ncols + m].x;
              jy =zeta[xx * (ncols * nrows) + yy * ncols + m].iy;
              www=zeta[xx * (ncols * nrows) + yy * ncols + m].w;
              l_Aij[iy+ny*(jy-iy+2*osample)]+=sP[xxx]*sP[x]*www*ww*mask[yy * ncols + xx];
            }
            l_bj[iy]+=im[yy * ncols + xx]*mask[yy * ncols + xx]*sP[x]*ww;
          }
        }
        n=2;
        ww=xi[x * (ncols * nrows) + iy * ncols + n].w;
        if(ww>0)
        {
          xx=xi[x * (ncols * nrows) + iy * ncols + n].x;
          yy=xi[x * (ncols * nrows) + iy * ncols + n].y;
          if(m_zeta[xx *  ncols  + yy]>0 && xx>=0 && xx<ncols && yy>=0 && yy<nrows)
          {
            for(m=0; m<m_zeta[xx *  ncols  + yy]; m++)
            {
              xxx=zeta[xx * (ncols * nrows) + yy * ncols + m].x;
              jy =zeta[xx * (ncols * nrows) + yy * ncols + m].iy;
              www=zeta[xx * (ncols * nrows) + yy * ncols + m].w;
              l_Aij[iy+ny*(jy-iy+2*osample)]+=sP[xxx]*sP[x]*www*ww*mask[yy * ncols + xx];
            }
            l_bj[iy]+=im[yy * ncols + xx]*mask[yy * ncols + xx]*sP[x]*ww;
          }
        }
        n=3;
        ww=xi[x * (ncols * nrows) + iy * ncols + n].w;
        if(ww>0)
        {
          xx=xi[x * (ncols * nrows) + iy * ncols + n].x;
          yy=xi[x * (ncols * nrows) + iy * ncols + n].y;
          if(m_zeta[xx *  ncols  + yy]>0 && xx>=0 && xx<ncols && yy>=0 && yy<nrows)
          {
            for(m=0; m<m_zeta[xx *  ncols  + yy]; m++)
            {
              xxx=zeta[xx * (ncols * nrows) + yy * ncols + m].x;
              jy =zeta[xx * (ncols * nrows) + yy * ncols + m].iy;
              www=zeta[xx * (ncols * nrows) + yy * ncols + m].w;
              l_Aij[iy+ny*(jy-iy+2*osample)]+=sP[xxx]*sP[x]*www*ww*mask[yy * ncols + xx];
            }
            l_bj[iy]+=im[yy * ncols + xx]*mask[yy * ncols + xx]*sP[x]*ww;
          }
        }
      }
      diag_tot+=l_Aij[iy+ny*2*osample];
    }

/* Scale regularization parameters */

    lambda=lambda_sL*diag_tot/ny;

/* Add regularization parts for the SLE matrix */

    l_Aij[ny*2*osample]    +=lambda;           /* Main diagonal  */
    l_Aij[ny*(2*osample+1)]-=lambda;           /* Upper diagonal */
    for(iy=1; iy<ny-1; iy++)
    {
      l_Aij[iy+ny*(2*osample-1)]-=lambda;      /* Lower diagonal */
      l_Aij[iy+ny*2*osample    ]+=lambda*2.e0; /* Main diagonal  */
      l_Aij[iy+ny*(2*osample+1)]-=lambda;      /* Upper diagonal */
    }
    l_Aij[ny-1+ny*(2*osample-1)]-=lambda;      /* Lower diagonal */
    l_Aij[ny-1+ny*2*osample]    +=lambda;      /* Main diagonal  */

/* Solve the system of equations */

    info=cr2res_extract_slitdec_bandsol(l_Aij, l_bj, ny, 4*osample+1);
    if(info) printf("info(sL)=%d\n", info);

/* Normalize the slit function */

    norm=0.e0;
    for(iy=0; iy<ny; iy++)
    {
      sL[iy]=l_bj[iy];
      norm+=sL[iy];
    }
    norm/=osample;
    for(iy=0; iy<ny; iy++) sL[iy]/=norm;
/*
  ====================================================================
  Compute spectrum sP
*/
    for(x=0; x<ncols; x++)
    {
      for(xx=0; xx<5; xx++) p_Aij[xx*ncols+x]=0.;
      p_bj[x]=0;
    }

    for(x=0; x<ncols; x++)
    {
      for(iy=0; iy<ny; iy++)
      {
        n=0;
        ww=xi[x * (ncols * nrows) + iy * ncols + n].w;
        if(ww>0)
        {
          xx=xi[x * (ncols * nrows) + iy * ncols + n].x;
          yy=xi[x * (ncols * nrows) + iy * ncols + n].y;
          if(m_zeta[xx *  ncols  + yy]>0 && xx>=0 && xx<ncols && yy>=0 && yy<nrows)
          {
            for(m=0; m<m_zeta[xx *  ncols  + yy]; m++)
            {
              xxx=zeta[xx * (ncols * nrows) + yy * ncols + m].x;
              jy =zeta[xx * (ncols * nrows) + yy * ncols + m].iy;
              www=zeta[xx * (ncols * nrows) + yy * ncols + m].w;
              p_Aij[x+ncols*(xxx-x+2)]+=sL[jy]*sL[iy]*www*ww*mask[yy * ncols + xx];
            }
            p_bj[x]+=im[yy * ncols + xx]*mask[yy * ncols + xx]*sL[iy]*ww;
          }
        }
        n=1;
        ww=xi[x * (ncols * nrows) + iy * ncols + n].w;
        if(ww>0)
        {
          xx=xi[x * (ncols * nrows) + iy * ncols + n].x;
          yy=xi[x * (ncols * nrows) + iy * ncols + n].y;
          if(m_zeta[xx *  ncols  + yy]>0 && xx>=0 && xx<ncols && yy>=0 && yy<nrows)
          {
            for(m=0; m<m_zeta[xx *  ncols  + yy]; m++)
            {
              xxx=zeta[xx * (ncols * nrows) + yy * ncols + m].x;
              jy =zeta[xx * nrows *  ncols  + yy * ncols + m].iy;
              www=zeta[xx * nrows *  ncols  + yy * ncols + m].w;
              p_Aij[x+ncols*(xxx-x+2)]+=sL[jy]*sL[iy]*www*ww*mask[yy * ncols + xx];
            }
            p_bj[x]+=im[yy * ncols + xx]*mask[yy * ncols + xx]*sL[iy]*ww;
          }
        }
        n=2;
        ww=xi[x * (ncols * nrows) + iy * ncols + n].w;
        if(ww>0)
        {
          xx=xi[x * (ncols * nrows) + iy * ncols + n].x;
          yy=xi[x * (ncols * nrows) + iy * ncols + n].y;
          if(m_zeta[xx *  ncols  + yy]>0 && xx>=0 && xx<ncols && yy>=0 && yy<nrows)
          {
            for(m=0; m<m_zeta[xx *  ncols  + yy]; m++)
            {
              xxx=zeta[xx * nrows *  ncols  + yy * ncols + m].x;
              jy =zeta[xx * nrows *  ncols  + yy * ncols + m].iy;
              www=zeta[xx * nrows *  ncols  + yy * ncols + m].w;
              p_Aij[x+ncols*(xxx-x+2)]+=sL[jy]*sL[iy]*www*ww*mask[yy * ncols + xx];
            }
            p_bj[x]+=im[yy * ncols + xx]*mask[yy * ncols + xx]*sL[iy]*ww;
          }
        }
        n=3;
        ww=xi[x * (ncols * nrows) + iy * ncols + n].w;
        if(ww>0)
        {
          xx=xi[x * (ncols * nrows) + iy * ncols + n].x;
          yy=xi[x * (ncols * nrows) + iy * ncols + n].y;
          if(m_zeta[xx *  ncols  + yy]>0 && xx>=0 && xx<ncols && yy>=0 && yy<nrows)
          {
            for(m=0; m<m_zeta[xx *  ncols  + yy]; m++)
            {
              xxx=zeta[xx * (ncols * nrows) + yy * ncols + m].x;
              jy =zeta[xx * (ncols * nrows) + yy * ncols + m].iy;
              www=zeta[xx * (ncols * nrows) + yy * ncols + m].w;
              p_Aij[x+ncols*(xxx-x+2)]+=sL[jy]*sL[iy]*www*ww*mask[yy * ncols + xx];
            }
            p_bj[x]+=im[yy * ncols + xx]*mask[yy * ncols + xx]*sL[iy]*ww;
          }
        }
      }
    }

    for(x=0; x<ncols; x++) sP_old[x]=sP[x];

    if(lambda_sP>0.e0)
    {
      norm=0.e0;
      for(x=0; x<ncols; x++)
      {
        norm+=sP[x];
      }
      norm/=ncols;
      lambda=lambda_sP*norm;           /* Scale regularization parameter */

      p_Aij[ncols*2]+=lambda;          /* Main diagonal  */
      p_Aij[ncols*3]-=lambda;          /* Upper diagonal */
      for(x=1; x<ncols-1; x++)
      {
        p_Aij[x+ncols  ]-=lambda;      /* Lower diagonal */
        p_Aij[x+ncols*2]+=lambda*2.e0; /* Main diagonal  */
        p_Aij[x+ncols*3]-=lambda;      /* Upper diagonal */
      }
      p_Aij[ncols-1+ncols  ]-=lambda;  /* Lower diagonal */
      p_Aij[ncols-1+ncols*2]+=lambda;  /* Main diagonal  */
    }

/* Solve the system of equations */

    info=cr2res_extract_slitdec_bandsol(p_Aij, p_bj, ncols, 5);
    if(info) printf("info(sP)=%d\n", info);

    for(x=0; x<ncols; x++) sP[x]=p_bj[x];

/* Compute the model */

    for(y=0; y<nrows; y++)
    {
      for(x=0; x<ncols; x++)
      {
        model[y * ncols + x]=0.;
      }
    }

  	for(y=0; y<nrows; y++)
  	{
      for(x=0; x<ncols; x++)
      {
        for(m=0; m<m_zeta[x *  ncols  + y]; m++)
        {
          xx=zeta[x * (ncols * nrows) + y * ncols + m].x;
          iy=zeta[x * (ncols * nrows) + y * ncols + m].iy;
          ww=zeta[x * (ncols * nrows) + y * ncols + m].w;
          model[y * ncols + x]+=sP[xx]*sL[iy]*ww;
        }
      }
    }

/* Compare model and data */

    sum=0.e0;
    isum=0;
    for(y=0; y<nrows; y++)
    {
     	for(x=delta_x; x<ncols-delta_x; x++)
     	{
        sum+=mask[y * ncols + x]*(model[y * ncols + x]-im[y * ncols + x])*(model[y * ncols + x]-im[y * ncols + x]);
        isum+=mask[y * ncols + x];
     	}
    }
    dev=sqrt(sum/isum);

/* Adjust the mask marking outlyers */

    for(y=0; y<nrows; y++)
    {
     	for(x=delta_x; x<ncols-delta_x; x++)
     	{
        if(fabs(model[y * ncols + x]-im[y * ncols + x])>6.*dev) mask[y * ncols + x]=0; else mask[y * ncols + x]=1;
     	}
    }
    //printf("iter=%d, dev=%g\n", iter, dev);

/* Compute the change in the spectrum */

    sP_change=0.e0;
    sP_max=1.e0;
    for(x=0; x<ncols; x++)
    {
      if(sP[x]>sP_max) sP_max=sP[x];
      if(fabs(sP[x]-sP_old[x])>sP_change) sP_change=fabs(sP[x]-sP_old[x]);
    }

/* Check for convergence */

  } while(iter++ < maxiter && sP_change > sP_stop*sP_max);

/* Uncertainty estimate */

  for(x=0; x<ncols; x++)
  {
    unc[x]=0.;
    p_bj[x]=0.;
  }

  for(y=0; y<nrows; y++)
  {
    for(x=0; x<ncols; x++)
    {
      for(m=0; m<m_zeta[x *  ncols  + y]; m++) // Loop through all pixels contributing to x,y
      {
        xx=zeta[x * (ncols * nrows) + y * ncols + m].x;
        iy=zeta[x * (ncols * nrows) + y * ncols + m].iy;
        ww=zeta[x * (ncols * nrows) + y * ncols + m].w;
        unc[xx]+=(im[y * ncols + x]-model[y * ncols + x])*(im[y * ncols + x]-model[y * ncols + x])*
                  ww*mask[y * ncols + x];
        unc[xx]+=pix_unc[y * ncols + x]*pix_unc[y * ncols + x]*ww*mask[y * ncols + x];
        p_bj[xx]+=ww*mask[y * ncols + x];   // Norm
      }
    }
  }

  for(x=0; x<ncols; x++)
  {
    unc[x]=sqrt(unc[x]/p_bj[x]*nrows);
  }

  cpl_free(l_Aij);
  cpl_free(l_bj);
  cpl_free(p_Aij);
  cpl_free(p_bj);
  cpl_free(sP_old);

	return 0;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Solve a sparse system of linear equations
  @param    a   2D array [n,nd]i
  @param    r   array of RHS of size n
  @param    n   number of equations
  @param    nd  width of the band (3 for tri-diagonal system)
  @return   0 on success, -1 on incorrect size of "a" and -4 on
            degenerate matrix

  Solve a sparse system of linear equations with band-diagonal matrix.
  Band is assumed to be symmetrix relative to the main diaginal.

  nd must be an odd number. The main diagonal should be in a(*,nd/2)
  The first lower subdiagonal should be in a(1:n-1,nd/2-1), the first
  upper subdiagonal is in a(0:n-2,nd/2+1) etc. For example:
                    / 0 0 X X X \
                    | 0 X X X X |
                    | X X X X X |
                    | X X X X X |
              A =   | X X X X X |
                    | X X X X X |
                    | X X X X X |
                    | X X X X 0 |
                    \ X X X 0 0 /
 */
/*----------------------------------------------------------------------------*/
static int cr2res_extract_slitdec_bandsol(
        double  *   a,
        double  *   r,
        int         n,
        int         nd)
{
    double aa;
    int i, j, k;

    //  if(mod(nd,2)==0) return -1;

    /* Forward sweep */
    for(i=0; i<n-1; i++)
    {
        aa=a[i+n*(nd/2)];
        //    if(aa==0.e0) return -3;
        r[i]/=aa;
        for(j=0; j<nd; j++) a[i+j*n]/=aa;
        for(j=1; j<min(nd/2+1,n-i); j++)
        {
            aa=a[i+j+n*(nd/2-j)];
            //      if(aa==0.e0) return -j;
            r[i+j]-=r[i]*aa;
            for(k=0; k<n*(nd-j); k+=n) a[i+j+k]-=a[i+k+n*j]*aa;
        }
    }

    /* Backward sweep */
    r[n-1]/=a[n-1+n*(nd/2)];
    for(i=n-1; i>0; i--)
    {
        for(j=1; j<=min(nd/2,i); j++) r[i-j]-=r[i]*a[i-j+n*(nd/2+j)];
        //    if(a[i-1+n*(nd/2)]==0.e0) return -5;
        r[i-1]/=a[i-1+n*(nd/2)];
    }

    //  if(a[n*(nd/2)]==0.e0) return -6;
    r[0]/=a[n*(nd/2)];
    return 0;
}


/** @} */
