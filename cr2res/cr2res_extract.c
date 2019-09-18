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

#pragma GCC optimize("O0")

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
#include "cr2res_wave.h"

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

typedef unsigned char byte;
#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))
#define signum(a) (((a)>0)?1:((a)<0)?-1:0)
#define zeta_index(x, y, z) (z * ncols * nrows) + (y * ncols) + x
#define mzeta_index(x, y) (y * ncols) + x
#define xi_index(x, y, z) (z * ncols * ny) + (y * ncols) + x

typedef struct {
    int     x ;
    int     y ;     /* Coordinates of target pixel x,y  */
    double  w ;     /* Contribution weight <= 1/osample */
} xi_ref;

typedef struct {
    int     x ;
    int     iy ;    /* Contributing subpixel  x,iy      */
    double  w;      /* Contribution weight <= 1/osample */
} zeta_ref;

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_extract_slit_func_vert(
        int         ncols,
        int         nrows,
        int         osample,
        double  *   im,
        double  *   pix_unc,
        int     *   mask,
        double  *   ycen,
        double  *   sL,
        double  *   sP,
        double  *   unc,
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
        //double  *   PSF_curve,
        cpl_polynomial ** slitcurves,
        int         delta_x,
        double  *   sL,
        double  *   sP,
        double  *   model,
        double  *   unc,
        double      lambda_sP,
        double      lambda_sL,
        double      sP_stop,
        int         maxiter) ;

static int cr2res_extract_xi_zeta_tensors(
        int         ncols,
        int         nrows,
        int         ny,
        double  *   ycen,
        int     *   ycen_offset,
        int         y_lower_lim,
        int         osample,
        cpl_polynomial ** slitcurves,
        xi_ref   *  xi,
        zeta_ref *  zeta,
        int      *  m_zeta) ;

static int cr2res_extract_slitdec_bandsol(double *, double *, int, int) ;

static int cr2res_extract_slitdec_adjust_swath(
        int             sw,
        int             nx,
        cpl_vector  *   bins_begin,
        cpl_vector  *   bins_end) ;

/*----------------------------------------------------------------------------*/
/**
  @defgroup cr2res_extract  Extraction routines (Slit Decomposition,...)
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Extracts all the passed traces at once
  @param    img             Full detector image
  @param    traces          The traces table
  @param    reduce_order    The order to extract (-1 for all)
  @param    reduce_trace    The Trace to extract (-1 for all)
  @param    extr_method     The wished extraction method
  @param    extr_height     number of pix above and below mid-line or -1
  @param    swath_width     width per swath
  @param    oversample      factor for oversampling
  @param    smooth_slit     
  @param    extracted       [out] the extracted spectra 
  @param    slit_func       [out] the slit functions
  @param    model_master    [out] the model
  @return   0 if ok, -1 otherwise

  This func takes a single image (contining many orders), and a traces table.
 */
/*----------------------------------------------------------------------------*/
int cr2res_extract_traces(
        const hdrl_image    *   img,
        const cpl_table     *   traces,
        int                     reduce_order,
        int                     reduce_trace,
        cr2res_extr_method      extr_method,
        int                     extr_height,
        int                     swath_width,
        int                     oversample,
        double                  smooth_slit,
        cpl_table           **  extracted,
        cpl_table           **  slit_func,
        hdrl_image          **  model_master)
{
    cpl_bivector        **  spectrum ;
    cpl_vector          **  slit_func_vec ;
    cpl_table           *   slit_func_loc ;
    cpl_table           *   extract_loc ;
    hdrl_image          *   model_loc ;
    hdrl_image          *   model_loc_one ;
    int                     nb_traces, i, order, trace_id ;

    /* Check Entries */
    if (img == NULL || traces == NULL) return -1 ;

    /* Initialise */
    nb_traces = cpl_table_get_nrow(traces) ;

    /* Allocate Data containers */
    spectrum = cpl_malloc(nb_traces * sizeof(cpl_bivector *)) ;
    slit_func_vec = cpl_malloc(nb_traces * sizeof(cpl_vector *)) ;
    model_loc = hdrl_image_duplicate(img) ;
    hdrl_image_mul_scalar(model_loc, (hdrl_value){0.0, 0.0}) ;

    /* Loop over the traces and extract them */
    for (i=0 ; i<nb_traces ; i++) {
        /* Initialise */
        slit_func_vec[i] = NULL ;
        spectrum[i] = NULL ;
        model_loc_one = NULL ;

        /* Get Order and trace id */
        order = cpl_table_get(traces, CR2RES_COL_ORDER, i, NULL) ;
        trace_id = cpl_table_get(traces, CR2RES_COL_TRACENB, i, NULL) ;

        /* Check if this order needs to be skipped */
        if (reduce_order > -1 && order != reduce_order) continue ;

        /* Check if this trace needs to be skipped */
        if (reduce_trace > -1 && trace_id != reduce_trace) continue ;

        cpl_msg_info(__func__, "Process Order %d/Trace %d",order,trace_id) ;
        cpl_msg_indent_more() ;

        /* Call the Extraction */
        if (extr_method == CR2RES_EXTR_SUM) {
            /* Call the SUM ONLY extraction */
            if (cr2res_extract_sum_vert(img, traces, order,
                        trace_id, extr_height, &(slit_func_vec[i]),
                        &(spectrum[i]), &model_loc_one) != 0) {
                cpl_msg_error(__func__, "Cannot (sum-)extract the trace") ;
                slit_func_vec[i] = NULL ;
                spectrum[i] = NULL ;
                model_loc_one = NULL ;
                cpl_error_reset() ;
                cpl_msg_indent_less() ;
                continue ;
            }
        } else if (extr_method == CR2RES_EXTR_OPT_VERT) {
            /* Call the vertical SLIT DECOMPOSITION */
            if (cr2res_extract_slitdec_vert(img, traces,
                        order, trace_id, extr_height, swath_width,
                        oversample, smooth_slit, &(slit_func_vec[i]),
                        &(spectrum[i]), &model_loc_one) != 0) {
                cpl_msg_error(__func__,
                        "Cannot (slitdec-vert-) extract the trace") ;
                slit_func_vec[i] = NULL ;
                spectrum[i] = NULL ;
                model_loc_one = NULL ;
                cpl_error_reset() ;
                cpl_msg_indent_less() ;
                continue ;
            }
        } else if (extr_method == CR2RES_EXTR_OPT_CURV) {
            /* Call the curved SLIT DECOMPOSITION */
            if (cr2res_extract_slitdec_curved(img, traces,
                        order, trace_id, extr_height, swath_width,
                        oversample, smooth_slit, &(slit_func_vec[i]),
                        &(spectrum[i]), &model_loc_one) != 0) {
                cpl_msg_error(__func__,
                        "Cannot (slitdec-curved-) extract the trace") ;
                slit_func_vec[i] = NULL ;
                spectrum[i] = NULL ;
                model_loc_one = NULL ;
                cpl_error_reset() ;
                cpl_msg_indent_less() ;
                continue ;
            }
        }

        /* Update the model global image */
        if (model_loc_one != NULL) {
            hdrl_image_add_image(model_loc, model_loc_one) ;
            hdrl_image_delete(model_loc_one) ;
        }
        cpl_msg_indent_less() ;
    }

    /* Create the slit_func_tab for the current detector */
    if ((slit_func_loc = cr2res_extract_SLITFUNC_create(slit_func_vec,
                    traces)) == NULL) {
        cpl_msg_error(__func__, "Cannot compute the slit function") ;
        for (i=0 ; i<nb_traces ; i++) {
            if (slit_func_vec[i] != NULL) cpl_vector_delete(slit_func_vec[i]) ;
            if (spectrum[i] != NULL) cpl_bivector_delete(spectrum[i]) ;
        }
        cpl_free(spectrum) ;
        cpl_free(slit_func_vec) ;
        hdrl_image_delete(model_loc) ;
        return -1; 
    }

    /* Create the extracted_tab for the current detector */
    extract_loc = cr2res_extract_EXTRACT1D_create(spectrum, traces) ;

    /* Deallocate Vectors */
    for (i=0 ; i<nb_traces ; i++) {
        if (slit_func_vec[i] != NULL) cpl_vector_delete(slit_func_vec[i]) ;
        if (spectrum[i] != NULL) cpl_bivector_delete(spectrum[i]) ;
    }
    cpl_free(spectrum) ;
    cpl_free(slit_func_vec) ;

    /* Return  */
    *extracted = extract_loc ;
    *slit_func = slit_func_loc ;
    *model_master = model_loc;
    return 0 ;
}

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
        const hdrl_image    *   hdrl_in,
        const cpl_table     *   trace_tab,
        int                     order,
        int                     trace_id,
        int                     height,
        cpl_vector          **  slit_func,
        cpl_bivector        **  spec,
        hdrl_image          **  model)
{
    int             *   ycen_int;
    cpl_vector      *   ycen ;
    cpl_image       *   img_tmp;
    cpl_image       *   img_1d;
    const cpl_image       *   img_in;
    const cpl_image       *   err_in;
    cpl_vector      *   spc;
    cpl_vector      *   slitfu;
    cpl_vector      *   sigma;
    cpl_size            lenx, leny;
    int                 i, j;
    int                 ymin, ymax;
    int                 empty_bottom = 0;
    cpl_type            imtyp;

    /* Check Entries */
    if (hdrl_in == NULL || trace_tab == NULL) return -1 ;

    img_in = hdrl_image_get_image_const(hdrl_in);
    err_in = hdrl_image_get_error_const(hdrl_in);

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

    img_tmp = cpl_image_multiply_create(err_in, err_in);
    img_1d = cpl_image_collapse_create(img_tmp, 0);
    sigma = cpl_vector_new_from_image_row(img_1d, 1);
    cpl_vector_sqrt(sigma);
    cpl_image_delete(img_tmp);
    cpl_image_delete(img_1d);

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
    *spec = cpl_bivector_wrap_vectors(spc, sigma);
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
        cpl_bivector    **  spectrum,
        const cpl_table *   trace_table)
{
    cpl_table       *   out ;
    char            *   col_name ;
    const double    *   pspec ;
    const double    *   perr ;
    cpl_vector      *   wave_vec ;
    const double    *   pwl ;
    int                 nrows, all_null, i, order, trace_id, nb_traces ;

    /* Check entries */
    if (spectrum == NULL || trace_table == NULL) return NULL ;

    /* Initialise */
    nb_traces = cpl_table_get_nrow(trace_table) ;

    /* Check if all vectorѕ are not null */
    all_null = 1 ;
    for (i=0 ; i<nb_traces ; i++)
        if (spectrum[i] != NULL) {
            nrows = cpl_bivector_get_size(spectrum[i]) ;
            all_null = 0 ;
        }
    if (all_null == 1) return NULL ;

    /* Check the sizes */
    for (i=0 ; i<nb_traces ; i++)
        if (spectrum[i] != NULL && cpl_bivector_get_size(spectrum[i]) != nrows)
            return NULL ;

    /* Create the table */
    out = cpl_table_new(nrows);
    for (i=0 ; i<nb_traces ; i++) {
        order = cpl_table_get(trace_table, CR2RES_COL_ORDER, i, NULL) ;
        trace_id = cpl_table_get(trace_table, CR2RES_COL_TRACENB, i, NULL) ;
        /* Create SPEC column */
        col_name = cr2res_dfs_SPEC_colname(order, trace_id) ;
        cpl_table_new_column(out, col_name, CPL_TYPE_DOUBLE);
        cpl_free(col_name) ;
        /* Create SPEC_ERR column */
        col_name = cr2res_dfs_SPEC_ERR_colname(order, trace_id) ;
        cpl_table_new_column(out, col_name, CPL_TYPE_DOUBLE);
        cpl_free(col_name) ;
        /* Create WAVELENGTH column */
        col_name = cr2res_dfs_WAVELENGTH_colname(order, trace_id) ;
        cpl_table_new_column(out, col_name, CPL_TYPE_DOUBLE);
        cpl_free(col_name) ;
    }

    /* Fill the table */
    for (i=0 ; i<nb_traces ; i++) {
        if (spectrum[i] != NULL) {
            order = cpl_table_get(trace_table, CR2RES_COL_ORDER, i, NULL) ;
            trace_id = cpl_table_get(trace_table, CR2RES_COL_TRACENB, i, NULL) ;
            pspec = cpl_bivector_get_x_data_const(spectrum[i]) ;
            perr = cpl_bivector_get_y_data_const(spectrum[i]);
            /* Fill SPEC column */
            col_name = cr2res_dfs_SPEC_colname(order, trace_id) ;
            cpl_table_copy_data_double(out, col_name, pspec) ;
            cpl_free(col_name) ;
            /* Fill SPEC_ERR column */
            col_name = cr2res_dfs_SPEC_ERR_colname(order, trace_id) ;
            cpl_table_copy_data_double(out, col_name, perr) ;
            cpl_free(col_name) ;

            /* Compute Wavelength column */
            wave_vec = cr2res_trace_get_wl(trace_table, order, trace_id,
                    CR2RES_DETECTOR_SIZE);
            if (wave_vec == NULL) 
                wave_vec = cpl_vector_new(CR2RES_DETECTOR_SIZE) ;
            pwl = cpl_vector_get_data(wave_vec) ;

            /* Fill WAVELENGTH column */
            col_name = cr2res_dfs_WAVELENGTH_colname(order, trace_id) ;
            cpl_table_copy_data_double(out, col_name, pwl) ;
            cpl_free(col_name) ;
            cpl_vector_delete(wave_vec) ;
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
        const cpl_table *   trace_table)
{
    cpl_table       *   out ;
    char            *   col_name ;
    const double    *   pslit ;
    int                 nrows, nrows_max, all_null, i, order, trace_id, 
                        nb_traces ;

    /* Check entries */
    if (slit_func == NULL || trace_table == NULL) return NULL ;

    /* Initialise */
    nrows_max = -1 ;
    nb_traces = cpl_table_get_nrow(trace_table) ;

    /* Check that all vectorѕ are not null */
    all_null = 1 ;
    for (i=0 ; i<nb_traces ; i++)
        if (slit_func[i] != NULL) {
            nrows = cpl_vector_get_size(slit_func[i]) ;
            if (nrows > nrows_max) nrows_max = nrows ;
            all_null = 0 ;
        }
    if (all_null == 1) return NULL ;

    /* Create the table */
    out = cpl_table_new(nrows_max);
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
  @brief    Get a Spectrum and its error from the EXTRACT_1D table
  @param    tab         the EXTRACT_1D table
  @param    order       the order
  @param    trace_nb    the wished trace
  @param    spec        [out] The wavelength/spectrum bivector
  @param    spec_err    [out] The wavelength/spectrum error bivector
  @return   O if ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
int cr2res_extract_EXTRACT1D_get_spectrum(
        const cpl_table     *   tab,
        int                     order,
        int                     trace_nb,
        cpl_bivector        **  spec,
        cpl_bivector        **  spec_err)
{
    char            *   spec_name ;
    char            *   wave_name ;
    char            *   spec_err_name ;
    const double    *   pspec ;
    const double    *   pwave ;
    const double    *   pspec_err ;
    double          *   pxspec ;
    double          *   pyspec ;
    double          *   pxspec_err ;
    double          *   pyspec_err ;
    int                 i, tab_size ;

    /* Check entries */
    if (tab == NULL || spec == NULL || spec_err == NULL) return -1 ;

    /* Get the Spectrum */
    spec_name = cr2res_dfs_SPEC_colname(order, trace_nb) ;
    if ((pspec = cpl_table_get_data_double_const(tab, spec_name)) == NULL) {
        cpl_msg_error(__func__, "Cannot find the spectrum") ;
        cpl_free(spec_name) ;
        return -1 ;
    }
    cpl_free(spec_name) ;

    /* Get the Wavelength */
    wave_name = cr2res_dfs_WAVELENGTH_colname(order, trace_nb) ;
    if ((pwave = cpl_table_get_data_double_const(tab, wave_name)) == NULL) {
        cpl_msg_error(__func__, "Cannot find the wavelength") ;
        cpl_free(wave_name) ;
        return -1 ;
    }
    cpl_free(wave_name) ;

    /* Get the Spectrum Error */
    spec_err_name = cr2res_dfs_SPEC_ERR_colname(order, trace_nb) ;
    if ((pspec_err = cpl_table_get_data_double_const(tab, spec_err_name)) 
            == NULL) {
        cpl_msg_error(__func__, "Cannot find the spectrum error") ;
        cpl_free(spec_err_name) ;
        return -1 ;
    }
    cpl_free(spec_err_name) ;

    /* Create the output */
    tab_size = cpl_table_get_nrow(tab) ;
    *spec = cpl_bivector_new(tab_size) ;
    *spec_err = cpl_bivector_new(tab_size) ;
    pxspec = cpl_bivector_get_x_data(*spec) ;
    pyspec = cpl_bivector_get_y_data(*spec) ;
    pxspec_err = cpl_bivector_get_x_data(*spec_err) ;
    pyspec_err = cpl_bivector_get_y_data(*spec_err) ;
    for (i=0 ; i<tab_size ; i++) {
        pxspec[i] = pwave[i] ;
        pyspec[i] = pspec[i] ;
        pxspec_err[i] = pwave[i] ;
        pyspec_err[i] = pspec_err[i] ;
    }
    return 0 ;
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
        const hdrl_image    *   img_hdrl,
        const cpl_table     *   trace_tab,
        int                     order,
        int                     trace_id,
        int                     height,
        int                     swath,
        int                     oversample,
        double                  smooth_slit,
        cpl_vector          **  slit_func,
        cpl_bivector        **  spec,
        hdrl_image          **  model)
{
    cpl_polynomial  **  traces ;
    int             *   ycen_int;
    double          *   ycen_rest;
    double          *   ycen_sw;
    double          *   img_sw_data;
    double          *   err_sw_data;
    double          *   spec_sw_data;
    double          *   slitfu_sw_data;
    double          *   model_sw;
    int             *   mask_sw;
    double          *   unc_sw_data;
    const cpl_image *   img_in;
    const cpl_image *   err_in;
    cpl_image       *   err_sw;
    cpl_image       *   img_sw;
    cpl_image       *   img_rect;
    cpl_image       *   err_rect;
    cpl_image       *   model_rect;
    cpl_vector      *   ycen ;
    cpl_image       *   img_tmp;
    cpl_image       *   img_out;
    cpl_vector      *   spec_sw;
    cpl_vector      *   slitfu_sw;
    cpl_vector      *   spc;
    cpl_vector      *   slitfu;
    cpl_vector      *   unc_sw;
    cpl_vector      *   weights_sw;
    cpl_vector      *   tmp_vec;
    cpl_vector      *   bins_begin;
    cpl_vector      *   bins_end;
    cpl_vector      *   unc_decomposition;
    cpl_size            lenx, leny;
    cpl_type            imtyp;
    double              pixval, errval, img_median, unc, model_unc, img_unc, 
                        norm;
    int                 i, j, k, nswaths, row, col, x, y, ny_os,
                        sw_start, sw_end, badpix;

    /* Check Entries */
    if (img_hdrl == NULL || trace_tab == NULL) return -1 ;

    /* Initialise */
    img_in = hdrl_image_get_image_const(img_hdrl);
    err_in = hdrl_image_get_error_const(img_hdrl);
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
    if (height > leny) {
        height = leny;
        cpl_msg_warning(__func__,
                "Given height larger than image, clipping height");
    }

    /* Get ycen */
    if ((ycen = cr2res_trace_get_ycen(trace_tab, order,
                    trace_id, lenx)) == NULL) {
        cpl_msg_error(__func__, "Cannot get ycen");
        return -1 ;
    }

    if (oversample <= 0) oversample = 1;

    /* Number of rows after oversampling */
    bins_begin = cpl_vector_new(1);
    bins_end = cpl_vector_new(1);
    ny_os = oversample*(height+1) +1;
    if ((swath = cr2res_extract_slitdec_adjust_swath(swath, lenx, bins_begin,
                    bins_end)) == -1){
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
        cpl_image_save(img_rect, "debug_rectorder_vert.fits", imtyp,
                NULL, CPL_IO_CREATE);
    }
    err_rect = cr2res_image_cut_rectify(err_in, ycen, height);
    ycen_rest = cr2res_vector_get_rest(ycen);

    // Work vectors
    slitfu_sw = cpl_vector_new(ny_os);
    for (j=0; j < ny_os; j++) cpl_vector_set(slitfu_sw, j, 0);
    slitfu_sw_data = cpl_vector_get_data(slitfu_sw);

    /* Allocate */
    mask_sw = cpl_malloc(height*swath*sizeof(int));
    model_sw = cpl_malloc(height*swath*sizeof(double));
    ycen_sw = cpl_malloc(swath*sizeof(double));
    img_sw = cpl_image_new(swath, height, CPL_TYPE_DOUBLE);
    err_sw = cpl_image_new(swath, height, CPL_TYPE_DOUBLE);
    unc_sw = cpl_vector_new(swath);


    weights_sw = cpl_vector_new(swath);

    /* Pre-calculate the weights for overlapping swaths*/
    for (i=0; i < swath/2; i++) {
        cpl_vector_set(weights_sw, i, i + 1);
        cpl_vector_set(weights_sw, swath - i - 1, i+1);
    }
    // normalize such that max(w)=1
    cpl_vector_divide_scalar(weights_sw, i + 1);

    // Local versions of return data
    slitfu = cpl_vector_new(ny_os);
    spc = cpl_vector_new(lenx);
    unc_decomposition = cpl_vector_new(lenx);
    for (j=0; j<lenx ; j++){
        cpl_vector_set(spc, j, 0.);
        cpl_vector_set(unc_decomposition, j, 0.);
    }
    img_out = cpl_image_new(lenx, leny, CPL_TYPE_DOUBLE);
    model_rect = cpl_image_new(lenx, height, CPL_TYPE_DOUBLE);


    for (i=0;i<nswaths;i++){
        sw_start = cpl_vector_get(bins_begin, i);
        sw_end = cpl_vector_get(bins_end, i);


        // Copy swath image into seperate image
        for(col=1; col<=swath; col++){      // col is x-index in swath
            x = sw_start + col;          // coords in large image
            for(y=1;y<=height;y++){
                pixval = cpl_image_get(img_rect, x, y, &badpix);
                errval = cpl_image_get(err_rect, x, y, &badpix);
                cpl_image_set(img_sw, col, y, pixval);
                cpl_image_set(err_sw, col, y, errval);
                if(cpl_error_get_code() != CPL_ERROR_NONE)
                    cpl_msg_error(__func__, "%d %d %s",
                            x, y, cpl_error_get_where());
                // raw index for mask, start with 0!
                j = (y-1)*swath + (col-1) ;
                if (badpix == 0) mask_sw[j] = 1;
                else mask_sw[j] = 0;
            }
        }

        img_sw_data = cpl_image_get_data_double(img_sw);
        err_sw_data = cpl_image_get_data_double(err_sw);
        img_tmp = cpl_image_collapse_median_create(img_sw, 0, 0, 0);
        spec_sw = cpl_vector_new_from_image_row(img_tmp,1);
        unc_sw_data = cpl_vector_get_data(unc_sw);

        cpl_image_delete(img_tmp);
        spec_sw_data = cpl_vector_get_data(spec_sw);
        for (j=sw_start;j<sw_end;j++) ycen_sw[j-sw_start] = ycen_rest[j];


        /* Finally ready to call the slit-decomp */
        cr2res_extract_slit_func_vert(swath, height, oversample, img_sw_data, 
                err_sw_data, mask_sw, ycen_sw, slitfu_sw_data, spec_sw_data,
                model_sw, unc_sw_data, 0.0, smooth_slit, 1.0e-5, 20);

        for(col=1; col<=swath; col++){   // col is x-index in cut-out
            x = sw_start + col;          // coords in large image
            for(y=1;y<=height;y++){
                // raw index for mask, start with 0!
                j = (y-1)*swath + (col-1) ;
                cpl_image_set(model_rect, x, y, model_sw[j]);
                if(cpl_error_get_code() != CPL_ERROR_NONE)
                    cpl_msg_error(__func__, "%d %d %s", x, y,
                            cpl_error_get_where());
            }
        }

        // add up slit-functions, divide by nswaths below to get average
        if (i==0) cpl_vector_copy(slitfu, slitfu_sw);
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
            cpl_image_save(img_tmp, "debug_model_sw.fits", CPL_TYPE_DOUBLE,
                    NULL, CPL_IO_CREATE);
            cpl_image_unwrap(img_tmp);
            cpl_image_save(img_sw, "debug_img_sw.fits", CPL_TYPE_DOUBLE, NULL,
                    CPL_IO_CREATE);
        }


        // The last bins are shifted, overwriting the first k values
        // this is the same amount the bin was shifted to the front before
        // (when creating the bins)
        // The duplicate values in the vector will not matter as they are
        // not used below
        if (i != 0)
            k = cpl_vector_get(bins_end, i-1) -
                cpl_vector_get(bins_begin, i) - swath / 2;
        else k = 0;

        if (k != 0) {
            for (j= 0; j < swath - k; j++){
                cpl_vector_set(spec_sw, j, cpl_vector_get(spec_sw, j + k));
                cpl_vector_set(unc_sw, j, cpl_vector_get(unc_sw, j + k));
            }
            sw_start = cpl_vector_get(bins_begin, i-1) + swath / 2;
            cpl_vector_set(bins_begin, i, sw_start);
            // for the following k's
            cpl_vector_set(bins_end, i, sw_start + swath);
        }

        // first and last half swath are not weighted
        if (i==0){
            for (j = swath/2; j < swath; j++) {
                cpl_vector_set(spec_sw, j,
                    cpl_vector_get(spec_sw,j) * cpl_vector_get(weights_sw,j));
                cpl_vector_set(unc_sw, j,
                    cpl_vector_get(unc_sw,j) * cpl_vector_get(weights_sw,j));
            }
        }
        else if (i == nswaths - 1) {
            for (j = 0; j < swath / 2; j++) {
                cpl_vector_set(spec_sw, j,
                cpl_vector_get(spec_sw,j) * cpl_vector_get(weights_sw,j));

                cpl_vector_set(unc_sw, j,
                cpl_vector_get(unc_sw, j) * cpl_vector_get(weights_sw, j));
            }
        }
        else {
            /* Multiply by weights and add to output array */
            cpl_vector_multiply(spec_sw, weights_sw);
            cpl_vector_multiply(unc_sw, weights_sw);
        }

        // Save swath to output vector
        for (j=sw_start;j<sw_end;j++) {
            cpl_vector_set(spc, j,
                cpl_vector_get(spec_sw, j-sw_start) + cpl_vector_get(spc, j) );
            cpl_vector_set(unc_decomposition, j, 
                cpl_vector_get(unc_sw, j - sw_start) + cpl_vector_get(unc_decomposition, j));
        }
        cpl_vector_delete(spec_sw);

    } // End loop over swaths
    cpl_vector_delete(slitfu_sw);
    cpl_vector_delete(weights_sw);
    cpl_vector_delete(bins_begin);
    cpl_vector_delete(bins_end);
    cpl_vector_delete(unc_sw);

    cpl_free(mask_sw);
    cpl_free(model_sw);
    cpl_free(ycen_sw);
    cpl_image_delete(img_sw);
    cpl_image_delete(err_sw);

    // insert model_rect into large frame
    cr2res_image_insert_rect(model_rect, ycen, img_out);

    // divide by nswaths to make the slitfu into the average over all swaths.
    cpl_vector_divide_scalar(slitfu, nswaths);

    // TODO: Deallocate return arrays in case of error, return -1
    cpl_image_delete(img_rect);
    cpl_image_delete(model_rect);
    cpl_image_delete(err_rect);
    cpl_vector_delete(ycen);
    cpl_free(ycen_rest);

    *slit_func = slitfu;
    *spec = cpl_bivector_wrap_vectors(spc, unc_decomposition);
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
        const hdrl_image    *   img_hdrl,
        const cpl_table     *   trace_tab,
        int                     order,
        int                     trace_id,
        int                     height,
        int                     swath,
        int                     oversample,
        double                  smooth_slit,
        cpl_vector          **  slit_func,
        cpl_bivector        **  spec,
        hdrl_image          **  model)
{
    double          *   ycen_rest;
    double          *   ycen_sw;
    int             *   ycen_offset_sw;
    double          *   img_sw_data;
    double          *   err_sw_data;
    double          *   spec_sw_data;
    double          *   slitfu_sw_data;
    double          *   model_sw;
    double          *   unc_sw_data;
    int             *   mask_sw;
    const cpl_image *   img_in;
    const cpl_image *   err_in;
    cpl_image       *   img_sw;
    cpl_image       *   err_sw;
    cpl_image       *   img_rect;
    cpl_image       *   err_rect;
    cpl_image       *   model_rect;
    cpl_vector      *   ycen ;
    cpl_image       *   img_tmp;
    cpl_image       *   img_out;
    cpl_vector      *   spec_sw;
    cpl_vector      *   slitfu_sw;
    cpl_vector      *   unc_sw;
    cpl_vector      *   spc;
    cpl_vector      *   slitfu;
    cpl_vector      *   weights_sw;
    cpl_vector      *   tmp_vec;
    cpl_vector      *   bins_begin;
    cpl_vector      *   bins_end;
    cpl_vector      *   unc_decomposition;
    cpl_size            lenx, leny, pow;
    cpl_type            imtyp;
    cpl_polynomial      *slitcurve_A, *slitcurve_B, *slitcurve_C;
    cpl_polynomial  **  slitcurves_sw;

    double              pixval, errval, img_median, norm, model_unc, img_unc,
                        unc, delta_tmp, a, b, c, yc;
    int                 i, j, k, nswaths, halfswath, row, col, x, y, ny_os,
                        sw_start, sw_end, badpix, y_lower_limit, y_upper_limit,
                        delta_x;

    /* Check Entries */
    if (img_hdrl == NULL || trace_tab == NULL) return -1 ;

    img_in = hdrl_image_get_image_const(img_hdrl);
    err_in = hdrl_image_get_error_const(img_hdrl);

    /* Initialise */
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
    if (height > leny){
        height = leny;
        cpl_msg_warning(__func__,
                "Given height larger than image, clipping height");
    }

    /* Get ycen */
    if ((ycen = cr2res_trace_get_ycen(trace_tab, order,
                    trace_id, lenx)) == NULL) {
        cpl_msg_error(__func__, "Cannot get ycen");
        return -1 ;
    }

    /* Number of rows after oversampling */
    bins_begin = cpl_vector_new(1);
    bins_end = cpl_vector_new(1);
    ny_os = oversample*(height+1) +1;
    if ((swath = cr2res_extract_slitdec_adjust_swath(swath, lenx,
                    bins_begin, bins_end)) == -1){
        cpl_msg_error(__func__, "Cannot calculate swath size");
        cpl_vector_delete(ycen);
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
        cpl_image_save(img_rect, "debug_rectorder_curved.fits", imtyp,
                NULL, CPL_IO_CREATE);
    }
    err_rect = cr2res_image_cut_rectify(err_in, ycen, height);
    ycen_rest = cr2res_vector_get_rest(ycen);

    /* Retrieve the polynomials that describe the slit tilt and curvature*/
    slitcurve_A = cr2res_get_trace_wave_poly(trace_tab, CR2RES_COL_SLIT_CURV_A,
                    order, trace_id);
    slitcurve_B = cr2res_get_trace_wave_poly(trace_tab, CR2RES_COL_SLIT_CURV_B,
                    order, trace_id);
    slitcurve_C = cr2res_get_trace_wave_poly(trace_tab, CR2RES_COL_SLIT_CURV_C,
                    order, trace_id);
    if ((slitcurve_A == NULL) || (slitcurve_B == NULL) || (slitcurve_C == NULL)){
        cpl_msg_error(__func__, "No (or incomplete) slitcurve data found in trace table");
        cpl_vector_delete(ycen);
        cpl_vector_delete(bins_begin);
        cpl_vector_delete(bins_end);
        return -1;
    }

    /* Maximum horizontal shift in detector pixels due to slit image curv. */
    delta_x=0;
    for (i=1; i<=lenx; i+=swath/2){
        /* Do a coarse sweep through the order and evaluate the slitcurve */
        /* polynomials at  +- height/2, update the value. */
        /* Note: The index i is subtracted from a because the polys have */
        /* their origin at the edge of the full frame */
        a = cpl_polynomial_eval_1d(slitcurve_A, i, NULL);
        b =    cpl_polynomial_eval_1d(slitcurve_B, i, NULL);
        c =    cpl_polynomial_eval_1d(slitcurve_C, i, NULL);
        yc = cpl_vector_get(ycen, i-1);

        // shift polynomial to local frame
        a = a - i + yc * b + yc * yc * c; // a = 0
        b += 2 * yc * c;
        // cpl_msg_debug(__func__, "a: %f", a);
        // cpl_msg_debug(__func__, "b: %f", b);
        // cpl_msg_debug(__func__, "c: %f", c);
        // cpl_msg_debug(__func__, "h: %d", height);
        // cpl_msg_debug(__func__, "d: %f", b * height * 0.5);


        delta_tmp = max( fabs(a + (c*height/2. + b)*height/2.),
                fabs(a + (c*height/-2. + b)*height/-2.));
        if (delta_tmp > delta_x) delta_x = (int)ceil(delta_tmp);
    }
    cpl_msg_debug(__func__, "Max delta_x from slit curvature is %d pix.", delta_x);

    /* Allocate */
    mask_sw = cpl_malloc(height * swath*sizeof(int));
    model_sw = cpl_malloc(height * swath*sizeof(double));
    unc_sw = cpl_vector_new(swath);
    img_sw = cpl_image_new(swath, height, CPL_TYPE_DOUBLE);
    err_sw = cpl_image_new(swath, height, CPL_TYPE_DOUBLE);
    ycen_sw = cpl_malloc(swath*sizeof(double));
    ycen_offset_sw = cpl_malloc(swath * sizeof(int));

    slitcurves_sw = cpl_malloc(swath * sizeof(cpl_polynomial*));
    for (i=0; i<swath; i++) slitcurves_sw[i]= cpl_polynomial_new(1);

    // Local versions of return data
    slitfu = cpl_vector_new(ny_os);
    spc = cpl_vector_new(lenx);
    unc_decomposition = cpl_vector_new(lenx);
    for (j=0; j<lenx ; j++){
        cpl_vector_set(spc, j, 0.);
        cpl_vector_set(unc_decomposition, j, 0.);
    }
    img_out = cpl_image_new(lenx, leny, CPL_TYPE_DOUBLE);
    model_rect = cpl_image_new(lenx, height, CPL_TYPE_DOUBLE);

    // Work vectors
    slitfu_sw = cpl_vector_new(ny_os);
    for (j=0; j < ny_os; j++) cpl_vector_set(slitfu_sw, j, 0);
    slitfu_sw_data = cpl_vector_get_data(slitfu_sw);
    weights_sw = cpl_vector_new(swath);

    /* Pre-calculate the weights for overlapping swaths*/
    for (i=0; i < swath/2; i++) {
        cpl_vector_set(weights_sw, i, i + 1);
        cpl_vector_set(weights_sw, swath - i - 1, i+1);
    }
    // normalize such that max(w)=1
    cpl_vector_divide_scalar(weights_sw,i+1);

    // remove delta_x outermost weights and move them to the center 
    for (i=0; i < delta_x; i++){
        cpl_vector_set(weights_sw, i, 0);
        cpl_vector_set(weights_sw, swath - i - 1, 0);
        cpl_vector_set(weights_sw, swath / 2 + i, 1);
        cpl_vector_set(weights_sw, swath / 2 - i - 1, 1);
    }

    for (i=0;i<nswaths;i++){
        sw_start = cpl_vector_get(bins_begin, i);
        sw_end = cpl_vector_get(bins_end, i);

        /* Prepare swath cut-outs and auxiliary data */
        for(col=1; col<=swath; col++){   // col is x-index in swath
            x = sw_start + col;          // coords in large image

            /* prepare signal, error and mask */
            for(y=1;y<=height;y++){
                errval = cpl_image_get(err_rect, x, y, &badpix);
                pixval = cpl_image_get(img_rect, x, y, &badpix);
                if (cpl_error_get_code() != CPL_ERROR_NONE){
                    cpl_msg_error(__func__, "%d %d %s",
                            x, y, cpl_error_get_where());
                    cpl_error_reset();
                }
                cpl_image_set(img_sw, col, y, pixval);
                cpl_image_set(err_sw, col, y, errval);
                // raw index for mask, start with 0!
                j = (y-1)*swath + (col-1) ;
                if (badpix == 0) mask_sw[j] = 1;
                else mask_sw[j] = 0;
            }

            /* set slit curvature polynomials */
            /* subtract col because we want origin relative to here */
            pow = 2;
            cpl_polynomial_set_coeff(slitcurves_sw[col-1], &pow,
                cpl_polynomial_eval_1d(slitcurve_C, x, NULL));
            pow = 1;
            cpl_polynomial_set_coeff(slitcurves_sw[col-1], &pow,
                cpl_polynomial_eval_1d(slitcurve_B, x, NULL));
            pow = 0;
            cpl_polynomial_set_coeff(slitcurves_sw[col-1], &pow,
                cpl_polynomial_eval_1d(slitcurve_A, x, NULL) - x);
            /* and shift the polynomial to the central line */
            cpl_polynomial_shift_1d(slitcurves_sw[col-1], 0, cpl_vector_get(ycen, x-1));
        }

        img_sw_data = cpl_image_get_data_double(img_sw);
        err_sw_data = cpl_image_get_data_double(err_sw);
        unc_sw_data = cpl_vector_get_data(unc_sw);
        img_tmp = cpl_image_collapse_median_create(img_sw, 0, 0, 0);
        spec_sw = cpl_vector_new_from_image_row(img_tmp,1);
        cpl_image_delete(img_tmp);
        spec_sw_data = cpl_vector_get_data(spec_sw);

        for (j=sw_start;j<sw_end;j++){
            ycen_sw[j-sw_start] = ycen_rest[j];
        }
        y_lower_limit = (int) ycen_sw[0];
        for (j=0; j < swath; j++) {
            ycen_offset_sw[j] = (int) ycen_sw[j];
            y_lower_limit = min(y_lower_limit, (int) ycen_sw[j]);
        }
        y_lower_limit = height / 2 + y_lower_limit;

        if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
            img_tmp = cpl_image_wrap_int(swath, height, mask_sw);
            cpl_image_save(img_tmp, "debug_mask_before_sw.fits", CPL_TYPE_INT, NULL,
                    CPL_IO_CREATE);
            cpl_image_unwrap(img_tmp);
        }
        /* Finally ready to call the slit-decomp */
        cr2res_extract_slit_func_curved(swath, height, oversample, img_sw_data,
                err_sw_data, mask_sw, ycen_sw, ycen_offset_sw, y_lower_limit,
                slitcurves_sw, delta_x,
                slitfu_sw_data, spec_sw_data, model_sw, unc_sw_data, 0.,
                smooth_slit, 1e-7, 20);

        for (col=1; col<=swath; col++) {        // col is x-index in cut-out
            x = sw_start + col;                 // coords in large image
            for(y=1;y<=height;y++){
                j = (y-1)*swath + (col-1) ;     // raw index for mask,
                                                // start with 0!
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
            cpl_vector_save(weights_sw, "debug_weights.fits", CPL_TYPE_DOUBLE,
                    NULL, CPL_IO_CREATE);
            cpl_vector_save(slitfu_sw, "debug_slitfu.fits", CPL_TYPE_DOUBLE,
                    NULL, CPL_IO_CREATE);
            img_tmp = cpl_image_wrap_double(swath, height, model_sw);
            cpl_image_save(img_tmp, "debug_model_sw.fits", CPL_TYPE_DOUBLE,
                    NULL, CPL_IO_CREATE);
            cpl_image_unwrap(img_tmp);
            cpl_image_save(img_sw, "debug_img_sw.fits", CPL_TYPE_DOUBLE, NULL,
                    CPL_IO_CREATE);
            //cpl_image_unwrap(img_tmp);
        }

        // The last bins are shifted, overwriting the first k values
        // this is the same amount the bin was shifted to the front before
        // (when creating the bins)
        // The duplicate values in the vector will not matter as they are
        // not used below
        if (i != 0)
            k = cpl_vector_get(bins_end, i-1) -
                cpl_vector_get(bins_begin, i) - swath / 2;
        else k = 0;

        if (k != 0) {
            for (j= 0; j < swath - k; j++){
                cpl_vector_set(spec_sw, j, cpl_vector_get(spec_sw, j + k));
                cpl_vector_set(unc_sw, j, cpl_vector_get(unc_sw, j + k));
            }
            sw_start = cpl_vector_get(bins_begin, i-1) + swath / 2;
            cpl_vector_set(bins_begin, i, sw_start);
            // for the following k's
            cpl_vector_set(bins_end, i, sw_start + swath);
        }

        // first and last half swath are not weighted
        if (i==0){
            for (j = 0; j < delta_x; j++)
            {
                cpl_vector_set(spec_sw, j, 0);
                cpl_vector_set(unc_sw, j, 0.);
            }
            for (j = swath/2; j < swath; j++) {
                cpl_vector_set(spec_sw, j,
                    cpl_vector_get(spec_sw,j) * cpl_vector_get(weights_sw,j));
                
                cpl_vector_set(unc_sw, j,
                    cpl_vector_get(unc_sw, j) * cpl_vector_get(weights_sw, j));
            }
        } else if (i == nswaths - 1) {
            for (j = sw_end - sw_start - 1; j >= sw_end - sw_start - delta_x - 1; j--)
            {
                cpl_vector_set(spec_sw, j, 0);
                cpl_vector_set(unc_sw, j, 0);
            }
            for (j = 0; j < swath / 2; j++) {
                cpl_vector_set(spec_sw, j,
                cpl_vector_get(spec_sw,j) * cpl_vector_get(weights_sw,j));

                cpl_vector_set(unc_sw, j,
                cpl_vector_get(unc_sw,j) * cpl_vector_get(weights_sw,j));
            }
        } else {
            /* Multiply by weights and add to output array */
            cpl_vector_multiply(spec_sw, weights_sw);
            cpl_vector_multiply(unc_sw, weights_sw);
        }

        // Save swath to output vector
        for (j=sw_start;j<sw_end;j++) {
            cpl_vector_set(spc, j,
                cpl_vector_get(spec_sw, j-sw_start) + cpl_vector_get(spc, j));
            // just add weighted errors (instead of squared sum), as they are not independant
            cpl_vector_set(unc_decomposition, j,
                cpl_vector_get(unc_sw, j - sw_start) + cpl_vector_get(unc_decomposition, j));
        }

        cpl_vector_delete(spec_sw);
    } // End loop over swaths

    // insert model_rect into large frame
    cr2res_image_insert_rect(model_rect, ycen, img_out);

    // divide by nswaths to make the slitfu into the average over all swaths.
    cpl_vector_divide_scalar(slitfu, nswaths);

    // TODO: Update BPM in img_out


    // TODO: Deallocate return arrays in case of error, return -1
    cpl_image_delete(img_rect);
    cpl_image_delete(model_rect);
    cpl_image_delete(err_rect);
    cpl_image_delete(img_sw);
    cpl_image_delete(err_sw);

    cpl_free(mask_sw);
    cpl_free(model_sw);
    cpl_vector_delete(unc_sw);
    cpl_free(ycen_rest);
    cpl_free(ycen_sw);
    cpl_free(ycen_offset_sw);

    cpl_vector_delete(bins_begin);
    cpl_vector_delete(bins_end);
    cpl_vector_delete(slitfu_sw);
    cpl_vector_delete(weights_sw);
    cpl_vector_delete(ycen);

    cpl_polynomial_delete(slitcurve_A);
    cpl_polynomial_delete(slitcurve_B);
    cpl_polynomial_delete(slitcurve_C);
    for (i=0; i<swath; i++) cpl_polynomial_delete(slitcurves_sw[i]);
    cpl_free(slitcurves_sw);

    *slit_func = slitfu;
    *spec = cpl_bivector_wrap_vectors(spc, unc_decomposition);
    *model = hdrl_image_create(img_out, NULL);
    cpl_image_delete(img_out);
    return 0;
}


/*----------------------------------------------------------------------------*/
/*--------------------         EXTRACT 2d      -------------------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Extract2d all the passed traces at once
  @param    img             Full detector image
  @param    traces          The traces table
  @param    reduce_order    The order to extract (-1 for all)
  @param    reduce_trace    The Trace to extract (-1 for all)
  @param    smooth_slit     
  @param    extracted       [out] the extracted spectra 
  @return   0 if ok, -1 otherwise

  This func takes a single image (contining many orders), and a traces table.
 */
/*----------------------------------------------------------------------------*/
int cr2res_extract2d_traces(
        const hdrl_image    *   img,
        const cpl_table     *   traces,
        int                     reduce_order,
        int                     reduce_trace,
        cpl_table           **  extracted)
{
    cpl_bivector        **  spectrum ;
    cpl_bivector        **  position ;
    cpl_vector          **  wavelength ;
    cpl_vector          **  slit_fraction ;
    cpl_table           *   extract_loc ;
    cpl_image           *   wavemap;
    cpl_image           *   slitmap;
    int                     nb_traces, i, order, trace_id, npoints ;

    /* Check Entries */
    if (img == NULL || traces == NULL) return -1 ;

    /* Initialise */
    nb_traces = cpl_table_get_nrow(traces) ;
    npoints = CR2RES_DETECTOR_SIZE * CR2RES_DETECTOR_SIZE / nb_traces;

    /* Allocate Data containers */
    spectrum = cpl_malloc(nb_traces * sizeof(cpl_bivector *)) ;
    position = cpl_malloc(nb_traces * sizeof(cpl_bivector *)) ;
    wavelength = cpl_malloc(nb_traces * sizeof(cpl_vector *)) ;
    slit_fraction = cpl_malloc(nb_traces * sizeof(cpl_vector *)) ;

    // Calculate wavelength and slitfunction map once
    wavemap = cpl_image_new(CR2RES_DETECTOR_SIZE, CR2RES_DETECTOR_SIZE, CPL_TYPE_DOUBLE);
    slitmap = cpl_image_new(CR2RES_DETECTOR_SIZE, CR2RES_DETECTOR_SIZE, CPL_TYPE_DOUBLE);
    if (cr2res_slit_pos_image(traces, &slitmap, &wavemap) != 0)
    {
        cpl_msg_error(__func__, "Could not create wavelength / slit_fraction image");
        cpl_free(spectrum);
        cpl_free(position);
        cpl_free(wavelength);
        cpl_free(slit_fraction);
        return -1;
    }

    /* Loop over the traces and extract them */
    for (i=0 ; i<nb_traces ; i++) {
        /* Initialise */
        spectrum[i] = NULL ;
        position[i] = NULL ;
        wavelength[i] = NULL ;
        slit_fraction[i] = NULL ;

        /* Get Order and trace id */
        order = cpl_table_get(traces, CR2RES_COL_ORDER, i, NULL) ;
        trace_id = cpl_table_get(traces, CR2RES_COL_TRACENB, i, NULL) ;

        /* Check if this order needs to be skipped */
        if (reduce_order > -1 && order != reduce_order) continue ;

        /* Check if this trace needs to be skipped */
        if (reduce_trace > -1 && trace_id != reduce_trace) continue ;

        cpl_msg_info(__func__, "Process Order %d/Trace %d",order,trace_id) ;
        cpl_msg_indent_more() ;

        /* Call the Extraction */
        if (cr2res_extract2d_trace(img, traces, order, trace_id,
                    npoints, wavemap, slitmap, 
                    &(spectrum[i]), &(position[i]), &(wavelength[i]),
                    &(slit_fraction[i])) != 0) {
            cpl_msg_error(__func__, "Cannot extract2d the trace") ;
            spectrum[i] = NULL ;
            position[i] = NULL ;
            wavelength[i] = NULL ;
            slit_fraction[i] = NULL ;
            cpl_error_reset() ;
            cpl_msg_indent_less() ;
            continue ;
        }
        cpl_msg_indent_less() ;
    }

    /* Create the extracted_tab for the current detector */
    extract_loc = cr2res_extract_EXTRACT2D_create(spectrum, position,
            wavelength, slit_fraction, traces) ;

    /* Deallocate Vectors */
    for (i=0 ; i<nb_traces ; i++) {
        if (spectrum[i] != NULL) cpl_bivector_delete(spectrum[i]) ;
        if (position[i] != NULL) cpl_bivector_delete(position[i]) ;
        if (wavelength[i] != NULL) cpl_vector_delete(wavelength[i]) ;
        if (slit_fraction[i] != NULL) cpl_vector_delete(slit_fraction[i]) ;
    }
    cpl_free(spectrum) ;
    cpl_free(position) ;
    cpl_free(wavelength) ;
    cpl_free(slit_fraction) ;
    cpl_image_delete(wavemap);
    cpl_image_delete(slitmap);

    /* Return  */
    *extracted = extract_loc ;
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Extraction2d function
  @param    img_in          full detector image
  @param    trace_tab       The traces table
  @param    order           The order to extract
  @param    trace_id        The Trace to extract
  @param    wavemap         Map of the wavelength for each pixel
  @param    slitmap         Map of the slit fraction for each pixel
  @param    spectrum        [out] the spectrum and error
  @param    position        [out] the x/y positions
  @param    wavelength      [out] the wavelength values
  @param    slit_fraction   [out] the slit_fraction values
  @return   0 if ok, -1 otherwise

  Return the position, value, wavelength, and slitfraction of each pixel inside a trace

 */
/*----------------------------------------------------------------------------*/
int cr2res_extract2d_trace(
        const hdrl_image    *   in,
        const cpl_table     *   trace_tab,
        int                     order,
        int                     trace_id,
        int                     npoints,
        const cpl_image     *   wavemap,
        const cpl_image     *   slitmap,
        cpl_bivector        **  spectrum,
        cpl_bivector        **  position,
        cpl_vector          **  wavelength,
        cpl_vector          **  slit_fraction)
{
    cpl_bivector    *   spectrum_local ;
    cpl_vector      *   spectrum_flux ;
    cpl_vector      *   spectrum_error ;
    cpl_bivector    *   position_local ;
    cpl_vector      *   position_x;
    cpl_vector      *   position_y;
    cpl_vector      *   wavelength_local ;
    cpl_vector      *   slit_fraction_local ;
    const cpl_array *   lower_array;
    const cpl_array *   upper_array;
    cpl_polynomial  *   lower_poly;
    cpl_polynomial  *   upper_poly;
    cpl_vector      *   lower;
    cpl_vector      *   upper;
    cpl_vector      *   x;

    int k, bad_pix;
    cpl_size npix_max;
    cpl_size i, j, row;
    double flux, err;

    /* Check Entries */
    if (in==NULL || trace_tab==NULL || spectrum==NULL || position==NULL
            || wavelength==NULL || slit_fraction==NULL) return -1 ;

    if ((k = cr2res_get_trace_table_index(trace_tab, order, trace_id)) == -1)
    {
        cpl_msg_error(__func__, "Order and/or Trace not found in trace table");
        return -1;
    }

    // Step 0: Initialise output arrays
    npix_max = CR2RES_DETECTOR_SIZE * CR2RES_DETECTOR_SIZE;
    wavelength_local = cpl_vector_new(npoints);
    slit_fraction_local = cpl_vector_new(npoints);
    position_x = cpl_vector_new(npoints);
    position_y = cpl_vector_new(npoints);
    spectrum_flux = cpl_vector_new(npoints);
    spectrum_error = cpl_vector_new(npoints);

    // Step 1: Figure out pixels in the current trace
    // i.e. everything between upper and lower in trace_wave

    // The x coordinate of the detector
    x = cpl_vector_new(CR2RES_DETECTOR_SIZE);
    for (i = 0; i < CR2RES_DETECTOR_SIZE; i++) cpl_vector_set(x, i, i + 1);

    lower_array = cpl_table_get_array(trace_tab, CR2RES_COL_LOWER, k);
    upper_array = cpl_table_get_array(trace_tab, CR2RES_COL_UPPER, k);
    lower_poly = cr2res_convert_array_to_poly(lower_array);
    upper_poly = cr2res_convert_array_to_poly(upper_array);
    lower = cr2res_polynomial_eval_vector(lower_poly, x);
    upper = cr2res_polynomial_eval_vector(upper_poly, x);

    // Step 2: Iterate over pixels in the given trace
    // and fill the vectors
    row = -1;
    for (i = 0; i < CR2RES_DETECTOR_SIZE; i++)
    {
        for (j = cpl_vector_get(lower, i); j < cpl_vector_get(upper, i); j++)
        {
            row++;
            cpl_vector_set(position_x, row, i +1);
            cpl_vector_set(position_y, row, j);
            flux = cpl_image_get(hdrl_image_get_image_const(in), i + 1, j, &bad_pix);
            err = cpl_image_get(hdrl_image_get_error_const(in), i + 1, j, &bad_pix);
            cpl_vector_set(spectrum_flux, row, flux);
            cpl_vector_set(spectrum_error, row, err);

            // Set wavelength
            cpl_vector_set(wavelength_local, row, cpl_image_get(
                wavemap, i + 1, j, &bad_pix));
            // Set Slitfraction
            cpl_vector_set(slit_fraction_local, row, cpl_image_get(
                slitmap, i + 1, j, &bad_pix));
        }
    }

    // Step 3: resize output
    for (i = row; i < npoints; i++){
        cpl_vector_set(position_x, i, NAN);
        cpl_vector_set(position_y, i, NAN);
        cpl_vector_set(spectrum_flux, i, NAN);
        cpl_vector_set(spectrum_error, i, NAN);
        cpl_vector_set(wavelength_local, i, NAN);
        cpl_vector_set(slit_fraction_local, i, NAN);
    }

    position_local = cpl_bivector_wrap_vectors(position_x, position_y);
    spectrum_local = cpl_bivector_wrap_vectors(spectrum_flux, spectrum_error);

    cpl_polynomial_delete(upper_poly);
    cpl_polynomial_delete(lower_poly);
    cpl_vector_delete(upper);
    cpl_vector_delete(lower);
    cpl_vector_delete(x);

    *spectrum = spectrum_local ;
    *position = position_local ;
    *wavelength = wavelength_local ;
    *slit_fraction = slit_fraction_local ;
    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Create the extract 2D table to be saved
  @param    spectrum        value and error of the spectrum
  @param    position        x and y positions
  @param    wavelength      WL value
  @param    slit_fraction   slit fraction
  @param    traces          Trace wave file used for extraction
  @return   the extract_2D table or NULL
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_extract_EXTRACT2D_create(
        cpl_bivector    **  spectrum,
        cpl_bivector    **  position,
        cpl_vector      **  wavelength,
        cpl_vector      **  slit_fraction,
        const cpl_table *   trace_table)
{
    cpl_table       *   out ;
    char            *   col_name ;
    const double    *   pspec ;
    const double    *   perr ;
    const double    *   pxposition ;
    const double    *   pyposition ;
    const double    *   pwave ;
    const double    *   pslit_frac ;
    int                 nrows, all_null, i, order, trace_id, nb_traces ;

    /* Check entries */
    if (spectrum==NULL || trace_table==NULL || position==NULL ||
            wavelength==NULL || slit_fraction==NULL || trace_table==NULL) 
        return NULL ;

    /* Initialise */
    nb_traces = cpl_table_get_nrow(trace_table) ;

    /* Check if all vectorѕ are not null */
    all_null = 1 ;
    for (i=0 ; i<nb_traces ; i++)
        if (spectrum[i] != NULL) {
            nrows = cpl_bivector_get_size(spectrum[i]) ;
            all_null = 0 ;
        }
    if (all_null == 1) return NULL ;

    /* Check the sizes */
    for (i=0 ; i<nb_traces ; i++)
        if (spectrum[i] != NULL && cpl_bivector_get_size(spectrum[i]) != nrows)
            return NULL ;

    /* Create the table */
    out = cpl_table_new(nrows);
    for (i=0 ; i<nb_traces ; i++) {
        order = cpl_table_get(trace_table, CR2RES_COL_ORDER, i, NULL) ;
        trace_id = cpl_table_get(trace_table, CR2RES_COL_TRACENB, i, NULL) ;
        /* Create SPEC column */
        col_name = cr2res_dfs_SPEC_colname(order, trace_id) ;
        cpl_table_new_column(out, col_name, CPL_TYPE_DOUBLE);
        cpl_free(col_name) ;
        /* Create SPEC_ERR column */
        col_name = cr2res_dfs_SPEC_ERR_colname(order, trace_id) ;
        cpl_table_new_column(out, col_name, CPL_TYPE_DOUBLE);
        cpl_free(col_name) ;
        /* Create WAVELENGTH column */
        col_name = cr2res_dfs_WAVELENGTH_colname(order, trace_id) ;
        cpl_table_new_column(out, col_name, CPL_TYPE_DOUBLE);
        cpl_free(col_name) ;
        /* Create POSITIONX column */
        col_name = cr2res_dfs_POSITIONX_colname(order, trace_id) ;
        cpl_table_new_column(out, col_name, CPL_TYPE_DOUBLE);
        cpl_free(col_name) ;
        /* Create POSITIONY column */
        col_name = cr2res_dfs_POSITIONY_colname(order, trace_id) ;
        cpl_table_new_column(out, col_name, CPL_TYPE_DOUBLE);
        cpl_free(col_name) ;
        /* Create SLIT_FRACTION column */
        col_name = cr2res_dfs_SLIT_FRACTION_colname(order, trace_id) ;
        cpl_table_new_column(out, col_name, CPL_TYPE_DOUBLE);
        cpl_free(col_name) ;
    }

    /* Fill the table */
    for (i=0 ; i<nb_traces ; i++) {
        if (spectrum[i]!=NULL && position[i]!=NULL &&
                wavelength[i]!=NULL && slit_fraction[i]!=NULL) {
            order = cpl_table_get(trace_table, CR2RES_COL_ORDER, i, NULL) ;
            trace_id = cpl_table_get(trace_table, CR2RES_COL_TRACENB, i, NULL) ;
            pspec = cpl_bivector_get_x_data_const(spectrum[i]) ;
            perr = cpl_bivector_get_y_data_const(spectrum[i]);
            pxposition = cpl_bivector_get_x_data_const(position[i]) ;
            pyposition = cpl_bivector_get_y_data_const(position[i]) ;
            pwave = cpl_vector_get_data_const(wavelength[i]) ;
            pslit_frac = cpl_vector_get_data_const(slit_fraction[i]) ;
            /* Fill SPEC column */
            col_name = cr2res_dfs_SPEC_colname(order, trace_id) ;
            cpl_table_copy_data_double(out, col_name, pspec) ;
            cpl_free(col_name) ;
            /* Fill SPEC_ERR column */
            col_name = cr2res_dfs_SPEC_ERR_colname(order, trace_id) ;
            cpl_table_copy_data_double(out, col_name, perr) ;
            cpl_free(col_name) ;
            /* Fill WAVELENGTH column */
            col_name = cr2res_dfs_WAVELENGTH_colname(order, trace_id) ;
            cpl_table_copy_data_double(out, col_name, pwave) ;
            cpl_free(col_name) ;
            /* Fill POSITIONX column */
            col_name = cr2res_dfs_POSITIONX_colname(order, trace_id) ;
            cpl_table_copy_data_double(out, col_name, pxposition) ;
            cpl_free(col_name) ;
            /* Fill POSITIONY column */
            col_name = cr2res_dfs_POSITIONY_colname(order, trace_id) ;
            cpl_table_copy_data_double(out, col_name, pyposition) ;
            cpl_free(col_name) ;
            /* Fill SLIT_FRACTION column */
            col_name = cr2res_dfs_SLIT_FRACTION_colname(order, trace_id) ;
            cpl_table_copy_data_double(out, col_name, pslit_frac) ;
            cpl_free(col_name) ;
        }
    }
    return out ;
}

/** @} */

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
        double  *   pix_unc,
        int     *   mask,
        double  *   ycen,
        double  *   sL,
        double  *   sP,
        double  *   model,
        double  *   unc,
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
    double * Aij = cpl_malloc(ny*(2 * osample + 1)*sizeof(double)); // double Aij[ny*ny];
    double * bj = cpl_malloc(ny*sizeof(double)); // double bj[ny];
    // double Adiag[ncols*3];
    double * Adiag = cpl_malloc(ncols*3*sizeof(double));
    // double omega[ny][nrows][ncols];
    double * omega = cpl_malloc(ny*nrows*ncols*sizeof(double));
    // index as: [iy+(y*ny)+(x*ny*nrows)]
    double *p_bj   = cpl_malloc(ncols * sizeof(double));

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
        iy2 = osample - floor(ycen[x] / step) - 1;
        iy1 = iy2 - osample;

        d1 = fmod(ycen[x], step);
        if (d1 == 0)
            d1 = step;
        d2 = step - d1;

        for(y=0; y<nrows; y++) {
            iy1+=osample;
            iy2+=osample;
            for(iy=0; iy<ny; iy++) {
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
    do {
        /* Compute slit function sL */

        /* Fill in SLE arrays */
        diag_tot=0.e0;
        for(iy=0; iy<ny; iy++) {
            bj[iy]=0.e0;
            for(jy=max(iy-osample,0); jy<=min(iy+osample,ny-1); jy++) {
                Aij[iy+ny*(jy-iy+osample)]=0.e0;
                for(x=0; x<ncols; x++) {
                    sum=0.e0;
                    for(y=0; y<nrows; y++)
                        sum+=omega[iy+(y*ny)+(x*ny*nrows)] *
                            omega[jy+(y*ny)+(x*ny*nrows)]*mask[y*ncols+x];
                    Aij[iy+ny*(jy-iy+osample)]+=sum*sP[x]*sP[x];
                }
            }
            for(x=0; x<ncols; x++) {
                sum=0.e0;
                for(y=0; y<nrows; y++)
                    sum+=omega[iy+(y*ny)+(x*ny*nrows)]*
                        mask[y*ncols+x]*im[y*ncols+x];
                bj[iy]+=sum*sP[x];
            }
            diag_tot+=Aij[iy+ny*osample];
        }

        /* Scale regularization parameters */
        lambda=lambda_sL*diag_tot/ny;

        /* Add regularization parts for the slit function */
        Aij[ny*osample]    +=lambda;           /* Main diagonal  */
        Aij[ny*(osample+1)]-=lambda;           /* Upper diagonal */
        for(iy=1; iy<ny-1; iy++) {
            Aij[iy+ny*(osample-1)]-=lambda;      /* Lower diagonal */
            Aij[iy+ny*osample    ]+=lambda*2.e0; /* Main diagonal  */
            Aij[iy+ny*(osample+1)]-=lambda;      /* Upper diagonal */
        }
        Aij[ny-1+ny*(osample-1)]-=lambda;      /* Lower diagonal */
        Aij[ny-1+ny*osample]    +=lambda;      /* Main diagonal  */

        /* Solve the system of equations */
        info=cr2res_extract_slitdec_bandsol(Aij, bj, ny, nd);
        if(info) cpl_msg_warning(__func__, "Bandsol exited with %d", info);

        /* Normalize the slit function */
        norm=0.e0;
        for(iy=0; iy<ny; iy++) {
            sL[iy]=bj[iy];
            norm+=sL[iy];
        }
        norm/=osample;
        for(iy=0; iy<ny; iy++) sL[iy]/=norm;

        /* Compute spectrum sP */
        for(x=0; x<ncols; x++) {
            Adiag[x]=0.e0;
            Adiag[x+ncols]=0.e0;
            Adiag[x+2 *ncols]=0.e0;

            E[x]=0.e0;
            for(y=0; y<nrows; y++) {
                sum=0.e0;
                for(iy=0; iy<ny; iy++) {
                    sum+=omega[iy+(y*ny)+(x*ny*nrows)]*sL[iy];
                }

                Adiag[x+ncols]+=sum*sum*mask[y*ncols+x];
                E[x]+=sum*im[y*ncols+x]*mask[y*ncols+x];
            }
        }
        if(lambda_sP>0.e0) {
            norm=0.e0;
            for(x=0; x<ncols; x++) {
                sP_old[x]=sP[x];
                norm+=sP[x];
            }
            norm/=ncols;
            lambda=lambda_sP*norm;
            Adiag[0        ]  = 0.e0;
            Adiag[0+ncols  ] += lambda;
            Adiag[0+ncols*2] -= lambda;
            for(x=1; x<ncols-1; x++) {
                Adiag[x]=-lambda;
                Adiag[x+ncols  ]+= 2.e0*lambda;
                Adiag[x+ncols*2] =-lambda;
            }
            Adiag[ncols - 1            ] -= lambda;
            Adiag[ncols - 1 +     ncols] += lambda;
            Adiag[ncols - 1 + 2 * ncols] = 0.e0;

            info=cr2res_extract_slitdec_bandsol(Adiag, E, ncols, 3);
            for(x=0; x<ncols; x++) sP[x]=E[x];
        } else {
            for(x=0; x<ncols; x++) {
                sP_old[x]=sP[x];
                sP[x]=E[x]/Adiag[x+ncols];
            }
        }

        /* Compute the model */
        for(y=0; y<nrows; y++) {
            for(x=0; x<ncols; x++) {
                sum=0.e0;
                for(iy=0; iy<ny; iy++)
                    sum+=omega[iy+(y*ny)+(x*ny*nrows)]*sL[iy];
                model[y*ncols+x]=sum*sP[x];
            }
        }

        /* Compare model and data */
        sum=0.e0;
        isum=0;
        for(y=0; y<nrows; y++) {
            for(x=0;x<ncols; x++) {
                sum+=mask[y*ncols+x]*(model[y*ncols+x]-im[y*ncols+x]) *
                                (model[y*ncols+x]-im[y*ncols+x]);
                isum+=mask[y*ncols+x];
            }
        }
        dev=sqrt(sum/isum);

        /* Adjust the mask marking outlyers */
        for(y=0; y<nrows; y++) {
            for(x=0;x<ncols; x++) {
                if(fabs(model[y*ncols+x]-im[y*ncols+x])>6.*dev)
                    mask[y*ncols+x]=0;
                else mask[y*ncols+x]=1;
            }
        }

        /* Compute the change in the spectrum */
        sP_change=0.e0;
        sP_max=1.e0;
        for(x=0; x<ncols; x++) {
            if(sP[x]>sP_max) sP_max=sP[x];
            if(fabs(sP[x]-sP_old[x])>sP_change) sP_change=fabs(sP[x]-sP_old[x]);
        }
        /* Check the convergence */
    } while(iter++ < maxiter && sP_change > sP_stop*sP_max);

    /* Uncertainty estimate */
    for (x = 0; x < ncols; x++) {
        unc[x] = 0.;
        p_bj[x] = 0.;
    }
    for (x = 0; x < ncols; x++) {
        // Loop through all pixels contributing to x
        for (y = 0; y < nrows; y++) {
            unc[x] += (im[y * ncols + x] - model[y * ncols + x]) *
                (im[y * ncols + x] - model[y * ncols + x]) *
                sL[y] * mask[y * ncols + x];
            unc[x] += pix_unc[y * ncols + x] * pix_unc[y * ncols + x] *
                sL[y] * mask[y * ncols + x];
            // Norm
            p_bj[x] += sL[y] * mask[y * ncols + x];
        }
        
    }
    for (x = 0; x < ncols; x++) {
        unc[x] = sqrt(unc[x] / p_bj[x] * nrows);
    }


    cpl_free(E);
    cpl_free(sP_old);
    cpl_free(omega);
    cpl_free(Aij);
    cpl_free(bj);
    cpl_free(Adiag);
    cpl_free(p_bj);

    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Helper function for cr2res_extract_slit_func_curved
  @param ncols          Swath width in pixels
  @param nrows          Extraction slit height in pixels
  @param ny             Size of the slit function array: ny=osample(nrows+1)+1
  @param ycen           Order centre line offset from pixel row boundary
  @param ycen_offset    Order image column shift
  @param y_lower_lim    Number of detector pixels below the pixel
                        containing the central line yc
  @param osample        Subpixel ovsersampling factor
  @param PSF_curve      Parabolic fit to the slit image curvature
                        For column d_x = PSF_curve[ncols][0] +
                                        PSF_curve[ncols][1] *d_y +
                                        PSF_curve[ncols][2] *d_y^2
                        where d_y is the offset from the central line ycen.
                        Thus central subpixel of omega[x][y'][delta_x][iy']
                        does not stick out of column x
  @param xi[ncols][ny][4]   Convolution tensor telling the coordinates
                            of detector pixels on which {x, iy} element
                            falls and the corresponding projections
  @param zeta[ncols][nrows][3 * (osample + 1)]
                        Convolution tensor telling the coordinates
                        of subpixels {x, iy} contributing to detector pixel
                        {x, y}
  @param m_zeta[ncols][nrows]
                        The actual number of contributing elements in zeta

  @return
 */
/*----------------------------------------------------------------------------*/
static int cr2res_extract_xi_zeta_tensors(
        int         ncols,
        int         nrows,
        int         ny,
        double  *   ycen,
        int     *   ycen_offset,
        int         y_lower_lim,
        int         osample,
        cpl_polynomial ** slitcurves,
        xi_ref   *  xi,
        zeta_ref *  zeta,
        int      *  m_zeta)
{
    int     x, xx, y, yy, ix, ix1, ix2, iy, iy1, iy2;
    double  step, delta, dy, w, d1, d2, offset;
    const cpl_size power = 1;

    /* Initialise */
    step = 1.e0 / osample;

    /* Clean xi   */
    for (x = 0; x < ncols; x++) {
        for (iy = 0; iy < ny; iy++) {
            xi[xi_index(x,iy,0)].x=xi[xi_index(x,iy,1)].x = xi[xi_index(x,iy,2)].x = xi[xi_index(x,iy,3)].x=0;
            xi[xi_index(x,iy,0)].y=xi[xi_index(x,iy,1)].y = xi[xi_index(x,iy,2)].y = xi[xi_index(x,iy,3)].y=0;
            xi[xi_index(x,iy,0)].w=xi[xi_index(x,iy,1)].w = xi[xi_index(x,iy,2)].w = xi[xi_index(x,iy,3)].w=0.;
        }
    }

    /* Clean zeta */
    for (x = 0; x < ncols; x++) {
        for (y = 0; y < nrows; y++) {
            m_zeta[mzeta_index(x,y)] = 0;
            for (ix = 0; ix < 3 * (osample + 1); ix++) {
                zeta[zeta_index(x,y,ix)].x = 0;
                zeta[zeta_index(x,y,ix)].iy = 0;
                zeta[zeta_index(x,y,ix)].w = 0.;
            }
        }
    }

    /*
    printf("%g %g %g; %g %g %g; %g %g %g\n",
            PSF_curve[313][0],PSF_curve[313][1],PSF_curve[313][2],
            PSF_curve[314][0],PSF_curve[314][1],PSF_curve[314][2],
            PSF_curve[315][0],PSF_curve[315][1],PSF_curve[315][2]);
    */

    /*
       Construct the xi and zeta tensors. They contain pixel references and
       contribution.
       values going from a given subpixel to other pixels (xi) and coming from
       other subpixels to a given detector pixel (zeta).
       Note, that xi and zeta are used in the equations for sL, sP and for
       the model but they do not involve the data, only the geometry.
       Thus it can be pre-computed once.
     */
    for (x = 0; x < ncols; x++) {
        /*
           I promised to reconsider the initial offset. Here it is. For the
           original layout (no column shifts and discontinuities in ycen)
           there is pixel y that contains the central line yc. There are two
           options here (by construction of ycen that can be 0 but cannot be
           1): (1) yc is inside pixel y and (2) yc falls at the boundary
           between pixels y and y-1. yc cannot be at the foundary of pixels
           y+1 and y because we would select y+1 to be pixel y in that case.

           Next we need to define starting and ending indices iy for sL
           subpixels that contribute to pixel y. I call them iy1 and iy2.
           For both cases we assume osample+1 subpixels covering pixel y
           (weird). So for case 1 iy1 will be (y-1)*osample and
           iy2 == y*osample. Special treatment of the boundary subpixels
           will compensate for introducing extra subpixel in case 1.
           In case 2 things are more logical: iy1=(yc-y)*osample+(y-1)*osample;
           iy2=(y+1-yc)*osample)+(y-1)*osample. ycen is yc-y making things
           simpler. Note also that the same pattern repeates for all rows:
           we only need to initialize iy1 and iy2 and keep incrementing them
           by osample.
         */
        iy2 = osample - floor(ycen[x] / step) - 1;
        iy1 = iy2 - osample;

        /*
           Handling partial subpixels cut by detector pixel rows is again
           tricky.
           Here we have three cases (mostly because of the decision to assume
           that we always have osample+1 subpixels per one detector pixel).
           Here d1 is the fraction of the subpixel iy1 inside detector pixel y.
           d2 is then the fraction of subpixel iy2 inside detector pixel y.
           By definition d1+d2==step.
           Case 1: ycen falls on the top boundary of each detector pixel
           (ycen == 1). Here we conclude that the first subpixel is fully
           contained inside pixel y and d1 is set to step.
           Case 2: ycen falls on the bottom boundary of each detector pixel
           (ycen == 0). Here we conclude that the first subpixel is totally
           outside of pixel y and d1 is set to 0.
           Case 3: ycen falls inside of each pixel (0>ycen>1). In this case d1
           is set to the fraction of the first step contained inside of each
           pixel.  And BTW, this also means that central line coinsides with
           the upper boundary of subpixel iy2 when the y loop reaches pixel
           y_lower_lim. In other words:
           dy=(iy-(y_lower_lim+ycen[x])*osample)*step-0.5*step
         */

        d1 = fmod(ycen[x], step);
        if (d1 == 0)
            d1 = step;
        d2 = step - d1;

        /*
           The final hurdle for 2D slit decomposition is to construct two
           3D reference tensors. We proceed similar to 1D case except that
           now each iy subpixel can be shifted left or right following the
           curvature of the slit image on the detector. We assume for now that
           each subpixel is exactly 1 detector pixel wide. This may not be
           exactly true if the curvature changes accross the focal plane but
           will deal with it when the necessity will become apparent. For now
           we just assume that a shift delta the weight w assigned to subpixel
           iy is divided between ix1=int(delta) and
           ix2=int(delta)+signum(delta) as (1-|delta-ix1|)*w and |delta-ix1|*w.
           The curvature is given by a quadratic polynomial evaluated from
           an approximation for column x:
           delta = PSF_curve[x][0] + PSF_curve[x][1] * (y-yc[x]) +
               PSF_curve[x][2] * (y-yc[x])^2.
           It looks easy except that y and yc are set in the global detector
           coordinate system rather than in the shifted and cropped swath
           passed to slit_func_2d. One possible solution I will try here is
           to modify PSF_curve before the call such as:
           delta = PSF_curve'[x][0] + PSF_curve'[x][1] * (y'-ycen[x]) +
               PSF_curve'[x][2] * (y'-ycen[x])^2
            where y' = y - floor(yc).
         */

        //dy=-ceil((ycen[x]+y_lower_lim)*osample)*step-step*0.5;
        //printf("dy=%g,", dy);

        dy = -(y_lower_lim * osample + floor(ycen[x] / step) + 0.5) * step;

        /* Define initial distance from ycen       */
        /* ie the center of the first subpixel falling into pixel y_lower_lim */

        /* printf("x=%d, dy=%g, step=%g, iy1=%d, iy2=%d, ycen=%g\n",  */
                /* x, dy, step, iy1, iy2, ycen[x]); */

        /*
           Now we go detector pixels x and y incrementing subpixels looking
           for their controibutions to the current and adjacent pixels.
           Note that the curvature/tilt of the projected slit image could be
           so large that subpixel iy may no contribute to column x at all.
           On the other hand, subpixels around ycen by definition must
           contribute to pixel x,y.
           3rd index in xi refers corners of pixel xx,y: 0:LL, 1:LR, 2:UL, 3:UR.
         */
        for (y = 0; y < nrows; y++) {
            iy1 += osample; // Bottom subpixel falling in row y
            iy2 += osample; // Top subpixel falling in row y
            dy -= step;
            for (iy = iy1; iy <= iy2; iy++) {
                if (iy == iy1)      w = d1;
                else if (iy == iy2) w = d2;
                else                w = step;
                dy += step;
                //delta = (PSF_curve[1 + 3 * x] + PSF_curve[2 + 3 * x] * dy) * dy;
                //offset = - (0.5 - step) * cpl_polynomial_get_coeff(slitcurves[x], &power);
                delta = cpl_polynomial_eval_1d(slitcurves[x], dy, NULL);
                ix1 = delta;
                ix2 = ix1 + signum(delta);

                /* Three cases: bottom boundary of row y, intermediate i
                   subpixels and top boundary */

                if (iy == iy1) {
                    /* Subpixel iy is entering detector row y        */
                    if (ix1 < ix2) {
                        /* Subpixel iy shifts to the right from column x */
                        if (x + ix1 >= 0 && x + ix2 < ncols) {
                            xx = x + ix1;
                            yy = y + ycen_offset[x] - ycen_offset[xx];
                            xi[xi_index(x,iy,1)].x = xx;
                            xi[xi_index(x,iy,1)].y = yy;
                            xi[xi_index(x,iy,1)].w = w - fabs(delta - ix1) * w;
                            if (xx >= 0 && xx < ncols && yy >= 0 &&
                                    yy < nrows && xi[xi_index(x,iy,1)].w > 0) {
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].x = x;
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].iy = iy;
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].w = xi[xi_index(x,iy,1)].w;
                                m_zeta[mzeta_index(xx,yy)]++;
                            }
                            xx = x + ix2;
                            yy = y + ycen_offset[x] - ycen_offset[xx];
                            xi[xi_index(x,iy,0)].x = xx;
                            xi[xi_index(x,iy,0)].y = yy;
                            xi[xi_index(x,iy,0)].w = fabs(delta - ix1) * w;
                            if (xx >= 0 && xx < ncols && yy >= 0 &&
                                    yy < nrows && xi[xi_index(x,iy,0)].w > 0) {
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].x = x;
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].iy = iy;
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].w = xi[xi_index(x,iy,0)].w;
                                m_zeta[mzeta_index(xx,yy)]++;
                            }
                        }
                    } else if (ix1 > ix2) {
                        /* Subpixel iy shifts to the left from column x */
                        if (x + ix2 >= 0 && x + ix1 < ncols) {
                            xx = x + ix2;
                            yy = y + ycen_offset[x] - ycen_offset[xx];
                            xi[xi_index(x,iy,1)].x = xx;
                            xi[xi_index(x,iy,1)].y = yy;
                            xi[xi_index(x,iy,1)].w = fabs(delta - ix1) * w;
                            if (xx >= 0 && xx < ncols && yy >= 0 &&
                                    yy < nrows && xi[xi_index(x,iy,1)].w > 0) {
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].x = x;
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].iy = iy;
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].w = xi[xi_index(x,iy,1)].w;
                                m_zeta[mzeta_index(xx,yy)]++;
                            }
                            xx = x + ix1;
                            yy = y + ycen_offset[x] - ycen_offset[xx];
                            xi[xi_index(x,iy,0)].x = xx;
                            xi[xi_index(x,iy,0)].y = yy;
                            xi[xi_index(x,iy,0)].w = w - fabs(delta - ix1) * w;

                            if (xx >= 0 && xx < ncols && yy >= 0 &&
                                    yy < nrows && xi[xi_index(x,iy,0)].w > 0) {
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].x = x;
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].iy = iy;
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].w = xi[xi_index(x,iy,0)].w;
                                m_zeta[mzeta_index(xx,yy)]++;
                            }
                        }
                    } else {
                        /* Subpixel iy stays inside column x */
                        xx = x + ix1;
                        yy = y + ycen_offset[x] - ycen_offset[xx];
                        xi[xi_index(x,iy,0)].x = xx;
                        xi[xi_index(x,iy,0)].y = yy;
                        xi[xi_index(x,iy,0)].w = w;
                        if (xx >= 0 && xx < ncols && yy >= 0 && yy < nrows &&
                                w > 0) {
                            zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].x = x;
                            zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].iy = iy;
                            zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].w = w;
                            m_zeta[mzeta_index(xx,yy)]++;
                        }
                    }
                } else if (iy == iy2) {
                    /* Subpixel iy is leaving detector row y    */
                    if (ix1 < ix2) {
                        /* Subpixel iy shifts to the right from column x */
                        if (x + ix1 >= 0 && x + ix2 < ncols) {
                            xx = x + ix1;
                            yy = y + ycen_offset[x] - ycen_offset[xx];
                            xi[xi_index(x,iy,3)].x = xx;
                            xi[xi_index(x,iy,3)].y = yy;
                            xi[xi_index(x,iy,3)].w = w - fabs(delta - ix1) * w;
                            if (xx >= 0 && xx < ncols && yy >= 0 &&
                                    yy < nrows && xi[xi_index(x,iy,3)].w > 0) {
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].x = x;
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].iy = iy;
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].w = xi[xi_index(x,iy,3)].w;
                                m_zeta[mzeta_index(xx,yy)]++;
                            }
                            xx = x + ix2;
                            yy = y + ycen_offset[x] - ycen_offset[xx];
                            xi[xi_index(x,iy,2)].x = xx;
                            xi[xi_index(x,iy,2)].y = yy;
                            xi[xi_index(x,iy,2)].w = fabs(delta - ix1) * w;
                            if (xx >= 0 && xx < ncols && yy >= 0 &&
                                    yy < nrows && xi[xi_index(x,iy,2)].w > 0) {
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].x = x;
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].iy = iy;
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].w = xi[xi_index(x,iy,2)].w;
                                m_zeta[mzeta_index(xx,yy)]++;
                            }
                        }
                    } else if (ix1 > ix2) {
                        /* Subpixel iy shifts to the left from column x */
                        if (x + ix2 >= 0 && x + ix1 < ncols) {
                            xx = x + ix2;
                            yy = y + ycen_offset[x] - ycen_offset[xx];
                            xi[xi_index(x,iy,3)].x = xx;
                            xi[xi_index(x,iy,3)].y = yy;
                            xi[xi_index(x,iy,3)].w = fabs(delta - ix1) * w;
                            if (xx >= 0 && xx < ncols && yy >= 0 &&
                                    yy < nrows && xi[xi_index(x,iy,3)].w > 0) {
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].x = x;
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].iy = iy;
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].w = xi[xi_index(x,iy,3)].w;
                                m_zeta[mzeta_index(xx,yy)]++;
                            }
                            xx = x + ix1;
                            yy = y + ycen_offset[x] - ycen_offset[xx];
                            xi[xi_index(x,iy,2)].x = xx;
                            xi[xi_index(x,iy,2)].y = yy;
                            xi[xi_index(x,iy,2)].w = w - fabs(delta - ix1) * w;
                            if (xx >= 0 && xx < ncols && yy >= 0 &&
                                    yy < nrows && xi[xi_index(x,iy,2)].w > 0) {
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].x = x;
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].iy = iy;
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].w = xi[xi_index(x,iy,2)].w;
                                m_zeta[mzeta_index(xx,yy)]++;
                            }
                        }
                    } else {
                        /* Subpixel iy stays inside column x        */
                        xx = x + ix1;
                        yy = y + ycen_offset[x] - ycen_offset[xx];
                        xi[xi_index(x,iy,2)].x = xx;
                        xi[xi_index(x,iy,2)].y = yy;
                        xi[xi_index(x,iy,2)].w = w;
                        if (xx >= 0 && xx < ncols && yy >= 0 &&
                                yy < nrows && w > 0) {
                            zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].x = x;
                            zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].iy = iy;
                            zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].w = w;
                            m_zeta[mzeta_index(xx,yy)]++;
                        }
                    }
                } else {
                    /* Subpixel iy is fully inside detector row y */
                    if (ix1 < ix2) {
                        /* Subpixel iy shifts to the right from column x   */
                        if (x + ix1 >= 0 && x + ix2 < ncols) {
                            xx = x + ix1;
                            yy = y + ycen_offset[x] - ycen_offset[xx];
                            xi[xi_index(x,iy,1)].x = xx;
                            xi[xi_index(x,iy,1)].y = yy;
                            xi[xi_index(x,iy,1)].w = w - fabs(delta - ix1) * w;
                            if (xx >= 0 && xx < ncols && yy >= 0 &&
                                    yy < nrows && xi[xi_index(x,iy,1)].w > 0) {
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].x = x;
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].iy = iy;
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].w = xi[xi_index(x,iy,1)].w;
                                m_zeta[mzeta_index(xx,yy)]++;
                            }
                            xx = x + ix2;
                            yy = y + ycen_offset[x] - ycen_offset[xx];
                            xi[xi_index(x,iy,0)].x = xx;
                            xi[xi_index(x,iy,0)].y = yy;
                            xi[xi_index(x,iy,0)].w = fabs(delta - ix1) * w;
                            if (xx >= 0 && xx < ncols && yy >= 0 &&
                                    yy < nrows && xi[xi_index(x,iy,0)].w > 0) {
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].x = x;
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].iy = iy;
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].w = xi[xi_index(x,iy,0)].w;
                                m_zeta[mzeta_index(xx,yy)]++;
                            }
                        }
                    } else if (ix1 > ix2) {
                        /* Subpixel iy shifts to the left from column x */
                        if (x + ix2 >= 0 && x + ix1 < ncols) {
                            xx = x + ix2;
                            yy = y + ycen_offset[x] - ycen_offset[xx];
                            xi[xi_index(x,iy,1)].x = xx;
                            xi[xi_index(x,iy,1)].y = yy;
                            xi[xi_index(x,iy,1)].w = fabs(delta - ix1) * w;
                            if (xx >= 0 && xx < ncols && yy >= 0 &&
                                    yy < nrows && xi[xi_index(x,iy,1)].w > 0) {
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].x = x;
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].iy = iy;
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].w = xi[xi_index(x,iy,1)].w;
                                m_zeta[mzeta_index(xx,yy)]++;
                            }
                            xx = x + ix1;
                            yy = y + ycen_offset[x] - ycen_offset[xx];
                            xi[xi_index(x,iy,0)].x = xx;
                            xi[xi_index(x,iy,0)].y = yy;
                            xi[xi_index(x,iy,0)].w = w - fabs(delta - ix1) * w;
                            if (xx >= 0 && xx < ncols && yy >= 0 &&
                                    yy < nrows && xi[xi_index(x,iy,0)].w > 0) {
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].x = x;
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].iy = iy;
                                zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].w = xi[xi_index(x,iy,0)].w;
                                m_zeta[mzeta_index(xx,yy)]++;
                            }
                        }
                    } else {
                        /* Subpixel iy stays inside column x */
                        xx = x + ix2;
                        yy = y + ycen_offset[x] - ycen_offset[xx];
                        xi[xi_index(x,iy,0)].x = xx;
                        xi[xi_index(x,iy,0)].y = yy;
                        xi[xi_index(x,iy,0)].w = w;
                        if (xx >= 0 && xx < ncols && yy >= 0 &&
                                yy < nrows && w > 0) {
                            zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].x = x;
                            zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].iy = iy;
                            zeta[zeta_index(xx,yy,m_zeta[mzeta_index(xx,yy)])].w = w;
                            m_zeta[mzeta_index(xx,yy)]++;
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
  @param ncols      Swath width in pixels
  @param nrows      Extraction slit height in pixels
  @param osample    Subpixel ovsersampling factor
  @param im         Image to be decomposed [nrows][ncols]
  @param pix_unc
  @param mask       Initial and final mask for the swath [nrows][ncols]
  @param ycen       Order centre line offset from pixel row boundary [ncols]
  @param ycen_offset    Order image column shift     [ncols]
  @param y_lower_lim    Number of detector pixels below the pixel containing
                        the central line yc
  @param PSF_curve  Slit curvature
  @param delta_x    Maximum horizontal shift in detector pixels due to slit
                    image curvature
  @param sL         Slit function resulting from decomposition    [ny]
  @param sP         Spectrum resulting from decomposition      [ncols]
  @param model      Model constructed from sp and sf
  @param unc        Spectrum uncertainties based on data - model [ncols]
  @param lambda_sP  Smoothing parameter for the spectrum, could be zero
  @param lambda_sL  Smoothing parameter for the slit function, usually>0
  @param sP_stop
  @param maxiter
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
        cpl_polynomial  ** slitcurves,
        int         delta_x,
        double  *   sL,
        double  *   sP,
        double  *   model,
        double  *   unc,
        double      lambda_sP,
        double      lambda_sL,
        double      sP_stop,
        int         maxiter)
{
    int         x, xx, xxx, y, yy, iy, jy, n, m, ny, y_upper_lim, i;
    double      sum, norm, dev, lambda, diag_tot, ww, www, sP_change, sP_max;
    int         info, iter, isum;

    /* The size of the sL array. */
    /* Extra osample is because ycen can be between 0 and 1. */
    ny = osample * (nrows + 1) + 1;

    y_upper_lim = nrows - 1 - y_lower_lim;
    double *sP_old = cpl_malloc(ncols * sizeof(double));
    double *l_Aij  = cpl_malloc(ny * (4*osample+1) * sizeof(double));
    double *p_Aij  = cpl_malloc(ncols * 5 * sizeof(double));
    double *l_bj   = cpl_malloc(ny * sizeof(double));
    double *p_bj   = cpl_malloc(ncols * sizeof(double));

    /*
       Convolution tensor telling the coordinates of detector pixels on which
       {x, iy} element falls and the corresponding projections. [ncols][ny][4]
     */
    xi_ref * xi = cpl_malloc(ncols * ny * 4 * sizeof(xi_ref));

    /* Convolution tensor telling the coordinates of subpixels {x, iy}
       contributing to detector pixel {x, y}. [ncols][nrows][3*(osample+1)]
     */
    zeta_ref * zeta = cpl_malloc(ncols * nrows * 3 * (osample + 1) * sizeof(zeta_ref));

    /* The actual number of contributing elements in zeta  [ncols][nrows]  */
    int * m_zeta = cpl_malloc(ncols * nrows * sizeof(int));

    i = cr2res_extract_xi_zeta_tensors(ncols, nrows, ny, ycen, ycen_offset,
            y_lower_lim, osample, slitcurves, xi, zeta, m_zeta);

    /* Loop through sL , sP reconstruction until convergence is reached */
    iter = 0;
    do {
        /* Compute slit function sL */
        /* Prepare the RHS and the matrix */
        for (iy = 0; iy < ny; iy++) {
            l_bj[iy] = 0.e0;
            /* Clean RHS                */
            for (jy = 0; jy <= 4 * osample; jy++)
                l_Aij[iy + ny * jy] = 0.e0;
        }
        /* Fill in SLE arrays for slit function */
        diag_tot = 0.e0;
        for (iy = 0; iy < ny; iy++) {
            for (x = 0; x < ncols; x++) {
                for (n = 0; n < 4; n++) {
                    ww = xi[xi_index(x,iy,n)].w;
                    if (ww > 0) {
                        xx = xi[xi_index(x,iy,n)].x;
                        yy = xi[xi_index(x,iy,n)].y;
                        if (m_zeta[mzeta_index(xx,yy)] > 0 && xx >= 0 && xx < ncols &&
                                yy >= 0 && yy < nrows) {
                            for (m = 0; m < m_zeta[mzeta_index(xx,yy)]; m++) {
                                xxx = zeta[zeta_index(xx,yy,m)].x;
                                jy = zeta[zeta_index(xx,yy,m)].iy;
                                www = zeta[zeta_index(xx,yy,m)].w;
                                l_Aij[iy + ny * (jy - iy + 2 * osample)] +=
                                    sP[xxx] * sP[x] * www * ww * mask[yy *
                                    ncols + xx];
                            }
                            l_bj[iy] += im[yy * ncols + xx] * mask[yy * ncols +
                                xx] * sP[x] * ww;
                        }
                    }
                }
            }
            diag_tot += l_Aij[iy + ny * 2 * osample];
        }
        /* Scale regularization parameters */
        lambda = lambda_sL * diag_tot / ny;
        /* Add regularization parts for the SLE matrix */
        /* Main diagonal  */
        l_Aij[ny * 2 * osample] += lambda;
        /* Upper diagonal */
        l_Aij[ny * (2 * osample + 1)] -= lambda;
        for (iy = 1; iy < ny - 1; iy++) {
            /* Lower diagonal */
            l_Aij[iy + ny * (2 * osample - 1)] -= lambda;
            /* Main diagonal  */
            l_Aij[iy + ny * 2 * osample] += lambda * 2.e0;
            /* Upper diagonal */
            l_Aij[iy + ny * (2 * osample + 1)] -= lambda;
        }
        /* Lower diagonal */
        l_Aij[ny - 1 + ny * (2 * osample - 1)] -= lambda;
        /* Main diagonal  */
        l_Aij[ny - 1 + ny * 2 * osample] += lambda;

        /* Solve the system of equations */
        info = cr2res_extract_slitdec_bandsol(l_Aij, l_bj, ny, 4 * osample + 1);
        if (info) cpl_msg_info(__func__, "info(sL)=%d\n", info);

        /* Normalize the slit function */
        norm = 0.e0;
        for (iy = 0; iy < ny; iy++) {
            sL[iy] = l_bj[iy];
            norm += sL[iy];
        }
        norm /= osample;
        for (iy = 0; iy < ny; iy++) sL[iy] /= norm;

        /*  Compute spectrum sP */
        for (x = 0; x < ncols; x++) {
            for (xx = 0; xx < 5; xx++) p_Aij[xx * ncols + x] = 0.;
            p_bj[x] = 0;
        }
        for (x = 0; x < ncols; x++) {
            for (iy = 0; iy < ny; iy++) {
                for (n=0; n < 4; n++) {
                    ww = xi[xi_index(x,iy,n)].w;
                    if (ww > 0) {
                        xx = xi[xi_index(x,iy,n)].x;
                        yy = xi[xi_index(x,iy,n)].y;
                        if (m_zeta[mzeta_index(xx,yy)] > 0 && xx >= 0 && xx < ncols &&
                                yy >= 0 && yy < nrows) {
                            for (m = 0; m < m_zeta[mzeta_index(xx,yy)]; m++) {
                                xxx = zeta[zeta_index(xx,yy,m)].x;
                                jy = zeta[zeta_index(xx,yy,m)].iy;
                                www = zeta[zeta_index(xx,yy,m)].w;
                                p_Aij[x + ncols * (xxx - x + 2)] += sL[jy] *
                                    sL[iy] * www * ww * mask[yy * ncols + xx];
                            }
                            p_bj[x] += im[yy * ncols + xx] * mask[yy * ncols +
                                xx] * sL[iy] * ww;
                        }
                    }
                }
            }
        }

        for (x = 0; x < ncols; x++) sP_old[x] = sP[x];
        if (lambda_sP > 0.e0) {
            norm = 0.e0;
            for (x = 0; x < ncols; x++) {
                norm += sP[x];
            }
            norm /= ncols;
            lambda = lambda_sP * norm; /* Scale regularization parameter */
            p_Aij[ncols * 2] += lambda; /* Main diagonal  */
            p_Aij[ncols * 3] -= lambda; /* Upper diagonal */
            for (x = 1; x < ncols - 1; x++) {
                p_Aij[x + ncols] -= lambda;            /* Lower diagonal */
                p_Aij[x + ncols * 2] += lambda * 2.e0; /* Main diagonal  */
                p_Aij[x + ncols * 3] -= lambda;        /* Upper diagonal */
            }
            p_Aij[ncols - 1 + ncols] -= lambda;     /* Lower diagonal */
            p_Aij[ncols - 1 + ncols * 2] += lambda; /* Main diagonal  */
        }

        /* Solve the system of equations */
        info = cr2res_extract_slitdec_bandsol(p_Aij, p_bj, ncols, 5);
        if (info) cpl_msg_info(__func__, "info(sP)=%d\n", info);
        for (x = 0; x < ncols; x++) sP[x] = p_bj[x];

        /* Compute the model */
        for (y = 0; y < nrows; y++) {
            for (x = 0; x < ncols; x++) {
                model[y * ncols + x] = 0.;
            }
        }
        for (y = 0; y < nrows; y++) {
            for (x = 0; x < ncols; x++) {
                for (m = 0; m < m_zeta[mzeta_index(x,y)]; m++) {
                    xx = zeta[zeta_index(x,y,m)].x;
                    iy = zeta[zeta_index(x,y,m)].iy;
                    ww = zeta[zeta_index(x,y,m)].w;
                    model[y * ncols + x] += sP[xx] * sL[iy] * ww;
                }
            }
        }
        /* Compare model and data */
        sum = 0.e0;
        isum = 0;
        for (y = 0; y < nrows; y++) {
            for (x = delta_x; x < ncols - delta_x; x++) {
                sum += mask[y * ncols + x] * (model[y * ncols + x] -
                        im[y * ncols + x]) * (model[y * ncols + x] -
                            im[y * ncols + x]);
                isum += mask[y * ncols + x];
            }
        }
        dev = sqrt(sum / isum);
        /* Adjust the mask marking outlyers */
        for (y = 0; y < nrows; y++) {
            for (x = delta_x; x < ncols - delta_x; x++) {
                if (fabs(model[y * ncols + x] - im[y * ncols + x]) > 6. * dev)
                    mask[y * ncols + x] = 0;
                else
                    mask[y * ncols + x] = 1;
            }
        }

        /* Compute the change in the spectrum */
        sP_change = 0.e0;
        sP_max = 1.e0;
        for (x = 0; x < ncols; x++) {
            if (sP[x] > sP_max)
                sP_max = sP[x];
            if (fabs(sP[x] - sP_old[x]) > sP_change)
                sP_change = fabs(sP[x] - sP_old[x]);
        }
        /* Check for convergence */
    } while (iter++ < maxiter && sP_change > sP_stop * sP_max);

    /* Uncertainty estimate */
    for (x = 0; x < ncols; x++) {
        unc[x] = 0.;
        p_bj[x] = 0.;
    }
    for (y = 0; y < nrows; y++) {
        for (x = 0; x < ncols; x++) {
            // Loop through all pixels contributing to x,y
            for (m = 0; m < m_zeta[mzeta_index(x,y)]; m++) {
                xx = zeta[zeta_index(x,y,m)].x;
                iy = zeta[zeta_index(x,y,m)].iy;
                ww = zeta[zeta_index(x,y,m)].w;
                unc[xx] += (im[y * ncols + x] - model[y * ncols + x]) *
                    (im[y * ncols + x] - model[y * ncols + x]) *
                    ww * mask[y * ncols + x];
                unc[xx] += pix_unc[y * ncols + x] * pix_unc[y * ncols + x] *
                    ww * mask[y * ncols + x];
                // Norm
                p_bj[xx] += ww * mask[y * ncols + x];
            }
        }
    }
    for (x = 0; x < ncols; x++) {
        unc[x] = sqrt(unc[x] / p_bj[x] * nrows);
    }
    cpl_free(sP_old);
    cpl_free(l_Aij);
    cpl_free(p_Aij);
    cpl_free(l_bj);
    cpl_free(p_bj);

    cpl_free(xi);
    cpl_free(zeta);
    cpl_free(m_zeta);
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

    //if(fmod(nd,2)==0) return -1;

    /* Forward sweep */
    for(i=0; i<n-1; i++)
    {
        aa=a[i+n*(nd/2)];
        //if(aa==0.e0) return -3;
        r[i]/=aa;
        for(j=0; j<nd; j++) a[i+j*n]/=aa;
        for(j=1; j<min(nd/2+1,n-i); j++)
        {
            aa=a[i+j+n*(nd/2-j)];
            //if(aa==0.e0) return -j;
            r[i+j]-=r[i]*aa;
            for(k=0; k<n*(nd-j); k+=n) a[i+j+k]-=a[i+k+n*j]*aa;
        }
    }

    /* Backward sweep */
    r[n-1]/=a[n-1+n*(nd/2)];
    for(i=n-1; i>0; i--)
    {
        for(j=1; j<=min(nd/2,i); j++) r[i-j]-=r[i]*a[i-j+n*(nd/2+j)];
        //if(a[i-1+n*(nd/2)]==0.e0) return -5;
        r[i-1]/=a[i-1+n*(nd/2)];
    }

    //if(a[n*(nd/2)]==0.e0) return -6;
    r[0]/=a[n*(nd/2)];
    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief   Adjust the swath width to match the length of detector
  @param sw     Swath width to start from
  @param nx     number of pixel columns to match
  @return   The new swath width and swath edges. All bins end up with
            the same even swath size
    Note that the last swath shifted forward to have the same length as the
    others, therefore the overlap will be larger
 */
/*----------------------------------------------------------------------------*/
static int cr2res_extract_slitdec_adjust_swath(
        int             sw,
        int             nx,
        cpl_vector  *   bins_begin,
        cpl_vector  *   bins_end)
{
    if (sw <= 0  || nx <= 0) return -1;
    if (bins_begin == NULL || bins_end == NULL) return -1;

    if (sw > nx) sw = nx-1;
    if (sw % 2 == 1) sw += 1;

    int nbin, i = 0;
    double step = 0;

    // Calculate number of bins
    nbin = 2 * (nx / sw);
    if (nx % sw != 0) nbin++;
    if (nbin < 1) nbin = 1;

    // Step / 2, to get half width swaths
    step = sw / 2;
    int bin=0;
    cpl_vector_set_size(bins_begin, nbin);
    cpl_vector_set_size(bins_end, nbin);

    // boundaries of bins
    for(i = 2; i < nbin + 2; i++)
    {
        bin = min(i * step, nx);
        cpl_vector_set(bins_end, i-2, bin);
        cpl_vector_set(bins_begin, i-2, bin - sw);
    }
    return sw;
}

