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
#include "cr2res_calib.h"
#include "cr2res_bpm.h"
#include "cr2res_pfits.h"
#include "cr2res_io.h"
#include "cr2res_utils.h"

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

#define pow2(x) (x)*(x)
#define pow3(x) (x)*(x)*(x)


/*----------------------------------------------------------------------------*/
/**
  @defgroup cr2res_calib
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    The images calibration routine for a given chip on a list
  @param    in          the input hdrl image list
  @param    chip        the chip to calibrate (1 to CR2RES_NB_DETECTORS)
  @param    cosmics_corr    Flag to correct for cosmics
  @param    flat        the flat frame or NULL
  @param    dark        the dark frame or NULL
  @param    bpm         the bpm frame or NULL
  @param    detlin      the detlin frame or NULL
  @param    dits        the DITs of the images for the dark correction
  The flat, dark and bpm must have the same size as the input in.
  In the case of detlin, data are only taken in normal mode.
  @return   the newly allocated imagelist or NULL in error case
 */
/*----------------------------------------------------------------------------*/
hdrl_imagelist * cr2res_calib_imagelist(
        const hdrl_imagelist    *   in,
        int                         chip,
        int                         cosmics_corr,
        const cpl_frame         *   flat,
        const cpl_frame         *   dark,
        const cpl_frame         *   bpm,
        const cpl_frame         *   detlin,
        const cpl_vector        *   dits)
{
    hdrl_imagelist      *   out ;
    const hdrl_image    *   cur_ima ;
    hdrl_image          *   cur_ima_calib ;
    double                  dit ;
    cpl_size                i ;

    /* Check Inputs */
    if (in == NULL) return NULL ;

    /* Initialise */
    dit = 0.0 ;

    /* Create calibrated image list */
    out = hdrl_imagelist_new() ;

    /* Loop on the images */
    for (i=0 ; i<hdrl_imagelist_get_size(in) ; i++) {
        cur_ima = hdrl_imagelist_get(in, i) ;
        if (dark != NULL) dit = cpl_vector_get(dits, i) ;

        /* Calibrate */
        if ((cur_ima_calib = cr2res_calib_image(cur_ima, chip, cosmics_corr, 
                        flat, dark, bpm, detlin, dit)) == NULL) {
            cpl_msg_error(__func__, "Failed to Calibrate the Data") ;
            hdrl_imagelist_delete(out) ;
            return NULL ;
        } else {
            /* All the calibrated image in the list */
            hdrl_imagelist_set(out, cur_ima_calib, i);
        }
    }
    return out ;
}
 
/*----------------------------------------------------------------------------*/
/**
  @brief    The images calibration routine for a given chip
  @param    in          the input hdrl image
  @param    chip        the chip to calibrate (1 to CR2RES_NB_DETECTORS)
  @param    cosmics_corr    Flag to correct for cosmics
  @param    flat        the flat frame or NULL
  @param    dark        the dark frame or NULL
  @param    bpm         the bpm frame or NULL
  @param    detlin      the detlin frame or NULL
  @param    dit         the DIT for the dark correction
  The flat, dark and bpm must have the same size as the input in.
  In the case of detlin, data are only taken in normal mode.
  @return   the newly allocated image or NULL in error case
 */
/*----------------------------------------------------------------------------*/
hdrl_image * cr2res_calib_image(
        const hdrl_image    *   in,
        int                     chip,
        int                     cosmics_corr,
        const cpl_frame     *   flat,
        const cpl_frame     *   dark,
        const cpl_frame     *   bpm,
        const cpl_frame     *   detlin,
        double                  dit)
{
    hdrl_image          *   out ;
    hdrl_image          *   calib ;
    hdrl_imagelist      *   calib_list ;
    cpl_propertylist    *   plist ;
    double                  dark_dit ;

    /* Test entries */
    if (in == NULL) return NULL ;
    if (chip < 1 || chip > CR2RES_NB_DETECTORS) return NULL ;

    /* Create out image */
    out = hdrl_image_duplicate(in) ;

    /* Clean the bad pixels */
    if (bpm != NULL) {
        cpl_msg_info(__func__, "Correct the bad pixels") ;
        if (cr2res_bpm_set_and_correct_image(hdrl_image_get_image(out),
                    cpl_frame_get_filename(bpm), chip, 0) != 0) {
            cpl_msg_error(__func__, "Cannot clean the bad pixels");
            hdrl_image_delete(out);
            return NULL ;
        }
    }

    /* Apply the dark */
    if (dark != NULL) {
        cpl_msg_info(__func__, "Correct for the dark") ;

        /* Load the dark */
        if ((calib = cr2res_io_load_MASTER_DARK(cpl_frame_get_filename(dark), 
                        chip)) == NULL) {
            cpl_msg_error(__func__, "Cannot load the dark") ;
            hdrl_image_delete(out);
            return NULL ;
        }

        /* Get the dark DIT */
        plist = cpl_propertylist_load(cpl_frame_get_filename(dark), 0);
        dark_dit = cr2res_pfits_get_dit(plist) ;
        cpl_propertylist_delete(plist) ;

        /* Multiply the dark by dit/dark_dit */
        hdrl_value hdrl_dit_corr = {dit/dark_dit, 0.0};
        hdrl_image_mul_scalar(calib, hdrl_dit_corr) ;

        /* Subtract the dark */
        if (hdrl_image_sub_image(out, calib) != CPL_ERROR_NONE) {
            cpl_msg_error(__func__, "Cannot apply the dark") ;
            hdrl_image_delete(calib) ;
            hdrl_image_delete(out);
            return NULL ;
        }
        hdrl_image_delete(calib) ;
    }

    /* Apply the non linearity correction */
    if (detlin != NULL) {
        /* Load the detlin coeffs */
        cpl_msg_info(__func__, "Load the Non-Linearity coefficients") ;
        if ((calib_list = cr2res_io_load_DETLIN_COEFFS(
                        cpl_frame_get_filename(detlin), chip)) == NULL) {
            cpl_msg_error(__func__, "Cannot load the detlin") ;
            hdrl_image_delete(out);
            return NULL ;
        }

        /* Detlin correction */
        cpl_msg_info(__func__, "Correct for the Non-Linearity") ;
        if (cr2res_detlin_correct(out, calib_list)) {
            hdrl_imagelist_delete(calib_list) ;
            hdrl_image_delete(out);
            cpl_msg_error(__func__, "Cannot correct for the Non-Linearity") ;
            return NULL ;
        }
        hdrl_imagelist_delete(calib_list) ;
    }

    /* Apply the flatfield */
    if (flat != NULL) {
        /* Load the flat */
        cpl_msg_info(__func__, "Load the flat field") ;
        if ((calib = cr2res_io_load_MASTER_FLAT(
                        cpl_frame_get_filename(flat), chip)) == NULL) {
            cpl_msg_error(__func__, "Cannot load the flat field") ;
            hdrl_image_delete(out);
            return NULL ;
        }
        
        /* Divide */
        cpl_msg_info(__func__, "Correct for the flat field") ;
        if (hdrl_image_div_image(out, calib) != CPL_ERROR_NONE) {
            cpl_msg_error(__func__, "Cannot apply the flat field") ;
            hdrl_image_delete(calib) ;
            hdrl_image_delete(out);
            return NULL ;
        }
        hdrl_image_delete(calib) ;
    }

    /* Comics correction */
    if (cosmics_corr) {
        cpl_msg_info(__func__, "Apply the cosmics corrections") ;
        /* TODO */
        cpl_msg_info(__func__, "NOT YET IMPLEMENTED") ;
    }
    return out ;
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Apply the detector linearity correction
  @param    in      the input image 
  @param    detlin  the detlin coeffs
  @return   0 if everything is ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
int cr2res_detlin_correct(
        hdrl_image              *   in,
        const hdrl_imagelist    *   detlin)
{
    const cpl_image     *   ima ;
    const cpl_image     *   erra ;
    const cpl_image     *   imb ;
    const cpl_image     *   errb ;
    const cpl_image     *   imc ;
    const cpl_image     *   errc ;
    const double        *   pima ;
    const double        *   perra ;
    const double        *   pimb ;
    const double        *   perrb ;
    const double        *   pimc ;
    const double        *   perrc ;
    cpl_image           *   cur_ima ;
    double              *   pdata ;
    double              *   perr ;
    int                     nx, ny ;
    double                  correction_factor;
    int                     i, j ;

    /* Test entries */
    if (!in || !detlin) return -1 ;

    /* Initialise */
    pdata = cpl_image_get_data_double(hdrl_image_get_image(in)) ;
    perr = cpl_image_get_data_double(hdrl_image_get_error(in)) ;


    /* Load the 3 coeffs images */
    ima = hdrl_image_get_image_const(hdrl_imagelist_get_const(detlin, 0)) ;
    erra = hdrl_image_get_error_const(hdrl_imagelist_get_const(detlin, 0)) ;
    imb = hdrl_image_get_image_const(hdrl_imagelist_get_const(detlin, 1)) ;
    errb = hdrl_image_get_error_const(hdrl_imagelist_get_const(detlin, 1)) ;
    imc = hdrl_image_get_image_const(hdrl_imagelist_get_const(detlin, 2)) ;
    errc = hdrl_image_get_error_const(hdrl_imagelist_get_const(detlin, 2)) ;
    
    if (!ima || !imb || !imc) {
        cpl_msg_error(cpl_func, "Cannot access the detlin images") ;
        return -1 ;
    }
    pima = cpl_image_get_data_double_const(ima) ;
    pimb = cpl_image_get_data_double_const(imb) ;
    pimc = cpl_image_get_data_double_const(imc) ;
    perra = cpl_image_get_data_double_const(erra);
    perrb = cpl_image_get_data_double_const(errb);
    perrc = cpl_image_get_data_double_const(errc);

    /* Test sizes */
    cur_ima = hdrl_image_get_image(in) ;
    nx = cpl_image_get_size_x(cur_ima) ;
    ny = cpl_image_get_size_y(cur_ima) ;
    if ((cpl_image_get_size_x(ima) != nx) ||
            (cpl_image_get_size_x(imb) != nx) ||
            (cpl_image_get_size_x(imc) != nx) ||
            (cpl_image_get_size_y(ima) != ny) ||
            (cpl_image_get_size_y(imb) != ny) ||
            (cpl_image_get_size_y(imc) != ny)) {
        cpl_msg_error(cpl_func, "Incompatible sizes") ;
        return -1 ;
    }

    /* Loop on pixels */
    for (i=0 ; i<nx*ny ; i++) {
        // for each pixel p' = a + b * p + c * p * p

        perr[i] = pow2(perra[i] * pdata[i]) + pow2(perrb[i] * pow2(pdata[i])) + pow2(perrc[i] * pow3(pdata[i])) 
                + pow2(perr[i] * (pima[i] + 2. * pimb[i] * pdata[i] + 3. * pimc[i] * pow2(pdata[i])));
        perr[i] = sqrt(perr[i]);

        correction_factor = pima[i] + (pimb[i] + pimc[i] * pdata[i]) * pdata[i];
        pdata[i] = pdata[i] * correction_factor;
    }
    /* return */
    return 0 ;
}



