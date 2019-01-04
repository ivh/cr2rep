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

static int cr2res_detlin_correct(
        cpl_imagelist       *   ilist,
        const cpl_imagelist *   detlin_coeffs) ;

/*----------------------------------------------------------------------------*/
/**
  @defgroup cr2res_calib
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    The images calibration routine for a given chip
  @param    ilist       the list of frames to calibrate 
  @param    chip        the chip to calibrate (1 to CR2RES_NB_DETECTORS)
  @param    cosmics_corr    Flag to correct for cosmics
  @param    flat        the flat frame or NULL
  @param    dark        the dark frame or NULL
  @param    bpm         the bpm frame or NULL
  @param    detlin      the detlin frame or NULL
  @param    dit         the DIT for the dark correction
  The flat, dark and bpm must have the same size as the input ilist.
  In the case of detlin, data are only taken in normal mode.
  @return   0 if ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
int cr2res_calib_chip_list(
        cpl_imagelist       *   ilist,
        int                     chip,
        int                     cosmics_corr,
        const cpl_frame     *   flat,
        const cpl_frame     *   dark,
        const cpl_frame     *   bpm,
        const cpl_frame     *   detlin,
        double                  dit)
{
    cpl_image           *   calib ;
    cpl_imagelist       *   calib_list ;
    cpl_propertylist    *   plist ;
    double                  dark_dit ;
    int                     i ;

    /* Test entries */
    if (ilist == NULL) return -1 ;
    if (chip < 1 || chip > CR2RES_NB_DETECTORS) return -1 ;

    /* Clean the bad pixels */
    if (bpm != NULL) {
        cpl_msg_info(__func__, "Correct the bad pixels") ;
        for (i=0 ; i<cpl_imagelist_get_size(ilist) ; i++) {
            if (cr2res_bpm_correct_image(cpl_imagelist_get(ilist, i),
                        cpl_frame_get_filename(bpm), chip) != 0) {
                cpl_msg_error(__func__,
                        "Cannot clean the bad pixels in obj %d", i+1);
                return -1 ;
            }
        }
    }

    /* Apply the dark */
    if (dark != NULL) {
        cpl_msg_info(__func__, "Correct for the dark") ;
        for (i=0 ; i<cpl_imagelist_get_size(ilist) ; i++) {

            /* Load the dark */
            if ((calib = cr2res_io_load_MASTER_DARK(
                            cpl_frame_get_filename(dark), chip, 1)) == NULL) {
                cpl_msg_error(__func__, "Cannot load the dark") ;
                return -1 ;
            }

            /* Get the dark DIT */
            plist = cpl_propertylist_load(cpl_frame_get_filename(dark), 0);
            dark_dit = cr2res_pfits_get_dit(plist) ;
            cpl_propertylist_delete(plist) ;

            /* Multiply the dark by dit/dark_dit */
            cpl_image_multiply_scalar(calib, dit/dark_dit) ;

            /* Subtract the dark */
            if (cpl_image_subtract(cpl_imagelist_get(ilist, i), 
                        calib) != CPL_ERROR_NONE) {
                cpl_msg_error(__func__, "Cannot apply the dark") ;
                cpl_image_delete(calib) ;
                return -1 ;
            }
            cpl_image_delete(calib) ;
        }
    }

    /* Apply the non linearity correction */
    if (detlin != NULL) {
        /* Load the detlin coeffs */
        cpl_msg_info(__func__, "Load the Non-Linearity coefficients") ;
        if ((calib_list = cr2res_io_load_DETLIN_COEFFS(
                        cpl_frame_get_filename(detlin), chip, 1)) == NULL) {
            cpl_msg_error(__func__, "Cannot load the detlin") ;
            return -1 ;
        }

        /* Detlin correction */
        cpl_msg_info(__func__, "Correct for the Non-Linearity") ;
        if (cr2res_detlin_correct(ilist, calib_list)) {
            cpl_imagelist_delete(calib_list) ;
            cpl_msg_error(__func__, "Cannot correct for the Non-Linearity") ;
            return -1 ;
        }
        cpl_imagelist_delete(calib_list) ;
    }

    /* Apply the flatfield */
    if (flat != NULL) {
        /* Load the flat */
        cpl_msg_info(__func__, "Load the flat field") ;
        if ((calib = cr2res_io_load_MASTER_FLAT(
                        cpl_frame_get_filename(flat), chip, 1)) == NULL) {
            cpl_msg_error(__func__, "Cannot load the flat field") ;
            return -1 ;
        }
        
        /* Divide */
        cpl_msg_info(__func__, "Correct for the flat field") ;
        if (cpl_imagelist_divide_image(ilist, calib) != CPL_ERROR_NONE) {
            cpl_msg_error(__func__, "Cannot apply the flat field") ;
            cpl_image_delete(calib) ;
            return -1 ;
        }
        cpl_image_delete(calib) ;
    }

    /* Comics correction */
    if (cosmics_corr) {
        cpl_msg_info(__func__, "Apply the cosmics corrections") ;
        /* TODO */
        cpl_msg_info(__func__, "STILL TO BE DONE") ;
    }
    return 0 ;
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Apply the detector linearity correction
  @param    ilist   the input image list
  @param    detlin      the detlin coeffs
  @return   0 if everything is ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
static int cr2res_detlin_correct(
        cpl_imagelist       *   ilist,
        const cpl_imagelist *   detlin_coeffs)
{
    const cpl_image     *   ima ;
    const cpl_image     *   imb ;
    const cpl_image     *   imc ;
    const float         *   pima ;
    const float         *   pimb ;
    const float         *   pimc ;
    float               *   pdata ;
    int                     nx, ny, ni ;
    float                   val, val2, val3 ;
    int                     i, j ;

    /* Test entries */
    if (!ilist || !detlin_coeffs) return -1 ;

    /* Load the 3 coeffs images */
    ima = cpl_imagelist_get_const(detlin_coeffs, 0) ;
    imb = cpl_imagelist_get_const(detlin_coeffs, 1) ;
    imc = cpl_imagelist_get_const(detlin_coeffs, 2) ;
    if (!ima || !imb || !imc) {
        cpl_msg_error(cpl_func, "Cannot access the detlin images") ;
        return -1 ;
    }
    pima = cpl_image_get_data_float_const(ima) ;
    pimb = cpl_image_get_data_float_const(imb) ;
    pimc = cpl_image_get_data_float_const(imc) ;

    /* Test sizes */
    nx = cpl_image_get_size_x(cpl_imagelist_get(ilist, 0)) ;
    ny = cpl_image_get_size_y(cpl_imagelist_get(ilist, 0)) ;
    ni = cpl_imagelist_get_size(ilist) ;
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
        if (fabs(pimc[i]) < 1e-5) {
            /* Correct this pixel in each plane */
            for (j=0 ; j<ni ; j++) {
                pdata = cpl_image_get_data_float(cpl_imagelist_get(ilist, j)) ;
                val = pdata[i] ;
                pdata[i] = val-pima[i] ;
            }
        } else if (fabs(pimb[i]) < 1e-3) {
            for (j=0 ; j<ni ; j++) {
                pdata = cpl_image_get_data_float(cpl_imagelist_get(ilist, j)) ;
                pdata[i] = 0.0 ;
            }
        } else {
            val2 = 2 * pimc[i] / (pimb[i] * pimb[i]) ;
            /* Correct this pixel in each plane */
            for (j=0 ; j<ni ; j++) {
                pdata = cpl_image_get_data_float(cpl_imagelist_get(ilist, j)) ;
                val = pdata[i] ;
                val3 = 1-2*val2*(pima[i]-val) ;
                if (val3 < 0.0) {
                    pdata[i] = val-pima[i] ;
                } else {
                    pdata[i]=((float)sqrt(val3)-1) / val2 ;
                }
            }
        }
    }
    /* return */
    return 0 ;
}



