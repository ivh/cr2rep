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
#include "cr2res_detlin.h"
#include "cr2res_utils.h"

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

int cr2res_add_shotnoise(hdrl_image * in, int ndit, int chip);
int cr2res_calib_subtract_interorder_column(hdrl_image * in,
    const hdrl_image * flat, cpl_size degree);

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
  @param    clean_bad   Flag to activate the cleaning of the bad pixels 
  @param    subtract_nolight_rows   
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
        int                         clean_bad,
        int                         subtract_nolight_rows,
        int                         subtract_interorder_column,
        int                         cosmics_corr,
        const cpl_frame         *   flat,
        const cpl_frame         *   dark,
        const cpl_frame         *   bpm,
        const cpl_frame         *   detlin,
        const cpl_vector        *   dits,
        const cpl_vector        *   ndits)
{
    hdrl_imagelist      *   out ;
    const hdrl_image    *   cur_ima ;
    hdrl_image          *   cur_ima_calib ;
    double                  dit ;
    int                     ndit ;
    cpl_size                i ;

    /* Check Inputs */
    if (in == NULL) return NULL ;

    /* Initialise */
    dit = 0.0 ;
    ndit = 1 ;

    /* Create calibrated image list */
    out = hdrl_imagelist_new() ;

    /* Loop on the images */
    for (i=0 ; i<hdrl_imagelist_get_size(in) ; i++) {
        cur_ima = hdrl_imagelist_get(in, i) ;
        if (dark != NULL)  dit = cpl_vector_get(dits, i) ;
        if (ndits != NULL) ndit = (int)cpl_vector_get(ndits, i) ;

        /* Calibrate */
        cpl_msg_info(__func__, 
                "Apply calibrations for image #%"CPL_SIZE_FORMAT, i+1) ;
        cpl_msg_indent_more() ;
        if ((cur_ima_calib = cr2res_calib_image(cur_ima, chip, clean_bad, 
                        subtract_nolight_rows, subtract_interorder_column,
                        cosmics_corr, flat, dark, bpm, 
                        detlin, dit, ndit)) == NULL) {
            cpl_msg_error(__func__, "Failed to Calibrate the Data") ;
            hdrl_imagelist_delete(out) ;
            cpl_msg_indent_less() ;
            return NULL ;
        } else {
            /* All the calibrated image in the list */
            hdrl_imagelist_set(out, cur_ima_calib, i);
        }
        cpl_msg_indent_less() ;
    }
    return out ;
}
 
/*----------------------------------------------------------------------------*/
/**
  @brief    The images calibration routine for a given chip
  @param    in          the input hdrl image
  @param    chip        the chip to calibrate (1 to CR2RES_NB_DETECTORS)
  @param    clean_bad   Flag to activate the cleaning of the bad pixels 
  @param    subtract_nolight_rows   
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
        int                     clean_bad,
        int                     subtract_nolight_rows,
        int                     subtract_interorder_column,
        int                     cosmics_corr,
        const cpl_frame     *   flat,
        const cpl_frame     *   dark,
        const cpl_frame     *   bpm,
        const cpl_frame     *   detlin,
        double                  dit,
        int                     ndit)
{
    hdrl_image          *   out ;
    hdrl_image          *   calib ;
    hdrl_imagelist      *   calib_list ;
    cpl_propertylist    *   plist ;
    double                  dark_dit ;
    cpl_image           *   img_tmp ;
    cpl_size                i;

    /* Test entries */
    if (in == NULL) return NULL ;
    if (chip < 1 || chip > CR2RES_NB_DETECTORS) return NULL ;

    /* Create out image */
    out = hdrl_image_duplicate(in) ;

    /* Clean the bad pixels */
    if (bpm != NULL) {
        cpl_msg_info(__func__, "Correct the bad pixels") ;
        if (cr2res_bpm_set_and_correct_image(hdrl_image_get_image(out),
                    cpl_frame_get_filename(bpm), chip, clean_bad) != 0) {
            cpl_msg_error(__func__, "Cannot clean the bad pixels");
            hdrl_image_delete(out);
            return NULL ;
        }
    }

    /* Apply the non linearity correction */
    if (detlin != NULL) {
        cpl_msg_info(__func__, "Correct for the Non-Linearity") ;
        if ((calib_list = cr2res_io_load_DETLIN_COEFFS(
                        cpl_frame_get_filename(detlin), chip)) == NULL) {
            cpl_msg_error(__func__, "Cannot load the detlin") ;
            hdrl_image_delete(out);
            return NULL ;
        }
        /* Detlin correction */
        if (cr2res_detlin_correct(out, calib_list)) {
            hdrl_imagelist_delete(calib_list) ;
            hdrl_image_delete(out);
            cpl_msg_error(__func__, "Cannot correct for the Non-Linearity") ;
            return NULL ;
        }
    }

    /* Add shot-noise */
    cpl_msg_info(__func__, "Add shot-noise") ;
    if (cr2res_add_shotnoise(out, ndit, chip)){
        cpl_msg_error(__func__, "Cannot add shot-noise") ;
        hdrl_imagelist_delete(calib_list) ;
        hdrl_image_delete(out);
        return NULL ;
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
        if (detlin != NULL) {
            cpl_msg_indent_more() ;
            cpl_msg_info(__func__, "Correct DARK for Non-Linearity") ;
            if (cr2res_detlin_correct(calib, calib_list)) {
                hdrl_imagelist_delete(calib_list) ;
                hdrl_image_delete(calib) ;
                hdrl_image_delete(out);
                cpl_msg_error(__func__,"Cannot correct DARK for Non-Linearity");
                cpl_msg_indent_less() ;
                return NULL ;
            }
            cpl_msg_indent_less() ;
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
    if (detlin != NULL) {
       hdrl_imagelist_delete(calib_list) ;
    }

    /* Subtract residual bias/dark from vignetted rows at bottom */
    if (subtract_nolight_rows) {
        cpl_msg_info(__func__, "Subtract median of bottom rows");

        img_tmp = cpl_image_collapse_median_create(
                hdrl_image_get_image(out), 0, CR2RES_NB_BPM_EDGEPIX, 
                CR2RES_DETECTOR_SIZE-CR2RES_NB_BPM_VIGN_BOTTOM);
        calib = hdrl_image_new(CR2RES_DETECTOR_SIZE,CR2RES_DETECTOR_SIZE);
        for (i=1; i<=CR2RES_DETECTOR_SIZE; i++) {
            hdrl_image_insert(calib, img_tmp, NULL, 1, i);
        }
        hdrl_image_sub_image(out, calib);
        hdrl_image_delete(calib) ;
        cpl_image_delete(img_tmp);
    }

    /* Apply the interorder correction */
    if (subtract_interorder_column) {
        if (flat == NULL) {
            cpl_msg_info(__func__, 
                "Skip subtracting fit to inter-order pixels (no flat-field)");
        } else {
            cpl_msg_info(__func__, "Subtract fit to inter-order pixels");
            if ((calib = cr2res_io_load_MASTER_FLAT(
                            cpl_frame_get_filename(flat), chip)) == NULL) {
                cpl_msg_error(__func__, "Cannot load the flat field") ;
                hdrl_image_delete(out);
                return NULL ;
            }
            if (cr2res_calib_subtract_interorder_column(out, calib, 0)){
                cpl_msg_error(__func__, "Could not subtract inter-order fit");
            }
            hdrl_image_delete(calib);
        }
    }

    /* Apply the flatfield */
    if (flat != NULL) {
        cpl_msg_info(__func__, "Correct for the flat field") ;
        /* Load the flat */
        if ((calib = cr2res_io_load_MASTER_FLAT(
                        cpl_frame_get_filename(flat), chip)) == NULL) {
            cpl_msg_error(__func__, "Cannot load the flat field") ;
            hdrl_image_delete(out);
            return NULL ;
        }

        /* Divide */
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
        hdrl_parameter * cospar = hdrl_lacosmic_parameter_create(12.0,2.5,1);
        if (hdrl_lacosmic_parameter_verify(cospar)!=CPL_ERROR_NONE) 
            cpl_msg_warning(__func__,"invalid cospar");
        cpl_mask* cosmask = hdrl_lacosmic_edgedetect(out,cospar);
        cpl_mask_or(hdrl_image_get_mask(out), cosmask);
        cpl_mask_delete(cosmask);
        hdrl_parameter_delete(cospar);
    }

    return out ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Add shot-noise to errors in HDRL-image
  @param    in          the input hdrl image, gets modified
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_add_shotnoise(hdrl_image * in, int ndit, int chip){
    double gain_sqrt;
    cpl_image * error = hdrl_image_get_error(in);
    cpl_image * adu  = hdrl_image_get_image(in);
    cpl_image * tmp_im;

    if (adu==NULL){
        cpl_msg_error(__func__,"Broken input image");
        return -1;
    }
    if (error==NULL){
        cpl_msg_error(__func__,"Error input null not supported");
    }

    if      (chip == 1) gain_sqrt = sqrt(CR2RES_GAIN_CHIP1);
    else if (chip == 2) gain_sqrt = sqrt(CR2RES_GAIN_CHIP2);
    else if (chip == 3) gain_sqrt = sqrt(CR2RES_GAIN_CHIP3);
    else {
        cpl_msg_error(__func__,"Unknown detector");
        return -1;
    }

    cpl_msg_debug(__func__, "chip:%d, sqrtgain:%g, ndit:%d",
                            chip, gain_sqrt, ndit);
    
    if ( (tmp_im=cpl_image_abs_create(adu)) == NULL){
        cpl_msg_error(__func__,"Abs failed");
        return -1;
    }
    if ( (cpl_image_power(tmp_im, 0.5)) != CPL_ERROR_NONE){
        cpl_msg_error(__func__,"Sqrt failed");
        return -1;
    }
    cpl_image_divide_scalar(tmp_im, gain_sqrt);
    cpl_image_divide_scalar(tmp_im, sqrt((float)ndit));

    /* The BPM should not be stored in the error image */
    /* Therefore it is removed from tmp_im before the addition to error */
    cpl_image_accept_all(tmp_im) ;

    cpl_image_add(error, tmp_im);
    cpl_image_delete(tmp_im);
    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Remove background column by column
  @param    in          the input hdrl image, gets modified
  @return   0 if ok, -1 in error case

  Fits the background in each column, using a 1d polynomial with given degree 
  to the out of order pixels.

  Here the orders are defined by the flat-field.

 */
/*----------------------------------------------------------------------------*/
int cr2res_calib_subtract_interorder_column(hdrl_image * in, 
                                    const hdrl_image * flat, cpl_size degree)
{
    // Define Parameters
    cpl_image * img;
    const cpl_image * flat_img;
    cpl_polynomial * poly;
    cpl_matrix * px;
    cpl_vector * py;
    cpl_vector * tmp;
    cpl_size ncolumns, nrow, npixel;
    cpl_size i, j;
    int badpix;
    int badcols = 0;
    double value, median, mad;

    // Check Input 
    if (in == NULL || flat == NULL) return -1;
    if (hdrl_image_get_size_x(in) != hdrl_image_get_size_x(flat)) return -1;
    if (hdrl_image_get_size_y(in) != hdrl_image_get_size_y(flat)) return -1;

    // Loop over the columns
    img = hdrl_image_get_image(in);
    flat_img = hdrl_image_get_image_const(flat);
    ncolumns = cpl_image_get_size_x(img);
    nrow = cpl_image_get_size_y(img);


    for (i = CR2RES_NB_BPM_EDGEPIX + 1; 
                i < ncolumns - CR2RES_NB_BPM_EDGEPIX; i++)
    {
        // Fill vectors and matrix
        px = cpl_matrix_new(1, nrow);
        py = cpl_vector_new(nrow);
        npixel = 0;
        for ( j = CR2RES_NB_BPM_VIGN_BOTTOM + 5; j < nrow -20 ; j++)
        {
            // TODO: check that the pixels are actually rejected
            // Filter pixels that are rejected in the flat (i.e. between orders)
            // but not rejected by the image (i.e. not bad pixels)
            if (cpl_image_is_rejected(flat_img, i, j) 
                    && !cpl_image_is_rejected(img, i, j)){
                cpl_matrix_set(px, 0, npixel, j);
                cpl_vector_set(py, npixel, cpl_image_get(img, i, j, &badpix));
                npixel++;
            }
        }
        cpl_msg_debug(__func__, "Column %lld has %lld pixels for background",
                i,npixel);
        if (npixel == 0){
            badcols += 1;
            cpl_msg_debug(__func__,
              "Could not find any valid points for column %"CPL_SIZE_FORMAT, i);
            cpl_matrix_delete(px);
            cpl_vector_delete(py);
            continue;
        }
        cpl_matrix_set_size(px, 1, npixel);
        cpl_vector_set_size(py, npixel);

        // Reject outliers and filter
        median = cr2res_vector_get_mad(py, &mad);
        for (j=0;j<npixel;j++){
            if (cpl_vector_get(py,j) > median*2)
                cpl_vector_set(py,j,median);
        }
        //tmp = cpl_vector_filter_median_create(py,50);
        tmp = cpl_vector_filter_lowpass_create(py,CPL_LOWPASS_LINEAR,100);
        cpl_vector_delete(py);
        py=tmp;
        if (cpl_msg_get_level() == CPL_MSG_DEBUG){
            tmp = cpl_vector_wrap(npixel, cpl_matrix_get_data(px));
            cpl_vector_save(tmp, "debug_background.fits",
                         CPL_TYPE_DOUBLE, NULL, CPL_IO_CREATE);
            cpl_vector_save(py, "debug_background.fits",
                         CPL_TYPE_DOUBLE, NULL, CPL_IO_EXTEND);
            cpl_vector_unwrap(tmp);
        }

        // Fit polynomial
        poly = cpl_polynomial_new(1);
        if (cpl_polynomial_fit(poly, px, NULL, py, NULL, CPL_FALSE, NULL,
                                                                    &degree)
                != CPL_ERROR_NONE) {
            badcols += 1;
            cpl_msg_warning(__func__,
            "Could not fit the background of column %"CPL_SIZE_FORMAT, i);
            cpl_error_reset();
            cpl_matrix_delete(px);
            cpl_vector_delete(py);
            cpl_polynomial_delete(poly);
            continue;
        }

        // Subtract column
        for ( j = 1; j < nrow + 1; j++)
        {
            value = cpl_image_get(img, i, j, &badpix);
            // Just skip this pixel if it is bad
            if (badpix) continue; 
            value -= cpl_polynomial_eval_1d(poly, j, NULL);
            cpl_image_set(img, i, j, value);
        }

        // Clean up
        cpl_matrix_delete(px);
        cpl_vector_delete(py);
        cpl_polynomial_delete(poly);
    }

    if (badcols){
        cpl_msg_warning(__func__, 
        "%d columns could not be background-ubtracted from inter-order gaps",
            badcols);
    }
    return 0;
}

/**@}*/

