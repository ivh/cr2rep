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
#include "cr2res_bpm.h"
#include "cr2res_utils.h"
#include "cr2res_io.h"

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static cpl_mask * cr2res_bpm_compute_global_stats(
        cpl_image   *   img, 
        double          kappa) ;
static cpl_mask * cr2res_bpm_compute_local_stats(
        cpl_image   *   img, 
        double          kappa, 
        int             size) ;
static cpl_mask * cr2res_bpm_compute_running_filter(
        cpl_image   *   img, 
        double          kappa, 
        int             size) ;

/*----------------------------------------------------------------------------*/
/**
  @defgroup cr2res_bpm
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    The BPM computation with min/max threshold
  @param    in              The input image
  @param    method          CR2RES_BPM_GLOBAL_STATS
                            CR2RES_BPM_LOCAL_STATS
                            CR2RES_BPM_RUNNING_FILTER
  @param    kappa           The kappa value
  @param    clean_flag      Clean the image using the computed BPM 
  @param    lines_ratio     The maximum ratio of bad pixels per line
  @return   the BPM
 */
/*----------------------------------------------------------------------------*/
cpl_mask * cr2res_bpm_compute(
        cpl_image           *   in,
        cr2res_bpm_method       method,
        double                  kappa,
        double                  lines_ratio,
        int                     clean_flag) 
{
    cpl_mask    *   bpm ;
    cpl_binary  *   pmask_cur ;
    int             nx, ny, cur_bp_nb, i, j, k ;

    /* Test entries */
    if (in == NULL) return NULL ;
    if (method != CR2RES_BPM_GLOBAL_STATS &&
            method != CR2RES_BPM_LOCAL_STATS &&
            method != CR2RES_BPM_RUNNING_FILTER) return NULL ;

    /* This first method works for SW, not L and M */
    if (method == CR2RES_BPM_GLOBAL_STATS) {
        if ((bpm = cr2res_bpm_compute_global_stats(in, kappa)) == NULL) {
            cpl_msg_error(__func__, "Cannot create bad pixels map") ;
            return NULL ;
        }
    } else if (method == CR2RES_BPM_LOCAL_STATS) {
        if ((bpm = cr2res_bpm_compute_local_stats(in, kappa, 10)) == NULL) {
            cpl_msg_error(__func__, "Cannot create bad pixels map") ;
            return NULL ;
        }
    } else if (method == CR2RES_BPM_RUNNING_FILTER) {
        if ((bpm = cr2res_bpm_compute_running_filter(in, kappa, 10)) == NULL) {
            cpl_msg_error(__func__, "Cannot create bad pixels map") ;
            return NULL ;
        }
    }

    /*
        Post processing : Big zones of bad pixels are not considered as
        bad pixels. Each line containing more than lines_ratio percent bad 
        pixels is reset to contain only good pixels.
    */
    nx = cpl_mask_get_size_x(bpm) ;
    ny = cpl_mask_get_size_y(bpm) ;
    pmask_cur = cpl_mask_get_data(bpm) ;
    for (j=0 ; j<ny ; j++) {
        cur_bp_nb = cpl_mask_count_window(bpm, 1, j+1, nx, j+1) ;
        /* Check if the line has too many bad pixels */
        if (cur_bp_nb > lines_ratio * nx) {
            /* Reset the bad pixels on the current line */
            for (k=0 ; k<nx ; k++) {
                pmask_cur[k+j*nx] = CPL_BINARY_0 ;
            }
        }
    }

    /* Clean the image using the computed BPM */
    if (clean_flag) {
        cpl_image_reject_from_mask(in, bpm) ;
        cpl_detector_interpolate_rejected(in) ;
    }

    return bpm ;
}
/*----------------------------------------------------------------------------*/
/**
  @brief    Count BPM of a given type
  @param    bpm         the BPM
  @param    type        the bad pixel type
  @return   the count or -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_bpm_count(
        cpl_image       *   bpm,
        cr2res_bpm_type     type)
{
    cpl_image   *   tmp ;
    int             count ;
 
    /* Check Entries */
    if (bpm == NULL) return -1 ;

    tmp = cpl_image_duplicate(bpm);
    cpl_image_xor_scalar(tmp, NULL, type);
    cpl_image_reject_value(tmp, CPL_VALUE_ZERO);
    count = cpl_image_count_rejected(tmp);
    cpl_image_delete(tmp);

    return count ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Create a BPM from a mask
  @param    mask        the input mask
  @param    type        the bad pixel type
  @return   the BPM or NULL in error case
 */
/*----------------------------------------------------------------------------*/
cpl_image * cr2res_bpm_from_mask(
        cpl_mask        *   mask,
        cr2res_bpm_type     type)
{
    cpl_image   *   bpm_ima ;
 
    /* Check Entries */
    if (mask == NULL) return NULL ;

    bpm_ima = cpl_image_new_from_mask(mask) ;
    cpl_image_multiply_scalar(bpm_ima, type) ;
    return bpm_ima ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Set the BPM and optionally apply the correction to an image
  @param    in          the input image
  @param    bpm         the BPM
  @param    chip        the chip number (1 to CR2RES_NB_DETECTORS)
  @param    correction  flag to apply the correction
  @return   0 if everything is ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
int cr2res_bpm_set_and_correct_image(
        cpl_image           *   in,
        const char          *   bpm,
        int                     chip,
        int                     correct)
{
    cpl_image      		*   bpm_im ;
    cpl_mask            *   bpm_im_bin ;

    /* Check entries */
    if (in == NULL || bpm == NULL) return -1 ;
    if (chip < 1 || chip > CR2RES_NB_DETECTORS) return -1 ;

    /* Load the bpm */
    if ((bpm_im = cr2res_io_load_BPM(bpm, chip, 1)) == NULL) {
        cpl_msg_error(__func__, "Cannot load the bpm") ;
        return -1 ;
    }
    /* Convert the map to binary */
    bpm_im_bin = cpl_mask_threshold_image_create(bpm_im, -0.5, 0.5) ;
    cpl_mask_not(bpm_im_bin) ;
    cpl_image_delete(bpm_im) ;

    /* Set the Bad pixels */
    cpl_image_reject_from_mask(in, bpm_im_bin);
    cpl_mask_delete(bpm_im_bin) ;

    /* Apply the bad pixels cleaning */
    if (correct) {
        if (cpl_detector_interpolate_rejected(in) != CPL_ERROR_NONE) {
            cpl_error_reset();
            cpl_msg_error(__func__, "Cannot clean the BPM") ;
            return -1 ;
        }
    }

    return 0 ;
}


/*----------------------------------------------------------------------------*/
/** 
  @brief    Extract a mask from a BPM image
  @param    bpm_ima     BPM image to extract from
  @param    bpm_type    Can use a bitwise combination
  @return   The newly allocated mask
 */
/*----------------------------------------------------------------------------*/
cpl_mask * cr2res_bpm_extract_mask(
        const cpl_image     *   bpm_ima,
        cr2res_bpm_type         bpm_type)
{
    const int   *   pbpm_ima ;
    cpl_mask    *   bpm ;
    cpl_binary  *   pbpm ;
    cpl_size        i, j, nx, ny, idx ;

    /* Check entries */
    if (bpm_ima == NULL) return NULL ;

    /* Initialize */
    pbpm_ima = cpl_image_get_data_int_const(bpm_ima) ;
    nx = cpl_image_get_size_x(bpm_ima) ;
    ny = cpl_image_get_size_y(bpm_ima) ;

    /* Create output mask */
    bpm = cpl_mask_new(nx, ny) ;
    pbpm = cpl_mask_get_data(bpm) ;

    /* Loop on the pixels */
    for (j=0 ; j<ny ; j++) {
        for (i=0 ; i<nx ; i++) {
            idx = i+j*ny ;
            if (pbpm_ima[idx] & bpm_type) pbpm[idx] = CPL_BINARY_1 ;
        }
    }
    return bpm ;
} 

/*----------------------------------------------------------------------------*/
/**
  @brief  Add a mask to a BPM image with a dedicated code   
  @param    bpm_ima     The input BPM image
  @param    bpm         The mask to add
  @param    bpm_code    CR2RES_BPM_DARK, CR2RES_BPM_FLAT,...
  @return   0 if ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
int cr2res_bpm_add_mask(
        cpl_image   *   bpm_ima,
        cpl_mask    *   bpm,
        int             bpm_code)
{
    int         *   pbpm_ima ;
    cpl_binary  *   pbpm ;
    cpl_size        i, j, nx, ny, idx ;

    /* Check entries */
    if (bpm_ima == NULL || bpm == NULL) return -1 ;
    if (cpl_image_get_type(bpm_ima) != CPL_TYPE_INT) return -1 ;
    
    /* Initialise */
    nx = cpl_image_get_size_x(bpm_ima) ;
    ny = cpl_image_get_size_y(bpm_ima) ;
    if (nx != cpl_mask_get_size_x(bpm) || ny != cpl_mask_get_size_y(bpm)) 
        return -1 ;
    pbpm_ima = cpl_image_get_data_int(bpm_ima) ;
    pbpm = cpl_mask_get_data(bpm) ;

    /* Loop on the pixels */
    for (j=0 ; j<ny ; j++) {
        for (i=0 ; i<nx ; i++) {
            idx = i+j*ny ;
            if (pbpm[idx]) pbpm_ima[idx] = pbpm_ima[idx] | bpm_code ;
        }
    }
    return 0 ;
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Find BPM based on global statistics
  @param    img    input image of type CPL_TYPE_DOUBLE
  @param    kappa  Multiplier for threshold
  @return   A newly allocated BPM mask
 */
/*----------------------------------------------------------------------------*/
static cpl_mask * cr2res_bpm_compute_global_stats(
        cpl_image   *   img, 
        double          kappa)
{
    cpl_mask    *   out ;
    double          low, high, med, sigma ;

    /* Compute Thresholds */
    med = cpl_image_get_median_dev(img, &sigma) ;
    if (cpl_error_get_code()) {
        cpl_msg_error(__func__, "Cannot compute statistics") ;
        return NULL ;
    }
    low = med - kappa * sigma ;
    high = med + kappa * sigma ;

    cpl_msg_debug(__func__, 
            "Median %.1f, Sigma %.1f" "BPM_low %.1f, BPM_hi %.1f", 
            med, sigma, low, high);

    /* Threshold to get the BPMs */
    if ((out = cpl_mask_threshold_image_create(img, low, high)) == NULL) {
        cpl_msg_error(__func__, "Cannot create bad pixels map") ;
        return NULL ;
    }
    cpl_mask_not(out) ;

    return out ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Find BPM based on the median of the surounding pixels
  @param    img    input image of type CPL_TYPE_DOUBLE
  @param    kappa  Multiplier for threshold
  @param    size   Half-size of the statistics box
  @return   A newly allocated BPM mask
 */
/*----------------------------------------------------------------------------*/
static cpl_mask * cr2res_bpm_compute_local_stats(
        cpl_image   *   img, 
        double          kappa, 
        int             size)
{
    cpl_mask        *   out ;
    cpl_binary      *   pout ;
    double              threshold, med, sigma ;
    int                 badpix;
    cpl_size            i, j, nx, ny, llx, lly, urx, ury ;

    /* Check Entries */
    if (img == NULL) return NULL ;
    nx = cpl_image_get_size_x(img) ;
    ny = cpl_image_get_size_y(img) ;
    if (sigma < 0.0) return NULL ;

    /* Create Output mask */
    out = cpl_mask_new(nx, ny) ;
    pout = cpl_mask_get_data(out) ;

    /* Loop on the pixels */
    for (i = 0; i < nx ; i++) {
        for (j = 0; j < ny ; j++) {
            llx = i+1-size ;
            lly = j+1-size ;
            urx = i+1+size ;
            ury = j+1+size ;
            if (llx < 1) llx = 1 ;
            if (lly < 1) lly = 1 ;
            if (urx > nx) urx = nx ;
            if (ury > ny) ury = ny ;

            /* Local statistics */
            med = cpl_image_get_median_dev_window(img, llx, lly, urx,
                    ury, &sigma); 

            /* Compute Threshold */
            threshold = med + kappa*sigma ;

            /* Set Bad Pixel */
            if (fabs(cpl_image_get(img, i+1, j+1, &badpix)) > threshold) 
                pout[i + j*nx] = CPL_BINARY_1 ;
        }   
    }
    return out;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Find BPM based on the median of the surounding pixels
  @author   Ansgar
  @param    img    input image
  @param    kappa  Multiplier for threshold, equivalent to the gaussian sigma
  @param    size   Size of the Median filter applied to the image, must be odd
  @return   A newly allocated BPM mask

  Ansgar function
  This function applies a median filter to the image and uses the MAD of the 
  difference between the input image and the smoothed version to determine a
  threshold. Pixels that deviate by more than the threshold are discarded.
 */
/*----------------------------------------------------------------------------*/
static cpl_mask * cr2res_bpm_compute_running_filter(
        cpl_image   *   img, 
        double          kappa, 
        int             size)
{
    cpl_mask    *   out ;
    cpl_binary  *   pout ;
    cpl_image   *   copy ;
    cpl_mask    *   kernel ;
    double          mad, median;
    int             badpix;
    double          threshold;
    cpl_size        i, j, nx, ny ;

    /* Check Entries */
    if (img == NULL) return NULL ;
    nx = cpl_image_get_size_x(img) ;
    ny = cpl_image_get_size_y(img) ;
    if (nx < size || ny < size) return NULL ;
    if (size % 2 != 1) return NULL ;

    /* Local copy */
    copy = cpl_image_duplicate(img);

    /* Filtering Kernel */
    kernel = cpl_mask_new(size, 1);
    cpl_mask_not(kernel);

    /* Apply Filter */
    cpl_image_filter_mask(copy, img, kernel, CPL_FILTER_MEDIAN, 
            CPL_BORDER_FILTER);
    cpl_mask_delete(kernel);

    /* Reject 0 values */
    for (i = 1; i <= nx ; i++) {
        for (j = 1; j <= ny ; j++) {
            if (cpl_image_get(copy, i, j, &badpix) == 0) 
                cpl_image_reject(copy, i, j);
        }   
    }

    /* Subtract smoothed image */
    cpl_image_subtract(copy, img);

    /* Apply running thresholding */
    median = cpl_image_get_mad(copy, &mad);
    mad *= CPL_MATH_STD_MAD;
    cpl_image_subtract_scalar(copy, median);
    cpl_image_abs(copy);
    
    threshold = kappa * mad;
    out = cpl_mask_new(nx, ny) ;
    pout = cpl_mask_get_data(out) ;

    for (i = 0; i < nx ; i++) {
        for (j = 0; j < ny ; j++) {
            if (cpl_image_get(copy, i+1, j+1, &badpix) > threshold) 
                pout[i + j*nx] = CPL_BINARY_1 ;
        }   
    }
    cpl_image_delete(copy);
    return out;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Find BPM based on the median of the surounding pixels
  @author   Thomas
  @param    bpm    input BPM as cpl_image
  @return   modified copy of BPM cpl_image

    Set the CR2RES_NB_BPM_EDGEPIX pixels around the
    detector edges to CR2RES_BPM_EDGEPIX
 */
/*----------------------------------------------------------------------------*/
cpl_image * cr2res_bpm_mask_edgepix(cpl_image * bpm){
    cpl_size    i,j;
    cpl_image * out;
    out = cpl_image_duplicate(bpm);
    for (i=1; i<=CR2RES_DETECTOR_SIZE; i++){
        for (j=1; j<=CR2RES_NB_BPM_EDGEPIX; j++){
            cpl_image_set(out,i,j,CR2RES_BPM_EDGEPIX);
            cpl_image_set(out,i,CR2RES_DETECTOR_SIZE-j+1,CR2RES_BPM_EDGEPIX);
            cpl_image_set(out,CR2RES_DETECTOR_SIZE-j+1,i,CR2RES_BPM_EDGEPIX);
            cpl_image_set(out,j,i,CR2RES_BPM_EDGEPIX);
        }
    }
    return out;
};