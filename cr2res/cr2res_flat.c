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
#include "cr2res_flat.h"

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static cpl_mask * cr2res_bpm_from_master_flat(
        const hdrl_image    *   master_flat,
        double                  low,
        double                  high,
        double                  bad_per_line_limit) ;

/*----------------------------------------------------------------------------*/
/**
  @defgroup cr2res_flat
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Calibrate and Collapse
  @param    imlist      input image list 
  @param  
  @return The newly allocated calibrated image
 */
/*----------------------------------------------------------------------------*/
hdrl_image * cr2res_calib_collapse(
        const cpl_imagelist     *   imlist,
        const cpl_image         *   master_dark,
        const cpl_imagelist     *   detlin_coeffs,
        int                         cosmics_corr)
{
    cpl_imagelist   *   imlist_loc ;
    hdrl_image      *   out ;
    cpl_image       *   collapsed ;
    
    /* Check Entries */
    if (imlist == NULL) return NULL ;

    /* Duplicate the input list - Not optimal but more readable */
    imlist_loc = cpl_imagelist_duplicate(imlist) ;

    /* Correct the Cosmics */
    if (cosmics_corr) {
        cpl_msg_info(__func__, "Correct the COSMICS") ;
        /* TODO */
    }

    /* Correct the DETLIN */
    if (detlin_coeffs != NULL) {
        cpl_msg_info(__func__, "Correct the DETLIN") ;
        if (cr2res_detlin_correct(imlist_loc, detlin_coeffs)) {
            cpl_msg_error(__func__, "Failed to correct the detlin") ;
            cpl_imagelist_delete(imlist_loc) ;
            return NULL ;
        }
    }

    /* Correct the dark */
    if (master_dark != NULL) {
        cpl_msg_info(__func__, "Correct the DARK") ;
        if (cpl_imagelist_subtract_image(imlist_loc, 
                    master_dark) != CPL_ERROR_NONE) {
            cpl_msg_error(__func__, "Failed to correct the dark") ;
            cpl_imagelist_delete(imlist_loc) ;
            return NULL ;
        }
    }
    
    /* Collapse the calibrated frames */
    cpl_msg_info(__func__, "Collapse the calibrated frames") ;
    if ((collapsed = cpl_imagelist_collapse_create(imlist_loc)) == NULL) {
        cpl_msg_error(__func__, "Failed to collapse") ;
        cpl_imagelist_delete(imlist_loc) ;
        return NULL ;
    }
    cpl_imagelist_delete(imlist_loc) ;
    out = hdrl_image_create(collapsed, NULL) ;
    cpl_image_delete(collapsed) ;
    return out ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Compute the Master Flat
  @param    imlist      input image list 
  @param  
  @return The newly allocated calibrated image
 */
/*----------------------------------------------------------------------------*/
hdrl_image * cr2res_master_flat(
        const hdrl_image    *   collapsed,
        const hdrl_image    *   model_master,
        double                  low,
        double                  high,
        double                  bad_per_line_limit,
        cpl_mask            **  bpm)
{
    hdrl_image      *   master_flat ;
    cpl_mask        *   bpm_loc ;
    
    /* Check Entries */
    if (collapsed == NULL) return NULL ;
    if (model_master == NULL) return NULL ;

    /* Compute the master flat */
    if ((master_flat = hdrl_image_div_image_create(collapsed, 
                    model_master)) == NULL) {
        cpl_msg_error(__func__, "Failed to divide by the model") ;
        return NULL ;
    }

    /* Compute BPM */
    if ((bpm_loc = cr2res_bpm_from_master_flat(master_flat, low, high,
                    bad_per_line_limit)) == NULL) {

        cpl_msg_error(__func__, "Failed to compute the BPM") ;
        hdrl_image_delete(master_flat) ;
        return NULL ;
    }

    *bpm = bpm_loc ;
    return master_flat ;
} 

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the BPM from the master flat
  @param    master_flat The master flat image
  @return The newly allocated BPM
 */
/*----------------------------------------------------------------------------*/
static cpl_mask * cr2res_bpm_from_master_flat(
        const hdrl_image    *   master_flat,
        double                  low,
        double                  high,
        double                  bad_per_line_limit)
{
    cpl_mask        *   mask ;
    cpl_binary      *   pmask ;
    int                 nx, ny, cur_bp_nb ;
    int                 j, k ;

    /* Test entries */
    if (master_flat == NULL) return NULL ;
        
    /* Threshold to get the BPMs */
    if ((mask = cpl_mask_threshold_image_create(
                    hdrl_image_get_image_const(master_flat), low, high))==NULL){
        cpl_msg_error(__func__, "Cannot create bad pixels map") ;
        return NULL ;
    }
    cpl_mask_not(mask) ;

    /*
        Post processing : Big zones of bad pixels are not considered as
        bad pixels. Each line containing more than
        100*bad_per_line_limit percent bad pixels is reset to contain
        anly good pixels.
    */
    nx = cpl_mask_get_size_x(mask) ;
    ny = cpl_mask_get_size_y(mask) ;
    pmask = cpl_mask_get_data(mask) ;
    for (j=0 ; j<ny ; j++) {
        cur_bp_nb = cpl_mask_count_window(mask, 1, j+1, nx, j+1) ;
        /* Check if the line has too many bad pixels */
        if (cur_bp_nb > bad_per_line_limit * nx) {
            /* Reset the bad pixels on the current line */
            for (k=0 ; k<nx ; k++) pmask[k+j*nx] = CPL_BINARY_0 ;
        }
    }

    return mask ;
}

