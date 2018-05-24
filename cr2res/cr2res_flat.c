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
    cpl_image       *   ima ;
    double          *   pima ;
    cpl_mask        *   mask ;
    cpl_binary      *   pmask ;
    int                 nx, ny, cur_bp_nb ;
    int                 i, j ;

    /* Test entries */
    if (master_flat == NULL) return NULL ;
        
    /* Initialise */
    ima = hdrl_image_get_image(master_flat) ;
    pima = cpl_image_get_data_double(ima) ;
    nx = cpl_image_get_size_x(ima) ;
    ny = cpl_image_get_size_y(ima) ;
    mask = cpl_mask_new(nx, ny) ; 
    pmask = cpl_mask_get_data(mask) ;
        
    /* Threshold to get the BPMs */
    for (j=0 ; j<ny ; j++) {
        for (i=0 ; i<nx ; i++) {
            if (isinf(pima[i+j*nx])) continue ;
            if (isnan(pima[i+j*nx])) continue ;
            if (pima[i+j*nx] < low || pima[i+j*nx] > high) {
                pmask[i+j*nx] = CPL_BINARY_1 ;
            }
        }
    }

    /*
        Post processing : Big zones of bad pixels are not considered as
        bad pixels. Each line containing more than
        100*bad_per_line_limit percent bad pixels is reset to contain
        anly good pixels.
    */
    for (j=0 ; j<ny ; j++) {
        cur_bp_nb = cpl_mask_count_window(mask, 1, j+1, nx, j+1) ;
        /* Check if the line has too many bad pixels */
        if (cur_bp_nb > bad_per_line_limit * nx) {
            /* Reset the bad pixels on the current line */
            for (i=0 ; i<nx ; i++) pmask[i+j*nx] = CPL_BINARY_0 ;
        }
    }
    return mask ;
}

