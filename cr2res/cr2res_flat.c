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
#include "cr2res_flat.h"

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static cpl_image * cr2res_bpm_from_master_flat(
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
  @param    collapsed
  @param    model_master
  @param    low
  @param    high
  @param    bad_per_line_limit
  @param    bpm                 [output] the computed BPM image
  @return The newly allocated master flat hdrl image
 */
/*----------------------------------------------------------------------------*/
hdrl_image * cr2res_master_flat(
        const hdrl_image    *   collapsed,
        const hdrl_image    *   model_master,
        double                  low,
        double                  high,
        double                  bad_per_line_limit,
        cpl_image           **  bpm)
{
    hdrl_image      *   master_flat ;
    cpl_image       *   bpm_loc ;
    
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

    if (bpm != NULL)    *bpm = bpm_loc ;
    else                cpl_image_delete(bpm_loc) ;
    return master_flat ;
} 

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the BPM from the master flat
  @param    master_flat The master flat image
  @return   The newly allocated BPM image
 */
/*----------------------------------------------------------------------------*/
static cpl_image * cr2res_bpm_from_master_flat(
        const hdrl_image    *   master_flat,
        double                  low,
        double                  high,
        double                  bad_per_line_limit)
{
    const cpl_image *   ima ;
    const double    *   pima ;
    cpl_mask        *   mask_flat ;
    cpl_mask        *   mask_inter_order ;
    cpl_binary      *   pmask_flat ;
    cpl_binary      *   pmask_inter_order ;
    cpl_image       *   bpm ;
    int                 nx, ny;
    int                 i, j ;

    /* Test entries */
    if (master_flat == NULL) return NULL ;
        
    /* Initialise */
    ima = hdrl_image_get_image_const(master_flat) ;
    pima = cpl_image_get_data_double_const(ima) ;
    nx = cpl_image_get_size_x(ima) ;
    ny = cpl_image_get_size_y(ima) ;
    mask_flat = cpl_mask_new(nx, ny) ; 
    pmask_flat = cpl_mask_get_data(mask_flat) ;
    mask_inter_order = cpl_mask_new(nx, ny) ; 
    pmask_inter_order = cpl_mask_get_data(mask_inter_order) ;
        
    /* Threshold to get the BPMs */
    for (j=0 ; j<ny ; j++) {
        for (i=0 ; i<nx ; i++) {
            if (isinf(pima[i+j*nx])) continue ;
            if (isnan(pima[i+j*nx])) {
                pmask_inter_order[i+j*nx] = CPL_BINARY_1 ;
            } else if (pima[i+j*nx] < low || pima[i+j*nx] > high) {
                pmask_flat[i+j*nx] = CPL_BINARY_1 ;
            }
        }
    }

    /*
        Post processing : Big zones of bad pixels are not considered as
        bad pixels. Each line containing more than
        100*bad_per_line_limit percent bad pixels is reset to contain
        only good pixels.
    */
    for (j=0 ; j<ny ; j++) {
        int cur_bp_nb;
        cur_bp_nb = cpl_mask_count_window(mask_flat, 1, j+1, nx, j+1) ;
        /* Check if the line has too many bad pixels */
        if (cur_bp_nb > bad_per_line_limit * nx) {
            /* Reset the bad pixels on the current line */
            for (i=0 ; i<nx ; i++) pmask_flat[i+j*nx] = CPL_BINARY_0 ;
        }
    }

    /* Create BPM */
    bpm = cpl_image_new(nx, ny, CPL_TYPE_INT) ;
    cr2res_bpm_add_mask(bpm, mask_flat, CR2RES_BPM_FLAT) ;
    cr2res_bpm_add_mask(bpm, mask_inter_order, CR2RES_BPM_OUTOFORDER) ;

    cpl_mask_delete(mask_flat) ;
    cpl_mask_delete(mask_inter_order) ;
    return bpm ;
}

