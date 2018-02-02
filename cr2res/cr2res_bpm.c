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

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @defgroup cr2res_bpm
 */
/*----------------------------------------------------------------------------*/

/**@{*/

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


