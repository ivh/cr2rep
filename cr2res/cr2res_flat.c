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


/*----------------------------------------------------------------------------*/
/**
  @defgroup cr2res_flat
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief  Main function for the flat computation
  @param    imlist      input image list 
  @param  
  @return The newly allocated master flat
 */
/*----------------------------------------------------------------------------*/
hdrl_image * cr2res_flat(
        const cpl_imagelist *   imlist,
		const cpl_image		*	master_dark,
		const cpl_imagelist	*	detlin_coeffs)
{
    cpl_imagelist   *   imlist_loc ;
    hdrl_image      *   master_flat ;
    cpl_image       *   collapsed ;
	
    /* Check Entries */
    if (imlist == NULL) return NULL ;

    /* Duplicate the input list - Not optimal but more readable */
    imlist_loc = cpl_imagelist_duplicate(imlist) ;

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
    
    /* Compute the MASTER FLAT */
    collapsed = cpl_imagelist_collapse_create(imlist_loc) ;
    cpl_imagelist_delete(imlist_loc) ;
    master_flat = hdrl_image_create(collapsed, NULL) ;
    cpl_image_delete(collapsed) ;
    return master_flat ;
}

/**@}*/

