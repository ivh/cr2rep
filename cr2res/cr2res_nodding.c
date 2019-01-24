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

#include <string.h>
#include <math.h>
#include <cpl.h>

#include "cr2res_nodding.h"
#include "cr2res_pfits.h"

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_utils       Nodding  Utilities
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the nodding positions from a frame set
  @param    set     Input frame set
  @return   the NODDING positions or NULL in error case
 */
/*----------------------------------------------------------------------------*/
cr2res_nodding_pos * cr2res_nodding_read_positions(const cpl_frameset * in)
{
    cr2res_nodding_pos  *   out ;
    cpl_propertylist    *   plist ;
    const char          *   fname ;
    cpl_size                nframes, i ;

    /* Check entries */
    if (in == NULL) return NULL ;

    /* Initialise */
    nframes = cpl_frameset_get_size(in) ;

    /* Allocate the vector */
    out = cpl_malloc(nframes * sizeof(cr2res_nodding_pos)) ;

    /* Loop on the frames */
    for (i=0 ; i< nframes ; i++) {
        plist = cpl_propertylist_load(cpl_frame_get_filename(
                    cpl_frameset_get_position_const(in, i)), 0) ;
        out[i] = cr2res_pfits_get_nodding_pos(plist) ;
        cpl_propertylist_delete(plist) ;
    }
    return out ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the nodding position character for display
  @param    pos     the nodding position
  @return   the character for display
 */
/*----------------------------------------------------------------------------*/
char cr2res_nodding_position_char(cr2res_nodding_pos pos) 
{
    if (pos == CR2RES_NODDING_A) return 'A' ;
    if (pos == CR2RES_NODDING_B) return 'B' ;

    return '-' ;
}

/**@}*/
