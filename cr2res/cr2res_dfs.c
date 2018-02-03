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

#include "cr2res_dfs.h"
#include "cr2res_utils.h"

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_dfs  DFS related functions
 *
 * TBD
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Set the group as RAW or CALIB in a frameset
  @param    set     the input frameset
  @return   CPL_ERROR_NONE iff OK
 */
/*----------------------------------------------------------------------------*/
cpl_error_code cr2res_dfs_set_groups(cpl_frameset * set)
{

    cpl_frame   *   frame ;
    const char  *   tag ; 
    int             nframes, i ;
    
    /* Check entries */
    if (set == NULL) return -1 ;
    
    /* Initialize */
    nframes = cpl_frameset_get_size(set) ;
    
    /* Loop on frames */
    for (i = 0 ; i < nframes ; i++) {
        frame = cpl_frameset_get_position(set, i);
        tag = cpl_frame_get_tag(frame);

        if (tag == NULL) {
            cpl_msg_warning(cpl_func, "Frame %d has no tag", i);
        } else if (!strcmp(tag, CR2RES_DETLIN_RAW) ||
                !strcmp(tag, CR2RES_DARK_RAW) ||
                !strcmp(tag, CR2RES_FLAT_RAW) ||
                !strcmp(tag, CR2RES_COMMAND_LINE)) {
            /* RAW frames */
            cpl_frame_set_group(frame, CPL_FRAME_GROUP_RAW);
        } else if (!strcmp(tag, CR2RES_DETLIN_BPM_PROCATG) ||
                !strcmp(tag, CR2RES_DETLIN_COEFFS_PROCATG) ||
                !strcmp(tag, CR2RES_MASTER_DARK_PROCATG) ||
                !strcmp(tag, CR2RES_DARK_BPM_PROCATG) ||
                !strcmp(tag, CR2RES_BLAZE_IMAGE_PROCATG) ||
                !strcmp(tag, CR2RES_BLAZE_PROCATG) ||
                !strcmp(tag, CR2RES_MASTER_FLAT_PROCATG) ||
                !strcmp(tag, CR2RES_SLIT_ILLUM_PROCATG) ||
                !strcmp(tag, CR2RES_FLAT_BPM_PROCATG) ||
                !strcmp(tag, CR2RES_EMISSION_LINES_PROCATG) ||
                !strcmp(tag, CR2RES_TILT_COEFFS_PROCATG)) {
            /* CALIB frames */
            cpl_frame_set_group(frame, CPL_FRAME_GROUP_CALIB);
        }
    }
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the TILT table column name for a given order
  @param    order       The order number (1->) 
  @return   the column name or NULL in error case
  The return string needs to be deallocated with cpl_free() 
 */
/*----------------------------------------------------------------------------*/
char * cr2res_dfs_TILT_colname(int order)
{
    int         order_loc ;
    if ((order_loc = cr2res_convert_order_to_idx(order)) < 0) return NULL ;
    return cpl_sprintf("%02d_%s", order_loc, CR2RES_COL_TILT_SUFFIX);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the SPEC table column name for a given order/trace
  @param    order       The order number (1->) 
  @param    trace       The trace number (1->)
  @return   the column name or NULL in error case
  The return string needs to be deallocated with cpl_free() 
 */
/*----------------------------------------------------------------------------*/
char * cr2res_dfs_SPEC_colname(int order, int trace)
{
    int         order_loc ;
    if ((order_loc = cr2res_convert_order_to_idx(order)) < 0) return NULL ;
    return cpl_sprintf("%02d_%02d_%s", order_loc, trace,CR2RES_COL_SPEC_SUFFIX);
}
 
/*----------------------------------------------------------------------------*/
/**
  @brief    Get the SLIT_FUNC table column name for a given order/trace
  @param    order       The order number (1->)
  @param    trace       The trace number (1->)
  @return   the column name or NULL in error case
  The return string needs to be deallocated with cpl_free() 
 */
/*----------------------------------------------------------------------------*/
char * cr2res_dfs_SLIT_FUNC_colname(int order, int trace)
{
    int         order_loc ;
    if ((order_loc = cr2res_convert_order_to_idx(order)) < 0) return NULL ;
    return cpl_sprintf("%02d_%02d_%s", order_loc,trace,
            CR2RES_COL_SLIT_FUNC_SUFFIX);
}

/**@}*/
