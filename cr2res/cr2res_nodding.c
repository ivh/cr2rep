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
    cpl_size                nframes, i ;

    /* Check entries */
    if (in == NULL) return NULL ;

    /* Initialise */
    nframes = cpl_frameset_get_size(in) ;

    /* Allocate the vector */
    out = cpl_malloc(nframes * sizeof(cr2res_nodding_pos)) ;
    for (i = 0; i < nframes; i++) out[i] = CR2RES_NODDING_NONE;

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

/*----------------------------------------------------------------------------*/
/**
  @brief  Split A/B positions in 2 framesets
  @param    in          Input frameset
  @param    positions   nodding positions
  @param    pos_a       [out] A position frames
  @param    pos_b       [out] B position frames
  @return   0 if ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
int cr2res_combine_nodding_split_frames(
        const cpl_frameset      *   in,
        cr2res_nodding_pos      *   positions,
        cpl_frameset            **  pos_a,
        cpl_frameset            **  pos_b)
{
    cpl_frameset    *   nod_a ;
    cpl_frameset    *   nod_b ;
    const cpl_frame *   cur_frame ;
    cpl_size            alist_idx, blist_idx, i, nframes ;

    /* Check entries */
    if (in==NULL || positions==NULL || pos_a==NULL || pos_b==NULL)
        return -1;

    /* Initialise */
    nframes = cpl_frameset_get_size(in) ;
    alist_idx = blist_idx = 0 ;

    /* Create A/B positions */
    nod_a = cpl_frameset_new() ;
    nod_b = cpl_frameset_new() ;

    /* Loop on the positions */
    for (i=0 ; i<nframes ; i++) {
        if (positions[i] == CR2RES_NODDING_A) {
            cur_frame = cpl_frameset_get_position_const(in, i) ;
            cpl_frameset_insert(nod_a, cpl_frame_duplicate(cur_frame)) ;
            alist_idx++ ;
            cur_frame = NULL ;
        }
        if (positions[i] == CR2RES_NODDING_B) {
            cur_frame = cpl_frameset_get_position_const(in, i) ;
            cpl_frameset_insert(nod_b, cpl_frame_duplicate(cur_frame)) ;
            blist_idx++ ;
            cur_frame = NULL ;
        }
    }
    *pos_a = nod_a ;
    *pos_b = nod_b ;
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief  Split A/B positions in 2 image lists
  @param    in          Input image list
  @param    positions   nodding positions
  @param    pos_a       [out] A position images
  @param    pos_b       [out] B position images
  @return   0 if ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
int cr2res_combine_nodding_split(
        const hdrl_imagelist    *   in,
        cr2res_nodding_pos      *   positions,
        hdrl_imagelist          **  pos_a,
        hdrl_imagelist          **  pos_b)
{
    hdrl_imagelist  *   nod_a ;
    hdrl_imagelist  *   nod_b ;
    hdrl_image      *   cur_ima ;
    cpl_size            alist_idx, blist_idx, i, nima ;

    /* Check entries */
    if (in==NULL || positions==NULL || pos_a==NULL || pos_b==NULL)
        return -1;

    /* Initialise */
    nima = hdrl_imagelist_get_size(in) ;
    alist_idx = blist_idx = 0 ;

    /* Create A/B positions */
    nod_a = hdrl_imagelist_new() ;
    nod_b = hdrl_imagelist_new() ;

    /* Loop on the positions */
    for (i=0 ; i<nima ; i++) {
        if (positions[i] == CR2RES_NODDING_A) {
            cur_ima = hdrl_image_duplicate(hdrl_imagelist_get(in, i)) ;
            hdrl_imagelist_set(nod_a, cur_ima, alist_idx) ;
            alist_idx++ ;
            cur_ima = NULL ;
        }
        if (positions[i] == CR2RES_NODDING_B) {
            cur_ima = hdrl_image_duplicate(hdrl_imagelist_get(in, i)) ;
            hdrl_imagelist_set(nod_b, cur_ima, blist_idx) ;
            blist_idx++ ;
            cur_ima = NULL ;
        }
    }
    *pos_a = nod_a ;
    *pos_b = nod_b ;
    return 0 ;
}

/**@}*/
