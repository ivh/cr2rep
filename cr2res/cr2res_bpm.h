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

#ifndef CR2RES_BPM_H
#define CR2RES_BPM_H

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>
#include "hdrl.h"

/*-----------------------------------------------------------------------------
                                    Define
 -----------------------------------------------------------------------------*/

typedef enum _cr2res_bpm_type_ {
    CR2RES_BPM_DARK  	= 1 << 0,
    CR2RES_BPM_FLAT     = 1 << 1
} cr2res_bpm_type ;

#define CR2RES_NB_BPM_TYPES     2 
static cr2res_bpm_type bpm_types[CR2RES_NB_BPM_TYPES] = {
    CR2RES_BPM_DARK, 
    CR2RES_BPM_FLAT};

/*-----------------------------------------------------------------------------
                                Prototypes
 -----------------------------------------------------------------------------*/

cpl_mask * cr2res_bpm_compute(
        cpl_image   *   in,
        double          low,
        double          high,
        double          lines_ratio,
        int             clean_flag) ;

int cr2res_bpm_count(
        cpl_image       *   bpm,
        cr2res_bpm_type     type) ;

cpl_image * cr2res_bpm_merge(
        cpl_image   *   bpm1,
        cpl_image   *   bpm2) ;

cpl_image * cr2res_bpm_from_mask(
        cpl_mask        *   mask,
        cr2res_bpm_type     type) ;

int cr2res_bpm_correct_image(
        cpl_image           *   in,
        const char          *   bpm,
        int                     chip) ;

cpl_mask * cr2res_bpm_extract_mask(
        const cpl_image     *   bpm_ima,
        cr2res_bpm_type         bpm_type) ;

int cr2res_bpm_add_mask(
        cpl_image   *   bpm_ima,
        cpl_mask    *   bpm,
        int             bpm_code) ;

#endif
