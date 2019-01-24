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

#ifndef CR2RES_NODDING_H
#define CR2RES_NODDING_H

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>
#include <hdrl.h>

/*-----------------------------------------------------------------------------
                                   Define
 -----------------------------------------------------------------------------*/

typedef enum {
    CR2RES_NODDING_A,
    CR2RES_NODDING_B,
    CR2RES_NODDING_NONE,
} cr2res_nodding_pos ; 

/*-----------------------------------------------------------------------------
                                       Prototypes
 -----------------------------------------------------------------------------*/

cr2res_nodding_pos * cr2res_nodding_read_positions( const cpl_frameset * in) ;
char cr2res_nodding_position_char(cr2res_nodding_pos pos) ;
int cr2res_combine_nodding_split(
        const hdrl_imagelist    *   in,
        cr2res_nodding_pos      *   positions,
        hdrl_imagelist          **  pos_a,
        hdrl_imagelist          **  pos_b) ;

#endif
