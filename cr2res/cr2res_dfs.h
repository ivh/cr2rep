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

#ifndef CR2RES_DFS_H
#define CR2RES_DFS_H

/*-----------------------------------------------------------------------------
                                   Define
 -----------------------------------------------------------------------------*/

/* Define here the PRO.CATG keywords */
#define CR2RES_DETLIN_BPM_PROCATG       "DETLIN_BPM" 
#define CR2RES_MASTER_DARK_PROCATG      "MASTER_DARK"
#define CR2RES_MASTER_BPM_PROCATG       "MASTER_BPM"
#define CR2RES_DARK_BPM_PROCATG         "DARK_BPM"

/* Define here the DO.CATG keywords */
#define CR2RES_DARK_RAW                 "DARK"

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

cpl_error_code cr2res_dfs_set_groups(cpl_frameset *);

#endif
