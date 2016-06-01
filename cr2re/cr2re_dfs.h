/* 
 * This file is part of the CR2RE Pipeline
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

#ifndef CR2RE_DFS_H
#define CR2RE_DFS_H

/*-----------------------------------------------------------------------------
                                   Define
 -----------------------------------------------------------------------------*/

/* Define here the PRO.CATG keywords */
#define CR2RE_BIAS_PROCATG              "CALPRO_BIAS"
#define CR2RE_TRACE_PROCATG             "CALPRO_TRACE"

/* Define here the DO.CATG keywords */
#define CR2RE_BIAS_RAW                  "BIAS"
#define CR2RE_TRACE_RAW                 "TRACE"
#define CR2RE_ETALON_RAW                "ETALON"
#define CR2RE_TEST_RAW                  "TEST"

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

cpl_error_code cr2re_dfs_set_groups(cpl_frameset *);

#endif
