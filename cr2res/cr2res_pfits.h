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

#ifndef CR2RES_PFITS_H
#define CR2RES_PFITS_H

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>
#include "cr2res_nodding.h"
#include "cr2res_utils.h"

/*-----------------------------------------------------------------------------
                                   Functions prototypes
 -----------------------------------------------------------------------------*/

cr2res_nodding_pos cr2res_pfits_get_nodding_pos(const cpl_propertylist * plist);

const char * cr2res_pfits_get_procatg(const cpl_propertylist *) ;
const char * cr2res_pfits_get_protype(const cpl_propertylist *) ;
const char * cr2res_pfits_get_wlen_id(const cpl_propertylist *) ;
const char * cr2res_pfits_get_arcfile(const cpl_propertylist *) ;

double cr2res_pfits_get_dit(const cpl_propertylist *) ;
double cr2res_pfits_get_wmin(const cpl_propertylist * plist, int order) ;
double cr2res_pfits_get_wmax(const cpl_propertylist * plist, int order) ;
double cr2res_pfits_get_wstrt(const cpl_propertylist * plist, int order) ;
double cr2res_pfits_get_wend(const cpl_propertylist * plist, int order) ;
double cr2res_pfits_get_ceny(const cpl_propertylist * plist, int order) ;

int cr2res_pfits_get_naxis1(const cpl_propertylist * plist) ;
int cr2res_pfits_get_naxis2(const cpl_propertylist * plist) ;
int cr2res_pfits_get_expno(const cpl_propertylist * plist) ;
int cr2res_pfits_get_ndit(const cpl_propertylist * plist) ;
int cr2res_pfits_get_order(const cpl_propertylist * plist, double yposition) ;
cr2res_decker cr2res_pfits_get_decker_position(const cpl_propertylist * plist) ;

#endif
