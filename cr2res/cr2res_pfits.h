/* $Id: cr2re_pfits.h,v 1.1.1.1 2013-03-21 13:53:45 yjung Exp $
 *
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

/*
 * $Author: yjung $
 * $Date: 2013-03-21 13:53:45 $
 * $Revision: 1.1.1.1 $
 * $Name: not supported by cvs2svn $
 */

#ifndef CR2RE_PFITS_H
#define CR2RE_PFITS_H

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>

/*-----------------------------------------------------------------------------
                                   Functions prototypes
 -----------------------------------------------------------------------------*/

const char * cr2re_pfits_get_arcfile(const cpl_propertylist *) ;
double cr2re_pfits_get_dit(const cpl_propertylist *) ;

#endif
