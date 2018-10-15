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

#ifndef CR2RES_SPLICE_H
#define CR2RES_SPLICE_H

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>
#include <hdrl.h>

/*-----------------------------------------------------------------------------
                                   Define
 -----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
                                       Prototypes
 -----------------------------------------------------------------------------*/

int cr2res_splice(
        cpl_table       **  extracted_1d,
        cpl_table       **  blaze,
        cpl_table       **  trace_wave,
        int                 ninputs,
        cpl_bivector    **  spliced,
        cpl_bivector    **  spliced_err) ;

int cr2res_splice_orders(
        cpl_table   *   trace_wave,
        cpl_table   *   spectra,
        int             trace) ;

cpl_table * cr2res_splice_SPLICED_create(
        cpl_bivector    *   spectrum,
        cpl_bivector    *   spectrum_error) ;

#endif
