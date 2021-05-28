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

#ifndef CR2RES_ETALON_H
#define CR2RES_ETALON_H

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>

/*-----------------------------------------------------------------------------
                                       Prototypes
 -----------------------------------------------------------------------------*/

cpl_image * cr2res_etalon_computation(const cpl_image * in) ;
cpl_vector * cr2res_etalon_get_maxpos(const cpl_vector *) ;
cpl_vector * cr2res_etalon_find_peaks(
    const cpl_vector * in, 
    double height, 
    double distance);
cpl_polynomial * cr2res_etalon_wave_2d(
    cpl_bivector        **  spectra,
    cpl_bivector        **  spectra_err,
    cpl_polynomial      **  wavesol_init,
    cpl_array           **  wavesol_init_err,
    int                 *   orders,
    int                 *   traces_nb,
    int                     ninputs,
    cpl_size                degree_x,
    cpl_size                degree_y,
    cpl_array           **  wavelength_error,
    cpl_table           **  line_diagnostics);
cpl_polynomial * cr2res_etalon_wave_2d_nikolai(
    cpl_bivector        **  spectra,
    cpl_bivector        **  spectra_err,
    cpl_polynomial      **  wavesol_init,
    cpl_array           **  wavesol_init_err,
    int                 *   orders,
    int                 *   traces_nb,
    int                     ninputs,
    cpl_size                degree_x,
    cpl_size                degree_y,
    cpl_array           **  wavelength_error,
    cpl_table           **  line_diagnostics);
#endif
