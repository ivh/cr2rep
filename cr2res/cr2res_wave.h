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

#ifndef CR2RES_WAVE_H
#define CR2RES_WAVE_H

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>
#include "hdrl.h"
#include "cr2res_utils.h"

/*-----------------------------------------------------------------------------
                                       Prototypes
 -----------------------------------------------------------------------------*/

cpl_polynomial * cr2res_wave(
        cpl_vector          *   spectrum,
        cpl_polynomial      *   initial_guess,
        cr2res_wavecal_type     wavecal_type,
        int                     line_fitting,
        const char          *   static_file,
        int                     degree,
        int                 	display) ;
cpl_polynomial * cr2res_wave_xcorr(
        cpl_vector      *   spectrum,
        cpl_polynomial  *   initial_guess,
        int                 wl_error,
        cpl_bivector    *   lines_list,
        int                 degree,
        int                 display) ;
cpl_polynomial * cr2res_wave_line_fitting(
        cpl_vector      *   spectrum,
        cpl_polynomial  *   initial_guess,
        cpl_table       *   catalog) ;

cpl_polynomial * cr2res_wave_etalon(
        cpl_vector      *   spectrum,
        cpl_polynomial  *   initial_guess) ;
cpl_vector * cr2res_wave_etalon_measure_fringes(cpl_vector * spectrum);
double cr2res_wave_etalon_fringe_stats(
        cpl_vector * peaks, 
        cpl_polynomial * initial_guess);
cpl_bivector * cr2res_wave_etalon_assign_fringes(
            const cpl_vector      * peaks_found,
            const cpl_vector      * peaks_should);

cpl_vector * cr2res_wave_line_detection(
        cpl_vector      *   spectrum) ;

cpl_bivector * cr2res_wave_gen_lines_spectrum(
        const char      *   catalog,
        cpl_polynomial  *   initial_guess,
        int                 wl_error) ;

cpl_array * cr2res_wave_get_estimate(
        const char  *   filename,
        int             detector,
        int             order) ;

hdrl_image * cr2res_wave_gen_wave_map(
        const cpl_table *   trace_wave) ;

#endif
