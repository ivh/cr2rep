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
                                    Define
 -----------------------------------------------------------------------------*/

#define CR2RES_WAVELENGTH_ERROR_DEFAULT     0.5

typedef enum {
    CR2RES_XCORR,
    CR2RES_LINE1D,
    CR2RES_LINE2D,
    CR2RES_ETALON,
    CR2RES_UNSPECIFIED
} cr2res_wavecal_type ;

/*-----------------------------------------------------------------------------
                                       Prototypes
 -----------------------------------------------------------------------------*/

cpl_polynomial * cr2res_wave_1d(
        cpl_bivector        *   spectrum,
        cpl_bivector        *   spectrum_err,
        cpl_polynomial      *   wavesol_init,
        const cpl_array     *   wave_error_init,
        int                     order,
        int                     trace_nb,
        cr2res_wavecal_type     wavecal_type,
        const char          *   catalog,
        int                     degree,
        int                     log_flag,
        int                     display,
        cpl_array           **  wavelength_error,
        cpl_table           **  lines_diagnostics) ;

cpl_polynomial * cr2res_wave_2d(
        cpl_bivector        **  spectra,
        cpl_bivector        **  spectra_err,
        cpl_polynomial      **  wavesol_init,
        cpl_array           **  wavesol_init_err,
        int                 *   orders,
        int                 *   traces_nb,
        int                     ninputs,
        const char          *   catalog,
        cpl_size                degree_x,
        cpl_size                degree_y,
        int                     display,
        cpl_array           **  wavelength_error,
        cpl_table           **  lines_diagnostics) ;

cpl_polynomial * cr2res_wave_xcorr(
        cpl_bivector    *   spectrum,
        cpl_polynomial  *   wavesol_init,
        double              wl_error,
        cpl_bivector    *   lines_list,
        int                 degree,
        int                 cleaning_filter_size,
        double              slit_width,
        double              fwhm,
        int                 display,
        cpl_array       **  wavelength_error) ;

cpl_polynomial * cr2res_wave_line_fitting(
        cpl_bivector    *   spectrum,
        cpl_bivector    *   spectrum_err,
        cpl_polynomial  *   wavesol_init,
        const cpl_array *   wave_error_init,
        int                 order,
        int                 trace_nb,
        cpl_bivector    *   lines_list,
        int                 degree,
        int                 display,
        cpl_vector      **  sigma_fit,
        cpl_array       **  wavelength_error,
        cpl_table       **  lines_diagnostics) ;

cpl_polynomial * cr2res_wave_etalon(
        cpl_bivector    *   spectrum,
        cpl_bivector    *   spectrum_err,
        cpl_polynomial  *   wavesol_init,
        int                 degree,
        cpl_array       **  wavelength_error) ;

cpl_vector * cr2res_wave_etalon_measure_fringes(cpl_vector * spectrum);

double cr2res_wave_etalon_get_x0(
        cpl_vector * xi, 
        cpl_polynomial * initial_guess);

double cr2res_wave_etalon_get_D(
        cpl_vector * peaks);

cpl_bivector * cr2res_wave_etalon_assign_fringes(
            const cpl_vector      * peaks_found,
            const cpl_vector      * peaks_should);

cpl_bivector * cr2res_wave_gen_lines_spectrum(
        const char      *   catalog,
        cpl_polynomial  *   wavesol_init,
        double              wl_error,
        double              max_intensity,
        int                 log_flag) ;

cpl_array * cr2res_wave_get_estimate(
        const char  *   filename,
        int             detector,
        int             order) ;

hdrl_image * cr2res_wave_gen_wave_map(
        const cpl_table *   trace_wave) ;

#endif
