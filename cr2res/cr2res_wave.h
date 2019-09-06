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
#define CR2RES_WAVELENGTH_MIN_FIT_PIX       20

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

int cr2res_wave_apply(
        cpl_table           *       tw_in,
        cpl_table           *       spectra_tab,
        const cpl_frame     *       catalog_frame,
        int                         reduce_order,
        int                         reduce_trace,
        cr2res_wavecal_type         wavecal_type,
        int                         degree,
        double                      wl_start,
        double                      wl_end,
        double                      wl_err_start,
        double                      wl_err_end,
        double                      wl_shift,
        int                         log_flag,
        int                         propagate_flag,
        int                         display,
        double                      display_wmin,
        double                      display_wmax,
        cpl_propertylist    **      qcs,
        cpl_table           **      lines_diagnostics,
        cpl_table           **      extracted_out,
        cpl_table           **      trace_wave_out) ;

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
        double                  display_wmin,
        double                  display_wmax,
        double              *       best_xcorr,
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
        double          *   best_xcorr,
        cpl_array       **  wavelength_error) ;

cpl_polynomial * cr2res_wave_etalon(
        cpl_bivector    *   spectrum,
        cpl_bivector    *   spectrum_err,
        cpl_polynomial  *   wavesol_init,
        int                 degree,
        cpl_array       **  wavelength_error) ;

cpl_polynomial * cr2res_wave_estimate_compute(
        double          wmin,
        double          wmax) ;

cpl_array * cr2res_wave_get_estimate(
        const char  *   filename,
        int             detector,
        int             order) ;

hdrl_image * cr2res_wave_gen_wave_map(
        const cpl_table *   trace_wave) ;

cpl_polynomial * cr2res_wave_polys_1d_to_2d(
        cpl_polynomial  **  poly_1ds,
        int             *   orders,
        int                 npolys,
        cpl_size            degree) ;

cpl_polynomial * cr2res_wave_poly_2d_to_1d(
        cpl_polynomial  *   poly_2d,
        int                 order) ;

char * cr2res_wave_method_print(cr2res_wavecal_type wavecal_type) ;

cr2res_wavecal_type cr2res_wave_guess_method(const cpl_frame * in) ;

#endif
