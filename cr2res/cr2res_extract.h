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

#ifndef CR2RES_EXTRACT_H
#define CR2RES_EXTRACT_H

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>
#include "hdrl.h"

/*-----------------------------------------------------------------------------
                                   Define
 -----------------------------------------------------------------------------*/

typedef enum {
    CR2RES_EXTR_SUM,
    CR2RES_EXTR_MEDIAN,
    CR2RES_EXTR_TILTSUM,
    CR2RES_EXTR_OPT_CURV,
} cr2res_extr_method ;

/*-----------------------------------------------------------------------------
                                       Prototypes
 -----------------------------------------------------------------------------*/

int cr2res_extract_traces(
        const hdrl_image    *   img,
        const cpl_table     *   traces,
        const cpl_table     *   slit_func_in,
        const cpl_table     *   blaze_table_in,
        int                     reduce_order,
        int                     reduce_trace,
        cr2res_extr_method      extr_method,
        int                     extr_height,
        int                     swath_width,
        int                     oversample,
        double                  smooth_slit,
        double                  smooth_spec,
        int                     niter,
        double                  kappa,
        double                  gain,
        int                     display,
        int                     disp_order_idx,
        int                     disp_trace,
        cpl_table           **  extracted,
        cpl_table           **  slit_func,
        hdrl_image          **  model_master) ;

int cr2res_extract_sum_vert(
        const hdrl_image    *   hdrl_in,
        const cpl_table     *   trace_tab,
        int                     order,
        int                     trace_id,
        int                     height,
        cpl_vector          **  slit_func,
        cpl_bivector        **  spec,
        hdrl_image          **  model) ;

int cr2res_extract_median(
        const hdrl_image    *   hdrl_in,
        const cpl_table     *   trace_tab,
        int                     order,
        int                     trace_id,
        int                     height,
        cpl_vector          **  slit_func,
        cpl_bivector        **  spec,
        hdrl_image          **  model) ;

int cr2res_extract_sum_tilt(
        const hdrl_image    *   hdrl_in,
        const cpl_table     *   trace_tab,
        int                     order,
        int                     trace_id,
        int                     height,
        cpl_vector          **  slit_func,
        cpl_bivector        **  spec,
        hdrl_image          **  model) ;

int cr2res_extract_slitdec_curved(
        const hdrl_image    *   img_hdrl,
        const cpl_table     *   trace_tab,
        const cpl_vector    *   slit_func_vec_in,
        int                     order,
        int                     trace_id,
        int                     height,
        int                     swath,
        int                     oversample,
        double                  smooth_slit,
        double                  smooth_spec,
        int                     niter,
        double                  kappa,
        double                  gain,
        cpl_vector          **  slit_func,
        cpl_bivector        **  spec,
        hdrl_image          **  model) ;

cpl_table * cr2res_extract_SLITFUNC_create(
        cpl_vector      **  slit_func,
        const cpl_table *   trace_table) ; 

cpl_table * cr2res_extract_EXTRACT1D_create(
        cpl_bivector    **  spectrum,
        const cpl_table *   trace_table) ;

int cr2res_extract_EXTRACT1D_get_spectrum(
        const cpl_table     *   tab,
        int                     order,
        int                     trace_nb,
        cpl_bivector        **  spec,
        cpl_bivector        **  spec_err) ;

int cr2res_extract_SLIT_FUNC_get_vector(
        const cpl_table     *   tab,
        int                     order,
        int                     trace_nb,
        cpl_vector          **  vec) ;

int cr2res_extract2d_traces(
        const hdrl_image    *   img,
        const cpl_table     *   traces,
        int                     reduce_order,
        int                     reduce_trace,
        cpl_table           **  extracted) ;

int cr2res_extract2d_trace(
        const hdrl_image    *   in,
        const cpl_table     *   trace_tab,
        int                     order,
        int                     trace_id,
        int                     npoints,
        const cpl_image     *   wavemap,
        const cpl_image     *   slitmap,
        cpl_bivector        **  spectrum,
        cpl_bivector        **  position,
        cpl_vector          **  wavelength,
        cpl_vector          **  slit_fraction) ;

cpl_table * cr2res_extract_EXTRACT2D_create(
        cpl_bivector    **  spectrum,
        cpl_bivector    **  position,
        cpl_vector      **  wavelength,
        cpl_vector      **  slit_fraction,
        const cpl_table *   trace_table) ;

int cr2res_extract_slitdec_bandsol(double *, double *, int, int, double) ;

#endif
