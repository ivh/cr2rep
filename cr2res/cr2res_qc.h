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

#ifndef CR2RES_QC_H
#define CR2RES_QC_H

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>

/*-----------------------------------------------------------------------------
                           Functions prototypes
 -----------------------------------------------------------------------------*/

/* DARK */
double cr2res_dark_qc_ron(
        const cpl_image     *   ima1,
        const cpl_image     *   ima2,
        int                     hsize,
        int                     nsamples,
        int                     ndit) ;

/* DETLIN */
double cr2res_qc_detlin_median(
        const cpl_imagelist     *   coeffs) ;
double cr2res_qc_detlin_gain(
        const cpl_imagelist     *   coeffs) ;
int cr2res_qc_detlin_min_max_level(
        const cpl_image     *   ima,
        double              *   min_level,
        double              *   max_level) ;

/* FLAT */
double cr2res_qc_flat_lamp_ints(
        const cpl_image     *   ima) ;
double cr2res_qc_flat_mean_level(
        const cpl_image     *   ima) ;
double cr2res_qc_flat_med_snr(
        const cpl_image     *   ima) ;
int cr2res_qc_flat_mean_med_flux(
        const cpl_image     *   ima,
        double              *   mean_flux,
        double              *   med_flux) ;
double cr2res_qc_flat_trace_center_y(
        const cpl_table     *   trace) ;
int cr2res_qc_flat_nb_overexposed(
        const cpl_image     *   ima) ;
int cr2res_qc_flat_order_positions(
        const cpl_table *   tw,
        int             **  order_nb,
        double          **  order_pos,
        int             *   nbvals) ;

/* WAVE */




/* OBS_1D */
double cr2res_qc_obs_nodding_signal(
        const cpl_table     *   extracted) ;
double cr2res_qc_obs_nodding_transmission(
        const cpl_table     *   extracted) ;
double cr2res_qc_obs_nodding_slit_psf(
        const cpl_table     *   slitfu);

#endif

