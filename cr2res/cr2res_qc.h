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
#include <hdrl.h>

#define CR2RES_NONLIN_LEVEL 20000
#define CR2RES_QC_ORDER 4
#define CR2RES_QC_TRACE 1
#define CR2RES_QC_SIZE  100
#define CR2RES_QC_WINDOW 20


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
double cr2res_qc_detlin(
        const hdrl_imagelist    *   hdrl_coeffs,
        double                      bpm_thresh,
        cpl_mask               **  outmask,
        double                  *   min_level,
        double                  *   max_level) ;
int cr2res_qc_detlin_stat(
        const hdrl_imagelist    *   hdrl_coeffs,
        double                  *   meda,
        double                  *   medb,
        double                  *   medc,
        double                  *   meda_err) ;

/* FLAT */
double cr2res_qc_flat_trace_center_y(
        const cpl_table     *   trace) ;
int cr2res_qc_flat_order_positions(
        const cpl_table *   tw,
        int             **  order_nb,
        double          **  order_pos,
        int             *   nbvals) ;
double cr2res_qc_flat_s2n(
        const cpl_table     *   extracted) ;

/* WAVE */
double cr2res_qc_wave_central(
        const cpl_table *   tw,
        int                 order) ;
double cr2res_qc_wave_disp(
        const cpl_table *   tw,
        int                 order) ;

cpl_vector * cr2res_qc_lines_collect(double wmin, double wmax) ;

double cr2res_qc_wave_lamp_effic(
        const cpl_bivector  *   spec) ;
cpl_bivector * cr2res_qc_lines_intens_bgd(
        const cpl_bivector  *   spec) ;

double cr2res_qc_wave_line_fwhm(
        const cpl_bivector  *   spec,
        double                  wl,
        double              *   peak_height) ;
double cr2res_qc_wave_resol_fwhm(
        const cpl_bivector  *   spec,
        double              *   wl) ;

/* OBS */
int cr2res_qc_numsat(const cpl_frameset * frameset) ;
double cr2res_qc_overexposed(
        const cpl_image *   ima,
        const cpl_table *   tw,
        int                 order_idx) ;
double cr2res_qc_obs_nodding_signal(
        const cpl_table     *   extracted) ;
double cr2res_qc_obs_nodding_standard_flux(
        const cpl_table     *   extracted,
        char                *   setting) ;
double cr2res_qc_obs_slit_psf(
        const cpl_table     *   slitfu,
        int                     order_idxp,
        int                     oversample) ;
double * cr2res_qc_snr(
    const cpl_table *   tw,
    const cpl_table *   extracted,
    int             **  out_order_idx_values,
    int             *   out_nb_order_idx_values) ;
double * cr2res_qc_der_snr(
    const cpl_table *   tw,
    const cpl_table *   extracted,
    int             **  out_order_idx_values,
    int             *   out_nb_order_idx_values) ;
double cr2res_qc_compute_snr(cpl_vector *spec,
                             cpl_vector *err);
double cr2res_qc_compute_der_snr(cpl_vector *spec,
                             cpl_vector *err);

int cr2res_qc_calculate_mean_and_rmsd(
        cpl_propertylist    ***      plists, 
        int                         size,
        const int                 *       sizes, 
        const char          *       ref_keyword,
        cpl_propertylist    *       qc_main,
        const char          *       result_avg_keyword,
        const char          *       result_rmsd_keyword);

int cr2res_qc_dup_mtrlgy_key(
        cpl_frameset * framelist, 
        cpl_propertylist *plist);

int cr2res_qc_dup_chip_idx(
        cpl_propertylist *plist);
#endif
