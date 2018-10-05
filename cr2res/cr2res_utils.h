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

#ifndef CR2RES_UTILS_H
#define CR2RES_UTILS_H

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>
#include <hdrl.h>

/*-----------------------------------------------------------------------------
                                   Define
 -----------------------------------------------------------------------------*/

#define CR2RES_NB_DETECTORS             3
#define CR2RES_NB_DECKER_POSITIONS      3
#define CR2RES_DETECTOR_SIZE            2048

typedef enum {
    CR2RES_LAMP,
    CR2RES_GAS,
    CR2RES_ETALON,
    CR2RES_UNSPECIFIED
} cr2res_wavecal_type ;

typedef enum {
    CR2RES_DECKER_INVALID,
    CR2RES_DECKER_NONE,
    CR2RES_DECKER_1_3,
    CR2RES_DECKER_2_4
} cr2res_decker ;

typedef enum {
    CR2RES_EXTR_SUM,
    CR2RES_EXTR_OPT_VERT,
    CR2RES_EXTR_OPT_CURV,
} cr2res_extr_method ; 

/*-----------------------------------------------------------------------------
                                       Prototypes
 -----------------------------------------------------------------------------*/

double * cr2res_vector_get_rest(
    const cpl_vector    * ycen);
int * cr2res_vector_get_int(
    const cpl_vector    * ycen);

cpl_image * cr2res_image_cut_rectify(
        const cpl_image     * img_in,
        const cpl_vector    * ycen,
        int                   height);
int cr2res_image_insert_rect(
        const cpl_image     * rect_in,
        const cpl_vector    * ycen,
        cpl_image           * img_out  );
cpl_vector * cr2res_polynomial_eval_vector(
        const cpl_polynomial * poly,
        const cpl_vector     * vec);

cpl_vector * cr2res_threshold_spec(const cpl_vector * invector, int smooth, 
                                double thresh) ;

char * cr2res_get_base_name(const char * filename) ;
char * cr2res_get_root_name(const char * filename) ;

const char * cr2res_extract_filename(const cpl_frameset *, const char *) ;
cpl_frameset * cr2res_extract_frameset(const cpl_frameset *, const char *) ;
cpl_frameset * cr2res_extract_decker_frameset(
        const cpl_frameset  *   in,
        const char          *   tag,
        cr2res_decker          	decker) ;

cpl_polynomial * cr2res_wlestimate_compute(
        double          wmin,
        double          wmax) ;

int cr2res_convert_order_to_idx(int order) ;
int cr2res_convert_idx_to_order(int order_idx) ;

cpl_polynomial * cr2res_convert_array_to_poly(const cpl_array * arr) ;
cpl_array * cr2res_convert_poly_to_array(
        const cpl_polynomial    *   poly,
        int                         size) ;

cpl_error_code cr2res_detector_shotnoise_model(
        const cpl_image *   ima_data,
        const double        gain,
        const double        ron,
        cpl_image       **  ima_errs) ;

cpl_table * cr2res_demod(hdrl_image *sum1, hdrl_image *sum2, cpl_table *trace_wave);

cpl_polynomial * cr2res_fit_noise(cpl_image *img, cpl_table *trace_wave, cpl_size order_x, cpl_size order_y);

int cr2res_slit_pos(cpl_table *tw_decker1, cpl_table *tw_decker2, cpl_polynomial **coef_slit, cpl_polynomial **coef_wave);
int cr2res_slit_pos_image(cpl_table *tw_decker1, cpl_table *tw_decker2, cpl_image *slitpos, cpl_image *wavelength);

int cr2res_splice_orders(cpl_table *trace_wave, cpl_table * spectra, int trace);

const char * cr2res_get_license(void) ;

#endif
