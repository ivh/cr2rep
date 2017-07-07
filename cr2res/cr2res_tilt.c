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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <math.h>
#include <string.h>

#include <cpl.h>

#include "cr2res_dfs.h"
#include "cr2res_io.h"
#include "cr2res_pfits.h"
#include "cr2res_utils.h"

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static cpl_polynomial * cr2res_tilt_compute_tilt_poly(
        cpl_polynomial  *   trace1_all_poly,
        cpl_polynomial  *   trace2_all_poly,
        cpl_polynomial  *   trace1_wave_poly,
        cpl_polynomial  *   trace2_wave_poly,
        int                 x_pos,
        double              slit_center_y) ;

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_tilt    Slit Tilt
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Compute the slit tilt polynomial for all x along the order
  @param    trace_wave  The trace wave table
  @param    order       The order that needѕ to be computed
  @param    slit_cen_y  The pixel Y position of the ѕlit center at x=1024
  @param    display     For Display
  @return   A list of CR2RES_DETECTOR_SIZE polynomials or NULL in error case

 */
/*----------------------------------------------------------------------------*/
cpl_polynomial ** cr2res_tilt(
        cpl_table           *   trace_wave,
        int                     order,
        double                  slit_cen_y,
        int                     display)
{
    cpl_polynomial  **  out_polys ;
    cpl_polynomial  *   all1 ;
    cpl_polynomial  *   all2 ;
    cpl_polynomial  *   wave1 ;
    cpl_polynomial  *   wave2 ;
    int                 trace_nb, trace_nb1, trace_nb2, nrows, cur_order ;
    int                 i ;

    /* Check Entries */
    if (trace_wave == NULL) return NULL ;

    /* Initialise */
    nrows = cpl_table_get_nrow(trace_wave) ;
    trace_nb1 = trace_nb2 = -1 ;

    /* Check that there are exactly 2 traces for this order */
    for (i=0 ; i<nrows ; i++) {
        cur_order = cpl_table_get(trace_wave, CR2RES_COL_ORDER, i, NULL) ;
        trace_nb = cpl_table_get(trace_wave, CR2RES_COL_TRACENB, i, NULL) ;
        if (order == cur_order) {
            /* Fill trace #1 */
            if (trace_nb1 < 0) {
                trace_nb1 = trace_nb ;
                continue ;
            }
            /* Fill trace #2 */
            if (trace_nb2 < 0) {
                trace_nb2 = trace_nb ;
                continue ;
            }
            /* We encounter here the 3rd trace - Exit for the moment */
            cpl_msg_error(__func__, "Third trace found for this Order"); 
            return NULL ;
        }
    }

    /* We need 2 traces to compute the tilt */
    if (trace_nb1 < 0 || trace_nb2 < 0) {
        cpl_msg_error(__func__, "Two traces are needed to compute the tilt"); 
        return NULL ;
    }

    /* Check that the wavelength polynomials are available */
    wave1 = cr2res_get_trace_wave_poly(trace_wave, CR2RES_COL_WAVELENGTH, order,
            trace_nb1) ;
    wave2 = cr2res_get_trace_wave_poly(trace_wave, CR2RES_COL_WAVELENGTH, order,
            trace_nb2) ;
    if (wave1 == NULL || wave2 == NULL) {
        cpl_msg_error(__func__, "The wavelength solution is missing"); 
        if (wave1 != NULL) cpl_polynomial_delete(wave1) ;
        if (wave2 != NULL) cpl_polynomial_delete(wave2) ;
        return NULL ;
    }
    all1 = cr2res_get_trace_wave_poly(trace_wave, CR2RES_COL_ALL, order,
            trace_nb1) ;
    all2 = cr2res_get_trace_wave_poly(trace_wave, CR2RES_COL_ALL, order,
            trace_nb2) ;
    if (all1 == NULL || all2 == NULL) {
        cpl_msg_error(__func__, "The trace is missing"); 
        if (all1 != NULL) cpl_polynomial_delete(all1) ;
        if (all2 != NULL) cpl_polynomial_delete(all2) ;
        cpl_polynomial_delete(wave1) ;
        cpl_polynomial_delete(wave2) ;
        return NULL ;
    }

    /* Create Output polynomials array */
    out_polys = cpl_calloc(CR2RES_DETECTOR_SIZE, sizeof(cpl_polynomial *)) ;

    /* Loop on the X coordinates and get the slit tilt for this X */
    for (i=0 ; i<CR2RES_DETECTOR_SIZE ; i++) {
        out_polys[i] = cr2res_tilt_compute_tilt_poly(all1, all2, wave1, wave2, 
                i+1, slit_cen_y);
    }

    /* Free and return */
    cpl_polynomial_delete(all1) ;
    cpl_polynomial_delete(all2) ;
    cpl_polynomial_delete(wave1) ;
    cpl_polynomial_delete(wave2) ;
    return out_polys ;
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Compute the slit tilt polynomial at a given position
  @param    y0_mid      Y pos of the slit at x=1024 (middle)
  @return   A slit tilt polynomial or NULL in error case
  Tilt Polynomial: Y=poly_xpos(X) = a.X+B
    X : pixel position along y (1 at the lower row)
    Y : pixel position along x (1 at the lower column)
 */
/*----------------------------------------------------------------------------*/
static cpl_polynomial * cr2res_tilt_compute_tilt_poly(
        cpl_polynomial  *   trace1_all_poly,
        cpl_polynomial  *   trace2_all_poly,
        cpl_polynomial  *   trace1_wave_poly,
        cpl_polynomial  *   trace2_wave_poly,
        int                 xpos,
        double              y0_mid)
{
    cpl_polynomial  *   poly_xpos ;
    cpl_polynomial  *   poly_tmp ;
    int                 x_mid ;
    double              y1_mid, y2_mid, y1_xpos, y2_xpos, factor, y0_xpos, 
                        w1_xpos, x2, y2, coeff, a, b ;
    cpl_size            degree ;
    int                 x_to_debug ;


    /* Check Entries */
    if (trace1_wave_poly == NULL || trace2_wave_poly == NULL) return NULL ;
    if (trace1_all_poly == NULL || trace2_all_poly == NULL) return NULL ;
    if (xpos < 1 || xpos > CR2RES_DETECTOR_SIZE) return NULL ;

    /* Initialise */
    x_mid = (int)(CR2RES_DETECTOR_SIZE/2.0) ;
    x_to_debug = -1 ;
 
    /* Get the Y slit center at the xpos position : y0_xpos */
    y1_mid = cpl_polynomial_eval_1d(trace1_all_poly, x_mid, NULL);
    y2_mid = cpl_polynomial_eval_1d(trace2_all_poly, x_mid, NULL);
    y1_xpos = cpl_polynomial_eval_1d(trace1_all_poly, xpos, NULL) ;
    y2_xpos = cpl_polynomial_eval_1d(trace2_all_poly, xpos, NULL) ;

    if (fabs(y2_mid-y0_mid) < 1e-3) {
        y0_xpos = y1_xpos - (y1_mid - y0_mid) ;
    } else {
        factor = (y1_mid-y0_mid)/(y2_mid-y0_mid) ;
        if (fabs(factor-1.0) < 1e-3) {
            y0_xpos = (y1_xpos + y2_xpos)/2.0 ;
        } else {
            y0_xpos = (y2_xpos*factor - y1_xpos)/(factor-1) ;
        }
    }
    if (xpos==x_to_debug) cpl_msg_debug(__func__, "y0_xpos=%g", y0_xpos) ;

    /* Get the Wavelength on trace1 : w1_xpos */
    w1_xpos = cpl_polynomial_eval_1d(trace1_wave_poly, xpos, NULL) ;
    if (xpos==x_to_debug) cpl_msg_debug(__func__, "w1_xpos=%g", w1_xpos) ;

    /* Get the Point on trace 2 with wavelength=w1_xpos : (x2, y2) */
    /* Solve trace2_wave_poly(x) = w1_xpos */
    /* poly_tmp = trace2_wave_poly() - w1_xpos -> Solve poly_tmp(x) = 0 */
    degree = 0 ;
    poly_tmp = cpl_polynomial_duplicate(trace2_wave_poly);
    coeff = cpl_polynomial_get_coeff(poly_tmp, &degree) - w1_xpos ;
    cpl_polynomial_set_coeff(poly_tmp, &degree, coeff) ;
    cpl_polynomial_solve_1d(poly_tmp, (double)xpos, &x2, 1) ; 
    cpl_polynomial_delete(poly_tmp) ;
    y2 = cpl_polynomial_eval_1d(trace2_all_poly, x2, NULL) ;
    if (xpos==x_to_debug) cpl_msg_debug(__func__, "(x2,y2)=(%g,%g)", x2, y2) ;
    
    /* Check the edges */
    if (x2 < 1 || x2 > CR2RES_DETECTOR_SIZE) return NULL ;

    /* Compute the Polynomial */
	/* (xpos, y1_xpos) and (x2, y2) for the slope a */
	/* poly_xpos(xpos) = y0_xpos for b */
    
    /* Slit is vertical */
    if (fabs(x2-xpos) < 1e-3) return NULL ;
    a = (y2-y1_xpos) / (x2-xpos) ;
    b = y0_xpos - a * xpos ;
     
    /* Create the polynomial */
    poly_xpos = cpl_polynomial_new(1) ;
    degree = 0 ;
    cpl_polynomial_set_coeff(poly_xpos, &degree, b) ;
    degree = 1 ;
    cpl_polynomial_set_coeff(poly_xpos, &degree, a) ;
    if (xpos==x_to_debug && cpl_msg_get_level() == CPL_MSG_DEBUG) 
        cpl_polynomial_dump(poly_xpos, stdout) ;

    /* Free and return */
    return poly_xpos ;
}

