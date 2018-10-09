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

#include "cr2res_slit_curv.h"
#include "cr2res_dfs.h"
#include "cr2res_io.h"
#include "cr2res_pfits.h"
#include "cr2res_utils.h"
#include "cr2res_trace.h"

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_slit_curv_get_position(
        cpl_polynomial  *   trace,
        cpl_polynomial  *   wave,
        double              ref_wl,
        double          *   xpos,
        double          *   ypos) ;

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_slit_curv   Slit Curvature
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief Compute the slit curvature along a specific trace
  @param trace_wave     The trace_wave table
  @param order          The order number
  @param trace_id       The trace_id number
  @param max_curv_degree    The maximum degree for the slit curv polynomial
  @param display        The flag used to display something
  @return   An array of polynomials or NULL in error case

  order and trace_id identify the trace along which the slit curve is computed.
  For each x along the trace, the slit curvature Px is computed.
  The wavelength ref_wl at the x position along the input trace is evaluated.
  All other traces in the current order are used.
  On each of these traces, the point (x,Y) on the trace that has the same 
    wavelength as ref_wl is identified. 
  All these points are along the slit and used to fit its shape.

  The returned array contains potentially CR2RES_DETECTOR_SIZE
  polynomials (1 for each x position in the detector). Each non-NULL 
  polynomial needs to be deallocated (cpl_polynomial_delete()) and the
  returned pointer needs to be freed (cpl_free()).

  The polynomial number x Px gives the slit curvature:
        X=Px(Y)
  where (X,Y) are the detector position, (1, 1) being the lower left pixel.
        Px is the polynomial to use along x
 */
/*----------------------------------------------------------------------------*/
cpl_polynomial ** cr2res_slit_curv_compute_order_trace(
        cpl_table           *   trace_wave,
        int                     order,
        int                     trace_id,
        int                     max_curv_degree,
        int                     display)
{
    cpl_polynomial  **  out_polys ;
    cpl_polynomial  *   trace_in ;
    cpl_polynomial  *   wave_poly_in ;
    cpl_polynomial  *   cur_trace_in ;
    cpl_polynomial  *   cur_wave_poly_in ;
    double              ref_wl, ref_y, cur_x, cur_y ;
    cpl_matrix      *   x_points ;
    cpl_vector      *   y_points ;
    cpl_polynomial  *   slit_curve ;
    int                 cur_order, cur_trace_id, abort_fit ;
    cpl_size            ntraces, ref_x, power, poly_degree, idx, i ;

    /* Check Entries */
    if (trace_wave == NULL) return NULL ;

    /* Get the number of traces */
    ntraces = cr2res_get_nb_traces_with_wavelength(trace_wave, order) ;
    if (ntraces < 2) return NULL ;

    /* Set the fitting polynomial degree */
    poly_degree = ntraces - 1 ;
    if (poly_degree > max_curv_degree) poly_degree = max_curv_degree ;

    /* Get the input trace wavelength */
    wave_poly_in = cr2res_get_trace_wave_poly(trace_wave, 
            CR2RES_COL_WAVELENGTH, order, trace_id) ;
    if (wave_poly_in == NULL) return NULL ;

    /* Get the input trace polynomial */
    trace_in = cr2res_get_trace_wave_poly(trace_wave, CR2RES_COL_ALL, order, 
            trace_id) ;

    /* Allocate Output array */
    out_polys = cpl_malloc(CR2RES_DETECTOR_SIZE * sizeof(cpl_polynomial*)) ;

    /* Loop on all x positions */
    for (ref_x=1 ; ref_x<=CR2RES_DETECTOR_SIZE ; ref_x++) {

        /* Get the Wavelength on the input trace */
        ref_wl = cpl_polynomial_eval_1d(wave_poly_in, (double)ref_x, NULL) ;

        /* Create Objects */
        x_points = cpl_matrix_new(1, ntraces) ;
        y_points = cpl_vector_new(ntraces) ;

		/* Store the reference point */
		idx = 0 ;
		ref_y = cpl_polynomial_eval_1d(trace_in, (double)ref_x, NULL) ;
        cpl_matrix_set(x_points, 0, idx, ref_y) ;
        cpl_vector_set(y_points, idx, (double)ref_x) ;

        /* Store the other traces points */
        abort_fit = 0 ;
        for (i=0 ; i<cpl_table_get_nrow(trace_wave) ; i++) {
            /* Find the other traces of the same order */
            cur_order = cpl_table_get(trace_wave, CR2RES_COL_ORDER, i, NULL) ;
            cur_trace_id = cpl_table_get(trace_wave, CR2RES_COL_TRACENB,i,NULL);
            cur_wave_poly_in = cr2res_get_trace_wave_poly(trace_wave, 
                    CR2RES_COL_WAVELENGTH, order, cur_trace_id) ;

            /* Search for the next trace in the current order */
            if (cur_order != order || cur_trace_id == trace_id ||
                    cur_wave_poly_in == NULL) {
                cpl_polynomial_delete(cur_wave_poly_in) ;
                continue ;
            }

            /* The current trace is used for the fit */
            idx++ ;

            /* Get the position on the current trace with a given wl */
            cur_trace_in=cr2res_get_trace_wave_poly(trace_wave,CR2RES_COL_ALL, 
                    order, cur_trace_id) ;
            cur_x = cur_y = 1.0 ;
            if (cr2res_slit_curv_get_position(cur_trace_in, cur_wave_poly_in, 
                    ref_wl, &cur_x, &cur_y) == -1) {
                /* Cannot get the position for this wavelength - abort */
                abort_fit = 1 ;
            }
            cpl_polynomial_delete(cur_trace_in) ;
            cpl_polynomial_delete(cur_wave_poly_in) ;
            cpl_matrix_set(x_points, 0, idx, cur_y) ;
            cpl_vector_set(y_points, idx, cur_x) ;
        }

        /* Compute the fit */
        slit_curve = cpl_polynomial_new(1) ;

        if (abort_fit) {
            cpl_msg_debug(__func__, "Abort the fit for x=%"CPL_SIZE_FORMAT, 
                    ref_x) ;
            /* Use the vertical Slit Polynomial */
            power = 0 ;
            cpl_polynomial_set_coeff(slit_curve, &power, (double)ref_x) ;
        } else if (cpl_polynomial_fit(slit_curve, x_points, NULL, y_points, 
                    NULL, CPL_FALSE, NULL, &poly_degree) != CPL_ERROR_NONE) {
            cpl_msg_warning(__func__, "Cannot fit the slit") ;
            cpl_polynomial_delete(slit_curve) ;

            /* Use the vertical Slit Polynomial */
            slit_curve = cpl_polynomial_new(1) ;
            power = 0 ;
            cpl_polynomial_set_coeff(slit_curve, &power, (double)ref_x) ;
        }
        /* Clean */
        cpl_matrix_delete(x_points) ;
        cpl_vector_delete(y_points) ;

        /* Store result */
        out_polys[ref_x-1] = slit_curve ;
    }
    cpl_polynomial_delete(trace_in) ;
    cpl_polynomial_delete(wave_poly_in) ;
    return out_polys ;
}

/* TODO */
/*----------------------------------------------------------------------------*/
/**
  @brief 
  @param 
  @return  
 */
/*----------------------------------------------------------------------------*/
int cr2res_slit_curv_fit_coefficients(
        cpl_polynomial  **  curvatures,
        int                 nb_polys,
        cpl_polynomial  **  slit_polya,
        cpl_polynomial  **  slit_polyb,
        cpl_polynomial  **  slit_polyc)
{
    cpl_polynomial  *   coeffs_fit ;
    cpl_size            poly_degree ;
    cpl_matrix      *   x_points ;
    cpl_vector      *   y_points_a ;
    cpl_vector      *   y_points_b ;
    cpl_vector      *   y_points_c ;
    double              coeff ;
    cpl_size            i, power ;

    /* Check Entries */
    if (curvatures ==  NULL || slit_polya == NULL || slit_polyb == NULL
            || slit_polyc == NULL) return -1 ;
    
    /* Initialise */
    poly_degree = 2 ;

    /* Create Objects */
    x_points = cpl_matrix_new(1, nb_polys) ;
    y_points_a = cpl_vector_new(nb_polys) ;
    y_points_b = cpl_vector_new(nb_polys) ;
    y_points_c = cpl_vector_new(nb_polys) ;

    /* Loop on all polynomials */
    for (i=0 ; i<nb_polys ; i++) {
        cpl_matrix_set(x_points, 0, i, (double)(i+1)) ;

        /* Retrieve the A coeeficient */
        power = 0 ;
        coeff = cpl_polynomial_get_coeff(curvatures[i], &power) ;
        cpl_vector_set(y_points_a, i, coeff) ;
        /* Retrieve the B coeeficient */
        power = 1 ;
        coeff = cpl_polynomial_get_coeff(curvatures[i], &power) ;
        cpl_vector_set(y_points_b, i, coeff) ;
        /* Retrieve the C coeeficient */
        power = 2 ;
        coeff = cpl_polynomial_get_coeff(curvatures[i], &power) ;
        cpl_vector_set(y_points_c, i, coeff) ;
    }

    /* Create Objets */
    *slit_polya = cpl_polynomial_new(1) ;
    power = 0 ; cpl_polynomial_set_coeff(*slit_polya, &power, 0.0) ;
    power = 1 ; cpl_polynomial_set_coeff(*slit_polya, &power, 1.0) ;
    power = 2 ; cpl_polynomial_set_coeff(*slit_polya, &power, 0.0) ;
    *slit_polyb = cpl_polynomial_duplicate(*slit_polya) ;
    *slit_polyc = cpl_polynomial_duplicate(*slit_polya) ;
    
    /* Compute the fit for a coefficients */
    cpl_polynomial_fit(*slit_polya, x_points, NULL, y_points_a, 
                    NULL, CPL_FALSE, NULL, &poly_degree) ;
    cpl_vector_delete(y_points_a) ;

    /* Compute the fit for b coefficients */
    cpl_polynomial_fit(*slit_polyb, x_points, NULL, y_points_b, 
                    NULL, CPL_FALSE, NULL, &poly_degree) ;
    cpl_vector_delete(y_points_b) ;
    
    /* Compute the fit for c coefficients */
    cpl_polynomial_fit(*slit_polyc, x_points, NULL, y_points_c, 
                    NULL, CPL_FALSE, NULL, &poly_degree) ;
    cpl_vector_delete(y_points_c) ;
    cpl_matrix_delete(x_points) ;
    return 0 ;
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief Get the (X,Y) image cordinate of a Wavelength on a trace
  @param trace  The trace center position polynomial on which we search
  @param wave   The wavelength polynomial of the trace
  @param ref_wl The wavelength value
  @param xpos   [out] The X position of the trace point with the ref_wl
  @param ypos   [out] The Y position of the trace point with the ref_wl
  @return   0 if ok, -1 in error case

  The exact X position on the trace for which the wl is ref_wl is found
  Y is computed using X and the trace polynomial.
 */
/*----------------------------------------------------------------------------*/
static int cr2res_slit_curv_get_position(
        cpl_polynomial  *   trace,
        cpl_polynomial  *   wave,
        double              ref_wl,
        double          *   xpos,
        double          *   ypos)
{
    double          cur_wl, tmp_wl ;
    cpl_size        x ;

    /* Check entries */
    if (trace == NULL || wave == NULL || xpos == NULL || ypos == NULL) 
        return -1 ;

    /* Loop on the x positions  */
    for (x=0 ; x<=CR2RES_DETECTOR_SIZE+1 ; x++) {
        cur_wl = cpl_polynomial_eval_1d(wave, (double)x, NULL) ;
        /* As soon as the WL is bigger than ref_wl, we keep X */
        if (cur_wl > ref_wl) break ;
    }

    tmp_wl = cpl_polynomial_eval_1d(wave, (double)(x-1), NULL) ;

    /* Linear interpolation */
    if (fabs(cur_wl-tmp_wl) < 1e-5)
        *xpos = (double)x ;
    else 
        *xpos = (double)(x -((cur_wl-ref_wl) / (cur_wl-tmp_wl))) ;

    *ypos = cpl_polynomial_eval_1d(trace, *xpos, NULL) ;
    /* Check results */
    if (*xpos < 1 || *xpos > CR2RES_DETECTOR_SIZE ||
            *ypos < 1 || *ypos > CR2RES_DETECTOR_SIZE) {
        *xpos = *ypos = -1 ;
        return -1 ;
    }
    return 0 ;
}

