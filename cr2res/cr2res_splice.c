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

#include <string.h>
#include <math.h>
#include <cpl.h>

#include "cr2res_splice.h"

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_splice      Splicing functions
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Main splicing function
  @param
  @param
  @param
  @param
  @return 
 */
/*----------------------------------------------------------------------------*/
int cr2res_splice_orders(
        cpl_table   *   trace_wave, 
        cpl_table   *   spectra, 
        int             trace)
{
    if (trace_wave == NULL | spectra == NULL) return -1;

    cpl_size i, j, k; 
    int nb_orders, *orders, central;
    int * loop0, *loop1;
    int iord0, iord1;
    int first0, first1, last0, last1, overlap0, overlap1;
    int n1 = CR2RES_DETECTOR_SIZE;
    double maximum = 0;
    double minW0, minW1, maxW0, maxW1;
    double wgt0, wgt1;
    // temporary work arrays
    char * colname;
    const cpl_array * wave1;
    cpl_vector ** spec, * s0, *s1, *tmpS0, *tmpS1;
    cpl_vector **uncs, *u0, *u1, *tmpU0, *tmpU1;
    cpl_vector **wave, *w0, *w1, *tmpW0, *tmpW1;
    cpl_vector **cont, *c0, *c1, *tmpC0, *tmpC1;

    double * tmp1;
    cpl_bivector * tmp2, * tmp3;
    cpl_vector * tmp4;
    double median;

    orders = cr2res_trace_get_order_numbers(trace_wave, &nb_orders);

    loop0 = cpl_malloc(nb_orders * sizeof(int));
    loop1 = cpl_malloc(nb_orders * sizeof(int));

    spec = cpl_malloc(nb_orders * sizeof(cpl_vector*));
    uncs = cpl_malloc(nb_orders * sizeof(cpl_vector*));
    wave = cpl_malloc(nb_orders * sizeof(cpl_vector*));
    cont = cpl_malloc(nb_orders * sizeof(cpl_vector*));

    // Load data into vector arrays
    for (i=0; i<nb_orders; i++){
        j  = cr2res_get_trace_table_index(trace_wave, orders[i], trace);
        // Get wavelength
        wave1 = cpl_table_get_array(trace_wave, CR2RES_COL_WAVELENGTH, j);
        tmp1 = cpl_array_get_data_double_const(wave1);
        wave[i] = cpl_vector_wrap(n1, tmp1);
        
        // Get Spectrum
        colname = cr2res_dfs_SPEC_colname(orders[i], trace);        
        tmp1 = cpl_table_get_data_double(spectra, colname);
        spec[i] = cpl_vector_wrap(n1, tmp1);
        cpl_free(colname);

        // Get Uncertainties
        colname = cr2res_dfs_SPEC_ERR_colname(orders[i], trace);
        tmp1 = cpl_table_get_data_double(spectra, colname);
        uncs[i] = cpl_vector_wrap(n1, tmp1);
        cpl_free(colname);

        // Get Continuum
        // ???

        // find order with largest signal
        median = cpl_vector_get_median(spec[i]);
        if (median > maximum){
            central = i;
            maximum = median;
        }

        // scale all orders to spec/cont = 1
        tmp4 = cpl_vector_new(n1);
        cpl_vector_copy(tmp4, spec[i]);
        cpl_vector_divide(tmp4, cont[i]);
        median = cpl_vector_get_median(tmp4);
        cpl_vector_multiply_scalar(cont[i], median);
    }

    // Loop over neighbouring orders, outwards from the central (stongest signal) order
    j = 0;
    for (i = central; i > 0; i--){
        loop0[j] = i;
        loop1[j] = i-1;
        j++;
    }
    for (i = central; i < nb_orders; i++){
        loop0[j] = i;
        loop1[j] = i+1;
        j++;
    }

    for (i=0; i < nb_orders; i++){
        // about the nomencalture
        // Order "0" is closer to central
        // Order "1" is the neighbour further away
        iord0 = loop0[i];
        iord1 = loop1[i];

        // Get relevant data vectors
        s0 = spec[iord0];
        s1 = spec[iord1];
        u0 = uncs[iord0];
        u1 = uncs[iord1];
        w0 = wave[iord0];
        w1 = wave[iord1];
        c0 = cont[iord0];
        c1 = cont[iord1];

        // Calculate overlap
        minW1 = cpl_vector_get_min(w1);
        minW0 = cpl_vector_get_min(w0);
        maxW1 = cpl_vector_get_max(w1);
        maxW0 = cpl_vector_get_max(w0);

        // indices of the first and last overlaping value
        first0 = 0;
        first1 = 0;
        last0 = 0;
        last1 = 0;

        // overlap0: number of overlap points of order "0"
        // overlap1: number of overlap points of order "1"

        overlap0 = 0;
        for (j = 0; j < cpl_vector_get_size(w0); j++){
            if ((cpl_vector_get(w0, j) >= minW1) & (cpl_vector_get(w0, j) <= maxW1)){
                // if its the first element
                if (first0 == 0) first0 = j;
                overlap0++;
                // update so that i1 will be the last element
                last0 = j;
            }
        }

        overlap1 = 0;
        for (j = 0; j < cpl_vector_get_size(w1); j++){
            if ((cpl_vector_get(w1, j) >= minW0) & (cpl_vector_get(w1, j) <= maxW0)){
                // if its the first element
                if (first1 == 0) first1 = j;
                overlap1++;
                // update so that i1 will be the last element
                last1 = j;
            }
        }

        // extract overlapping wavelength vectors
        // point to the latter parts of the wavelength vector, no need for copies
        tmpW0 = cpl_vector_wrap(last0 - first0, cpl_vector_get_data(w0) + first0);
        tmpW1 = cpl_vector_wrap(last1 - first1, cpl_vector_get_data(w1) + first1);
        // tmpW0 = cpl_vector_extract(w0, first0, last0, 1);
        // tmpW1 = cpl_vector_extract(w0, first0, last0, 1);

        // Prepare new vectors with interpolated data
        tmpS0 = cpl_vector_new(overlap0);
        tmpU0 = cpl_vector_new(overlap0);
        tmpC0 = cpl_vector_new(overlap0);

        tmpS1 = cpl_vector_new(overlap1);
        tmpU1 = cpl_vector_new(overlap1);
        tmpC1 = cpl_vector_new(overlap1);


        // Interpolate vectors on the other wavelength grid
        // Note that this is just linear interpolation
        tmp3 = cpl_bivector_wrap_vectors(s1, w1);
        tmp2 = cpl_bivector_wrap_vectors(tmpS0, tmpW0);
        cpl_bivector_interpolate_linear(tmp2, tmp3);
        cpl_bivector_unwrap_vectors(tmp3);
        cpl_bivector_unwrap_vectors(tmp2);

        tmp3 = cpl_bivector_wrap_vectors(s0, w0);
        tmp2 = cpl_bivector_wrap_vectors(tmpS1, tmpW1);
        cpl_bivector_interpolate_linear(tmp2, tmp3);
        cpl_bivector_unwrap_vectors(tmp3);
        cpl_bivector_unwrap_vectors(tmp2);

        tmp3 = cpl_bivector_wrap_vectors(u1, w1);
        tmp2 = cpl_bivector_wrap_vectors(tmpU0, tmpW0);
        cpl_bivector_interpolate_linear(tmp2, tmp3);
        cpl_bivector_unwrap_vectors(tmp3);
        cpl_bivector_unwrap_vectors(tmp2);

        tmp3 = cpl_bivector_wrap_vectors(u0, w0);
        tmp2 = cpl_bivector_wrap_vectors(tmpU1, tmpW1);
        cpl_bivector_interpolate_linear(tmp2, tmp3);
        cpl_bivector_unwrap_vectors(tmp3);
        cpl_bivector_unwrap_vectors(tmp2);

        tmp3 = cpl_bivector_wrap_vectors(c1, w1);
        tmp2 = cpl_bivector_wrap_vectors(tmpC0, tmpW0);
        cpl_bivector_interpolate_linear(tmp2, tmp3);
        cpl_bivector_unwrap_vectors(tmp3);
        cpl_bivector_unwrap_vectors(tmp2);

        tmp3 = cpl_bivector_wrap_vectors(c0, w0);
        tmp2 = cpl_bivector_wrap_vectors(tmpC1, tmpW1);
        cpl_bivector_interpolate_linear(tmp2, tmp3);
        cpl_bivector_unwrap_vectors(tmp3);
        cpl_bivector_unwrap_vectors(tmp2);

        // Combine Spectra with weighted average at each point
        // its important to remember the continuum in the uncertainties/weights
        k = 0;
        for (j=first0; j<=last0; j++){
            
            wgt0 = pow(cpl_vector_get(c0, j) / cpl_vector_get(u0, j), 2.);
            wgt1 = pow(cpl_vector_get(tmpC0, k) / cpl_vector_get(tmpU0, k), 2.);
            
            cpl_vector_set(s0, j, (cpl_vector_get(s0, j) * wgt0 
                    + cpl_vector_get(tmpS0, k) * wgt1) / (wgt0 + wgt1));

            cpl_vector_set(c0, j, (cpl_vector_get(c0, j) * wgt0 
                    + cpl_vector_get(tmpC0, k) * wgt1) / (wgt0 + wgt1));

            cpl_vector_set(u0, j, cpl_vector_get(c0, j) / sqrt(wgt0 + wgt1));
            
            k++;
        }

        k = 0;
        for (j=first1; j<=last1; j++){
            // Weighted average
            wgt0 = pow(cpl_vector_get(c1, j) / cpl_vector_get(u1, j), 2.);
            wgt1 = pow(cpl_vector_get(tmpC1, k) / cpl_vector_get(tmpU1, k), 2.);
            
            cpl_vector_set(s1, j, (cpl_vector_get(s1, j) * wgt0 
                    + cpl_vector_get(tmpS1, k) * wgt1) / (wgt0 + wgt1));

            cpl_vector_set(c1, j, (cpl_vector_get(c1, j) * wgt0 
                    + cpl_vector_get(tmpC1, k) * wgt1) / (wgt0 + wgt1));

            cpl_vector_set(u1, j, cpl_vector_get(c1, j) / sqrt(wgt0 + wgt1));
            
            k++;
        }

        // remove temporary vectors
        cpl_vector_delete(tmpS0);
        cpl_vector_delete(tmpU0);
        cpl_vector_delete(tmpC0);
        cpl_vector_unwrap(tmpW0);

        cpl_vector_delete(tmpS1);
        cpl_vector_delete(tmpU1);
        cpl_vector_delete(tmpC1);
        cpl_vector_unwrap(tmpW1);

    }
    
    // clean up pointers/data
    for (i=0; i< nb_orders; i++){
        cpl_vector_unwrap(spec[i]);
        cpl_vector_unwrap(uncs[i]);
        cpl_vector_unwrap(wave[i]);
        cpl_vector_unwrap(cont[i]);
    }

    cpl_free(spec);
    cpl_free(uncs);
    cpl_free(wave);
    cpl_free(cont);

    cpl_free(loop0);
    cpl_free(loop1);
    cpl_free(orders);
    return 0;
}

/**@}*/
