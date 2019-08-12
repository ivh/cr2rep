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

#include "cr2res_dfs.h"
#include "cr2res_utils.h"
#include "cr2res_trace.h"
#include "cr2res_splice.h"

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_splice      Splicing functions
 */
/*----------------------------------------------------------------------------*/

int cr2res_combine_spectra(cpl_bivector * spliced[],
    cpl_bivector * spliced_err[],
    cpl_vector * spectrum_order,
    cpl_bivector * first,
    cpl_bivector * last,
    int nspectra,
    cpl_bivector ** spectrum,
    cpl_bivector ** spectrum_err);

int cr2res_extract_data(
    cpl_table   *   trace_wave,
    cpl_table   *   blaze, 
    cpl_table   *   spectra, 
    int             order, 
    int             trace,
    cpl_vector **   wave,
    cpl_vector **   spec,
    cpl_vector **   uncs,
    cpl_vector **   cont);

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Main splicing function
  @param    extracted   List of extracted tables
  @param    blaze       List of blaze tables
  @param    trace_wave  List of trace_wave
  @param    ninputs     Size of the input lists
  @param    spliced     [out] Spliced spectrum
  @param    spliced_err [out] Spliced spectrum error
  @return   0 if ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
int cr2res_splice(
        cpl_table       **  extracted_1d,
        cpl_table       **  blaze,
        cpl_table       **  trace_wave,
        int                 ninputs,
        cpl_bivector    **  spliced,
        cpl_bivector    **  spliced_err)
{
    int i, j, k, trace, nspectra, nb_orders, nb_traces, sucess;
    int * orders, * traces;
    int count = 0;
    // trace determines which trace to get from the table
    // trace == -1, means all traces
    trace = -1;

    cpl_bivector ** sp;
    cpl_bivector ** sp_err;
    cpl_vector * sp_order;
    cpl_bivector * first;
    cpl_bivector * last;

    cpl_vector ** wave;
    cpl_vector ** spec;
    cpl_vector ** uncs;
    cpl_vector ** cont;


    // Count total number of spectra
    nspectra = 0;
    for (i = 0; i < ninputs; i++){
        if (trace == -1){
            // Count all traces
            nspectra += cpl_table_get_nrow(trace_wave[i]);
        }
        else{
            // Count number of orders with the requested trace
            orders = cr2res_trace_get_order_numbers(trace_wave[i], &nb_orders);
            for (j = 0; j < nb_orders; j++){
                if ((k = cr2res_get_trace_table_index(trace_wave[i], orders[j],
                                                        trace)) != -1) {
                    nspectra++;
                }
            }
            cpl_free(orders);
        }
    }

    cpl_msg_info(__func__, "Number of segments: %i", nspectra);


    // Prepare data vectors
    wave = cpl_malloc(nspectra * sizeof(cpl_vector *));
    spec = cpl_malloc(nspectra * sizeof(cpl_vector *));
    uncs = cpl_malloc(nspectra * sizeof(cpl_vector *));
    cont = cpl_malloc(nspectra * sizeof(cpl_vector *));

    // fill data vectors by reading table data
    count = 0;
    for (i = 0; i < ninputs; i++){
        orders = cr2res_trace_get_order_numbers(trace_wave[i], &nb_orders);
        for (j = 0; j < nb_orders; j++){
            if (trace == -1){
                traces = cr2res_get_trace_numbers(trace_wave[i], orders[j],
                                                    &nb_traces);
                for (k = 0; k < nb_traces; k++) {
                    sucess = cr2res_extract_data(trace_wave[i], blaze[i],
                        extracted_1d[i], orders[j], traces[k], &wave[count],
                        &spec[count], &uncs[count], &cont[count]);
                    if (sucess == 0){
                        count++;
                    }

                }
                cpl_free(traces);
            }
            else
            {
                // Skip orders without the requested trace
                if (cr2res_extract_data(trace_wave[i], blaze[i], extracted_1d[i], 
                    orders[j], trace, &wave[count], &spec[count], &uncs[count],
                    &cont[count]) != -1)
                {
                    count++;
                }
            }
        }
        cpl_free(orders);
    }

    nspectra = count;
    cpl_msg_info(__func__, "%i segments with data found", nspectra);

    // Splice orders, but keep them seperate
    cr2res_splice_orders(wave, spec, uncs, cont, nspectra, &sp, &sp_err,
        &sp_order, &first, &last);

    // Combine orders into one big spectrum
    cr2res_combine_spectra(sp, sp_err, sp_order, first, last, nspectra, spliced,
        spliced_err);

    cpl_bivector_delete(first);
    cpl_bivector_delete(last);
    cpl_vector_delete(sp_order);
    
    for (i = 0; i < nspectra; i++){
        cpl_bivector_delete(sp[i]);
        cpl_bivector_delete(sp_err[i]);
    }

    cpl_free(sp);
    cpl_free(sp_err);

    for (i = 0; i < nspectra; i++){
        // spec, uncs, cont were just wrapping table data
        // while wave, is its own data (the evaluated polynomial)
        cpl_vector_delete(wave[i]);
        cpl_vector_unwrap(spec[i]);
        cpl_vector_unwrap(uncs[i]);
        cpl_vector_unwrap(cont[i]);
    }
    cpl_free(wave);
    cpl_free(spec);
    cpl_free(uncs);
    cpl_free(cont);

    if (cpl_error_get_code() == CPL_ERROR_DIVISION_BY_ZERO){
        cpl_msg_warning(__func__,
            "Division by Zero encountered. Resetting Error State");
        cpl_error_reset();
    }
    return 0 ;
}


int cr2res_extract_data(
    cpl_table   *   trace_wave,
    cpl_table   *   blaze, 
    cpl_table   *   spectra, 
    int             order, 
    int             trace,
    cpl_vector **   wave,
    cpl_vector **   spec,
    cpl_vector **   uncs,
    cpl_vector **   cont)
{
    int j, k, n = CR2RES_DETECTOR_SIZE;
    const cpl_array * wave1;
    cpl_polynomial * wave_poly;
    double * tmp1;
    char * colname;

    j  = cr2res_get_trace_table_index(trace_wave, order, trace);
    if (j == -1){
        // trace or order not found
        return -1;
    }

    // Get wavelength
    wave1 = cpl_table_get_array(trace_wave, CR2RES_COL_WAVELENGTH, j);
    if (wave1 == NULL){
        cpl_msg_warning(__func__, "Order %i, Trace %i", order, trace);
        cpl_msg_warning(__func__, "Wavelength Array not set");
        return -1;
    }
    wave_poly = cr2res_convert_array_to_poly(wave1);
    *wave = cpl_vector_new(n);
    for (k = 0; k < n; k++)
        cpl_vector_set(*wave, k, cpl_polynomial_eval_1d(wave_poly, k, NULL));
    cpl_polynomial_delete(wave_poly);

    // Get Spectrum
    colname = cr2res_dfs_SPEC_colname(order, trace);
    tmp1 = cpl_table_get_data_double(spectra, colname);
    cpl_free(colname);

    if (tmp1 == NULL){
        cpl_msg_warning(__func__, "Order %i, Trace %i", order, trace);
        cpl_msg_warning(__func__, "Spectrum not set");
        cpl_vector_delete(*wave);
        return -1;
    }
    *spec = cpl_vector_wrap(n, tmp1);


    // Get Uncertainties
    colname = cr2res_dfs_SPEC_ERR_colname(order, trace);
    tmp1 = cpl_table_get_data_double(spectra, colname);
    cpl_free(colname);
    if (tmp1 == NULL){
        cpl_msg_warning(__func__, "Order %i, Trace %i", order, trace);
        cpl_msg_warning(__func__, "Error not set");
        cpl_vector_delete(*wave);
        cpl_vector_unwrap(*spec);
        return -1;
    }
    *uncs = cpl_vector_wrap(n, tmp1);

    // Get Continuum
    colname = cr2res_dfs_SPEC_colname(order, trace);
    tmp1 = cpl_table_get_data_double(blaze, colname);
    cpl_free(colname);
    if (tmp1 == NULL){
        cpl_msg_warning(__func__, "Order %i, Trace %i", order, trace);
        cpl_msg_warning(__func__, "Blaze not set");
        cpl_vector_delete(*wave);
        cpl_vector_unwrap(*spec);
        cpl_vector_unwrap(*uncs);
        return -1;
    }
    *cont = cpl_vector_wrap(n, tmp1);

    return 0;
}


int cr2res_combine_spectra(
    cpl_bivector * spliced[],
    cpl_bivector * spliced_err[],
    cpl_vector * spectrum_order,
    cpl_bivector * first,
    cpl_bivector * last,
    int nspectra,
    cpl_bivector ** spectrum,
    cpl_bivector ** spectrum_err)
{

    if (spliced == NULL | spliced_err == NULL | nspectra <= 0 |
        spectrum == NULL | spectrum_err == NULL | spectrum_order == NULL)
            return -1;

    int i = 0, j = 0, k = 0;
    int iord = 0, iord_p1 = 0, iord_m1 = 0;
    volatile int fx, lx, fy, ly;

    *spectrum = cpl_bivector_new(nspectra * CR2RES_DETECTOR_SIZE);
    *spectrum_err = cpl_bivector_new(nspectra * CR2RES_DETECTOR_SIZE);

    for (i = 0; i < nspectra; i++){
        iord = cpl_vector_get(spectrum_order, i);

        fx = cpl_bivector_get_x_data(first)[iord];
        lx = cpl_bivector_get_x_data(last)[iord] + 1;
        fy = cpl_bivector_get_y_data(first)[iord];
        ly = cpl_bivector_get_y_data(last)[iord] + 1;

        if ( i == 0 ){
            fx = 0;
            lx = 0;
        }
        if (i == nspectra-1){
            fy = CR2RES_DETECTOR_SIZE;            
            ly = CR2RES_DETECTOR_SIZE;
        }
        cpl_msg_debug(__func__, "lx: %i, ly: %i", lx, ly);

        // in the overlap region keep the data from the lower wavelength order
        for (j = lx; j < ly; j++){
            cpl_vector_set(cpl_bivector_get_x(*spectrum), k, 
                cpl_vector_get(cpl_bivector_get_x(spliced[iord]), j));
            cpl_vector_set(cpl_bivector_get_y(*spectrum), k, 
                cpl_vector_get(cpl_bivector_get_y(spliced[iord]), j));

            cpl_vector_set(cpl_bivector_get_x(*spectrum_err), k, 
                cpl_vector_get(cpl_bivector_get_x(spliced_err[iord]), j));
            cpl_vector_set(cpl_bivector_get_y(*spectrum_err), k, 
                cpl_vector_get(cpl_bivector_get_y(spliced_err[iord]), j));

            k++;
        }

    }

    cpl_vector_set_size(cpl_bivector_get_x(*spectrum), k-1);
    cpl_vector_set_size(cpl_bivector_get_y(*spectrum), k-1);

    cpl_vector_set_size(cpl_bivector_get_x(*spectrum_err), k-1);
    cpl_vector_set_size(cpl_bivector_get_y(*spectrum_err), k-1);


    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Splice any number of orders and traces into adjacent spectra
  @param    trace_wave     trace wave table
  @param    blaze          blaze table
  @param    spectra        spectrum and error table
  @param    trace          trace to merge, or -1 for all
  @param    spliced        [out] list of spliced spectra
  @param    spliced_err    [out] list of spliced error spectra
  @return   0 on success, -1 on failure

  Spectra are spliced by first interpolating the overlapping regions on the 
  wavelength grid of the other order, and then adding the spectra 
  weighted by the errors

 */
/*----------------------------------------------------------------------------*/
int cr2res_splice_orders(
        cpl_vector   **  wave,
        cpl_vector   **  spec,
        cpl_vector   **  uncs,
        cpl_vector   **  cont,
        int              nspectra,
        cpl_bivector **  spliced[],
        cpl_bivector **  spliced_err[],
        cpl_vector   **  spectrum_order,
        cpl_bivector **  first,
        cpl_bivector **  last)
{
    if (wave == NULL | spec == NULL | cont == NULL | uncs == NULL
        | spliced == NULL | spliced_err == NULL | nspectra <= 0
        | spectrum_order == NULL | first == NULL | last == NULL) return -1;

    cpl_size i, j, k; 
    int central = 0;
    int current_order = 0, current_trace = 0;
    int iorder = 0, itrace = 0;
    int * loop0, *loop1;
    int iord0, iord1;
    int first0, first1, last0, last1, overlap0, overlap1;
    int n = CR2RES_DETECTOR_SIZE;
    double maximum = 0;
    double minW0, minW1, maxW0, maxW1;
    cpl_size minW0pos, minW1pos, maxW0pos, maxW1pos;
    double wgt0, wgt1;
    // temporary work arrays
    char * colname;
    const cpl_array * wave1;
    cpl_polynomial * wave_poly;
    cpl_bivector * wave_center;
    cpl_vector *s0, *s1, *tmpS0, *tmpS1;
    cpl_vector *u0, *u1, *tmpU0, *tmpU1;
    cpl_vector *w0, *w1, *tmpW0, *tmpW1;
    cpl_vector *c0, *c1, *tmpC0, *tmpC1;

    double * tmp1;
    cpl_bivector * tmp2, * tmp3;
    cpl_vector * tmp4;
    double median;

    loop0 = cpl_malloc((nspectra-1) * sizeof(int));
    loop1 = cpl_malloc((nspectra-1) * sizeof(int));

    *spliced = cpl_malloc(nspectra * sizeof(cpl_bivector*));
    *spliced_err = cpl_malloc(nspectra * sizeof(cpl_bivector*));
    *first = cpl_bivector_new(nspectra);
    *last = cpl_bivector_new(nspectra);


    for (i = 0; i < nspectra; i++){
        (*spliced)[i] = cpl_bivector_new(n);
        (*spliced_err)[i] = cpl_bivector_new(n);
    }
    wave_center = cpl_bivector_new(nspectra);

    // Load data into vector arrays
    for (i=0; i<nspectra; i++){

        // Get Wavelength Center
        cpl_vector_set(cpl_bivector_get_x(wave_center), i, i);        
        cpl_vector_set(cpl_bivector_get_y(wave_center), i, 
            cpl_vector_get(wave[i], n/2));

        // scale all orders to spec/cont = 1
        // also replace nans with 0 in spec
        // and nan and 0 with 1 in cont
        tmp4 = cpl_vector_duplicate(spec[i]);
        for (j=0; j < n; j++)
        {
            if (isfinite(cpl_vector_get(spec[i], j)) == 0){
                cpl_vector_set(tmp4, j, 0);
                cpl_vector_set(spec[i], j, 0);
            }
            if (isfinite(cpl_vector_get(uncs[i], j)) == 0){
                cpl_vector_set(uncs[i], j, 1);
            }
            if ((cpl_vector_get(cont[i], j) == 0) | (isfinite(cpl_vector_get(cont[i], j)) == 0)){
                cpl_vector_set(cont[i], j, 1);
            }
            cpl_vector_set(tmp4, j, 
                    cpl_vector_get(tmp4, j) / cpl_vector_get(cont[i], j));
        }
        // cpl_vector_divide(tmp4, cont[i]);
        median = cpl_vector_get_median(tmp4);
        cpl_vector_multiply_scalar(cont[i], median);
        cpl_vector_delete(tmp4);
    }

    // Determine order of orders (wavelength sections)
    cpl_bivector_sort(wave_center, wave_center, CPL_SORT_ASCENDING, 
        CPL_SORT_BY_Y);

    // just loop from left to right?
    for (i = 0; i < nspectra-1; i++){
        loop0[i] = cpl_bivector_get_x_data(wave_center)[i];
        loop1[i] = cpl_bivector_get_x_data(wave_center)[i+1];
    }

    for (i=0; i < nspectra-1; i++){
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
        // If the extraction was curved, the border points should be 0
        // with the number of points depending on the curvature
        minW1pos = cpl_vector_get_minpos(w1);
        minW0pos = cpl_vector_get_minpos(w0);
        maxW1pos = cpl_vector_get_maxpos(w1);
        maxW0pos = cpl_vector_get_maxpos(w0);

        while (cpl_vector_get(s1, minW1pos) == 0)   minW1pos++;
        while (cpl_vector_get(s0, minW0pos) == 0)   minW0pos++;
        while (cpl_vector_get(s1, maxW1pos) == 0)   maxW1pos--;
        while (cpl_vector_get(s0, maxW0pos) == 0)   maxW0pos--;
        
        minW1 = cpl_vector_get(w1, minW1pos);
        minW0 = cpl_vector_get(w0, minW0pos);
        maxW1 = cpl_vector_get(w1, maxW1pos);
        maxW0 = cpl_vector_get(w0, maxW0pos);

        // indices of the first and last overlaping value
        first0 = -1;
        first1 = -1;
        last0 = -1;
        last1 = -1;

        // overlap0: number of overlap points of order "0"
        // overlap1: number of overlap points of order "1"

        overlap0 = 0;
        for (j = 0; j < cpl_vector_get_size(w0); j++){
            if ((cpl_vector_get(w0, j) >= minW1) & (cpl_vector_get(w0, j) <= maxW1)){
                // if its the first element
                if (first0 == -1) first0 = j;
                overlap0++;
                // update so that it will be the last element
                last0 = j;
            }
        }

        overlap1 = 0;
        for (j = 0; j < cpl_vector_get_size(w1); j++){
            if ((cpl_vector_get(w1, j) >= minW0) & (cpl_vector_get(w1, j) <= maxW0)){
                // if its the first element
                if (first1 == -1) first1 = j;
                overlap1++;
                // update so that it will be the last element
                last1 = j;
            }
        }
        cpl_msg_debug(__func__, "Overlap: %i, %i", overlap0, overlap1);
        if (overlap0 == 0 || overlap1 == 0){
            cpl_vector_set(cpl_bivector_get_y(*first), iord0, CR2RES_DETECTOR_SIZE);
            cpl_vector_set(cpl_bivector_get_y(*last), iord0, maxW0pos-1);
            cpl_vector_set(cpl_bivector_get_x(*first), iord1, 0);
            cpl_vector_set(cpl_bivector_get_x(*last), iord1, minW1pos-1);

            for (k = 0; k < CR2RES_DETECTOR_SIZE; k++){
                cpl_vector_set(cpl_bivector_get_x((*spliced)[iord0]), k, 
                    cpl_vector_get(w0, k));
                cpl_vector_set(cpl_bivector_get_x((*spliced)[iord1]), k, 
                    cpl_vector_get(w1, k));

                cpl_vector_set(cpl_bivector_get_x((*spliced_err)[iord0]), k, 
                    cpl_vector_get(w0, k));
                cpl_vector_set(cpl_bivector_get_x((*spliced_err)[iord1]), k, 
                    cpl_vector_get(w1, k));


                cpl_vector_set(cpl_bivector_get_y((*spliced)[iord0]), k, 
                    cpl_vector_get(s0, k));
                cpl_vector_set(cpl_bivector_get_y((*spliced)[iord1]), k, 
                    cpl_vector_get(s1, k));

                cpl_vector_set(cpl_bivector_get_y((*spliced_err)[iord0]), k, 
                    cpl_vector_get(u0, k));
                cpl_vector_set(cpl_bivector_get_y((*spliced_err)[iord1]), k, 
                    cpl_vector_get(u1, k));

            }
            
            continue;
        }

        // It should be:
        // first.x : first index of the overlap on the left (lower wavelength) side of the order
        // first.y : first index of the overlap on the right (higher wavelength) side of the order
        // last.x : last index of the overlap on the left (lower wavelength) side of the order
        // last.y : last index of the overlap on the right (higher wavelength) side of the order

        cpl_vector_set(cpl_bivector_get_y(*first), iord0, first0);
        cpl_vector_set(cpl_bivector_get_y(*last), iord0, last0);
        cpl_vector_set(cpl_bivector_get_x(*first), iord1, first1);
        cpl_vector_set(cpl_bivector_get_x(*last), iord1, last1);

        // extract overlapping wavelength vectors
        // point to the latter parts of the wavelength vector, no need for copies
        tmpW0 = cpl_vector_wrap(last0 - first0 + 1, cpl_vector_get_data(w0) + first0);
        tmpW1 = cpl_vector_wrap(last1 - first1 + 1, cpl_vector_get_data(w1) + first1);

        // Prepare new vectors with interpolated data
        tmpS0 = cpl_vector_new(overlap0);
        tmpU0 = cpl_vector_new(overlap0);
        tmpC0 = cpl_vector_new(overlap0);

        tmpS1 = cpl_vector_new(overlap1);
        tmpU1 = cpl_vector_new(overlap1);
        tmpC1 = cpl_vector_new(overlap1);

        // Interpolate vectors on the other wavelength grid
        // Note that this is just linear interpolation
        tmp3 = cpl_bivector_wrap_vectors(w1, s1);
        tmp2 = cpl_bivector_wrap_vectors(tmpW0, tmpS0);
        cpl_bivector_interpolate_linear(tmp2, tmp3);
        cpl_bivector_unwrap_vectors(tmp3);
        cpl_bivector_unwrap_vectors(tmp2);

        tmp3 = cpl_bivector_wrap_vectors(w0, s0);
        tmp2 = cpl_bivector_wrap_vectors(tmpW1, tmpS1);
        cpl_bivector_interpolate_linear(tmp2, tmp3);
        cpl_bivector_unwrap_vectors(tmp3);
        cpl_bivector_unwrap_vectors(tmp2);

        tmp3 = cpl_bivector_wrap_vectors(w1, u1);
        tmp2 = cpl_bivector_wrap_vectors(tmpW0, tmpU0);
        cpl_bivector_interpolate_linear(tmp2, tmp3);
        cpl_bivector_unwrap_vectors(tmp3);
        cpl_bivector_unwrap_vectors(tmp2);

        tmp3 = cpl_bivector_wrap_vectors(w0, u0);
        tmp2 = cpl_bivector_wrap_vectors(tmpW1, tmpU1);
        cpl_bivector_interpolate_linear(tmp2, tmp3);
        cpl_bivector_unwrap_vectors(tmp3);
        cpl_bivector_unwrap_vectors(tmp2);

        tmp3 = cpl_bivector_wrap_vectors(w1, c1);
        tmp2 = cpl_bivector_wrap_vectors(tmpW0, tmpC0);
        cpl_bivector_interpolate_linear(tmp2, tmp3);
        cpl_bivector_unwrap_vectors(tmp3);
        cpl_bivector_unwrap_vectors(tmp2);

        tmp3 = cpl_bivector_wrap_vectors(w0, c0);
        tmp2 = cpl_bivector_wrap_vectors(tmpW1, tmpC1);
        cpl_bivector_interpolate_linear(tmp2, tmp3);
        cpl_bivector_unwrap_vectors(tmp3);
        cpl_bivector_unwrap_vectors(tmp2);

        // Combine Spectra with weighted average at each point
        // its important to remember the continuum in the uncertainties/weights
        k = 0;
        for (j=first0; j<=last0; j++){
            
            wgt0 = cpl_vector_get(c0, j) / cpl_vector_get(u0, j);
            wgt0 = wgt0 * wgt0;
            wgt1 = cpl_vector_get(tmpC0, k) / cpl_vector_get(tmpU0, k);
            wgt1 = wgt1 * wgt1;
            
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
            wgt0 = cpl_vector_get(c1, j) / cpl_vector_get(u1, j);
            wgt0 = wgt0 * wgt0;
            wgt1 = cpl_vector_get(tmpC1, k) / cpl_vector_get(tmpU1, k);
            wgt1 = wgt1 * wgt1;
            
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

        // Set data in output vector
        // TODO, avoid duplicate copies

        for (k = 0; k < n; k++){
            cpl_vector_set(cpl_bivector_get_x((*spliced)[iord0]), k, 
                cpl_vector_get(w0, k));
            cpl_vector_set(cpl_bivector_get_x((*spliced)[iord1]), k, 
                cpl_vector_get(w1, k));

            cpl_vector_set(cpl_bivector_get_x((*spliced_err)[iord0]), k, 
                cpl_vector_get(w0, k));
            cpl_vector_set(cpl_bivector_get_x((*spliced_err)[iord1]), k, 
                cpl_vector_get(w1, k));


            cpl_vector_set(cpl_bivector_get_y((*spliced)[iord0]), k, 
                cpl_vector_get(s0, k));
            cpl_vector_set(cpl_bivector_get_y((*spliced)[iord1]), k, 
                cpl_vector_get(s1, k));

            cpl_vector_set(cpl_bivector_get_y((*spliced_err)[iord0]), k, 
                cpl_vector_get(u0, k));
            cpl_vector_set(cpl_bivector_get_y((*spliced_err)[iord1]), k, 
                cpl_vector_get(u1, k));

        }
    }
    
    for (i = 0; i < nspectra; i++)
    {
        cpl_vector_divide(cpl_bivector_get_y((*spliced)[i]), cont[i]);
        cpl_vector_divide(cpl_bivector_get_y((*spliced_err)[i]), cont[i]);
    }
    
    *spectrum_order = cpl_bivector_get_x(wave_center);
    cpl_vector_delete(cpl_bivector_get_y(wave_center));
    cpl_free(wave_center);

    cpl_free(loop0);
    cpl_free(loop1);
    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Create the SPLICED_1D table to be saved
  @param    spectrum        The extracted spectra of the different orders
  @param    spectrum_error  The spectrum error
  @return   the SPLICED_1D table or NULL
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_splice_SPLICED_1D_create(
        cpl_bivector    *   spectrum,
        cpl_bivector    *   spectrum_error)
{
    cpl_table       *   out ;
    const double    *   pspec ;
    const double    *   perr ;
    cpl_vector      *   wave_vec ;
    const double    *   pwl ;
    cpl_size            nbins, i;

    /* Check entries */
    if (spectrum == NULL || spectrum_error == NULL) return NULL ;

    /* Initialise */
    nbins = cpl_bivector_get_size(spectrum) ;
    if (cpl_bivector_get_size(spectrum_error) != nbins) return NULL ;

    /* Create the table */
    out = cpl_table_new(nbins);

    /* Create SPLICED_1D_SPEC columns */
    cpl_table_new_column(out, CR2RES_COL_SPLICED_1D_WL, CPL_TYPE_DOUBLE);
    cpl_table_new_column(out, CR2RES_COL_SPLICED_1D_SPEC, CPL_TYPE_DOUBLE);
    cpl_table_new_column(out, CR2RES_COL_SPLICED_1D_ERROR, CPL_TYPE_DOUBLE);

    /* Fill the table */
    pwl = cpl_bivector_get_x_data_const(spectrum) ;
    pspec = cpl_bivector_get_y_data_const(spectrum) ;
    perr = cpl_bivector_get_y_data_const(spectrum_error);

    cpl_table_copy_data_double(out, CR2RES_COL_SPLICED_1D_WL, pwl) ;
    cpl_table_copy_data_double(out, CR2RES_COL_SPLICED_1D_SPEC, pspec) ;
    cpl_table_copy_data_double(out, CR2RES_COL_SPLICED_1D_ERROR, perr) ;

    return out ;
}

/**@}*/
