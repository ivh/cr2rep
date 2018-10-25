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

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cpl.h>
#include <hdrl.h>
#include <cr2res_dfs.h>
#include <cr2res_splice.h>

#define CR2RES_DETECTOR_SIZE            2048

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static void cr2res_splice_orders_test(void);
static void cr2res_splice_test(void);


/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_splice-test    Unit test of cr2res_splice
 *
 */
/*----------------------------------------------------------------------------*/
/**@{*/

cpl_table * make_trace_wave(){
    int wave_poly_size = 3;
    int i;
    cpl_array * wave = cpl_array_new(wave_poly_size, CPL_TYPE_DOUBLE);
    // columns: order, trace, wavelength
    cpl_table * trace_wave = cpl_table_new(2);
    cpl_table_new_column(trace_wave, CR2RES_COL_ORDER, CPL_TYPE_INT);
    cpl_table_new_column(trace_wave, CR2RES_COL_TRACENB, CPL_TYPE_INT);
    cpl_table_new_column_array(trace_wave, CR2RES_COL_WAVELENGTH, CPL_TYPE_DOUBLE, wave_poly_size);

    cpl_array_set_double(wave, 0, 2000);
    cpl_array_set_double(wave, 1, 1);
    cpl_array_set_double(wave, 2, 0.);

    // first order
    i = 0;
    cpl_table_set_int(trace_wave, CR2RES_COL_ORDER, i, 1);
    cpl_table_set_int(trace_wave, CR2RES_COL_TRACENB, i, 1);
    cpl_table_set_array(trace_wave, CR2RES_COL_WAVELENGTH, i, wave);


    // second order
    i = 1;
    cpl_table_set_int(trace_wave, CR2RES_COL_ORDER, i, 2);
    cpl_table_set_int(trace_wave, CR2RES_COL_TRACENB, i, 1);

    cpl_array_set(wave, 0, 500.5);
    cpl_table_set_array(trace_wave, CR2RES_COL_WAVELENGTH, i, wave);

    cpl_array_delete(wave);

    return trace_wave;
}

cpl_table * make_spectra(cpl_table ** blaze){
    cpl_size i;
    char * cn1, * cn2, * cn3, * cn4; // column names
    cpl_table * spectra = cpl_table_new(CR2RES_DETECTOR_SIZE);
    *blaze = cpl_table_new(CR2RES_DETECTOR_SIZE);

    cn1 = cr2res_dfs_SPEC_colname(1, 1);
    cn2 = cr2res_dfs_SPEC_ERR_colname(1, 1);
    cn3 = cr2res_dfs_SPEC_colname(2, 1);
    cn4 = cr2res_dfs_SPEC_ERR_colname(2, 1);

    cpl_table_new_column(spectra, cn1, CPL_TYPE_DOUBLE);
    cpl_table_new_column(*blaze, cn1, CPL_TYPE_DOUBLE);

    cpl_table_new_column(spectra, cn2, CPL_TYPE_DOUBLE);

    cpl_table_new_column(spectra, cn3, CPL_TYPE_DOUBLE);
    cpl_table_new_column(*blaze, cn3, CPL_TYPE_DOUBLE);

    cpl_table_new_column(spectra, cn4, CPL_TYPE_DOUBLE);
    
    for (i = 0; i < CR2RES_DETECTOR_SIZE; i++){
        cpl_table_set_double(spectra, cn1, i, i + 2000);
        cpl_table_set_double(spectra, cn2, i, 0.1);

        cpl_table_set_double(spectra, cn3, i, i + 500.5);
        cpl_table_set_double(spectra, cn4, i, 0.1);

        cpl_table_set_double(*blaze, cn1, i, 1.);
        cpl_table_set_double(*blaze, cn3, i, 1.);
    }

    cpl_free(cn1);
    cpl_free(cn2);
    cpl_free(cn3);
    cpl_free(cn4);

    return spectra;
}

static void cr2res_splice_orders_test(){

    int i = 0, j = 0;
    int ntotal = 2;

    cpl_vector ** wave = cpl_malloc(ntotal * sizeof(cpl_vector*));
    cpl_vector ** spec = cpl_malloc(ntotal * sizeof(cpl_vector*));
    cpl_vector ** uncs = cpl_malloc(ntotal * sizeof(cpl_vector*));
    cpl_vector ** cont = cpl_malloc(ntotal * sizeof(cpl_vector*));

    for (j = 0; j < ntotal; j++){
        wave[j] = cpl_vector_new(CR2RES_DETECTOR_SIZE);
        spec[j] = cpl_vector_new(CR2RES_DETECTOR_SIZE);
        uncs[j] = cpl_vector_new(CR2RES_DETECTOR_SIZE);
        cont[j] = cpl_vector_new(CR2RES_DETECTOR_SIZE);
    }

    for (i = 0; i < CR2RES_DETECTOR_SIZE; i++){
        cpl_vector_set(wave[0], i, 2000 + i);
        cpl_vector_set(spec[0], i, 2000 + i);
        cpl_vector_set(uncs[0], i, 0.1);
        cpl_vector_set(cont[0], i, 1);


        cpl_vector_set(wave[1], i, 500.5 + i);
        cpl_vector_set(spec[1], i, 500.5 + i);
        cpl_vector_set(uncs[1], i, 0.1);
        cpl_vector_set(cont[1], i, 1);
    }

    cpl_bivector ** spliced, ** spliced_err;
    cpl_vector * spectrum_order;
    cpl_bivector * first, *last;

    int res;

    // Test NULL input
    cpl_test_eq(-1, cr2res_splice_orders(NULL, spec, uncs, cont, ntotal, &spliced, &spliced_err, &spectrum_order, &first, &last));
    cpl_test_eq(-1, cr2res_splice_orders(wave, NULL, uncs, cont, ntotal, &spliced, &spliced_err, &spectrum_order, &first, &last));
    cpl_test_eq(-1, cr2res_splice_orders(wave, spec, NULL, cont, ntotal, &spliced, &spliced_err, &spectrum_order, &first, &last));
    cpl_test_eq(-1, cr2res_splice_orders(wave, spec, uncs, NULL, ntotal, &spliced, &spliced_err, &spectrum_order, &first, &last));
    cpl_test_eq(-1, cr2res_splice_orders(wave, spec, uncs, cont, 0, &spliced, &spliced_err, &spectrum_order, &first, &last));
    cpl_test_eq(-1, cr2res_splice_orders(wave, spec, uncs, cont, ntotal, NULL, &spliced_err, &spectrum_order, &first, &last));
    cpl_test_eq(-1, cr2res_splice_orders(wave, spec, uncs, cont, ntotal, &spliced, NULL, &spectrum_order, &first, &last));
    cpl_test_eq(-1, cr2res_splice_orders(wave, spec, uncs, cont, ntotal, &spliced, &spliced_err, NULL, &first, &last));
    cpl_test_eq(-1, cr2res_splice_orders(wave, spec, uncs, cont, ntotal, &spliced, &spliced_err, &spectrum_order, NULL, &last));
    cpl_test_eq(-1, cr2res_splice_orders(wave, spec, uncs, cont, ntotal, &spliced, &spliced_err, &spectrum_order, &first, NULL));


    // Test proper values
    res = cr2res_splice_orders(wave, spec, uncs, cont, ntotal, &spliced, &spliced_err, &spectrum_order, &first, &last);
    cpl_test_eq(0, res);
    
    // Check results, this depends on the initial values of course
    // Here the spectrum was equal to the wavelength, this should still be the case after splicing
    for (i = 0; i < ntotal; i++){
        cpl_test_vector_abs(cpl_bivector_get_x(spliced[i]), cpl_bivector_get_y(spliced[i]), DBL_EPSILON);
    }

    // Delete all memory
    cpl_vector_delete(spectrum_order);
    cpl_bivector_delete(first);
    cpl_bivector_delete(last);

    for (int i = 0; i < ntotal; i++){
        cpl_bivector_delete(spliced[i]);
        cpl_bivector_delete(spliced_err[i]);
        cpl_vector_delete(wave[i]);
        cpl_vector_delete(spec[i]);
        cpl_vector_delete(uncs[i]);
        cpl_vector_delete(cont[i]);
    }

    cpl_free(wave);
    cpl_free(spec);
    cpl_free(uncs);
    cpl_free(cont);

    cpl_free(spliced);
    cpl_free(spliced_err);
    
    return;
}

static void cr2res_splice_test(){
    int i = 0, j = 0;
    cpl_table * trace_wave = make_trace_wave();
    cpl_table * blaze;
    cpl_table * spectra = make_spectra(&blaze);
    cpl_bivector * spliced, * spliced_err;

    int res;
    res = cr2res_splice(&spectra, &blaze, &trace_wave, 1, &spliced, &spliced_err);

    cpl_test_eq(0, res);
    for (i = 0; i < CR2RES_DETECTOR_SIZE - 1; i++){
        // Wavelength grid is increasing
        cpl_test(cpl_bivector_get_x_data(spliced)[i] < cpl_bivector_get_x_data(spliced)[i+1]);

        // Here: spectrum == wavelength grid, i.e. also increasing
        cpl_test(cpl_bivector_get_y_data(spliced)[i] < cpl_bivector_get_y_data(spliced)[i+1]);
        cpl_test_abs(cpl_bivector_get_x_data(spliced)[i], cpl_bivector_get_y_data(spliced)[i], DBL_EPSILON);
    
        // Here: Error <= 0.1, initial error
        cpl_test(cpl_bivector_get_y_data(spliced_err)[i] <= 0.1);
    }

    cpl_bivector_dump(spliced_err, NULL);

    cpl_table_delete(trace_wave);
    cpl_table_delete(blaze);
    cpl_table_delete(spectra);
    cpl_bivector_delete(spliced);
    cpl_bivector_delete(spliced_err);


}

/*----------------------------------------------------------------------------*/
/**
  @brief    Run the Unit tests
 */
/*----------------------------------------------------------------------------*/
int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_DEBUG);

    // cr2res_splice_orders_test();
    cr2res_splice_test();

    return cpl_test_end(0);
}

/**@}*/
