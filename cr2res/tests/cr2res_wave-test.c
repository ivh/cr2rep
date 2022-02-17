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
#include <cr2res_wave.h>
#include <cr2res_wave.c>
#include <cr2res_etalon.h>

#define CR2RES_DETECTOR_SIZE            2048

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static void test_cr2res_wave_1d(void);
static void test_cr2res_wave_2d(void);
static void test_cr2res_wave_line_fitting_2d_other(void);
static void test_cr2res_wave_etalon(void);
static void test_cr2res_wave_etalon_other(void);
static void test_cr2res_wave_polys_1d_to_2d(void);
static void test_cr2res_wave_poly_2d_to_1d(void);
static void test_cr2res_wave_estimate_compute(void);
static void test_cr2res_wave_clean_spectrum(void);
static void test_cr2res_wave_etalon_2d(void);

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_wave-test    Unit test of cr2res_wave
 *
 */
/*----------------------------------------------------------------------------*/
/**@{*/

// make a test line list catalog with just 2 lines
static cpl_table * make_test_catalog()
{
    cpl_table * catalog = cpl_table_new(2);
    cpl_table_new_column(catalog, CR2RES_COL_WAVELENGTH, CPL_TYPE_DOUBLE);
    cpl_table_new_column(catalog, CR2RES_COL_EMISSION, CPL_TYPE_DOUBLE);
    //cpl_table_new_column(catalog, CR2RES_COL_WIDTH, CPL_TYPE_DOUBLE);

    // from lines_thar.txt catalog
    // 2.551218908614062912e+03 9.935860000000000127e+02
    // 2.603302966671994909e+03 2.992800000000000082e+01

    cpl_table_set_double(catalog, CR2RES_COL_WAVELENGTH, 0, 2.551218908614062912e+03);
    cpl_table_set_double(catalog, CR2RES_COL_EMISSION, 0, 9.935860000000000127e+02);

    cpl_table_set_double(catalog, CR2RES_COL_WAVELENGTH, 1, 2.603302966671994909e+03);
    cpl_table_set_double(catalog, CR2RES_COL_EMISSION, 1, 2.992800000000000082e+01);

    return catalog;
}

static const char * save_catalog(cpl_table * catalog){
    const char * filename = "TEST_linelist.fits";
    cpl_propertylist * header = cpl_propertylist_new();
    cpl_propertylist_append_string(header, CPL_DFS_PRO_TYPE, CR2RES_PROTYPE_CATALOG);

    cpl_table_save(catalog, header, NULL, filename, CPL_IO_CREATE);
    cpl_propertylist_delete(header);
    return filename;
}

static const char * save_linelist(cpl_bivector * linelist){
    const char * filename = "TEST_linelist.fits";
    cpl_propertylist * header = cpl_propertylist_new();
    cpl_propertylist_append_string(header, CPL_DFS_PRO_TYPE, CR2RES_PROTYPE_CATALOG);

    cpl_table * table = cpl_table_new(cpl_bivector_get_size(linelist));
    cpl_table_new_column(table, CR2RES_COL_WAVELENGTH, CPL_TYPE_DOUBLE);
    cpl_table_new_column(table, CR2RES_COL_EMISSION, CPL_TYPE_DOUBLE);

    for (cpl_size i = 0; i < cpl_bivector_get_size(linelist); i++)
    {
        cpl_table_set_double(table, CR2RES_COL_WAVELENGTH, i, cpl_vector_get(cpl_bivector_get_x(linelist), i));
        cpl_table_set_double(table, CR2RES_COL_EMISSION, i, cpl_vector_get(cpl_bivector_get_y(linelist), i));
    }

    cpl_table_save(table, header, NULL, filename, CPL_IO_CREATE);

    cpl_propertylist_delete(header);
    cpl_table_delete(table);
    return filename;
}

// make a test spectrum based on a line list catalog
static cpl_bivector * make_test_spectrum(cpl_table * catalog, double wmin, double wmax, int size, cpl_bivector ** spectrum_err)
{
    double wl, mu, line_em, sig;
    double tmp;
    int i, j;
    cpl_bivector * spectrum = cpl_bivector_new(size);
    *spectrum_err = cpl_bivector_new(size);
    cpl_vector * wave1 = cpl_bivector_get_x(spectrum);
    cpl_vector * wave2 = cpl_bivector_get_x(*spectrum_err);

    cpl_vector * spec = cpl_bivector_get_y(spectrum);
    cpl_vector * unc = cpl_bivector_get_y(*spectrum_err);

    for (i = 0; i < size; i++){
        wl = wmin + i * (wmax - wmin) / (double)size;
        tmp = 0;
        for (j = 0; j < 2; j++){
            mu = cpl_table_get_double(catalog, CR2RES_COL_WAVELENGTH, j, NULL);
            line_em = cpl_table_get_double(catalog, CR2RES_COL_EMISSION, j, NULL);
            sig = 2;
            tmp += line_em * exp(- 1 * pow(wl - mu, 2) / (2 * pow(sig, 2)));
        }

        cpl_vector_set(spec, i, tmp);
        cpl_vector_set(unc, i, 0.01);
        cpl_vector_set(wave1, i, wl);
        cpl_vector_set(wave2, i, wl);
    }

    return spectrum;
}

static cpl_bivector * make_test_etalon_spectrum(int size, double freq, cpl_bivector ** spectrum_err)
{
    cpl_bivector * spectrum = cpl_bivector_new(size);
    *spectrum_err = cpl_bivector_new(size);

    cpl_vector * wave1 = cpl_bivector_get_x(spectrum);
    cpl_vector * wave2 = cpl_bivector_get_x(*spectrum_err);

    cpl_vector * spec = cpl_bivector_get_y(spectrum);
    cpl_vector * unc = cpl_bivector_get_y(*spectrum_err);


    for (int i = 0; i < size; i++){
        cpl_vector_set(spec, i, fabs(sin(i * freq)));
        cpl_vector_set(unc, i, 1);
        cpl_vector_set(wave1, i, i);
        cpl_vector_set(wave2, i, i);

    }
    return spectrum;

}

// make a simple linear polynomial from wmin to wmax
static cpl_polynomial * make_test_polynomial(double wmin, double wmax, int size)
{
    cpl_polynomial * poly = cpl_polynomial_new(1);
    cpl_size power = 0;

    power = 0;
    cpl_polynomial_set_coeff(poly, &power, wmin);
    power = 1;
    cpl_polynomial_set_coeff(poly, &power,(wmax - wmin)/size);

    return poly;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Make a sample wavecal spectrum with two lines, and check that we get the right linear fit
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_wave_1d()
{
    double wmin=2500, wmax=2650;
    int size = 200;
    cpl_table * catalog = make_test_catalog();
    cpl_table * diagnostics;
    cpl_bivector * spectrum_err;
    cpl_bivector * spectrum = make_test_spectrum(catalog, wmin, wmax, size, &spectrum_err);
    cpl_polynomial * initial_guess = make_test_polynomial(wmin, wmax, size);
    int window_size = 30;
    int degree = 1;
    int order = 0;
    int trace = 0;
    int log_flag = 0; // False
    int display = 0; // False
    int propagate_flag = 0; // False
    const char * catalog_name = save_catalog(catalog);
    cr2res_wavecal_type wavecal_type = CR2RES_LINE1D;
    cpl_array * wave_error_init = cpl_array_new(2, CPL_TYPE_DOUBLE);
    cpl_array_set_double(wave_error_init, 0, 3.1);
    cpl_array_set_double(wave_error_init, 1, 3.5);

    cpl_array * wavelength_error;
    cpl_polynomial * wavelength;
    cpl_size power;

    // bad inputs
    wavelength = cr2res_wave_1d(NULL, spectrum_err, initial_guess,
        wave_error_init, order, trace, wavecal_type, catalog_name,
        degree, 0, log_flag, propagate_flag, display, -1.0, -1.0, NULL,
        NULL, NULL, NULL, NULL, &wavelength_error, &diagnostics);
    cpl_test_null(wavelength);

    wavelength = cr2res_wave_1d(spectrum, NULL, initial_guess,
        wave_error_init, order, trace, wavecal_type, catalog_name, degree, 0,
        log_flag, propagate_flag, display, -1.0, -1.0, NULL, NULL, NULL, NULL, 
        NULL, &wavelength_error, &diagnostics);
    cpl_test_null(wavelength);

    wavelength = cr2res_wave_1d(spectrum, spectrum_err, NULL,
        wave_error_init, order, trace, wavecal_type, catalog_name, degree,
        0, log_flag, propagate_flag, display, -1.0, -1.0, NULL, NULL,
        NULL, NULL, NULL, &wavelength_error, &diagnostics);
    cpl_test_null(wavelength);

    wavelength = cr2res_wave_1d(spectrum, spectrum_err, initial_guess,
        wave_error_init, order, trace, wavecal_type, NULL, degree, 0,
        log_flag, propagate_flag, display, -1.0, -1.0, NULL, NULL, NULL, NULL, 
        NULL, &wavelength_error, &diagnostics);
    cpl_test_null(wavelength);

    wavelength = cr2res_wave_1d(spectrum, spectrum_err, initial_guess,
        wave_error_init, order, trace, wavecal_type, catalog_name,
        degree, 0, log_flag, propagate_flag, display, -1.0, -1.0, NULL, NULL, 
        NULL, NULL, NULL, NULL, &diagnostics);
    cpl_test_null(wavelength);

    wavelength = cr2res_wave_1d(spectrum, spectrum_err, initial_guess,
        wave_error_init, order, trace, wavecal_type, catalog_name,
        degree, 0, log_flag, propagate_flag, display, -1.0, -1.0,NULL, NULL, 
        NULL, NULL, NULL, &wavelength_error, NULL);
    cpl_test_null(wavelength);

    // // to many polynomial degrees
    wavelength = cr2res_wave_1d(spectrum, spectrum_err, initial_guess,
        wave_error_init, order, trace, wavecal_type, catalog_name, 5, 0,
        log_flag, propagate_flag, display, -1.0, -1.0, NULL, NULL, NULL, NULL, 
        NULL, &wavelength_error, &diagnostics);

    cpl_test_null(wavelength);
    cpl_test_null(wavelength_error);
    cpl_test_nonnull(diagnostics);
    cpl_table_delete(diagnostics);

    // regular run
    cpl_test(wavelength = cr2res_wave_1d(spectrum, spectrum_err, initial_guess,
                wave_error_init, order, trace, wavecal_type, catalog_name,
                degree, 0, log_flag, propagate_flag, display, -1.0, -1.0, NULL,
                NULL, NULL, NULL, NULL, &wavelength_error, &diagnostics));

    cpl_test_nonnull(wavelength);
    cpl_test_nonnull(wavelength_error);
    cpl_test_nonnull(diagnostics);

    // these values obviously need to be changed if the number of degrees is changed
    power = 0;
    cpl_test_abs(cpl_polynomial_get_coeff(wavelength, &power), wmin, 0.2);
    power = 1;
    cpl_test_abs(cpl_polynomial_get_coeff(wavelength, &power), (wmax-wmin)/(double)size, 0.01);

    // Fitting two points with a first order polynomial -> perfect fit
    cpl_test_abs(cpl_array_get_double(wavelength_error, 0, NULL), 0, DBL_EPSILON);
    cpl_test_abs(cpl_array_get_double(wavelength_error, 1, NULL), 0, DBL_EPSILON);


    cpl_table_delete(diagnostics);
    cpl_array_delete(wavelength_error);
    cpl_polynomial_delete(wavelength);
    cpl_table_delete(catalog);
    cpl_bivector_delete(spectrum);
    cpl_bivector_delete(spectrum_err);
    cpl_polynomial_delete(initial_guess);
    cpl_array_delete(wave_error_init);
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Use two identical orders, with two lines each, and check that the result is still linear
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_wave_2d()
{
    int i, norders = 5;
    double wmin=2500, wmax=2650;
    int size = 200;
    // Make test data
    cpl_table * catalog = make_test_catalog();
    cpl_bivector * spectrum_err;
    cpl_bivector * spectrum = make_test_spectrum(catalog, wmin, wmax, size, &spectrum_err);
    cpl_polynomial * initial_guess = make_test_polynomial(wmin, wmax, size);
    int * orders = cpl_malloc(norders * sizeof(int));
    int * traces = cpl_malloc(norders * sizeof(int));
    int display = FALSE; // False
    const char * catalog_name = save_catalog(catalog);
    cpl_array * wave_error_init = cpl_array_new(2, CPL_TYPE_DOUBLE);
    cpl_array_set_double(wave_error_init, 0, 3.1);
    cpl_array_set_double(wave_error_init, 1, 3.5);

    cpl_array * wavelength_error;
    cpl_table * diagnostics;
    cpl_polynomial * wavelength;
    cpl_size power;

    cpl_size degree_x = 1; // polynomial degree in wavelength direction
    cpl_size degree_y = 2; // polynomial degree in order direction

    // Make lists that contain the same spectrum several times
    cpl_bivector ** spec = cpl_malloc(norders * sizeof(cpl_bivector*));
    cpl_bivector ** spec_err = cpl_malloc(norders * sizeof(cpl_bivector*));
    cpl_polynomial ** guess = cpl_malloc(norders * sizeof(cpl_polynomial*));
    cpl_array ** init_error = cpl_malloc(norders * sizeof(cpl_array*));

    for (i = 0; i < norders; i++){
        orders[i] = i;
        traces[i] = 1;

        spec[i] = spectrum;
        spec_err[i] = spectrum_err;
        guess[i] = initial_guess;
        init_error[i] = wave_error_init;
    }

    // Run function
    cpl_test(wavelength = cr2res_wave_2d(spec, spec_err, guess, init_error,
            orders, traces, norders, catalog_name, degree_x, degree_y,
            0.0, 0, 0, display, &wavelength_error, &diagnostics));

    // Check output
    cpl_polynomial_dump(wavelength, stdout);
    // #----- 2 dimensional polynomial -----
    // 1.dim.power  2.dim.power  coefficient
    //     0            0      2500
    //     1            0      0.75
    //     0            1      -1.04049e-12
    //     0            2      2.48754e-13
    //     1            2      1.07764e-16
    //     1            1      -4.31057e-16
    // #------------------------------------

    // first two are the linear component in x direction
    cpl_size idx[2] = {0, 0};
    cpl_test_abs(wmin, cpl_polynomial_get_coeff(wavelength, idx), 1e-3);
    idx[0] = 1;
    cpl_test_abs((wmax-wmin)/(double)size, cpl_polynomial_get_coeff(wavelength, idx), 1e-5);
    // all others should be 0 (or close to it), as there is no y dependance
    idx[0] = 0;
    idx[1] = 1;
    cpl_test_abs(0, cpl_polynomial_get_coeff(wavelength, idx), 1e-10);
    idx[0] = 0;
    idx[1] = 2;
    cpl_test_abs(0, cpl_polynomial_get_coeff(wavelength, idx), 1e-10);
    idx[0] = 1;
    idx[1] = 2;
    cpl_test_abs(0, cpl_polynomial_get_coeff(wavelength, idx), 1e-10);
    idx[0] = 1;
    idx[1] = 1;
    cpl_test_abs(0, cpl_polynomial_get_coeff(wavelength, idx), 1e-10);

    // Free Memory
    cpl_free(orders);
    cpl_free(traces);
    cpl_array_delete(wavelength_error);
    cpl_polynomial_delete(wavelength);
    cpl_table_delete(diagnostics);
    cpl_table_delete(catalog);
    cpl_bivector_delete(spectrum);
    cpl_bivector_delete(spectrum_err);
    cpl_polynomial_delete(initial_guess);
    cpl_array_delete(wave_error_init);

    cpl_free(spec);
    cpl_free(spec_err);
    cpl_free(guess);
    cpl_free(init_error);
}

static void test_cr2res_wave_etalon(void){

    cpl_bivector * spectrum;
    cpl_bivector * spectrum_err;
    cpl_array * error;
    cpl_polynomial * initial;
    cpl_polynomial * result;
    cpl_size power;

    double wmin = 500;
    double wmax = 600;
    int size = 2000;
    int degree = 1;

    spectrum = make_test_etalon_spectrum(size, 0.1, &spectrum_err);

    // Remove one peak
    for(int i = 100; i < 150; i++)
    {
        cpl_vector_set(cpl_bivector_get_y(spectrum), i, 0);
    }

    // Add a peak
    for(int i = 1462; i < 1493 ; i++)
    {
        cpl_vector_set(cpl_bivector_get_y(spectrum), i, 1);
    }

    initial = make_test_polynomial(wmin, wmax, size);

    error = cpl_array_new(2, CPL_TYPE_DOUBLE);
    cpl_array_set_double(error, 0, 3.1);
    cpl_array_set_double(error, 1, 3.5);

    cpl_test(result = cr2res_wave_etalon(spectrum, spectrum_err, initial, degree, &error));

    // these values obviously need to be changed if the number of degrees is changed
    power = 0;
    cpl_test_abs(cpl_polynomial_get_coeff(result, &power), wmin, 0.2);
    power = 1;
    cpl_test_abs(cpl_polynomial_get_coeff(result, &power), (wmax-wmin)/(double)size, 0.01);


    cpl_bivector_delete(spectrum);
    cpl_bivector_delete(spectrum_err);
    cpl_array_delete(error);
    cpl_polynomial_delete(initial);
    cpl_polynomial_delete(result);
}

static void test_cr2res_wave_etalon_other(void){

    cpl_bivector * spectrum;
    cpl_bivector * spectrum_err;
    cpl_array * error;
    cpl_polynomial * initial;
    cpl_polynomial * result;
    cpl_size power;

    double wmin = 500;
    double wmax = 600;
    int size = 2000;
    int degree = 1;

    spectrum = make_test_etalon_spectrum(size, 0.1, &spectrum_err);
    initial = make_test_polynomial(wmin, wmax, size);

    error = cpl_array_new(2, CPL_TYPE_DOUBLE);
    cpl_array_set_double(error, 0, 3.1);
    cpl_array_set_double(error, 1, 3.5);

    result = cr2res_wave_etalon(spectrum, spectrum_err, initial, degree, &error);

    // these values obviously need to be changed if the number of degrees is changed
    power = 0;
    cpl_test_abs(cpl_polynomial_get_coeff(result, &power), wmin, 0.2);
    power = 1;
    cpl_test_abs(cpl_polynomial_get_coeff(result, &power), (wmax-wmin)/(double)size, 0.01);


    cpl_bivector_delete(spectrum);
    cpl_bivector_delete(spectrum_err);
    cpl_array_delete(error);
    cpl_polynomial_delete(initial);
    cpl_polynomial_delete(result);
}

static void test_cr2res_wave_polys_1d_to_2d(void)
{
    int                 npolys = 10;
    cpl_polynomial  **  poly_1ds = cpl_malloc(npolys * sizeof(cpl_polynomial*));
    int             *   orders = cpl_malloc(npolys * sizeof(int));
    cpl_size            degree = 2;
    cpl_size            coef_pos;
    cpl_size        *   power = cpl_malloc(2 * sizeof(cpl_size));
    cpl_polynomial  *   res;


    for (cpl_size i = 0; i < npolys; i++)
    {
        orders[i] = i;
        poly_1ds[i] = cpl_polynomial_new(1);
        coef_pos = 0;
        cpl_polynomial_set_coeff(poly_1ds[i], &coef_pos, i * 10);
        coef_pos = 1;
        cpl_polynomial_set_coeff(poly_1ds[i], &coef_pos, 1);
    }

    cpl_test(res = cr2res_wave_polys_1d_to_2d(poly_1ds, orders, npolys, degree));

    // Check results
    cpl_test_eq(2, cpl_polynomial_get_dimension(res));
    cpl_test_eq(3, cpl_polynomial_get_degree(res));

    power[0] = 0;
    power[1] = 0;
    cpl_test_abs(0, cpl_polynomial_get_coeff(res, power), FLT_EPSILON);

    power[0] = 1;
    power[1] = 0;
    cpl_test_abs(1, cpl_polynomial_get_coeff(res, power), FLT_EPSILON);

    power[0] = 0;
    power[1] = 1;
    cpl_test_abs(10, cpl_polynomial_get_coeff(res, power), FLT_EPSILON);

    power[0] = 1;
    power[1] = 1;
    cpl_test_abs(0, cpl_polynomial_get_coeff(res, power), FLT_EPSILON);

    power[0] = 0;
    power[1] = 2;
    cpl_test_abs(0, cpl_polynomial_get_coeff(res, power), FLT_EPSILON);

    power[0] = 1;
    power[1] = 2;
    cpl_test_abs(0, cpl_polynomial_get_coeff(res, power), FLT_EPSILON);


    cpl_polynomial_delete(res);
    for (size_t i = 0; i < npolys; i++)
    {
        cpl_polynomial_delete(poly_1ds[i]);
    }
    cpl_free(poly_1ds);
    cpl_free(orders);
    cpl_free(power);
}

static void test_cr2res_wave_poly_2d_to_1d()
{
    cpl_polynomial  *   poly_2d = cpl_polynomial_new(2);
    int                 order = 1;
    cpl_size        *   power = cpl_malloc(2 * sizeof(cpl_size));
    cpl_size            pow_1d;
    cpl_polynomial  *   res;

    power[0] = 0;
    power[1] = 0;
    cpl_polynomial_set_coeff(poly_2d, power, 10);

    power[0] = 0;
    power[1] = 1;
    cpl_polynomial_set_coeff(poly_2d, power, 1);

    power[0] = 1;
    power[1] = 0;
    cpl_polynomial_set_coeff(poly_2d, power, 1);

    cpl_test(res = cr2res_wave_poly_2d_to_1d(poly_2d, order));

    cpl_test_eq(1, cpl_polynomial_get_dimension(res));
    pow_1d = 0;
    cpl_test_abs(11, cpl_polynomial_get_coeff(res, &pow_1d), DBL_EPSILON);
    pow_1d = 1;
    cpl_test_abs(1, cpl_polynomial_get_coeff(res, &pow_1d), DBL_EPSILON);

    cpl_polynomial_delete(poly_2d);
    cpl_polynomial_delete(res);
    cpl_free(power);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Use a simple wavelength range to check the estimate, which is just a linear polynomial
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_wave_estimate_compute(void)
{
    //define input
    // these values return "simple" results
    double wmin = 2000;
    double wmax = 4047;
    cpl_polynomial *res;

    //run test
    cpl_test_null(cr2res_wave_estimate_compute(-1, 1));
    cpl_test_null(cr2res_wave_estimate_compute(5, -1));
    cpl_test_null(cr2res_wave_estimate_compute(5, 1));

    cpl_test(res = cr2res_wave_estimate_compute(wmin, wmax));
    //test output
    cpl_size power = 0;
    cpl_test_abs(1999.0, cpl_polynomial_get_coeff(res, &power), DBL_EPSILON);
    power = 1;
    cpl_test_abs(1.0, cpl_polynomial_get_coeff(res, &power), DBL_EPSILON);

    //deallocate memory
    cpl_polynomial_delete(res);

    // Test invalid wavelength ranges
    wmin = 5000;
    wmax = 4047;
    cpl_test_null(cr2res_wave_estimate_compute(wmin, wmax));

    wmin = -10;
    wmax = 0.11;
    cpl_test_null(cr2res_wave_estimate_compute(wmin, wmax));
}

static void test_cr2res_wave_clean_spectrum(void)
{
    int ncols = 100;
    cpl_size pow = 0;
    cpl_vector        *   spec_intens = cpl_vector_new(ncols);
    cpl_polynomial    *   wl_poly = cpl_polynomial_new(1);
    cpl_bivector      *   catalog = cpl_bivector_new(5);
    double wl_error = 5.;

    cpl_vector * res = NULL;

    // Flat Spectrum (doesn't actually matter as long as its not NaN)
    for (int i = 0; i < ncols; i++){
        cpl_vector_set(spec_intens, i, 1);
    }

    // Add two lines to the catalog
    cpl_vector_set(cpl_bivector_get_x(catalog), 0, 20);
    cpl_vector_set(cpl_bivector_get_y(catalog), 0, 1);
    cpl_vector_set(cpl_bivector_get_x(catalog), 1, 23);
    cpl_vector_set(cpl_bivector_get_y(catalog), 1, 1);
    cpl_vector_set(cpl_bivector_get_x(catalog), 2, 23.2);
    cpl_vector_set(cpl_bivector_get_y(catalog), 2, 1);
    cpl_vector_set(cpl_bivector_get_x(catalog), 3, 23.3);
    cpl_vector_set(cpl_bivector_get_y(catalog), 3, 1);
    cpl_vector_set(cpl_bivector_get_x(catalog), 4, 50);
    cpl_vector_set(cpl_bivector_get_y(catalog), 4, 1);

    // Just linear spectrum
    pow = 1;
    cpl_polynomial_set_coeff(wl_poly, &pow, 1);

    cpl_test(res = cr2res_wave_clean_spectrum(spec_intens, wl_poly, catalog, wl_error));

    for (int i = 0; i < 15 ; i++){
        cpl_test(isnan(cpl_vector_get(res, i)));
    }
    for (int i = 15; i < 29 ; i++){
        cpl_test_abs(1, cpl_vector_get(res, i), DBL_EPSILON);
    }
    for (int i = 29; i < 45 ; i++){
        cpl_test(isnan(cpl_vector_get(res, i)));
    }
    for (int i = 45; i < 56 ; i++){
        cpl_test_abs(1, cpl_vector_get(res, i), DBL_EPSILON);
    }
    for (int i = 56; i < ncols ; i++){
        cpl_test(isnan(cpl_vector_get(res, i)));
    }

    // FILE * f = fopen("TEST_spec_clean.txt", "w");
    // cpl_vector_dump(res, f);
    // fclose(f);

    cpl_vector_delete(spec_intens);
    cpl_polynomial_delete(wl_poly);
    cpl_bivector_delete(catalog);
    cpl_vector_delete(res);

}

/*----------------------------------------------------------------------------*/
/**
  @brief    Create the 2d wavecal fit using etalon peaks
  @param    spectra         List of extracted spectra
  @param    spectra_err     List of extracted spectra errors
  @param    wavesol_init    List of Initial wavelength solutions
  @param    wavesol_init_err List of Initial wavelength error (can be NULL)
  @param    orders          List of orders of the various spectra
  @param    ninputs         Number of entries in the previous parameters
  @param    degree_x        The polynomial degree in x
  @param    degree_y        The polynomial degree in y
  @return   Wavelength solution, i.e. polynomial that translates pixel
            values to wavelength.

 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_wave_etalon_2d(){
    cpl_bivector        **  spectra = NULL;
    cpl_bivector        **  spectra_err = NULL;
    cpl_polynomial      **  wavesol_init = NULL;
    cpl_array           **  wavesol_init_err = NULL;
    int                 *   orders = NULL;
    int                 *   traces_nb = NULL;
    cpl_table           *  line_diagnostics = NULL;
    cpl_array           *  wavelength_error = NULL;
    int ninputs = 2;
    cpl_size degree_x = 1;
    cpl_size degree_y = 1;
    cpl_polynomial * result = NULL;
    cpl_size degree, degree2d[2];
    double c00, c01, c10, c11;
    double wave, freq;

    cpl_matrix * px;
    cpl_vector * py;
    px = cpl_matrix_new(1, CR2RES_DETECTOR_SIZE);
    py = cpl_vector_new(CR2RES_DETECTOR_SIZE);

    spectra = cpl_malloc(ninputs * sizeof(cpl_bivector*));
    spectra_err = cpl_malloc(ninputs * sizeof(cpl_bivector*));
    wavesol_init = cpl_malloc(ninputs* sizeof(cpl_polynomial*));
    wavesol_init_err = cpl_malloc(ninputs * sizeof(cpl_array*));
    orders = cpl_malloc(ninputs * sizeof(int));
    traces_nb = cpl_malloc(ninputs * sizeof(int));
    for (cpl_size i = 0; i < ninputs; i++)
    {
        orders[i] = i + 1;
        traces_nb[1] = 1;
        spectra[i] = cpl_bivector_new(CR2RES_DETECTOR_SIZE);
        spectra_err[i] = cpl_bivector_new(CR2RES_DETECTOR_SIZE);
        wavesol_init[i] = cpl_polynomial_new(1);
        wavesol_init_err[i] = cpl_array_new(2, CPL_TYPE_DOUBLE);

        for (cpl_size j = 0; j < CR2RES_DETECTOR_SIZE; j++)
        {
            cpl_vector_set(cpl_bivector_get_y(spectra[i]), j,
                1 + sin(j * 100 * CPL_MATH_PI / CR2RES_DETECTOR_SIZE));

            // Fit the wavelength guess
            freq = 60 - i * 5 - j * 0.001;
            wave = SPEED_OF_LIGHT / freq;
            cpl_matrix_set(px, 0, j, j);
            cpl_vector_set(py, j, wave);
        }

        degree = 3;
        cpl_polynomial_fit(wavesol_init[i], px, NULL, py, NULL, CPL_FALSE, NULL, &degree);
    }

    cpl_matrix_delete(px);
    cpl_vector_delete(py);

    result = cr2res_etalon_wave_2d(spectra, spectra_err, wavesol_init, wavesol_init_err,
                            orders, traces_nb, ninputs, degree_x, degree_y, 0, 0,
                            &wavelength_error, &line_diagnostics);

    cpl_test_nonnull(result);
    if (result != NULL){
        degree2d[0] = 0;
        degree2d[1] = 0;
        c00 = cpl_polynomial_get_coeff(result, degree2d);
        degree2d[0] = 0;
        degree2d[1] = 1;
        c01 = cpl_polynomial_get_coeff(result, degree2d);
        degree2d[0] = 1;
        degree2d[1] = 0;
        c10 = cpl_polynomial_get_coeff(result, degree2d);
        degree2d[0] = 1;
        degree2d[1] = 1;
        c11 = cpl_polynomial_get_coeff(result, degree2d);

        cpl_test_abs(c00 + c01, SPEED_OF_LIGHT / 60, 2);
        cpl_test_abs(c10, 0.1, 0.05);
        cpl_test_abs(c11, -0.01, 0.03);
    }

    for (cpl_size i = 0; i < ninputs; i++)
    {
        cpl_bivector_delete(spectra[i]);
        cpl_bivector_delete(spectra_err[i]);
        cpl_polynomial_delete(wavesol_init[i]);
        cpl_array_delete(wavesol_init_err[i]);
    }
    cpl_free(spectra);
    cpl_free(spectra_err);
    cpl_free(wavesol_init);
    cpl_free(wavesol_init_err);
    cpl_free(orders);
    cpl_free(traces_nb);

    if (wavelength_error != NULL) cpl_array_delete(wavelength_error);
    if (line_diagnostics != NULL) cpl_table_delete(line_diagnostics);
    if (result != NULL) cpl_polynomial_delete(result);

}

/*----------------------------------------------------------------------------*/
/**
  @brief    Run the Unit tests
 */
/*----------------------------------------------------------------------------*/
int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_DEBUG);

    test_cr2res_wave_1d();
    //test_cr2res_wave_2d();
    // test_cr2res_wave_etalon();
    test_cr2res_wave_polys_1d_to_2d();
    test_cr2res_wave_poly_2d_to_1d();
	test_cr2res_wave_estimate_compute();
    test_cr2res_wave_clean_spectrum();
    test_cr2res_wave_etalon_2d();
    return cpl_test_end(0);
}

/**@}*/
