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

#define CR2RES_DETECTOR_SIZE            2048

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static void test_cr2res_wave_line_fitting(void);
static void test_cr2res_wave_line_fitting_2d(void);
static void test_cr2res_wave_etalon(void);



/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_wave-test    Unit test of cr2res_wave
 *
 */
/*----------------------------------------------------------------------------*/
/**@{*/

static cpl_bivector * load_etalon()
{
    cpl_table * table;
    double * data;
    cpl_size n;
    cpl_bivector * result;
    cpl_vector * wave, *spec;

    if (access("./tests/etalon.fits", F_OK) != -1){
        table = cpl_table_load("./tests/etalon.fits", 1, 0);
    } else{
        table = cpl_table_load("./etalon.fits", 1, 0);
    }

    n = cpl_table_get_nrow(table);
    data = cpl_table_get_data_double(table, "99_04_SPEC");

    result = cpl_bivector_new(n);
    wave = cpl_bivector_get_x(result);    
    spec = cpl_bivector_get_y(result);

    for (cpl_size i = 0; i < n; i++){
        cpl_vector_set(spec, i, data[i]);
    }

    cpl_table_delete(table);

    return result;
}

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
static void test_cr2res_wave_line_fitting()
{
    double wmin=2500, wmax=2650;
    int size = 200;
    cpl_table * catalog = make_test_catalog();
    cpl_bivector * linelist;
    cpl_bivector * spectrum_err;
    cpl_bivector * spectrum = make_test_spectrum(catalog, wmin, wmax, size, &spectrum_err);
    cpl_polynomial * initial_guess = make_test_polynomial(wmin, wmax, size);
    int window_size = 30;
    int degree = 1;
    int display = 0; // False
    cpl_array * wave_error_init = cpl_array_new(2, CPL_TYPE_DOUBLE);
    cpl_array_set_double(wave_error_init, 0, 3.1);
    cpl_array_set_double(wave_error_init, 1, 3.5);

    cpl_vector * sigma_fit = cpl_vector_new(2);
    cpl_array * wavelength_error = cpl_array_new(2, CPL_TYPE_DOUBLE);
    cpl_polynomial * wavelength;
    cpl_size power;

    int len_linelist = cpl_table_get_nrow(catalog);
    cpl_vector * tmp_w = cpl_vector_wrap(len_linelist, cpl_table_get_data_double(catalog, CR2RES_COL_WAVELENGTH));
    cpl_vector * tmp_h = cpl_vector_wrap(len_linelist, cpl_table_get_data_double(catalog, CR2RES_COL_EMISSION));
    linelist = cpl_bivector_wrap_vectors(tmp_w, tmp_h);

    // bad inputs
    wavelength = cr2res_wave_line_fitting(NULL, spectrum_err, initial_guess, 
        wave_error_init, linelist, degree, display, &sigma_fit, &wavelength_error);
    cpl_test_null(wavelength);
    
    wavelength = cr2res_wave_line_fitting(spectrum, NULL, initial_guess,
        wave_error_init, linelist, degree, display, &sigma_fit, &wavelength_error);
    cpl_test_null(wavelength);

    wavelength = cr2res_wave_line_fitting(spectrum, spectrum_err, NULL,
        wave_error_init, linelist, degree, display, &sigma_fit, &wavelength_error);
    cpl_test_null(wavelength);

    wavelength = cr2res_wave_line_fitting(spectrum, spectrum_err, initial_guess,
        NULL, linelist, degree, display, &sigma_fit, &wavelength_error);
    cpl_test_null(wavelength);

    wavelength = cr2res_wave_line_fitting(spectrum, spectrum_err, initial_guess,
        wave_error_init, NULL, degree, display, &sigma_fit, &wavelength_error);
    cpl_test_null(wavelength);

    // optional NULL inputs
    wavelength = cr2res_wave_line_fitting(spectrum, spectrum_err, initial_guess,
        wave_error_init, linelist, degree, display, NULL, NULL);
    cpl_test_nonnull(wavelength);
    cpl_polynomial_delete(wavelength);

    // to many polynomial degrees
    wavelength = cr2res_wave_line_fitting(spectrum, spectrum_err, initial_guess,
        wave_error_init, linelist, 5, display, &sigma_fit, &wavelength_error);

    cpl_test_null(wavelength);
    cpl_test_nonnull(sigma_fit);
    cpl_test_nonnull(wavelength_error);

    // regular run
    wavelength = cr2res_wave_line_fitting(spectrum, spectrum_err, initial_guess, wave_error_init, linelist, degree, display, &sigma_fit, &wavelength_error);

    cpl_test_nonnull(wavelength);
    cpl_test_nonnull(sigma_fit);
    cpl_test_nonnull(wavelength_error);

    // these values obviously need to be changed if the number of degrees is changed
    power = 0;
    cpl_test_abs(cpl_polynomial_get_coeff(wavelength, &power), wmin, 0.2);
    power = 1;
    cpl_test_abs(cpl_polynomial_get_coeff(wavelength, &power), (wmax-wmin)/(double)size, 0.01);

    // Fitting two points with a first order polynomial -> perfect fit
    cpl_test_abs(cpl_array_get_double(wavelength_error, 0, NULL), 0, DBL_EPSILON);
    cpl_test_abs(cpl_array_get_double(wavelength_error, 1, NULL), 0, DBL_EPSILON);


    cpl_bivector_unwrap_vectors(linelist);
    cpl_vector_unwrap(tmp_h);
    cpl_vector_unwrap(tmp_w);
    cpl_array_delete(wavelength_error);
    cpl_polynomial_delete(wavelength);
    cpl_vector_delete(sigma_fit);
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
static void test_cr2res_wave_line_fitting_2d()
{   
    int i, norders = 5;
    double wmin=2500, wmax=2650;
    int size = 200;
    // Make test data
    cpl_table * catalog = make_test_catalog();
    cpl_bivector * linelist;
    cpl_bivector * spectrum_err;
    cpl_bivector * spectrum = make_test_spectrum(catalog, wmin, wmax, size, &spectrum_err);
    cpl_polynomial * initial_guess = make_test_polynomial(wmin, wmax, size);
    cpl_vector * orders = cpl_vector_new(norders);
    int display = FALSE; // False
    cpl_array * wave_error_init = cpl_array_new(2, CPL_TYPE_DOUBLE);
    cpl_array_set_double(wave_error_init, 0, 3.1);
    cpl_array_set_double(wave_error_init, 1, 3.5);

    cpl_vector * sigma_fit = cpl_vector_new(2);
    cpl_array * wavelength_error = cpl_array_new(2, CPL_TYPE_DOUBLE);
    cpl_polynomial * wavelength;
    cpl_size power;

    int len_linelist = cpl_table_get_nrow(catalog);
    cpl_vector * tmp_w = cpl_vector_wrap(len_linelist, 
        cpl_table_get_data_double(catalog, CR2RES_COL_WAVELENGTH));
    cpl_vector * tmp_h = cpl_vector_wrap(len_linelist, 
        cpl_table_get_data_double(catalog, CR2RES_COL_EMISSION));
    linelist = cpl_bivector_wrap_vectors(tmp_w, tmp_h);

    cpl_size * degree = cpl_malloc(2 * sizeof(cpl_size));
    degree[0] = 1; // polynomial degree in wavelength direction
    degree[1] = 2; // polynomial degree in order direction

    // Make lists that contain the same spectrum several times
    cpl_bivector ** spec = cpl_malloc(norders * sizeof(cpl_bivector*));
    cpl_bivector ** spec_err = cpl_malloc(norders * sizeof(cpl_bivector*));
    cpl_polynomial ** guess = cpl_malloc(norders * sizeof(cpl_polynomial*));
    const cpl_array ** init_error = cpl_malloc(norders * sizeof(cpl_array*));

    for (i = 0; i < norders; i++){
        cpl_vector_set(orders, i, i);

        spec[i] = spectrum;
        spec_err[i] = spectrum_err;
        guess[i] = initial_guess;
        init_error[i] = wave_error_init;
    }

    // Run function
    wavelength = cr2res_wave_line_fitting_2D(spec, spec_err, guess, init_error,
            linelist, orders, norders, degree, display, NULL, NULL);

    // Check output
    // cpl_polynomial_dump(wavelength, stdout);
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
    cpl_test_abs(wmin, cpl_polynomial_get_coeff(wavelength, idx), 1e-10);
    idx[0] = 1;
    cpl_test_abs((wmax-wmin)/(double)size, cpl_polynomial_get_coeff(wavelength, idx), 1e-10);
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
    cpl_free(degree);
    cpl_vector_delete(orders);
    cpl_bivector_unwrap_vectors(linelist);
    cpl_vector_unwrap(tmp_h);
    cpl_vector_unwrap(tmp_w);
    cpl_array_delete(wavelength_error);
    cpl_polynomial_delete(wavelength);
    cpl_vector_delete(sigma_fit);
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



/*----------------------------------------------------------------------------*/
/**
  @brief    Use two identical orders, with two lines each, and check that the result is still linear
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_wave_line_fitting_2d()
{   
    int i, norders = 5;
    double wmin=2500, wmax=2650;
    int size = 200;
    // Make test data
    cpl_table * catalog = make_test_catalog();
    cpl_bivector * linelist;
    cpl_bivector * spectrum_err;
    cpl_bivector * spectrum = make_test_spectrum(catalog, wmin, wmax, size, &spectrum_err);
    cpl_polynomial * initial_guess = make_test_polynomial(wmin, wmax, size);
    cpl_vector * orders = cpl_vector_new(norders);
    int display = FALSE; // False
    cpl_array * wave_error_init = cpl_array_new(2, CPL_TYPE_DOUBLE);
    cpl_array_set_double(wave_error_init, 0, 3.1);
    cpl_array_set_double(wave_error_init, 1, 3.5);

    cpl_vector * sigma_fit = cpl_vector_new(2);
    cpl_array * wavelength_error = cpl_array_new(2, CPL_TYPE_DOUBLE);
    cpl_polynomial * wavelength;
    cpl_size power;

    int len_linelist = cpl_table_get_nrow(catalog);
    cpl_vector * tmp_w = cpl_vector_wrap(len_linelist, 
        cpl_table_get_data_double(catalog, CR2RES_COL_WAVELENGTH));
    cpl_vector * tmp_h = cpl_vector_wrap(len_linelist, 
        cpl_table_get_data_double(catalog, CR2RES_COL_EMISSION));
    linelist = cpl_bivector_wrap_vectors(tmp_w, tmp_h);

    cpl_size * degree = cpl_malloc(2 * sizeof(cpl_size));
    degree[0] = 1; // polynomial degree in wavelength direction
    degree[1] = 2; // polynomial degree in order direction

    // Make lists that contain the same spectrum several times
    cpl_bivector ** spec = cpl_malloc(norders * sizeof(cpl_bivector*));
    cpl_bivector ** spec_err = cpl_malloc(norders * sizeof(cpl_bivector*));
    cpl_polynomial ** guess = cpl_malloc(norders * sizeof(cpl_polynomial*));
    const cpl_array ** init_error = cpl_malloc(norders * sizeof(cpl_array*));

    for (i = 0; i < norders; i++){
        cpl_vector_set(orders, i, i);

        spec[i] = spectrum;
        spec_err[i] = spectrum_err;
        guess[i] = initial_guess;
        init_error[i] = wave_error_init;
    }

    // Run function
    wavelength = cr2res_wave_line_fitting_2D(spec, spec_err, guess, init_error,
            linelist, orders, norders, degree, display, NULL, NULL);

    // Check output
    // cpl_polynomial_dump(wavelength, stdout);
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
    cpl_test_abs(wmin, cpl_polynomial_get_coeff(wavelength, idx), 1e-10);
    idx[0] = 1;
    cpl_test_abs((wmax-wmin)/(double)size, cpl_polynomial_get_coeff(wavelength, idx), 1e-10);
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
    cpl_free(degree);
    cpl_vector_delete(orders);
    cpl_bivector_unwrap_vectors(linelist);
    cpl_vector_unwrap(tmp_h);
    cpl_vector_unwrap(tmp_w);
    cpl_array_delete(wavelength_error);
    cpl_polynomial_delete(wavelength);
    cpl_vector_delete(sigma_fit);
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



/*----------------------------------------------------------------------------*/
/**
  @brief    Run the Unit tests
 */
/*----------------------------------------------------------------------------*/
int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_DEBUG);

    test_cr2res_wave_line_fitting();
    test_cr2res_wave_line_fitting_2d();
    test_cr2res_wave_etalon();

    return cpl_test_end(0);
}

/**@}*/
