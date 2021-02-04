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
#include <cr2res_pol.h>

#define CR2RES_DETECTOR_SIZE            2048

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/


static void test_cr2res_pol_demod_stokes(void);
static void test_cr2res_pol_demod_stokes_different_wavelengths(void);
static void test_cr2res_pol_demod_stokes_zeros(void);

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_pol-test    Unit test of cr2res_pol
 *
 */
/*----------------------------------------------------------------------------*/
/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Demodulate extracted spectra into Stokes parameter
  @param    speclist  Array of bivectors
  @param    n         Length of speclist, needs to be 8
  @return   cpl_bivector with Stokes parameter spectrum (P/I) and error, needs
            to be de-allocated by caller.

  The input list of spectra needs to come in this order:
    1u, 1d, 2u , 2d, 3u, 3d, 4u, 4d
  i.e. first exposure upper beam, then down, then second exposure etc.

  Demodulation formula is P/I = (R^1/4 - 1) / (R^1/4 + 1) with
    R = 1u/2u * 2d/1d * 3u/4u * 4d/3d
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_pol_demod_stokes(){
    int n = 8; // as it has to be
    cpl_vector ** wl = cpl_malloc(n * sizeof(cpl_vector*));
    cpl_vector ** intens = cpl_malloc(n * sizeof(cpl_vector*));
    cpl_vector ** errors = cpl_malloc(n * sizeof(cpl_vector*));
    cpl_bivector * pol;
    double value = 0;
    double error = 1. / sqrt(8);

    // spec = sigma = 1
    // -> pol = 0
    // -> sigma pol = 1 / sqrt(8)
    cpl_bivector * spec = cpl_bivector_new(CR2RES_DETECTOR_SIZE);
    for (cpl_size i = 0; i < CR2RES_DETECTOR_SIZE; i++)
    {
        cpl_vector_set(cpl_bivector_get_x(spec), i, 1);
        cpl_vector_set(cpl_bivector_get_y(spec), i, 1);
    }

    // Set the wavelength
    // for this test they all have the same wavelength scale
    double wmin = 2000;
    double wmax = 3000;
    double step = (wmax - wmin) / CR2RES_DETECTOR_SIZE;
    for (int i = 0; i < n; i++){
      wl[i] = cpl_vector_new(CR2RES_DETECTOR_SIZE);
      for (int j = 0; j < CR2RES_DETECTOR_SIZE; j++){
        cpl_vector_set(wl[i], j, wmin + j * step);
      }
    }
    
    // just use the same spectrum 8 times
    for (int i = 0; i < n; i++)
    { 
        intens[i] = cpl_bivector_get_x(spec);
        errors[i] = cpl_bivector_get_y(spec);
    }

    cpl_test_nonnull(pol = cr2res_pol_demod_stokes(intens, wl, errors, n));

    // Check that the results are as expected    
    for (cpl_size i = 0; i < CR2RES_DETECTOR_SIZE; i++)
    {
        cpl_test_abs(value, cpl_vector_get(cpl_bivector_get_x(pol), i),
                    DBL_EPSILON);
        cpl_test_abs(error, cpl_vector_get(cpl_bivector_get_y(pol), i),
                    DBL_EPSILON);
    }

    // Clean memory
    for (int i = 0; i < n; i++){
      cpl_vector_delete(wl[i]);
    }
    cpl_bivector_delete(pol);
    cpl_bivector_delete(spec);
    cpl_free(intens);
    cpl_free(wl);
    cpl_free(errors);
}


static void test_cr2res_pol_demod_stokes_different_wavelengths(){
    int n = 8; // as it has to be
    cpl_vector ** wl = cpl_malloc(n * sizeof(cpl_vector*));
    cpl_vector ** intens = cpl_malloc(n * sizeof(cpl_vector*));
    cpl_vector ** errors = cpl_malloc(n * sizeof(cpl_vector*));
    cpl_bivector * pol;
    double value = 0;
    double error = 1. / sqrt(8);

    cpl_size i, j;

    // spec = sigma = 1
    // -> pol = 0
    // -> sigma pol = 1 / sqrt(8)
    cpl_bivector * spec = cpl_bivector_new(CR2RES_DETECTOR_SIZE);
    for (cpl_size i = 0; i < CR2RES_DETECTOR_SIZE; i++)
    {
        cpl_vector_set(cpl_bivector_get_x(spec), i, 1);
        cpl_vector_set(cpl_bivector_get_y(spec), i, 1);
    }

    // Set the wavelength
    // each spectra has a slightly shifted scale
    double wmin = 2000;
    double wmax = 3000;
    double step = (wmax - wmin) / CR2RES_DETECTOR_SIZE;
    double diff_between_spectra = 20;
    for (i = 0; i < n; i++){
      wl[i] = cpl_vector_new(CR2RES_DETECTOR_SIZE);
      for (j = 0; j < CR2RES_DETECTOR_SIZE; j++){
        cpl_vector_set(wl[i], j, wmin + j * step + i * diff_between_spectra);
      }
    }
    
    // just use the same spectrum 8 times
    // demod should not modify them (it makes copies)
    for (int i = 0; i < n; i++)
    { 
        intens[i] = cpl_bivector_get_x(spec);
        errors[i] = cpl_bivector_get_y(spec);
    }

    cpl_test_nonnull(pol = cr2res_pol_demod_stokes(intens, wl, errors, n));


    // These values depend on the chosen wavelength grid
    // make sure to adjust them appropiately when you change those values
    cpl_size xmin = 287;
    cpl_size xmax = 1760;
    // Check that the results are as expected
    // First we have some NANs
    i = 0;
    for (; i < xmin; i++){
      cpl_test(isnan(cpl_vector_get(cpl_bivector_get_x(pol), i)));
      cpl_test(isnan(cpl_vector_get(cpl_bivector_get_y(pol), i)));
    }
    // Then 0s as in the same wavelengths case
    for (; i < xmax + 1; i++)
    {
        cpl_test_abs(value, cpl_vector_get(cpl_bivector_get_x(pol), i), 
                    DBL_EPSILON);
        cpl_test_abs(error, cpl_vector_get(cpl_bivector_get_y(pol), i), 
                    DBL_EPSILON);
    }
    // Followed by more NANs
    for (; i < CR2RES_DETECTOR_SIZE; i++){
      cpl_test(isnan(cpl_vector_get(cpl_bivector_get_x(pol), i)));
      cpl_test(isnan(cpl_vector_get(cpl_bivector_get_y(pol), i)));
    }

    // Clean memory
    for (j = 0; j < n; j++){
      cpl_vector_delete(wl[j]);
    }
    cpl_bivector_delete(pol);
    cpl_bivector_delete(spec);
    cpl_free(intens);
    cpl_free(wl);
    cpl_free(errors);
}

static void test_cr2res_pol_demod_stokes_zeros(){
    int n = 8; // as it has to be
    cpl_vector ** wl = cpl_malloc(n * sizeof(cpl_vector*));
    cpl_vector ** intens = cpl_malloc(n * sizeof(cpl_vector*));
    cpl_vector ** errors = cpl_malloc(n * sizeof(cpl_vector*));
    cpl_bivector * pol;

    // spec = sigma = 1
    // -> pol = 0
    // -> sigma pol = 1 / sqrt(8)
    cpl_bivector * spec = cpl_bivector_new(CR2RES_DETECTOR_SIZE);
    for (cpl_size i = 0; i < CR2RES_DETECTOR_SIZE; i++)
    {
        cpl_vector_set(cpl_bivector_get_x(spec), i, 0);
        cpl_vector_set(cpl_bivector_get_y(spec), i, 0);
    }

    // Set the wavelength
    // for this test they all have the same wavelength scale
    double wmin = 2000;
    double wmax = 3000;
    double step = (wmax - wmin) / CR2RES_DETECTOR_SIZE;
    for (int i = 0; i < n; i++){
      wl[i] = cpl_vector_new(CR2RES_DETECTOR_SIZE);
      for (int j = 0; j < CR2RES_DETECTOR_SIZE; j++){
        cpl_vector_set(wl[i], j, wmin + j * step);
      }
    }
    
    // just use the same spectrum 8 times
    for (int i = 0; i < n; i++)
    { 
        intens[i] = cpl_bivector_get_x(spec);
        errors[i] = cpl_bivector_get_y(spec);
    }

    cpl_test_nonnull(pol = cr2res_pol_demod_stokes(intens, wl, errors, n));

    // Check that the results are as expected    
    for (cpl_size i = 0; i < CR2RES_DETECTOR_SIZE; i++)
    {
        cpl_test(isnan(cpl_vector_get(cpl_bivector_get_x(pol), i))); 
        cpl_test(isnan(cpl_vector_get(cpl_bivector_get_y(pol), i)));
    }

    // Clean memory
    for (int i = 0; i < n; i++){
      cpl_vector_delete(wl[i]);
    }
    cpl_bivector_delete(pol);
    cpl_bivector_delete(spec);
    cpl_free(intens);
    cpl_free(wl);
    cpl_free(errors);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Run the Unit tests
 */
/*----------------------------------------------------------------------------*/
int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_DEBUG);

    test_cr2res_pol_demod_stokes();
    test_cr2res_pol_demod_stokes_different_wavelengths();
    test_cr2res_pol_demod_stokes_zeros();

    return cpl_test_end(0);
}

/**@}*/
