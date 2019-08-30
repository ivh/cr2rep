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
#include <cpl.h>
#include <cr2res_qc.h>

#include "cr2res_trace.h"
#include "cr2res_dfs.h"
#include "cr2res_pfits.h"
#include "cr2res_io.h"


#define WLEN_BEGIN(i) ({char s[20]; sprintf(s, CR2RES_HEADER_WLEN_BEGIN, i); s;})
#define WLEN_END(i) ({char s[20]; sprintf(s, CR2RES_HEADER_WLEN_END, i); s;})
#define WLEN_CENY(i) ({char s[20]; sprintf(s, CR2RES_HEADER_WLEN_CENY, i); s;})

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static void test_cr2res_dark_qc_ron(void);
static void test_cr2res_qc_flat_lamp_ints(void);
static void test_cr2res_qc_flat_mean_level(void);
static void test_cr2res_qc_flat_mean_med_flux(void);
static void test_cr2res_qc_obs_nodding_slit_psf(void);

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_qc-test    Unit test of cr2res_qc
 *
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief Create a default table with test data
  @param
  @return test table with some data, needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
static cpl_table *create_test_table()
{
    cpl_table *traces;
    cpl_array *array, *slit_fraction, *wave, *wave_err, *slit_a, *slit_b, 
              *slit_c;
    int poly_order, norders;
    cpl_propertylist *hdr = cpl_propertylist_new();
    cpl_propertylist *main_header = cpl_propertylist_new();
    char * extname;

    /* Initialise */
    poly_order = 2;
    norders = 9;

    /* NULL Input */
    traces = cpl_table_new(norders);
    cpl_table_new_column_array(traces, CR2RES_COL_ALL, CPL_TYPE_DOUBLE, 
            poly_order);
    cpl_table_new_column_array(traces, CR2RES_COL_UPPER, CPL_TYPE_DOUBLE, 
            poly_order);
    cpl_table_new_column_array(traces, CR2RES_COL_LOWER, CPL_TYPE_DOUBLE, 
            poly_order);
    cpl_table_new_column(traces, CR2RES_COL_ORDER, CPL_TYPE_INT);
    cpl_table_new_column(traces, CR2RES_COL_TRACENB, CPL_TYPE_INT);

    cpl_table_new_column_array(traces, CR2RES_COL_WAVELENGTH, CPL_TYPE_DOUBLE, 
            2);
    cpl_table_new_column_array(traces, CR2RES_COL_WAVELENGTH_ERROR, 
            CPL_TYPE_DOUBLE, 2);
    cpl_table_new_column_array(traces, CR2RES_COL_SLIT_CURV_A, CPL_TYPE_DOUBLE,
            3);
    cpl_table_new_column_array(traces, CR2RES_COL_SLIT_CURV_B, CPL_TYPE_DOUBLE,
            3);
    cpl_table_new_column_array(traces, CR2RES_COL_SLIT_CURV_C, CPL_TYPE_DOUBLE,
            3);
    cpl_table_new_column_array(traces, CR2RES_COL_SLIT_FRACTION, 
            CPL_TYPE_DOUBLE, 3);

    /*
                 All|               Upper|               Lower|  Order
 54.3148, 0.00738623|  110.379, 0.0159501|1.63328, -8.95809e-05|      1
  226.289, 0.0169145|  311.885, 0.0167957|   139.106, 0.017396|      2
  437.881, 0.0172448|  524.126, 0.0171958|  350.398, 0.0170009|      3
  660.954, 0.0178211|  747.931, 0.0179547|  572.986, 0.0177137|      4
  897.266, 0.0185319|   985.06, 0.0186544|  808.823, 0.0185517|      5
   1148.31, 0.019388|  1237.01, 0.0193813|  1059.39, 0.0194215|      6
  1415.88, 0.0202877|  1505.34, 0.0202763|  1326.43, 0.0202534|      7
   1701.85, 0.021292|  1792.03, 0.0213662|  1611.65, 0.0212178|      8
  1982.59, 0.0111388|2047.88, 8.66878e-05|   1917.3, 0.0221835|      9

  */
    double all_1[] = {54.3148, 226.289, 437.881, 660.954, 897.266,
                      1148.31, 1415.88, 1701.85, 1982.59};
    double all_2[] = {0.00738623, 0.0169145, 0.0172448, 0.0178211,
                      0.0185319, 0.019388, 0.0202877, 0.021292, 0.0111388};
    double upper_1[] = {110.379, 311.885, 524.126, 747.931, 985.06,
                        1237.01, 1505.34, 1792.03, 2047.88};
    double upper_2[] = {0.0159501, 0.0167957, 0.0171958, 0.0179547,
                        0.0186544, 0.0193813, 0.0202763, 0.0213662, 8.66878e-05};
    double lower_1[] = {1.63328, 139.106, 350.398, 572.986, 808.823,
                        1059.39, 1326.43, 1611.65, 1917.3};
    double lower_2[] = {-8.95809e-05, 0.017396, 0.0170009, 0.0177137,
                        0.0185517, 0.0194215, 0.0202534, 0.0212178, 0.0221835};
    array = cpl_array_new(poly_order, CPL_TYPE_DOUBLE);
    slit_fraction = cpl_array_new(3, CPL_TYPE_DOUBLE);
    cpl_array_set_double(slit_fraction, 0, 0);
    cpl_array_set_double(slit_fraction, 1, 0.5);
    cpl_array_set_double(slit_fraction, 2, 1);
    wave = cpl_array_new(2, CPL_TYPE_DOUBLE);
    cpl_array_set_double(wave, 0, 9.45e2);
    cpl_array_set_double(wave, 1, 3.13e-3);
    wave_err = cpl_array_new(2, CPL_TYPE_DOUBLE);
    cpl_array_set_double(wave_err, 0, 5e-2);
    cpl_array_set_double(wave_err, 1, 5e-2);

    slit_a = cpl_array_new(3, CPL_TYPE_DOUBLE);
    cpl_array_set(slit_a, 0, 0);
    cpl_array_set(slit_a, 1, 1);
    cpl_array_set(slit_a, 2, 0);
    slit_b = cpl_array_new(3, CPL_TYPE_DOUBLE);
    cpl_array_set(slit_b, 0, 0);
    cpl_array_set(slit_b, 1, 0);
    cpl_array_set(slit_b, 2, 0);
    slit_c = cpl_array_new(3, CPL_TYPE_DOUBLE);
    cpl_array_set(slit_c, 0, 0);
    cpl_array_set(slit_c, 1, 0);
    cpl_array_set(slit_c, 2, 0);

    for (int i = 0; i < norders; i++)
    {
        cpl_array_set(array, 0, all_1[i]);
        cpl_array_set(array, 1, all_2[i]);
        cpl_table_set_array(traces, CR2RES_COL_ALL, i, array);
        cpl_array_set(array, 0, upper_1[i]);
        cpl_array_set(array, 1, upper_2[i]);
        cpl_table_set_array(traces, CR2RES_COL_UPPER, i, array);
        cpl_array_set(array, 0, lower_1[i]);
        cpl_array_set(array, 1, lower_2[i]);
        cpl_table_set_array(traces, CR2RES_COL_LOWER, i, array);
        cpl_table_set(traces, CR2RES_COL_ORDER, i, cr2res_io_convert_idx_to_order(i + 1));
        cpl_table_set(traces, CR2RES_COL_TRACENB, i, 1);
    
        cpl_table_set_array(traces, CR2RES_COL_SLIT_FRACTION, i, slit_fraction);
        cpl_table_set_array(traces, CR2RES_COL_WAVELENGTH, i, wave);
        cpl_table_set_array(traces, CR2RES_COL_WAVELENGTH_ERROR, i, wave_err);
        cpl_table_set_array(traces, CR2RES_COL_SLIT_CURV_A, i, slit_a);
        cpl_table_set_array(traces, CR2RES_COL_SLIT_CURV_B, i, slit_b);
        cpl_table_set_array(traces, CR2RES_COL_SLIT_CURV_C, i, slit_c);
    }

    extname = cr2res_io_create_extname(1, 1);
    cpl_propertylist_append_string(hdr, CR2RES_HEADER_EXTNAME, extname);

    double ceny[] = {1994.0945859223, 1723.67027599362, 1436.61298619847, 
                     1168.0222016174, 915.8934665223831, 678.542785839296,
                     454.468576982434, 242.388497032926, 63.5899165277783};
    double begin[] = {1756.78720770673, 1703.55123171562, 1653.44678372399,
                      1606.20544704616, 1561.58862907265, 1519.38353098961,
                      1479.3997538583, 1441.46642683629, -1};
    double end[] = {1768.81709603003, 1715.21657796851, 1664.76903155768,
                    1617.2042020846, 1572.2818631378, 1529.78775872867,
                    1489.53018613055, 1451.3371044349, -1};

    for (int i = 0; i < 9; i++)
    {
        cpl_propertylist_append_double(hdr, WLEN_CENY(i), ceny[i]);
        cpl_propertylist_append_double(hdr, WLEN_BEGIN(i), begin[i]);
        cpl_propertylist_append_double(hdr, WLEN_END(i), end[i]);
    }

    
    cpl_propertylist_append_int(main_header, CR2RES_HEADER_DECKER_POS, CR2RES_DECKER_2_4);

    cpl_table_save(traces, main_header, hdr, "test_table.fits", CPL_IO_CREATE);

    cpl_array_delete(array);
    cpl_array_delete(slit_fraction);
    cpl_array_delete(wave);
    cpl_array_delete(wave_err);
    cpl_array_delete(slit_a);
    cpl_array_delete(slit_b);
    cpl_array_delete(slit_c);

    cpl_propertylist_delete(hdr);
    cpl_propertylist_delete(main_header);
    cpl_free(extname);
    return traces;
}


/*----------------------------------------------------------------------------*/
/**
  @brief Create a default image to test against based on the default table data
  @param
  @return default test image(2048x2048), needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
static cpl_image *create_test_image(void)
{
    cpl_table *traces;
    cpl_image *trace_ima;
    //cpl_image *extract;

    traces = create_test_table();
    trace_ima = cr2res_trace_gen_image(traces, 2048, 2048);
    cpl_table_delete(traces);

    return trace_ima;
}

static void test_cr2res_dark_qc_ron(void)
{
    cpl_image *ima1, *ima2;
    int hsize, nsamples, ndit;
    double res;

    ima1 = cpl_image_new(100, 100, CPL_TYPE_DOUBLE);
    ima2 = cpl_image_new(100, 100, CPL_TYPE_DOUBLE);

    hsize = -1;
    nsamples = -1;
    ndit = 2;

    for(cpl_size x = 1; x <= 100; x++)
    {
        for(cpl_size y = 1; y <= 100; y++)
        {
            cpl_image_set(ima1, x, y, 1);
            cpl_image_set(ima2, x, y, 1);
        }
    }
    
    res = cr2res_dark_qc_ron(ima1, ima2, hsize, nsamples, ndit);
    cpl_test_error(CPL_ERROR_NONE);

    cpl_image_delete(ima1);
    cpl_image_delete(ima2);
}

static void test_cr2res_qc_flat_lamp_ints()
{
    cpl_image * ima;
    double res;

    ima = cpl_image_new(10, 10, CPL_TYPE_DOUBLE);
    for(cpl_size x = 1; x <= 10; x++)
    {
        for(cpl_size y = 1; y <= 10; y++)
        {
            cpl_image_set(ima, x, y, 1);
        }
    }

    cpl_test_eq(-1, cr2res_qc_flat_lamp_ints(NULL));

    cpl_test(res = cr2res_qc_flat_lamp_ints(ima));
    cpl_test_abs(res, 10*10, DBL_EPSILON);

    cpl_image_delete(ima);
}

static void test_cr2res_qc_flat_mean_level()
{
    cpl_image * ima;
    double res;

    ima = cpl_image_new(10, 10, CPL_TYPE_DOUBLE);
    for(cpl_size x = 1; x <= 10; x++)
    {
        for(cpl_size y = 1; y <= 10; y++)
        {
            cpl_image_set(ima, x, y, 1);
        }
    }

    cpl_test_eq(-1, cr2res_qc_flat_mean_level(NULL));

    cpl_test(res = cr2res_qc_flat_mean_level(ima));
    cpl_test_abs(res, 1, DBL_EPSILON);

    cpl_image_delete(ima);
}

static void test_cr2res_qc_flat_mean_med_flux()
{
    cpl_image * ima;
    double mean, median;
    double res;

    ima = cpl_image_new(10, 10, CPL_TYPE_DOUBLE);
    for(cpl_size x = 1; x <= 10; x++)
    {
        for(cpl_size y = 1; y <= 10; y++)
        {
            cpl_image_set(ima, x, y, 1);
        }
    }

    cpl_test_eq(-1, cr2res_qc_flat_mean_med_flux(NULL, &mean, &median));

    cpl_test_eq(0, cr2res_qc_flat_mean_med_flux(ima, &mean, &median));
    cpl_test_abs(mean, 1, DBL_EPSILON);
    cpl_test_abs(median, 1, DBL_EPSILON);


    cpl_image_delete(ima);
}

static void test_cr2res_qc_obs_nodding_slit_psf()
{
    int nrow = 100;
    // values for creating the default data
    double x0 = nrow / 2;
    double A = 1; // total area of the slitfunc should always be 1
    double offset = 0;
    double sigma = nrow / 20;
    double value = 0;
    double fwhm = -1;

    char col1[] = "01_01_SLIT_FUNC";
    char col2[] = "02_01_SLIT_FUNC";

    cpl_table * slitfu = cpl_table_new(nrow);
    cpl_table_new_column(slitfu, col1, CPL_TYPE_DOUBLE);
    cpl_table_new_column(slitfu, col2, CPL_TYPE_DOUBLE);

    for (cpl_size i = 0; i < nrow; i++)
    {
        // area / sqrt(2 pi sigma^2) * exp( -(x - x0)^2/(2 sigma^2)) + offset
        value = A / sqrt(CPL_MATH_2_PI * sigma * sigma) * exp( -(i - x0) * (i - x0) / (2. * sigma * sigma)) + offset;
        cpl_table_set_double(slitfu, col1, i, value);
        cpl_table_set_double(slitfu, col2, i, value);
    }

    cpl_test_abs(-1, cr2res_qc_obs_nodding_slit_psf(NULL), DBL_EPSILON);
    cpl_test(fwhm = cr2res_qc_obs_nodding_slit_psf(slitfu));
    cpl_test_abs(fwhm, 2.355 * sigma, FLT_EPSILON);

    cpl_table_delete(slitfu);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Run the Unit tests
 */
/*----------------------------------------------------------------------------*/
int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_DEBUG);

    test_cr2res_dark_qc_ron();
    test_cr2res_qc_flat_lamp_ints();
    test_cr2res_qc_flat_mean_level();
    test_cr2res_qc_flat_mean_med_flux();
    test_cr2res_qc_obs_nodding_slit_psf();

    return cpl_test_end(0);
}

/**@}*/
