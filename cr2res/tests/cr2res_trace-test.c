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
#include <cr2res_dfs.h>
#include <cr2res_pfits.h>
#include <cr2res_trace.h>
#include <cr2res_trace.c>
#include "cr2res_io.h"




/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static void test_cr2res_trace(void);
static void test_cr2res_trace_clean(void);
static void test_cr2res_trace_gen_image(void);
static void test_cr2res_trace_get_order_idx_values(void);
static void test_cr2res_trace_get_ycen(void);
static void test_cr2res_trace_get_height(void);
static void test_cr2res_trace_compute_middle(void);
static void test_cr2res_trace_compute_height(void);
static void test_cr2res_trace_get_trace_ypos(void);
static void test_cr2res_trace_signal_detect(void);
static void test_cr2res_trace_fit_traces(void);
static void test_cr2res_trace_fit_trace(void);
static void test_cr2res_trace_convert_cluster_to_labels(void);
static void test_cr2res_trace_convert_labels_to_cluster(void);
static void test_cr2res_trace_clean_blobs(void);
static void test_cr2res_trace_extract_edges(void);
static void test_cr2res_trace_new_slit_fraction(void);
static void test_cr2res_trace_add_extra_columns(void);
static void test_cr2res_get_trace_table_index(void);
static void test_cr2res_get_trace_wave_poly(void);

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_trace-test    Unit test of cr2res_trace
 *
 */
/*----------------------------------------------------------------------------*/

/**@{*/
#define WLEN_BEGIN(i) ({char s[20]; sprintf(s, CR2RES_HEADER_WLEN_BEGIN, i); s;})
#define WLEN_END(i) ({char s[20]; sprintf(s, CR2RES_HEADER_WLEN_END, i); s;})
#define WLEN_CENY(i) ({char s[20]; sprintf(s, CR2RES_HEADER_WLEN_CENY, i); s;})


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
        cpl_table_set(traces, CR2RES_COL_ORDER, i,
                cr2res_io_convert_order_idx_to_idxp(i + 1));
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

    cpl_table_save(traces, main_header, hdr, "TEST_table.fits", CPL_IO_CREATE);

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
    traces = create_test_table();
    trace_ima = cr2res_trace_gen_image(traces, 2048, 2048);
    cpl_image_add_scalar(trace_ima, 1.);
    cpl_image_multiply_scalar(trace_ima, 10.);
    cpl_table_delete(traces);
    cpl_image_save(trace_ima, "TEST_trace.fits", CPL_TYPE_INT, NULL, CPL_IO_CREATE);
    return trace_ima;
}

static cpl_table *create_cluster_table(void)
{
    // for each pixel in a cluster
    // XS, YS - coordinates
    // Cluster - which cluster that pixel belongs to
    int data[] = {2, 2, 0, 0, 0,
                  0, 0, 0, 1, 1,
                  0, 1, 1, 1, 1,
                  1, 1, 1, 1, 1,
                  1, 1, 1, 1, 0};
    // flip data, so the image will be the same
    int data_inverse[] = {1, 1, 1, 1, 0,
                          1, 1, 1, 1, 1,
                          0, 1, 1, 1, 1,
                          0, 0, 0, 1, 1,
                          2, 2, 0, 0, 0};

    int xs[] = {4, 5, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 1, 2};
    int ys[] = {4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 5, 5};
    int clusters[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2};

    cpl_table *cluster = cpl_table_new(17);
    cpl_table_wrap_int(cluster, xs, CR2RES_COL_XS);
    cpl_table_wrap_int(cluster, ys, CR2RES_COL_YS);
    cpl_table_wrap_int(cluster, clusters, CR2RES_COL_CLUSTERS);
    return cluster;
}

/*----------------------------------------------------------------------------*/
/**
  @brief  Test missing input data cases and simple regular case for order tracing with trace

        cpl_image       *   ima,
        int                 smooth_x,
        int                 smooth_y,
        double              threshold,
        int                 opening,
        int                 degree,
        int                 min_cluster)
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_trace(void)
{
    cpl_image *trace_ima = create_test_image();
    cpl_table *out;
    const cpl_array *all;

    double threshold = 5;

    cpl_test_null(cr2res_trace(NULL, 1.0, 1.0, threshold, 1, 2, 10));
    cpl_test_null(cr2res_trace(trace_ima, -1.0, 1.0, threshold, 1, 2, 10));
    cpl_test_null(cr2res_trace(trace_ima, 1.0, 1.0, threshold, -1, 2, 10));
    cpl_test_null(cr2res_trace(trace_ima, 1.0, 1.0, threshold, 1, 2, -10));

    cpl_test(out = cr2res_trace(trace_ima, 1.0, 1.0, threshold, 1, 2, 10));
    
    all = cpl_table_get_array(out, CR2RES_COL_ALL, 0);

    cpl_table_save(out, NULL, NULL, "TEST_table2.fits", CPL_IO_CREATE);
    cpl_table_delete(out);
    cpl_image_delete(trace_ima);
    return;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Test the removal of small blobs from 6x6 sample grid
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_trace_clean(void)
{

    /* test_cr2res_trace() ; */
    /* test_cr2res_trace_clean() ; */
    cpl_binary data[] = {1, 1, 0, 0, 0, 0,
                         1, 1, 0, 0, 0, 1,
                         0, 0, 0, 0, 0, 0,
                         0, 1, 0, 1, 1, 0,
                         0, 0, 0, 0, 0, 0,
                         1, 0, 0, 0, 1, 1};
    cpl_mask *mask = cpl_mask_wrap(6, 6, data);
    int opening = 0;
    int min_cluster = 3;
    cpl_mask *res;
    // without opening:
    cpl_binary data_cmp[] = {1, 1, 0, 0, 0, 0,
                             1, 1, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0};

    // with opening
    // cpl_binary data_cmp[] = {1, 1, 0, 0, 0, 0,
    //                          1, 1, 0, 0, 0, 0,
    //                          0, 0, 0, 0, 0, 0,
    //                          0, 1, 1, 1, 1, 0,
    //                          0, 0, 0, 0, 0, 0,
    //                          0, 0, 0, 0, 0, 0};
    cpl_mask *cmp = cpl_mask_wrap(6, 6, data_cmp);

    cpl_test_null(cr2res_trace_clean(NULL, opening, min_cluster));
    cpl_test(res = cr2res_trace_clean(mask, opening, min_cluster));
    cpl_mask_save(res, "TEST_res.fits", NULL, CPL_IO_CREATE);
    cpl_test_eq_mask(cmp, res);
    cpl_mask_unwrap(mask);
    cpl_mask_unwrap(cmp);
    cpl_mask_delete(res);
}

/**
  @brief    Test the generated image by comparing a small 10x10 patch to the expected result
 */
static void test_cr2res_trace_gen_image(void)
{
    int nx = 2048;
    int ny = 2048;
    cpl_image *res;
    cpl_image *extract;

    // order is inversed, first row is lowest in FITS
    int data[10 * 10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                         1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                         1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                         1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                         1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                         1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                         -1, -1, -1, -1, -1, -1, -1, 1, 1, 1,
                         -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                         -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
                         -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

    cpl_image *cmp = cpl_image_wrap(10, 10, CPL_TYPE_INT, data);
    cpl_table *trace = create_test_table();

    //run test
    // Null tests
    cpl_test_null(cr2res_trace_gen_image(NULL, nx, ny));
    cpl_test_null(cr2res_trace_gen_image(trace, 0, ny));
    cpl_test_null(cr2res_trace_gen_image(trace, nx, 0));
    // normal run
    cpl_test(res = cr2res_trace_gen_image(trace, nx, ny));
    //test output ?
    // This is more of a regression test than anything, but the image looks good
    extract = cpl_image_extract(res, 32, 105, 41, 114);
    cpl_test_image_abs(extract, cmp, DBL_EPSILON);

    //deallocate memory
    cpl_image_delete(res);
    cpl_image_delete(extract);
    cpl_image_unwrap(cmp);
    cpl_table_delete(trace);
}

/*----------------------------------------------------------------------------*/
/**
  @brief   Extracted order numbers are compared with known input table
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_trace_get_order_idx_values(void)
{
    //define input
    cpl_table *trace = create_test_table();
    int nb_orders;
    int *res;
    /* test_cr2res_trace_compute_height() ; */
    cpl_test_null(cr2res_trace_get_order_idx_values(NULL, &nb_orders));
    cpl_test_null(cr2res_trace_get_order_idx_values(trace, NULL));
    cpl_test(res = cr2res_trace_get_order_idx_values(trace, &nb_orders));
    //test output
    cpl_test_eq(nb_orders, 9);
    for (int i = 0; i < 9; i++)
    {
        cpl_test_eq(res[i], i + 1);
    }
    //deallocate memory
    cpl_table_delete(trace);
    cpl_free(res);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Run trace_get_ycen with linear input data and check that results are also linear
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_trace_get_ycen(void)
{
    //define input
    cpl_table *trace = create_test_table();
    cpl_size order_nb = 3;
    cpl_size trace_nb = 1;
    int size = CR2RES_DETECTOR_SIZE;
    cpl_vector *cmp;
    cpl_vector *res;
    //run test
    cpl_test_null(cr2res_trace_get_ycen(NULL, order_nb, trace_nb, size));
    cpl_test_null(cr2res_trace_get_ycen(trace, 100, trace_nb, size));
    cpl_test_null(cr2res_trace_get_ycen(trace, order_nb, 5, size));
    cpl_test_null(cr2res_trace_get_ycen(trace, order_nb, trace_nb, -1));

    // assemble comparison vector
    double data[size];
    for (int i = 0; i < size; i++)
    {
        // values from test table
        data[i] = 437.881 + (i + 1) * 0.0172448;
    }
    cmp = cpl_vector_wrap(size, data);

    // Run that sould not fail, compare output
    cpl_test(res = cr2res_trace_get_ycen(trace, order_nb, trace_nb, size));
    cpl_test_vector_abs(cmp, res, size * DBL_EPSILON);

    //deallocate memory
    cpl_vector_delete(res);
    cpl_vector_unwrap(cmp);

    cpl_table_delete(trace);
}

/*----------------------------------------------------------------------------*/
/**
  @brief     Run trace_get_height and compare with expected result
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_trace_get_height(void)
{
    //define input
    cpl_table *trace = create_test_table();
    cpl_size order_nb = 3;
    cpl_size trace_nb = 1;
    int res;

    //run test
    cpl_test_eq(cr2res_trace_get_height(NULL, order_nb, trace_nb), -1);
    cpl_test_eq(cr2res_trace_get_height(trace, 20, trace_nb), -1);
    cpl_test_eq(cr2res_trace_get_height(trace, order_nb, 5), -1);

    cpl_test(res = cr2res_trace_get_height(trace, order_nb, trace_nb));
    //test output
    cpl_test_eq(res, 174.1271552); // value analytically from test table

    //deallocate memory
    cpl_table_delete(trace);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Use two linear polynomials as input data and compare trace_compute_middle with expected result
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_trace_compute_middle(void)
{
    //define input
    cpl_polynomial *trace1 = cpl_polynomial_new(1);
    cpl_polynomial *trace2 = cpl_polynomial_new(1);

    cpl_size power = 0;
    cpl_polynomial_set_coeff(trace1, &power, 10.);
    cpl_polynomial_set_coeff(trace2, &power, 20.);

    power = 1;
    cpl_polynomial_set_coeff(trace1, &power, 1.);
    cpl_polynomial_set_coeff(trace2, &power, 3.);

    int vector_size = 10;
    cpl_vector *res;
    double data[] = {17., 19., 21., 23., 25., 27., 29., 31., 33., 35.};
    cpl_vector *cmp = cpl_vector_wrap(10, data);

    // run test
    cpl_test_null(cr2res_trace_compute_middle(NULL, trace2, vector_size));
    cpl_test_null(cr2res_trace_compute_middle(trace1, NULL, vector_size));
    cpl_test_null(cr2res_trace_compute_middle(trace1, trace2, -1));

    cpl_test(res = cr2res_trace_compute_middle(trace1, trace2, vector_size));
    //test output
    cpl_test_vector_abs(res, cmp, DBL_EPSILON);

    //deallocate memory
    cpl_vector_unwrap(cmp);
    cpl_vector_delete(res);
    cpl_polynomial_delete(trace1);
    cpl_polynomial_delete(trace2);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Use two linear polynomials as input data and compare trace_compute_height with expected result
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_trace_compute_height(void)
{
    //define input
    cpl_polynomial *trace1 = cpl_polynomial_new(1);
    cpl_polynomial *trace2 = cpl_polynomial_new(1);

    cpl_size power = 0;
    cpl_polynomial_set_coeff(trace1, &power, 10.);
    cpl_polynomial_set_coeff(trace2, &power, 20.);

    power = 1;
    cpl_polynomial_set_coeff(trace1, &power, 1.);
    cpl_polynomial_set_coeff(trace2, &power, 3.);

    int vector_size = 10;
    int res;
    int cmp = 21; // 10 + 10.2222

    //run test
    cpl_test_eq(-1, cr2res_trace_compute_height(NULL, trace2, vector_size));
    cpl_test_eq(-1, cr2res_trace_compute_height(trace1, NULL, vector_size));
    cpl_test_eq(-1, cr2res_trace_compute_height(trace1, trace2, 0));

    cpl_test(res = cr2res_trace_compute_height(trace1, trace2, vector_size));
    //test output
    cpl_test_eq(res, cmp);

    //deallocate memory
    cpl_polynomial_delete(trace1);
    cpl_polynomial_delete(trace2);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Check result of trace_get_trace_ypos against expected result from input data
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_trace_get_trace_ypos(void)
{
    //define input
    cpl_table *traces = create_test_table();
    int idx = 1;
    double res;

    //run test
    cpl_test_abs(-1.0, cr2res_trace_get_trace_ypos(NULL, idx), DBL_EPSILON);
    cpl_test_abs(-1.0, cr2res_trace_get_trace_ypos(NULL, -10), DBL_EPSILON);
    cpl_test_abs(-1.0, cr2res_trace_get_trace_ypos(NULL, 100), DBL_EPSILON);

    cpl_test(res = cr2res_trace_get_trace_ypos(traces, idx));
    //test output
    // 226.289 + 1024 * 0.0169145
    cpl_test_abs(res, 243.609448, DBL_EPSILON);
    //deallocate memory
    cpl_table_delete(traces);
}

/*----------------------------------------------------------------------------*/
/**
  @brief   Add ORDER, TRACE, WAVELENGTH columns to the plain trace table
  @param    traces          The plain traces table
  @param    file_for_wl     File used for WL information
  @param    det_nr          Detector
  @return   0 if ok
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_trace_add_extra_columns(void)
{
    // define input
    cpl_table *traces = create_test_table();
    // remove columns to fill
    cpl_table_erase_column(traces, CR2RES_COL_ORDER);
    cpl_table_erase_column(traces, CR2RES_COL_TRACENB);
    cpl_table_erase_column(traces, CR2RES_COL_WAVELENGTH);
    cpl_table_erase_column(traces, CR2RES_COL_WAVELENGTH_ERROR);
    cpl_table_erase_column(traces, CR2RES_COL_SLIT_FRACTION);
    cpl_table_erase_column(traces, CR2RES_COL_SLIT_CURV_A);
    cpl_table_erase_column(traces, CR2RES_COL_SLIT_CURV_B);
    cpl_table_erase_column(traces, CR2RES_COL_SLIT_CURV_C);

    cpl_table *tmp = cpl_table_duplicate(traces);
    char *file_for_wl = "TEST_table.fits";
    int det_nr = 1;
    int res;
    const cpl_array *wl;
    double cmp1[] = {-1, 1441.46160481499, 1479.3948049417, 1519.37844831851, 1561.58340521624, 1606.20007393671, 1653.44125258191, 1703.54553296318, 1756.78133086827};
    double cmp2[] = {0, 0.00482202129878357, 0.00494891659611621, 0.00508267109871028, 0.00522385640701021, 0.00537310944721049, 0.00553114207801171, 0.00569875244401075, 0.00587683845788956};

    //run test
    cpl_test_eq(-1, cr2res_trace_add_extra_columns(NULL, file_for_wl, det_nr));

    cpl_test_eq(-1, cr2res_trace_add_extra_columns(tmp, "invalid_path", det_nr));
    cpl_table_delete(tmp);
    tmp = cpl_table_duplicate(traces);

    cpl_test_eq(-1, cr2res_trace_add_extra_columns(tmp, file_for_wl, 10));
    cpl_table_delete(tmp);
    tmp = cpl_table_duplicate(traces);

    cpl_test_eq(-1, cr2res_trace_add_extra_columns(tmp, file_for_wl, -1));
    cpl_table_delete(tmp);
    tmp = cpl_table_duplicate(traces);

    cpl_test_eq(0, cr2res_trace_add_extra_columns(tmp, file_for_wl, det_nr));
    //test output
    cpl_table_save(tmp, NULL, NULL, "TEST_new_table.fits", CPL_IO_CREATE);

    // Check that all columns are there
    cpl_test(cpl_table_has_column(tmp, CR2RES_COL_ORDER));
    cpl_test(cpl_table_has_column(tmp, CR2RES_COL_TRACENB));
    cpl_test(cpl_table_has_column(tmp, CR2RES_COL_WAVELENGTH));
    cpl_test(cpl_table_has_column(tmp, CR2RES_COL_WAVELENGTH_ERROR));
    cpl_test(cpl_table_has_column(tmp, CR2RES_COL_SLIT_FRACTION));
    cpl_test(cpl_table_has_column(tmp, CR2RES_COL_SLIT_CURV_A));
    cpl_test(cpl_table_has_column(tmp, CR2RES_COL_SLIT_CURV_B));
    cpl_test(cpl_table_has_column(tmp, CR2RES_COL_SLIT_CURV_C));

    // Check wavelength
    for (int i = 0; i < 9; i++)
    {
        wl = cpl_table_get_array(tmp, CR2RES_COL_WAVELENGTH, i);
        cpl_test_abs(cpl_array_get(wl, 0, 0), cmp1[i], FLT_EPSILON);
        cpl_test_abs(cpl_array_get(wl, 1, 0), cmp2[i], FLT_EPSILON);
    }

    // TODO check that original data is the same (tmp == traces)
    // TODO check other new columns

    //deallocate memory
    cpl_table_delete(traces);
    cpl_table_delete(tmp);
}

/*----------------------------------------------------------------------------*/
/**
  @brief Test a small 10x10 patch of the result to expected result based on trace table polynomials
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_trace_new_slit_fraction(void)
{
    //define input
    cpl_table *trace_table = create_test_table();
    int norder = cpl_table_get_nrow(trace_table); // as defined in create test_table
    cpl_array * new_slit_fraction;
    cpl_table *res;
    cpl_array *array;
    const cpl_array *carray;

    // test for NULL input
    cpl_test_null(cr2res_trace_new_slit_fraction(NULL, new_slit_fraction));
    cpl_test_null(cr2res_trace_new_slit_fraction(trace_table, NULL));

    // test with too small new slit fraction, should get NULL
    new_slit_fraction = cpl_array_new(1, CPL_TYPE_DOUBLE);
    cpl_test_null(cr2res_trace_new_slit_fraction(trace_table, new_slit_fraction));
    cpl_array_delete(new_slit_fraction);


    // test with new slit fraction == old slit fractions
    new_slit_fraction = cpl_array_new(3, CPL_TYPE_DOUBLE);
    cpl_array_set_double(new_slit_fraction, 0, 0);
    cpl_array_set_double(new_slit_fraction, 1, 0.5);
    cpl_array_set_double(new_slit_fraction, 2, 1);

    cpl_test(res = cr2res_trace_new_slit_fraction(trace_table, new_slit_fraction));
    cpl_test_eq(cpl_table_get_nrow(res), norder);
    for(cpl_size i = 0; i < norder; i++)
    {
        cpl_test_array_abs(cpl_table_get_array(res, CR2RES_COL_SLIT_FRACTION, i), 
                cpl_table_get_array(trace_table, CR2RES_COL_SLIT_FRACTION, i), FLT_EPSILON);
    
        cpl_test_array_abs(cpl_table_get_array(res, CR2RES_COL_UPPER, i), 
                cpl_table_get_array(trace_table, CR2RES_COL_UPPER, i), FLT_EPSILON);

        cpl_test_array_abs(cpl_table_get_array(res, CR2RES_COL_ALL, i), 
            cpl_table_get_array(trace_table, CR2RES_COL_ALL, i), FLT_EPSILON);

        cpl_test_array_abs(cpl_table_get_array(res, CR2RES_COL_LOWER, i), 
                cpl_table_get_array(trace_table, CR2RES_COL_LOWER, i), FLT_EPSILON);

        cpl_test_array_abs(cpl_table_get_array(res, CR2RES_COL_WAVELENGTH, i), 
                cpl_table_get_array(trace_table, CR2RES_COL_WAVELENGTH, i), FLT_EPSILON);

        cpl_test_array_abs(cpl_table_get_array(res, CR2RES_COL_SLIT_CURV_A, i), 
                cpl_table_get_array(trace_table, CR2RES_COL_SLIT_CURV_A, i), FLT_EPSILON);

        cpl_test_array_abs(cpl_table_get_array(res, CR2RES_COL_SLIT_CURV_B, i), 
                cpl_table_get_array(trace_table, CR2RES_COL_SLIT_CURV_B, i), FLT_EPSILON);
        

        cpl_test_array_abs(cpl_table_get_array(res, CR2RES_COL_SLIT_CURV_C, i), 
                cpl_table_get_array(trace_table, CR2RES_COL_SLIT_CURV_C, i), FLT_EPSILON);
    }
    cpl_array_delete(new_slit_fraction);


    // test with slit fractions other than old ones
    new_slit_fraction = cpl_array_new(3, CPL_TYPE_DOUBLE);
    cpl_array_set_double(new_slit_fraction, 0, 0.2);
    cpl_array_set_double(new_slit_fraction, 1, 0.4);
    cpl_array_set_double(new_slit_fraction, 2, 0.6);

    // Set the first order of the table to some easy to estimate values, 1, 2, 3
    array = cpl_array_new(2, CPL_TYPE_DOUBLE);
    cpl_array_set_double(array, 0, 3);
    cpl_table_set_array(trace_table, CR2RES_COL_UPPER, 0, array);
    cpl_array_set_double(array, 0, 2);
    cpl_table_set_array(trace_table, CR2RES_COL_ALL, 0, array);
    cpl_array_set_double(array, 0, 1);
    cpl_table_set_array(trace_table, CR2RES_COL_LOWER, 0, array);
    cpl_array_delete(array);

    cpl_table_delete(res);
    cpl_test(res = cr2res_trace_new_slit_fraction(trace_table, new_slit_fraction));

    carray = cpl_table_get_array(res, CR2RES_COL_UPPER, 0);
    cpl_test_abs(cpl_array_get_double(carray, 0, NULL), 2.2, FLT_EPSILON);
    cpl_test_abs(cpl_array_get_double(carray, 1, NULL), 0, FLT_EPSILON);
    carray = cpl_table_get_array(res, CR2RES_COL_ALL, 0);
    cpl_test_abs(cpl_array_get_double(carray, 0, NULL), 1.8, FLT_EPSILON);
    cpl_test_abs(cpl_array_get_double(carray, 1, NULL), 0, FLT_EPSILON);
    carray = cpl_table_get_array(res, CR2RES_COL_LOWER, 0);
    cpl_test_abs(cpl_array_get_double(carray, 0, NULL), 1.4, FLT_EPSILON);
    cpl_test_abs(cpl_array_get_double(carray, 1, NULL), 0, FLT_EPSILON);

    // deallocate memory
    cpl_table_delete(res);
    cpl_table_delete(trace_table);
    cpl_array_delete(new_slit_fraction);
}



/*----------------------------------------------------------------------------*/
/**
  @brief    Test a small 10x10 patch of the result to expected result based on trace image
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_trace_signal_detect(void)
{
    //define input
    cpl_image *image = create_test_image();
    int trace_sep = 150;
    double smoothfactor = 1;
    double thresh = 0.5;
    cpl_mask *sub;

    cpl_binary data2[10 * 10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    cpl_mask *cmp = cpl_mask_wrap(10, 10, data2);

    cpl_mask *res;

    //run test
    cpl_test_null(cr2res_trace_signal_detect(NULL, trace_sep, smoothfactor, thresh, 0, 0));
    cpl_test_null(cr2res_trace_signal_detect(image, -10, smoothfactor, thresh, 0, 0));
    cpl_test_null(cr2res_trace_signal_detect(image, trace_sep, -1, thresh, 0, 0));
    //cpl_test_null(cr2res_trace_signal_detect(image, trace_sep, smoothfactor, 50000));

    cpl_test(res = cr2res_trace_signal_detect(image, trace_sep, smoothfactor, thresh, 0, 0));
    //test output
    sub = cpl_mask_extract(res, 32, 105, 41, 114);


/* TODO : Why is this failing ? */
    /* cpl_test_eq_mask(sub, cmp); */

    //deallocate memory
    cpl_image_delete(image);
    cpl_mask_delete(res);
    cpl_mask_unwrap(cmp);
    cpl_mask_delete(sub);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Compare fitted trace polynomial to expected result from simple input data with few points
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_trace_fit_traces(void)
{
    //define input
    int xs[] = {4, 5, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 1, 2};
    int ys[] = {4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 5, 5};
    int clusters[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2};

    cpl_table *cluster = cpl_table_new(17);
    cpl_table_wrap_int(cluster, xs, CR2RES_COL_XS);
    cpl_table_wrap_int(cluster, ys, CR2RES_COL_YS);
    cpl_table_wrap_int(cluster, clusters, CR2RES_COL_CLUSTERS);

    int degree = 1;
    cpl_table *res;
    const cpl_array *arr;

    //run test
    cpl_test_null(cr2res_trace_fit_traces(NULL, degree));
    cpl_test_null(cr2res_trace_fit_traces(cluster, -5));

    cpl_test(res = cr2res_trace_fit_traces(cluster, degree));
    //test output
    // use tolerance * 10, because reading the data in the array is not that precise ?

    arr = cpl_table_get_array(res, CR2RES_COL_UPPER, 0);
    cpl_test_abs(1.7, cpl_array_get(arr, 0, NULL), DBL_EPSILON * 10);
    cpl_test_abs(0.5, cpl_array_get(arr, 1, NULL), DBL_EPSILON * 10);
    arr = cpl_table_get_array(res, CR2RES_COL_UPPER, 1);
    cpl_test_abs(5, cpl_array_get(arr, 0, NULL), DBL_EPSILON * 10);
    cpl_test_abs(0, cpl_array_get(arr, 1, NULL), DBL_EPSILON * 10);

    arr = cpl_table_get_array(res, CR2RES_COL_LOWER, 0);
    cpl_test_abs(0.6, cpl_array_get(arr, 0, NULL), DBL_EPSILON * 10);
    cpl_test_abs(0.2, cpl_array_get(arr, 1, NULL), DBL_EPSILON * 10);
    arr = cpl_table_get_array(res, CR2RES_COL_LOWER, 1);
    cpl_test_abs(5, cpl_array_get(arr, 0, NULL), DBL_EPSILON * 10);
    cpl_test_abs(0, cpl_array_get(arr, 1, NULL), DBL_EPSILON * 10);

    arr = cpl_table_get_array(res, CR2RES_COL_ALL, 0);
    cpl_test_abs(1.15151515151515, cpl_array_get(arr, 0, NULL), DBL_EPSILON * 10);
    cpl_test_abs(0.348484848484849, cpl_array_get(arr, 1, NULL), DBL_EPSILON * 10);
    arr = cpl_table_get_array(res, CR2RES_COL_ALL, 1);
    cpl_test_abs(5, cpl_array_get(arr, 0, NULL), DBL_EPSILON * 10);
    cpl_test_abs(0, cpl_array_get(arr, 1, NULL), DBL_EPSILON * 10);

    //deallocate memory
    cpl_table_unwrap(cluster, CR2RES_COL_XS);
    cpl_table_unwrap(cluster, CR2RES_COL_YS);
    cpl_table_unwrap(cluster, CR2RES_COL_CLUSTERS);
    cpl_table_delete(cluster);
    cpl_table_delete(res);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Test fit of trace polynomial based on large number of points in a single trace and compare to expected result
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_trace_fit_trace(void)
{
    //define input
    // use only cluster 3 from test image
    cpl_image *test_image = create_test_image();
    cpl_table *all = cr2res_trace_convert_labels_to_cluster(test_image);
    cpl_table_and_selected_int(all, CR2RES_COL_CLUSTERS, CPL_EQUAL_TO, 30);
    cpl_table *table = cpl_table_extract_selected(all);

    int degree = 2;
    cpl_array *res;

    //run test
    cpl_test_null(cr2res_trace_fit_trace(NULL, degree));
    cpl_test_null(cr2res_trace_fit_trace(table, -10));

    cpl_test(res = cr2res_trace_fit_trace(table, degree));
    //test output
    cpl_array_dump(res, 0, 2, NULL);
    // input was 437.881, 0.0172448
    cpl_test_abs(cpl_array_get(res, 0, NULL), 225.007, 1.2);
    cpl_test_abs(cpl_array_get(res, 1, NULL), 0.0172448, 0.0002);

    //deallocate memory
    cpl_image_delete(test_image);
    cpl_table_delete(all);
    cpl_table_delete(table);
    cpl_array_delete(res);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Check that small 5x5 image is created properly based on simple input dataset
 */
/*----------------------------------------------------------------------------*/
// static void test_cr2res_trace_convert_cluster_to_labels(void)
// {
//     //define input
//     //cpl_table *cluster = create_cluster_table();

//     int xs[] = {4, 5, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 1, 2};
//     int ys[] = {4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 5, 5};
//     int clusters[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2};

//     cpl_table *cluster = cpl_table_new(17);
//     cpl_table_wrap_int(cluster, xs, CR2RES_COL_XS);
//     cpl_table_wrap_int(cluster, ys, CR2RES_COL_YS);
//     cpl_table_wrap_int(cluster, clusters, CR2RES_COL_CLUSTERS);

//     // flip data, so the image will be the same
//     int data_inverse[] = {1, 1, 1, 1, 0,
//                           1, 1, 1, 1, 1,
//                           0, 1, 1, 1, 1,
//                           0, 0, 0, 1, 1,
//                           2, 2, 0, 0, 0};

//     int nx = 5;
//     int ny = 5;
//     cpl_image *res;
//     cpl_image *cmp = cpl_image_wrap(nx, ny, CPL_TYPE_INT, data_inverse);

//     //run test
//     cpl_test_null(cr2res_trace_convert_cluster_to_labels(NULL, nx, ny));
//     cpl_test_null(cr2res_trace_convert_cluster_to_labels(cluster, 0, ny));
//     cpl_test_null(cr2res_trace_convert_cluster_to_labels(cluster, nx, 0));

//     cpl_test(res = cr2res_trace_convert_cluster_to_labels(cluster, nx, ny));
//     //test output
//     //cpl_image_save(res, "TEST_labels.fits", CPL_TYPE_INT, NULL, CPL_IO_CREATE);
//     cpl_test_image_abs(res, cmp, 0);

//     //deallocate memory
//     cpl_table_unwrap(cluster, CR2RES_COL_XS);
//     cpl_table_unwrap(cluster, CR2RES_COL_YS);
//     cpl_table_unwrap(cluster, CR2RES_COL_CLUSTERS);
//     cpl_table_delete(cluster);
//     cpl_image_delete(res);
//     cpl_image_unwrap(cmp);
// }

/*----------------------------------------------------------------------------*/
/**
  @brief    Check that small 5x5 image is converted to correct table of data points
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_trace_convert_labels_to_cluster(void)
{
    //define input
    int data_inverse[] = {1, 1, 1, 1, 0,
                          1, 1, 1, 1, 1,
                          0, 1, 1, 1, 1,
                          0, 0, 0, 1, 1,
                          2, 2, 0, 0, 0};
    cpl_image *labels = cpl_image_wrap(5, 5, CPL_TYPE_INT, data_inverse);
    cpl_table *res;

    int xs[] = {1, 2, 3, 4, 1, 2, 3, 4, 5, 2, 3, 4, 5, 4, 5, 1, 2};
    int ys[] = {1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5};
    int clusters[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2};

    //run test
    cpl_test_null(cr2res_trace_convert_labels_to_cluster(NULL));
    cpl_test(res = cr2res_trace_convert_labels_to_cluster(labels));
    //test output
    // cpl_table_save(res, NULL, NULL, "TEST_convert.fits", CPL_IO_CREATE);

    for (int i = 0; i < 17; i++)
    {
        cpl_test_eq(xs[i], cpl_table_get(res, CR2RES_COL_XS, i, NULL));
        cpl_test_eq(ys[i], cpl_table_get(res, CR2RES_COL_YS, i, NULL));
        cpl_test_eq(clusters[i], cpl_table_get(res, CR2RES_COL_CLUSTERS, i, NULL));
    }

    //deallocate memory
    cpl_image_unwrap(labels);
    cpl_table_delete(res);
}

/*----------------------------------------------------------------------------*/
/**
  @brief   Check the removal of small clusters in small 4x4 patch

 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_trace_clean_blobs(void)
{
    //define input
    cpl_binary data[] = {1, 1, 0, 0,
                         1, 1, 0, 0,
                         0, 0, 1, 0,
                         0, 0, 1, 0};
    cpl_mask *mask = cpl_mask_wrap(4, 4, data);
    int min_cluster = 3;
    cpl_binary data2[] = {1, 1, 0, 0,
                          1, 1, 0, 0,
                          0, 0, 0, 0,
                          0, 0, 0, 0};
    cpl_mask *cmp = cpl_mask_wrap(4, 4, data2);
    cpl_mask *res;

    //run test
    cpl_test_null(cr2res_trace_clean_blobs(NULL, min_cluster));
    cpl_test_null(cr2res_trace_clean_blobs(mask, -1));

    // if min_cluster <= 0 nothing changes
    cpl_test(res = cr2res_trace_clean_blobs(mask, 0));
    cpl_test_eq_mask(res, mask);
    cpl_mask_delete(res);

    cpl_test(res = cr2res_trace_clean_blobs(mask, min_cluster));
    //test output
    //small blob of size 2 removed
    cpl_test_eq_mask(res, cmp);

    //deallocate memory
    cpl_mask_unwrap(mask);
    cpl_mask_unwrap(cmp);
    cpl_mask_delete(res);
}

/*----------------------------------------------------------------------------*/
/**
  @brief  Test that edge pixels are identified correctly in small sample case
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_trace_extract_edges(void)
{
    //define input
    int xs[] = {4, 5, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4};
    int ys[] = {4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1};
    int clusters[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    cpl_table *pixels_table = cpl_table_new(15);
    cpl_table_wrap_int(pixels_table, xs, CR2RES_COL_XS);
    cpl_table_wrap_int(pixels_table, ys, CR2RES_COL_YS);
    cpl_table_wrap_int(pixels_table, clusters, CR2RES_COL_CLUSTERS);

    cpl_table *edge_lower_table;
    cpl_table *edge_upper_table;

    int cmp_xs_lower[] = {5, 1, 2, 3, 4};
    int cmp_ys_lower[] = {2, 1, 1, 1, 1};
    int cmp_cluster_lower[] = {1, 1, 1, 1, 1};

    int cmp_xs_upper[] = {4, 5, 2, 3, 1};
    int cmp_ys_upper[] = {4, 4, 3, 3, 2};
    int cmp_cluster_upper[] = {1, 1, 1, 1, 1};

    //run test
    cpl_test_eq(-1, cr2res_trace_extract_edges(NULL, &edge_lower_table, &edge_upper_table));
    cpl_test_eq(-1, cr2res_trace_extract_edges(pixels_table, NULL, &edge_upper_table));
    cpl_test_eq(-1, cr2res_trace_extract_edges(pixels_table, &edge_lower_table, NULL));

    cpl_test_eq(0, cr2res_trace_extract_edges(pixels_table, &edge_lower_table, &edge_upper_table));
    //test output
    for (int i = 0; i < 5; i++)
    {
        cpl_test_eq(cpl_table_get(edge_lower_table, CR2RES_COL_XS, i, NULL), cmp_xs_lower[i]);
        cpl_test_eq(cpl_table_get(edge_lower_table, CR2RES_COL_YS, i, NULL), cmp_ys_lower[i]);
        cpl_test_eq(cpl_table_get(edge_lower_table, CR2RES_COL_CLUSTERS, i, NULL), cmp_cluster_lower[i]);
        cpl_test_eq(cpl_table_get(edge_upper_table, CR2RES_COL_XS, i, NULL), cmp_xs_upper[i]);
        cpl_test_eq(cpl_table_get(edge_upper_table, CR2RES_COL_YS, i, NULL), cmp_ys_upper[i]);
        cpl_test_eq(cpl_table_get(edge_upper_table, CR2RES_COL_CLUSTERS, i, NULL), cmp_cluster_upper[i]);
    }

    //deallocate memory
    cpl_table_unwrap(pixels_table, CR2RES_COL_XS);
    cpl_table_unwrap(pixels_table, CR2RES_COL_YS);
    cpl_table_unwrap(pixels_table, CR2RES_COL_CLUSTERS);
    cpl_table_delete(pixels_table);

    cpl_table_delete(edge_lower_table);
    cpl_table_delete(edge_upper_table);
}

/*----------------------------------------------------------------------------*/
/**
   @brief   Create a table and check if the same index is recovered
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_get_trace_table_index(void)
{
    //define input
    int n = 10;
    int data[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    int data2[] = {1, 1, 1, 1, 1, 1, 2, 1, 1, 1};
    cpl_table *trace_wave = cpl_table_new(n);
    cpl_table_wrap_int(trace_wave, data, CR2RES_COL_ORDER); //what table do we need ?
    cpl_table_wrap_int(trace_wave, data2, CR2RES_COL_TRACENB);

    int order = 5; 
    int trace_nb = 1;
    cpl_size res;

    //run test
    cpl_test_eq(-1, cr2res_get_trace_table_index(NULL, order, trace_nb));
    cpl_test_eq(-1, cr2res_get_trace_table_index(trace_wave, -1, trace_nb));
    cpl_test_eq(-1, cr2res_get_trace_table_index(trace_wave, order, -50));

    cpl_test(res = cr2res_get_trace_table_index(trace_wave, order, trace_nb));
    //test output
    cpl_test_eq(res, 4);

    order = 7;
    // trace would be 2, but we just look for 1
    //run test
    cpl_test(res = cr2res_get_trace_table_index(trace_wave, order, trace_nb));
    //test output
    cpl_test_eq(res, -1);

    order = -10;
    // order does not exist
    //run test
    cpl_test(res = cr2res_get_trace_table_index(trace_wave, order, trace_nb));
    //test output
    cpl_test_eq(res, -1);

    //deallocate memory
    cpl_table_unwrap(trace_wave, CR2RES_COL_ORDER);
    cpl_table_unwrap(trace_wave, CR2RES_COL_TRACENB);
    cpl_table_delete(trace_wave);
}

/*----------------------------------------------------------------------------*/
/**
   @brief   Create a table and check if the same polynomial is recovered
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_get_trace_wave_poly(void)
{
    cpl_table *trace_wave = cpl_table_new(1);
    cpl_table_new_column_array(trace_wave, CR2RES_COL_WAVELENGTH, CPL_TYPE_DOUBLE, 3);
    cpl_table_new_column(trace_wave, CR2RES_COL_ORDER, CPL_TYPE_INT);
    cpl_table_new_column(trace_wave, CR2RES_COL_TRACENB, CPL_TYPE_INT);
    cpl_table_set(trace_wave, CR2RES_COL_ORDER, 0, 1);
    cpl_table_set(trace_wave, CR2RES_COL_TRACENB, 0, 1);
    double pdata[] = {1.1, 2.2, 3.3};
    cpl_array *parr = cpl_array_wrap_double(pdata, 3);
    cpl_table_set_array(trace_wave, CR2RES_COL_WAVELENGTH, 0, parr);

    //run test
    cpl_polynomial *res_poly;
    cpl_test_null(cr2res_get_trace_wave_poly(NULL, CR2RES_COL_WAVELENGTH, 1, 1));
    cpl_test_null(cr2res_get_trace_wave_poly(trace_wave, "blub", 1, 1));
    cpl_test_null(cr2res_get_trace_wave_poly(trace_wave, CR2RES_COL_WAVELENGTH, 20, 1));
    cpl_test_null(cr2res_get_trace_wave_poly(trace_wave, CR2RES_COL_WAVELENGTH, 1, -90));

    cpl_test(res_poly = cr2res_get_trace_wave_poly(trace_wave, CR2RES_COL_WAVELENGTH, 1, 1));
    //test output
    cpl_size power = 0;
    cpl_test_abs(1.1, cpl_polynomial_get_coeff(res_poly, &power), DBL_EPSILON);
    power = 1;
    cpl_test_abs(2.2, cpl_polynomial_get_coeff(res_poly, &power), DBL_EPSILON);
    power = 2;
    cpl_test_abs(3.3, cpl_polynomial_get_coeff(res_poly, &power), DBL_EPSILON);

    cpl_array_unwrap(parr);
    cpl_table_delete(trace_wave);
    cpl_polynomial_delete(res_poly);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Run the Unit tests
 */
/*----------------------------------------------------------------------------*/
int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_DEBUG);

    test_cr2res_trace();
    test_cr2res_trace_clean();
    test_cr2res_trace_gen_image();
    test_cr2res_trace_get_order_idx_values();
    test_cr2res_trace_get_ycen();
    test_cr2res_trace_get_height();
    test_cr2res_trace_compute_middle();
    test_cr2res_trace_compute_height();
    test_cr2res_trace_get_trace_ypos();
    test_cr2res_trace_add_extra_columns();
    test_cr2res_trace_signal_detect();
    test_cr2res_trace_fit_traces();
    test_cr2res_trace_fit_trace();
    test_cr2res_trace_convert_labels_to_cluster();
    test_cr2res_trace_clean_blobs();
    test_cr2res_trace_extract_edges();
    test_cr2res_trace_new_slit_fraction();
    test_cr2res_get_trace_table_index();
    test_cr2res_get_trace_wave_poly();
    return cpl_test_end(0);
}
/**@}*/
