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
#include <cr2res_trace.h>
#include <cr2res_trace.c>

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static void test_cr2res_trace(void);
static void test_cr2res_trace_clean(void);
static void test_cr2res_trace_gen_image(void);
static void test_cr2res_trace_get_order_numbers(void);
static void test_cr2res_trace_get_ycen(void);
static void test_cr2res_trace_get_height(void);
static void test_cr2res_trace_wave_get_polynomials(void);
static void test_cr2res_trace_compute_middle(void);
static void test_cr2res_trace_compute_height(void);
static void test_cr2res_trace_get_trace_ypos(void);
static void test_cr2res_trace_add_order_trace_wavelength_columns(void);
static void test_cr2res_trace_split_traces(void);
static void test_cr2res_trace_signal_detect(void);
static void test_cr2res_trace_fit_traces(void);
static void test_cr2res_trace_fit_trace(void);
static void test_cr2res_trace_convert_cluster_to_labels(void);
static void test_cr2res_trace_convert_labels_to_cluster(void);
static void test_cr2res_trace_clean_blobs(void);
static void test_cr2res_trace_extract_edges(void);

// these don't exist
// static void test_cr2res_trace_labelize(void);
// static void test_cr2res_trace_fit(void);
// static void test_cr2res_trace_compare(void);
// static void test_cr2res_trace_combine(void);
// static void test_cr2res_trace_open_get_polynomials(void);

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_trace-test    Unit test of cr2res_trace
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
    cpl_array *array;
    int poly_order, norders;
    cpl_propertylist *hdr = cpl_propertylist_new();

    /* Initialise */
    poly_order = 2;
    norders = 9;

    /* NULL Input */
    traces = cpl_table_new(norders);
    cpl_table_new_column_array(traces, CR2RES_COL_ALL, CPL_TYPE_DOUBLE, poly_order);
    cpl_table_new_column_array(traces, CR2RES_COL_UPPER, CPL_TYPE_DOUBLE, poly_order);
    cpl_table_new_column_array(traces, CR2RES_COL_LOWER, CPL_TYPE_DOUBLE, poly_order);
    cpl_table_new_column(traces, CR2RES_COL_ORDER, CPL_TYPE_INT);
    cpl_table_new_column(traces, CR2RES_COL_TRACENB, CPL_TYPE_INT);

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
        cpl_table_set(traces, CR2RES_COL_ORDER, i, i + 1);
        cpl_table_set(traces, CR2RES_COL_TRACENB, i, 1);
    }

    cpl_propertylist_append_string(hdr, "EXTNAME", "CHIP1");

    cpl_propertylist_append_double(hdr, "HIERARCH ESO INS WLEN CENY00", 1994.0945859223);
    cpl_propertylist_append_double(hdr, "HIERARCH ESO INS WLEN CENY01", 1723.67027599362);
    cpl_propertylist_append_double(hdr, "HIERARCH ESO INS WLEN CENY02", 1436.61298619847);
    cpl_propertylist_append_double(hdr, "HIERARCH ESO INS WLEN CENY03", 1168.0222016174);
    cpl_propertylist_append_double(hdr, "HIERARCH ESO INS WLEN CENY04", 915.8934665223831);
    cpl_propertylist_append_double(hdr, "HIERARCH ESO INS WLEN CENY05", 678.542785839296);
    cpl_propertylist_append_double(hdr, "HIERARCH ESO INS WLEN CENY06", 454.468576982434);
    cpl_propertylist_append_double(hdr, "HIERARCH ESO INS WLEN CENY07", 242.388497032926);
    cpl_propertylist_append_double(hdr, "HIERARCH ESO INS WLEN CENY08", 63.5899165277783);

    cpl_propertylist_append_double(hdr, "HIERARCH ESO INS WLEN STRT00", 1756.78720770673);
    cpl_propertylist_append_double(hdr, "HIERARCH ESO INS WLEN STRT01", 1703.55123171562);
    cpl_propertylist_append_double(hdr, "HIERARCH ESO INS WLEN STRT02", 1653.44678372399);
    cpl_propertylist_append_double(hdr, "HIERARCH ESO INS WLEN STRT03", 1606.20544704616);
    cpl_propertylist_append_double(hdr, "HIERARCH ESO INS WLEN STRT04", 1561.58862907265);
    cpl_propertylist_append_double(hdr, "HIERARCH ESO INS WLEN STRT05", 1519.38353098961);
    cpl_propertylist_append_double(hdr, "HIERARCH ESO INS WLEN STRT06", 1479.3997538583);
    cpl_propertylist_append_double(hdr, "HIERARCH ESO INS WLEN STRT07", 1441.46642683629);
    cpl_propertylist_append_double(hdr, "HIERARCH ESO INS WLEN STRT08", -1.);

    cpl_propertylist_append_double(hdr, "HIERARCH ESO INS WLEN END00", 1768.81709603003);
    cpl_propertylist_append_double(hdr, "HIERARCH ESO INS WLEN END01", 1715.21657796851);
    cpl_propertylist_append_double(hdr, "HIERARCH ESO INS WLEN END02", 1664.76903155768);
    cpl_propertylist_append_double(hdr, "HIERARCH ESO INS WLEN END03", 1617.2042020846);
    cpl_propertylist_append_double(hdr, "HIERARCH ESO INS WLEN END04", 1572.2818631378);
    cpl_propertylist_append_double(hdr, "HIERARCH ESO INS WLEN END05", 1529.78775872867);
    cpl_propertylist_append_double(hdr, "HIERARCH ESO INS WLEN END06", 1489.53018613055);
    cpl_propertylist_append_double(hdr, "HIERARCH ESO INS WLEN END07", 1451.3371044349);
    cpl_propertylist_append_double(hdr, "HIERARCH ESO INS WLEN END08", -1.);

    cpl_table_save(traces, NULL, hdr, "test_table.fits", CPL_IO_CREATE);

    cpl_array_delete(array);
    cpl_propertylist_delete(hdr);
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
    //extract = cpl_image_extract(trace_ima, 32, 105, 41, 114);

    //cpl_image_save(extract, "extract.fits", CPL_TYPE_INT, NULL, CPL_IO_CREATE);
    cpl_image_save(trace_ima, "TEST.fits", CPL_TYPE_INT, NULL, CPL_IO_CREATE);
    //cpl_image_delete(extract);

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

/* cpl_test_null(out); */
/**
  @brief  Main function for running all parts of the trace algorithm
  @param    ima             input image
  @param    smoothfactor    Used for detection
  @param    opening         Used for cleaning the mask
  @param    degree          Fitted polynomial degree
  @param    min_cluster     A trace must be bigger - discarded otherwise
  @param    split_single_trace_orders   Flag to split traces
  @return The newly allocated trace table or NULL in error case

  A detection is applied to create a mask. This one is labelised.
  The function converts the label image in the proper cluster table in
  trace to call the traces fitting function.
  The cluster table contains the label image information in the form of
  a table. One column per pixel. The columns are xs (pixel x position),
  ys (pixel y position) and cluster (label number).
  The returned table contains 1 line per trace. Each line has 3 polynomials
  (All, Upper and Lower).
    For example with degree 1 :
                 All|               Upper|               Lower|
  24.3593, 0.0161583|  34.6822, 0.0164165|  14.0261, 0.0159084|
  225.479, 0.0167469|  236.604, 0.0168986|  214.342, 0.0166058|
   436.94, 0.0173438|   448.436, 0.017493|   425.423, 0.017203|
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_trace(void)
{
    cpl_image *trace_ima = create_test_image();
    cpl_table *out;

    // run tests
    /* NULL Input */
    //cpl_test_null(cr2res_trace(NULL, 1.0, 1, 6, 500, 0));
    // regular run
    cpl_test(out = cr2res_trace(trace_ima, 1.0, 1, 2, 10, 0));
    // test results?

    cpl_table_save(out, NULL, NULL, "TEST2.fits", CPL_IO_CREATE);
    // debug output
    cpl_table_delete(out);
    cpl_image_delete(trace_ima);
    return;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Clean small blobs
  @param    mask        input mask with small blobs
  @param    opening     Flag to apply opening filtering to the traces
  @param    min_cluster Remove all clusters smaller than this
  @return A newly allocated mask or NULL in error case
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
    int opening = 1;
    int min_cluster = 3;
    cpl_mask *res;
    cpl_binary data_cmp[] = {1, 1, 0, 0, 0, 0,
                             1, 1, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0,
                             0, 1, 1, 1, 1, 0,
                             0, 0, 0, 0, 0, 0,
                             0, 0, 0, 0, 0, 0};
    cpl_mask *cmp = cpl_mask_wrap(6, 6, data_cmp);
    /* test_cr2res_trace_labelize() ; */

    cpl_test_null(cr2res_trace_clean(NULL, opening, min_cluster));
    cpl_test(res = cr2res_trace_clean(mask, opening, min_cluster));
    cpl_mask_save(res, "res.fits", NULL, CPL_IO_CREATE);
    /* test_cr2res_trace_fit() ; */
    cpl_test_eq_mask(cmp, res);
    /* test_cr2res_trace_compare() ; */
    cpl_mask_unwrap(mask);
    cpl_mask_unwrap(cmp);
    cpl_mask_delete(res);
}
/* test_cr2res_trace_combine() ; */
/* test_cr2res_trace_gen_image() ; */
/**
  @brief    Make an image out of the trace solution
  @param trace  The trace table
  @param nx     X size of the produced image
  @param ny     Y size of the produced image
  @return   A newly allocated image or NULL in error case
  The returned INT image is of size nx x ny, is filled with -1. The
  polynomials of the different trace edges are used to fill the traces with
  the value of the order.
 */
/* test_cr2res_trace_get_order_numbers() ; */
static void test_cr2res_trace_gen_image(void)
{
    /* test_cr2res_trace_open_get_polynomials() ; */
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

    /* test_cr2res_trace_compute_middle() ; */
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
  @brief    Count and return the order numbers in a trace table
  @param    trace       trace table
  @param    nb_orders   [output] number of orders
  @return   newly allocated int array

  The int array will need to be freed by the caller. Its size iѕ
  nb_orders. It contains the list of orders found in the trace table.
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_trace_get_order_numbers(void)
{
    //define input
    cpl_table *trace = create_test_table();
    int nb_orders;
    int *res;
    /* test_cr2res_trace_compute_height() ; */
    cpl_test_null(cr2res_trace_get_order_numbers(NULL, &nb_orders));
    cpl_test_null(cr2res_trace_get_order_numbers(trace, NULL));
    cpl_test(res = cr2res_trace_get_order_numbers(trace, &nb_orders));
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
  @brief    Retrieves the middle (All) polynomial from trace table and evaluates
  @param    trace   TRACE table
  @param    order_nb   Wished order
  @param    trace_nb   Wished trace
  @param    size    Output vector size
  @return
  The returned vector contains the poly evaluation reszult on vector from 1 to
  size. It needs to be destryed by the caller.
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
    cpl_test_null(cr2res_trace_get_ycen(trace, 100, trace_nb, size)); // XXX
    cpl_test_null(cr2res_trace_get_ycen(trace, order_nb, 5, size));
    cpl_test_null(cr2res_trace_get_ycen(trace, order_nb, trace_nb, -1));

    // assemble comparison vector
    double data[size];
    for (int i = 0; i < size; i++)
    {
        // values from test table
        data[i] = 437.881 + (i + 1) * 0.0172448;
    }
    cmp = cpl_vector_wrap(2048, data);

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
  @brief    Computes the average height (pix) of an order, from trace polys.
  @param    trace   TRACE table
  @param    order_nb   Wished order
  @param    trace_nb   Wished trace
  @return   height in pixels
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
  @brief    Select the upper and lower polynomials for the given order/trace
  @param trace      TRACE table
  @param order_nb   Wished order
  @param trace_nb   Wished trace
  @return   array of two polynomials or NULL in error case

  The polynomials will need to be destroyed by the caller:
  cpl_polynomial_delete(out[0]) ; -> Upper
  cpl_polynomial_delete(out[1]) ; -> Lower
  cpl_free(out) ;

 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_trace_wave_get_polynomials(void)
{
    //define input
    cpl_table *trace = create_test_table();
    cpl_size order_nb = 3;
    cpl_size trace_nb = 1;
    cpl_size power = 0;
    cpl_polynomial **res;

    //run test
    cpl_test_null(cr2res_trace_wave_get_polynomials(NULL, order_nb, trace_nb));
    cpl_test_null(cr2res_trace_wave_get_polynomials(trace, 100, trace_nb));
    cpl_test_null(cr2res_trace_wave_get_polynomials(trace, order_nb, 5));

    cpl_test(res = cr2res_trace_wave_get_polynomials(trace, order_nb, trace_nb));
    //test output
    // Upper | Lower
    // 524.126, 0.0171958|  350.398, 0.0170009
    power = 0;
    cpl_test_abs(cpl_polynomial_get_coeff(res[0], &power), 524.126, DBL_EPSILON);
    cpl_test_abs(cpl_polynomial_get_coeff(res[1], &power), 350.398, DBL_EPSILON);

    power = 1;
    cpl_test_abs(cpl_polynomial_get_coeff(res[0], &power), 0.0171958, DBL_EPSILON);
    cpl_test_abs(cpl_polynomial_get_coeff(res[1], &power), 0.0170009, DBL_EPSILON);

    //deallocate memory
    cpl_table_delete(trace);
    cpl_polynomial_delete(res[0]);
    cpl_polynomial_delete(res[1]);
    cpl_free(res);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Computes the positions between 2 trace polynomials
  @param    poly1   First trace
  @param    poly2   Second trace
  @param    size    Output vector size
  @return
  The returned vector contains the pixel positions of the middle of the
  2 traces.
  The nth vector value is trace1(n) + trace2(n) / 2
  n=1 for the first value
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
  @brief    Computes extraction height between 2 trace polynomials
  @param trace1         First trace
  @param trace2         Second trace
  @param vector_size    detector x size
  @return   The average height between 2 polynomials or -1 in error case

  The returned int is the rounded-up mean difference between the two
  input polynomials, evaluated on a vector from 1 to vector_size.
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
  @brief    Compute the y position of the trace
  @param traces     The traces table
  @param idx        The index of the trace row
  @return   The y position of the center of the trace
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
static void test_cr2res_trace_add_order_trace_wavelength_columns(void)
{
    //define input
    cpl_table *traces = create_test_table();
    cpl_table_erase_column(traces, CR2RES_COL_ORDER);
    cpl_table *tmp = cpl_table_duplicate(traces);
    char *file_for_wl = "test_table.fits";
    int det_nr = 1;
    int res;
    const cpl_array *wl;
    double cmp1[] = {0, 1441.46160481499, 1479.3948049417, 1519.37844831851, 1561.58340521624, 1606.20007393671, 1653.44125258191, 1703.54553296318, 1756.78133086827};
    double cmp2[] = {0, 0.00482202129878357, 0.00494891659611621, 0.00508267109871028, 0.00522385640701021, 0.00537310944721049, 0.00553114207801171, 0.00569875244401075, 0.00587683845788956};

    //run test
    cpl_test_eq(-1, cr2res_trace_add_order_trace_wavelength_columns(NULL, file_for_wl, det_nr));

    cpl_test_eq(-1, cr2res_trace_add_order_trace_wavelength_columns(tmp, "invalid_path", det_nr));
    cpl_table_delete(tmp);
    tmp = cpl_table_duplicate(traces);

    cpl_test_eq(-1, cr2res_trace_add_order_trace_wavelength_columns(tmp, file_for_wl, 10));
    cpl_table_delete(tmp);
    tmp = cpl_table_duplicate(traces);

    cpl_test_eq(-1, cr2res_trace_add_order_trace_wavelength_columns(tmp, file_for_wl, -1));
    cpl_table_delete(tmp);
    tmp = cpl_table_duplicate(traces);

    cpl_test_eq(0, cr2res_trace_add_order_trace_wavelength_columns(traces, file_for_wl, det_nr));
    //test output
    //cpl_table_save(traces, NULL, NULL, "new_table.fits", CPL_IO_CREATE);

    wl = cpl_table_get_array(traces, CR2RES_COL_WAVELENGTH, 0);
    cpl_test_abs(cpl_array_get(wl, 0, 0), 0, FLT_EPSILON);
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_abs(cpl_array_get(wl, 1, 0), 0, FLT_EPSILON);
    cpl_test_error(CPL_ERROR_NULL_INPUT);

    for (int i = 1; i < 9; i++)
    {
        wl = cpl_table_get_array(traces, CR2RES_COL_WAVELENGTH, i);
        cpl_test_abs(cpl_array_get(wl, 0, 0), cmp1[i], FLT_EPSILON);
        cpl_test_abs(cpl_array_get(wl, 1, 0), cmp2[i], FLT_EPSILON);
    }

    //deallocate memory
    cpl_table_delete(traces);
    cpl_table_delete(tmp);
}

/*----------------------------------------------------------------------------*/
/**
  @brief
  @param
  @return
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_trace_split_traces(void)
{
    //define input
    cpl_binary data[2048 * 2048];
    for (int i = 0; i < 2048 * 2048; i++)
        data[i] = 1;

    cpl_mask *mask = cpl_mask_wrap(2048, 2048, data);
    cpl_table *trace_table = create_test_table();
    cpl_mask *sub;

    cpl_binary data2[10 * 10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    cpl_mask *cmp = cpl_mask_wrap(10, 10, data2);
    cpl_mask *res;

    //run test
    cpl_test_null(cr2res_trace_split_traces(NULL, trace_table));
    cpl_test_null(cr2res_trace_split_traces(mask, NULL));

    cpl_test(res = cr2res_trace_split_traces(mask, trace_table));
    //test output
    sub = cpl_mask_extract(res, 36, 165, 45, 174);
    //cpl_mask_save(res, "res.fits", NULL, CPL_IO_CREATE);
    //cpl_mask_save(sub, "sub.fits", NULL, CPL_IO_CREATE);
    //cpl_mask_save(cmp, "cmp.fits", NULL, CPL_IO_CREATE);
    cpl_test_eq_mask(sub, cmp);

    //deallocate memory
    cpl_mask_unwrap(mask);
    cpl_mask_unwrap(cmp);
    cpl_mask_delete(sub);
    cpl_mask_delete(res);
    cpl_table_delete(trace_table);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Detect the Traces signal
  @param image          The input image with the traces
  @param trace_sep       The approximate number of pixels between 2 traces
  @param smoothfactor   Mult. factor for the low pass filter kernel size
  @param thresh         The threshold used for detection
  @return   A newly allocated mask or NULL in error case.

  The returned mask identifies the pixels belonging to a trace
  The input image is smoothed, subtracted to the result, and a simple
  thresholding is applied. The smoothing kernel size is
  trace_sep*smoothfactor x 1
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_trace_signal_detect(void)
{
    //define input
    cpl_image *image = create_test_image();
    int trace_sep = 150;
    double smoothfactor = 1;
    double thresh = 0;
    cpl_mask *sub;

    cpl_binary data2[10 * 10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    cpl_mask *cmp = cpl_mask_wrap(10, 10, data2);

    cpl_mask *res;

    //run test
    cpl_test_null(cr2res_trace_signal_detect(NULL, trace_sep, smoothfactor, thresh));
    cpl_test_null(cr2res_trace_signal_detect(image, -10, smoothfactor, thresh));
    cpl_test_null(cr2res_trace_signal_detect(image, trace_sep, -1, thresh));
    //cpl_test_null(cr2res_trace_signal_detect(image, trace_sep, smoothfactor, 50000));

    cpl_test(res = cr2res_trace_signal_detect(image, trace_sep, smoothfactor, thresh));
    //test output
    sub = cpl_mask_extract(res, 32, 105, 41, 114);
    cpl_test_eq_mask(sub, cmp);

    //deallocate memory
    cpl_image_delete(image);
    cpl_mask_delete(res);
    cpl_mask_unwrap(cmp);
    cpl_mask_delete(sub);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Fit a polynomial on the different traces (center and edges)
  @param clustertable   The table holding the traces pixels with their labels
  @param degree         Fitting polynomial degree
  @return   A newly allocated table or NULL in error case

  The function loops over the traces labels. For each of them, it
  identifies the upper and lower edges and fits a polynomial to them. It
  also fits a polynomial using all the pixels of the trace.

  The returned table contains 1 line per trace. Each line has 3 polynomials
  (All, Upper and Lower).
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
    //cpl_test_null(cr2res_trace_fit_traces(clustertable, -5));

    cpl_test(res = cr2res_trace_fit_traces(cluster, degree));
    //test output
    //cpl_table_save(res, NULL, NULL, "fit_traces.fits", CPL_IO_CREATE);
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
  @brief    Simple fitting
  @param table  Table containing all the pixels to fit
  @param degree Fitting polynomial degree
  @return   A newly allocated array or NULL in error case

  The pixels in the input table (columns are xs, ys, clusters) are all
  used for the fitting of a polynomial of degree degree. The clusters
  column is IGNORED.
  If the x range of pixels does not exceed 1500 pixels, a linear fit is
  applied.
  The polynomial coefficients are returned in an array.
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_trace_fit_trace(void)
{
    //define input
    // use only cluster 3 from test image
    cpl_image *test_image = create_test_image();
    cpl_table *all = cr2res_trace_convert_labels_to_cluster(test_image);
    cpl_table_and_selected_int(all, CR2RES_COL_CLUSTERS, CPL_EQUAL_TO, 3);
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
    cpl_test_abs(cpl_array_get(res, 0, NULL), 437.881, 1.2);
    cpl_test_abs(cpl_array_get(res, 1, NULL), 0.0172448, 0.0002);

    //deallocate memory
    cpl_image_delete(test_image);
    cpl_table_delete(all);
    cpl_table_delete(table);
    cpl_array_delete(res);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Convert cluster table to labels image
  @param cluster    A cluster table (xs, ys, cluster)
  @param nx         X size of the returned label image
  @param ny         Y size of the returned label image
  @return   A newly allocated INT image or NULL in error case

  A new label image is created from the cluster table. Each entry in the
  cluster table is used to set a pixel in the label image.
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_trace_convert_cluster_to_labels(void)
{
    //define input
    //cpl_table *cluster = create_cluster_table();

    int xs[] = {4, 5, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 1, 2};
    int ys[] = {4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 5, 5};
    int clusters[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2};

    cpl_table *cluster = cpl_table_new(17);
    cpl_table_wrap_int(cluster, xs, CR2RES_COL_XS);
    cpl_table_wrap_int(cluster, ys, CR2RES_COL_YS);
    cpl_table_wrap_int(cluster, clusters, CR2RES_COL_CLUSTERS);

    // flip data, so the image will be the same
    int data_inverse[] = {1, 1, 1, 1, 0,
                          1, 1, 1, 1, 1,
                          0, 1, 1, 1, 1,
                          0, 0, 0, 1, 1,
                          2, 2, 0, 0, 0};

    int nx = 5;
    int ny = 5;
    cpl_image *res;
    cpl_image *cmp = cpl_image_wrap(nx, ny, CPL_TYPE_INT, data_inverse);

    //run test
    cpl_test_null(cr2res_trace_convert_cluster_to_labels(NULL, nx, ny));
    cpl_test_null(cr2res_trace_convert_cluster_to_labels(cluster, 0, ny));
    cpl_test_null(cr2res_trace_convert_cluster_to_labels(cluster, nx, 0));

    cpl_test(res = cr2res_trace_convert_cluster_to_labels(cluster, nx, ny));
    //test output
    //cpl_image_save(res, "labels.fits", CPL_TYPE_INT, NULL, CPL_IO_CREATE);
    cpl_test_image_abs(res, cmp, 0);

    //deallocate memory
    cpl_table_unwrap(cluster, CR2RES_COL_XS);
    cpl_table_unwrap(cluster, CR2RES_COL_YS);
    cpl_table_unwrap(cluster, CR2RES_COL_CLUSTERS);
    cpl_table_delete(cluster);
    cpl_image_delete(res);
    cpl_image_unwrap(cmp);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Convert Labels image to the cluster table
  @param labels     The label image
  @return   A newly allocated cluster table or NULL in error case

  The cluster table contains the label image information in the form of
  a table. One column per pixel. The columns are xs (pixel x position),
  ys (pixel y position) and cluster (label number).
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
    // cpl_table_save(res, NULL, NULL, "convert.fits", CPL_IO_CREATE);

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
  @brief    Cleans small size group of pixels from a mask
  @param mask           Input mask
  @param min_cluster    Size of clusters under which they need to be removed
  @return   A newly allocated mask or NULL in error case
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
    // if min_cluster <= 0 nothing changes
    cpl_test(res = cr2res_trace_clean_blobs(mask, -1));
    cpl_test_eq_mask(res, mask);
    cpl_mask_delete(res);

    cpl_test(res = cr2res_trace_clean_blobs(mask, min_cluster));
    //test output
    cpl_test_eq_mask(res, cmp);

    //deallocate memory
    cpl_mask_unwrap(mask);
    cpl_mask_unwrap(cmp);
    cpl_mask_delete(res);
}

/*----------------------------------------------------------------------------*/
/**
  @brief  Extracts pixels on the upper and lower edges of a trace
  @param pixels_table       Input cluster table with a single trace
  @param edge_lower_table   [output] Lower edge pixels table
  @param edge_upper_table   [output] Upper edge pixels table
  @return   0 if ok, -1 otherwise

  For each found x in the input cluster table, the min and max y are
  used for adding an entry in the 2 output cluster tables.
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
  @brief    Run the Unit tests
 */
/*----------------------------------------------------------------------------*/
int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_DEBUG);

    //test_cr2res_trace();
    /* test_cr2res_trace_clean(); */
    /* test_cr2res_trace_gen_image(); */
    /* test_cr2res_trace_get_order_numbers(); */
    /* test_cr2res_trace_get_ycen(); */
    /* test_cr2res_trace_get_height(); */
    /* test_cr2res_trace_wave_get_polynomials(); */
    /* test_cr2res_trace_compute_middle(); */
    /* test_cr2res_trace_compute_height(); */
    /* test_cr2res_trace_get_trace_ypos(); */
    /* test_cr2res_trace_add_order_trace_wavelength_columns(); */
    /* test_cr2res_trace_split_traces(); */
    /* test_cr2res_trace_signal_detect(); */
    /* test_cr2res_trace_fit_traces(); */
    /* test_cr2res_trace_fit_trace(); */
    /* test_cr2res_trace_convert_cluster_to_labels(); */
    /* test_cr2res_trace_convert_labels_to_cluster(); */
    /* test_cr2res_trace_clean_blobs(); */
    /* test_cr2res_trace_extract_edges(); */

    return cpl_test_end(0);
}
/**@}*/
