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
static cpl_table *create_test_table(void)
{
    cpl_table *traces;
    cpl_array *array;
    int poly_order, norders;

    /* Initialise */
    poly_order = 2;
    norders = 9;

    /* Create a label image with 9 orders */
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
    cpl_array_delete(array);
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
    cpl_table_delete(traces);

    cpl_image_save(trace_ima, "TEST.fits", CPL_TYPE_INT, NULL, CPL_IO_CREATE);
    return trace_ima;
}

/*----------------------------------------------------------------------------*/
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
    cpl_image *trace_ima;
    cpl_table *out;
    // Get test image
    trace_ima = create_test_image();

    // run tests
    /* NULL Input */
    cpl_test_null(cr2res_trace(NULL, 1.0, 1, 6, 500, 0));
    // regular run
    cpl_test(out = cr2res_trace(trace_ima, 1.0, 1, 2, 10, 0));
    // test results?

    // debug output
    cpl_table_save(out, NULL, NULL, "TEST2.fits", CPL_IO_CREATE);

    // free memory
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
    //TODO what about opening ??

    //define input
    cpl_binary data[] = {1, 1, 0, 0,
                         1, 1, 0, 0,
                         0, 0, 1, 0,
                         0, 0, 1, 0};
    cpl_mask *mask = cpl_mask_wrap(4, 4, data);
    int opening = 0;
    int min_cluster = 3;
    cpl_mask *res;
    cpl_binary data_cmp[] = {1, 1, 0, 0,
                             1, 1, 0, 0,
                             0, 0, 0, 0,
                             0, 0, 0, 0};
    cpl_mask *cmp = cpl_mask_wrap(4, 4, data_cmp);
    //run test
    cpl_test_null(cr2res_trace_clean(NULL, opening, min_cluster));
    cpl_test(res = cr2res_trace_clean(mask, opening, min_cluster));
    //test output
    cpl_test_eq_mask(cmp, res);
    //deallocate memory
    cpl_mask_unwrap(mask);
    cpl_mask_unwrap(cmp);
    cpl_mask_delete(res);
}

/*----------------------------------------------------------------------------*/
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
/*----------------------------------------------------------------------------*/
static void test_cr2res_trace_gen_image(void)
{
    //define input
    int nx = 2048;
    int ny = 2048;
    cpl_image *res;
    cpl_image *cmp = cpl_image_load("test_cr2res_trace_gen_image.fits", CPL_TYPE_INT, 0, 0);
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
    cpl_test_image_abs(res, cmp, DBL_EPSILON);

    //deallocate memory
    cpl_image_delete(res);
    cpl_image_delete(cmp);
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
    //run test
    cpl_test_null(cr2res_trace_get_order_numbers(NULL, &nb_orders));
    cpl_test_null(cr2res_trace_get_order_numbers(trace, NULL));
    cpl_test(res = cr2res_trace_get_order_numbers(trace, &nb_orders));
    //test output
    cpl_test_eq(nb_orders, 9);
    for(int i = 0; i < 9; i++){
        cpl_test_eq(res[i], i+1);
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
    cpl_vector *res;
    //run test
    cpl_test_null(cr2res_trace_get_ycen(NULL, order_nb, trace_nb, size));
    cpl_test_null(cr2res_trace_get_ycen(trace, 100, trace_nb, size));
    cpl_test_null(cr2res_trace_get_ycen(trace, order_nb, 5, size));
    cpl_test_null(cr2res_trace_get_ycen(trace, order_nb, trace_nb, -1));

    cpl_test(res = cr2res_trace_get_ycen(trace, order_nb, trace_nb, size));
    //test output
    // is that really the expected behaviour ?
    double data[size];
    for(int i = 0; i < size; i++){
        // values from test table
        data[i] = 437.881 + (i+1) * 0.0172448;
    }
    cpl_vector * cmp = cpl_vector_wrap(2048, data);

    cpl_test_vector_abs(cmp, res, size * DBL_EPSILON); // comparing vector and array probably doesn't work like that

    //deallocate memory
    cpl_table_delete(trace);
    cpl_vector_delete(res);
    cpl_vector_unwrap(cmp);
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
    cpl_polynomial **res;

    //run test
    cpl_test_null(cr2res_trace_wave_get_polynomials(NULL, order_nb, trace_nb));
    cpl_test_null(cr2res_trace_wave_get_polynomials(trace, 100, trace_nb));
    cpl_test_null(cr2res_trace_wave_get_polynomials(trace, order_nb, 5));

    cpl_test(res = cr2res_trace_wave_get_polynomials(trace, order_nb, trace_nb));
    //test output
    // Upper | Lower
    // 524.126, 0.0171958|  350.398, 0.0170009
    cpl_size power = 0;
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
    int cmp = 23; // 10 + 12.2222

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
    cpl_test(res = cr2res_trace_get_trace_ypos(traces, idx));
    //test output

    //deallocate memory
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
    cpl_table *traces;
    char *file_for_wl;
    int det_nr;
    int res;

    //run test
    cpl_test(res = cr2res_trace_add_order_trace_wavelength_columns(traces, file_for_wl, det_nr));
    //test output

    //deallocate memory
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
    cpl_mask *mask;
    cpl_table *trace_table;
    cpl_mask *res;

    //run test
    cpl_test(res = cr2res_trace_split_traces(mask, trace_table));
    //test output

    //deallocate memory
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
    cpl_image *image;
    int trace_sep;
    double smoothfactor;
    double thresh;
    cpl_mask *res;

    //run test
    cpl_test(res = cr2res_trace_signal_detect(image, trace_sep, smoothfactor, thresh));
    //test output

    //deallocate memory
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
    cpl_table *clustertable;
    int degree;
    cpl_table *res;

    //run test
    cpl_test(res = cr2res_trace_fit_traces(clustertable, degree));
    //test output

    //deallocate memory
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
    cpl_table *table;
    int degree;
    cpl_array *res;

    //run test
    cpl_test(res = cr2res_trace_fit_trace(table, degree));
    //test output

    //deallocate memory
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
    cpl_table *cluster;
    int nx;
    int ny;
    cpl_image *res;

    //run test
    cpl_test(res = cr2res_trace_convert_cluster_to_labels(cluster, nx, ny));
    //test output

    //deallocate memory
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
    cpl_image *labels;
    cpl_table *res;

    //run test
    cpl_test(res = cr2res_trace_convert_labels_to_cluster(labels));
    //test output

    //deallocate memory
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
    cpl_mask *mask;
    int min_cluster;
    cpl_mask *res;

    //run test
    cpl_test(res = cr2res_trace_clean_blobs(mask, min_cluster));
    //test output

    //deallocate memory
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
    cpl_table *pixels_table;
    cpl_table **edge_lower_table;
    cpl_table **edge_upper_table;
    int res;

    //run test
    cpl_test(res = cr2res_trace_extract_edges(pixels_table, edge_lower_table, edge_upper_table));
    //test output

    //deallocate memory
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Run the Unit tests
 */
/*----------------------------------------------------------------------------*/
int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    test_cr2res_trace();
    test_cr2res_trace_clean();
    test_cr2res_trace_gen_image();
    test_cr2res_trace_get_order_numbers();
    test_cr2res_trace_get_ycen();
    test_cr2res_trace_get_height();
    test_cr2res_trace_wave_get_polynomials();
    test_cr2res_trace_compute_middle();
    test_cr2res_trace_compute_height();
    test_cr2res_trace_get_trace_ypos();
    test_cr2res_trace_add_order_trace_wavelength_columns();
    test_cr2res_trace_split_traces();
    test_cr2res_trace_signal_detect();
    test_cr2res_trace_fit_traces();
    test_cr2res_trace_fit_trace();
    test_cr2res_trace_convert_cluster_to_labels();
    test_cr2res_trace_convert_labels_to_cluster();
    test_cr2res_trace_clean_blobs();
    test_cr2res_trace_extract_edges();

    return cpl_test_end(0);
}

/**@}*/
