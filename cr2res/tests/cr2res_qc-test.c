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

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static void test_cr2res_qc_test1(void);
// static void test_cr2res_qc_count_badpix(void);
// static void test_cr2res_qc_read_out_noise(void);
// static void test_cr2res_qc_dark_mean(void);
// static void test_cr2res_qc_dark_median(void);
// static void test_cr2res_qc_dark_stddev(void);
// static void test_cr2res_qc_trace_count_orders(void);
// static void test_cr2res_qc_trace_count_traces(void);
static void test_cr2res_qc_trace_get_ypos(void);
static void test_cr2res_qc_wave_zeropoint(void);

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

/*---------------------------------------------------------------------------*/
/**
* @brief    count number of bad pixels in image
* @param    bpm    image with bad pixels
* @param    type   value to count as bad pixels
* @return   nbp    number of bad pixels
*/
/*---------------------------------------------------------------------------*/
// static void test_cr2res_qc_count_badpix(void)
// {
    // cpl_image *bpm = create_test_image();
    // cpl_image *tmp = cpl_image_duplicate(bpm);
    // int type = 2;
    // int res;
    // int cmp = 354639;

    // cpl_test_eq(-1, cr2res_qc_count_badpix(NULL, type));
    // cpl_test_eq(0, cr2res_qc_count_badpix(bpm, -2));

    // cpl_test(res = cr2res_qc_count_badpix(bpm, type));
    // cpl_test_eq(res, cmp);
    // cpl_test_image_abs(bpm, tmp, 0);

    // cpl_image_delete(bpm);
    // cpl_image_delete(tmp);
// }

/*---------------------------------------------------------------------------*/
/**
* @brief    get the read out noise of image1 and 2
* @param    im1      first image
* @param    im2      second image
* @return   rdnoise  read out noise
*/
/*---------------------------------------------------------------------------*/
// static void test_cr2res_qc_read_out_noise(void)
// {
    // cpl_image *im1 = create_test_image();
    // cpl_image *im2 = create_test_image();
    // cpl_image_fill_window(im2, 100, 100, 200, 200, 2);
    
    // double res;
    // double cmp = 0.0797785;

    // cpl_test_abs(-1, cr2res_qc_read_out_noise(NULL, im2), DBL_EPSILON);
    // cpl_test_abs(-1, cr2res_qc_read_out_noise(im1, NULL), DBL_EPSILON);

    // cpl_test(res = cr2res_qc_read_out_noise(im1, im2));
    // cpl_test_abs(res, cmp, 1e-5);

    // cpl_image_delete(im1);
    // cpl_image_delete(im2);
// }

// /*---------------------------------------------------------------------------*/
// /**
// * @brief    get the mean dark current of the image
// * @param    dark   dark frame image, i.e. with no light exposure
// * @return   mean   mean dark current
// */
// /*---------------------------------------------------------------------------*/
// static void test_cr2res_qc_dark_mean(void)
// {
//     cpl_image *dark = create_test_image();
//     double res;
//     double mean = 3.31055;

//     cpl_test_abs(-1, cr2res_qc_dark_mean(NULL), 0);
//     cpl_test(res = cr2res_qc_dark_mean(dark));
//     cpl_test_abs(res, mean, 1e-5);

//     cpl_image_delete(dark);
// }

// /*---------------------------------------------------------------------------*/
// /**
// * @brief    get the median of the dark current in the image
// * @param    dark   dark frame image, i.e. with no light exposure
// * @return   median median dark current
// */
// /*---------------------------------------------------------------------------*/
// static void test_cr2res_qc_dark_median(void)
// {
//     cpl_image *dark = create_test_image();
//     double res;
//     double median = 3;

//     cpl_test_abs(-1, cr2res_qc_dark_median(NULL), 0);
//     cpl_test(res = cr2res_qc_dark_median(dark));
//     cpl_test_abs(res, median, 1e-5);

//     cpl_image_delete(dark);
// }

// /*---------------------------------------------------------------------------*/
// /**
// * @brief    get the standard deviation of the dark current
// * @param    dark    dark frame image, i.e. with no light exposure
// * @return   stddev  standard deviation of dark current
// */
// /*---------------------------------------------------------------------------*/
// static void test_cr2res_qc_dark_stddev(void)
// {
//     cpl_image *dark = create_test_image();
//     double res;
//     double std = 3.38368;

//     cpl_test_abs(-1, cr2res_qc_dark_stddev(NULL), 0);
//     cpl_test(res = cr2res_qc_dark_stddev(dark));
//     cpl_test_abs(res, std, 1e-5);

//     cpl_image_delete(dark);
// }

// /*---------------------------------------------------------------------------*/
// /**
// * @brief    count the trace orders in the table
// * @param    tracewave    table with traces as polynomials
// * @return   nb_orders    number of individual orders
// */
// /*---------------------------------------------------------------------------*/
// static void test_cr2res_qc_trace_count_orders(void)
// {
//     cpl_table *tracewave = create_test_table();
//     int res;
//     int norders = 9;

//     cpl_test_eq(-1, cr2res_qc_trace_count_orders(NULL));
//     cpl_test(res = cr2res_qc_trace_count_orders(tracewave));
//     cpl_test_eq(res, norders);
//     cpl_table_delete(tracewave);
// }

// /*---------------------------------------------------------------------------*/
// /**
// * @brief    count the number of traces in the table
// * @param    tracewave    table with traces as polynomials
// * @return   ntraces      number of traces
// */
// /*---------------------------------------------------------------------------*/
// static void test_cr2res_qc_trace_count_traces(void)
// {
//     cpl_table *tracewave = create_test_table();
//     int res;
//     int ntraces = 9;

//     cpl_test_eq(-1, cr2res_qc_trace_count_traces(NULL));
//     cpl_test(res = cr2res_qc_trace_count_traces(tracewave));
//     cpl_test_eq(res, ntraces);

//     cpl_table_delete(tracewave);
// }

/*---------------------------------------------------------------------------*/
/**
* @brief    get the central y position of a given trace and order
* @param    tracewave    table with traces as polynomials
* @param    order        order to get values for
* @param    trace        trace of that order
* @return   ycen         y value of central pixel of trace and order
*/
/*---------------------------------------------------------------------------*/
static void test_cr2res_qc_trace_get_ypos(void)
{
    cpl_table *tracewave = create_test_table();
    int order = 4;
    int trace = 1;
    int res;
    int ycen = 679;

    cpl_test_eq(-1, cr2res_qc_trace_get_ypos(NULL, order, trace));
    cpl_test_eq(-1, cr2res_qc_trace_get_ypos(tracewave, -10, trace));
    cpl_test_eq(-1, cr2res_qc_trace_get_ypos(tracewave, order, -1));

    cpl_test(res = cr2res_qc_trace_get_ypos(tracewave, order, trace));
    cpl_test_eq(res, ycen);

    cpl_table_delete(tracewave);
}

/*---------------------------------------------------------------------------*/
/**
* @brief    get the zeropoint (i.e. y(x=0)) for a given order and trace 
* @param    tracewave    table with traces as polynomials
* @param    order        order to get values for
* @param    trace        trace of that order
* @return   y0           y position of the center leftmost pixel of the trace and order
*/
/*---------------------------------------------------------------------------*/
static void test_cr2res_qc_wave_zeropoint(void)
{
    cpl_table *tracewave = create_test_table();
    int order = 4;
    int trace = 1;
    int res;
    int y0 = 660;

    cpl_test_eq(-1, cr2res_qc_trace_get_ypos(NULL, order, trace));
    cpl_test_eq(-1, cr2res_qc_trace_get_ypos(tracewave, -10, trace));
    cpl_test_eq(-1, cr2res_qc_trace_get_ypos(tracewave, order, -1));

    cpl_test(res = cr2res_qc_wave_zeropoint(tracewave, order, trace));
    cpl_test_eq(res, y0);

    cpl_table_delete(tracewave);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Run the Unit tests
 */
/*----------------------------------------------------------------------------*/
int main(void)
{
    cpl_test_init("test@bugreport.se", CPL_MSG_DEBUG);

    // test_cr2res_qc_count_badpix();
    // test_cr2res_qc_read_out_noise();
    // test_cr2res_qc_dark_mean();
    // test_cr2res_qc_dark_median();
    // test_cr2res_qc_dark_stddev();
    // test_cr2res_qc_trace_count_orders();
    // test_cr2res_qc_trace_count_traces();
    // test_cr2res_qc_trace_get_ypos();
    // test_cr2res_qc_wave_zeropoint();

    return cpl_test_end(0);
}

/**@}*/
