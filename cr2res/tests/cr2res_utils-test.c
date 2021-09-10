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
#include <cr2res_utils.h>
#include <cr2res_dfs.h>
#include <cr2res_trace.h>
#include <cr2res_wave.h>
#include <cr2res_io.h>
#include <cr2res_slit_curv.h>


/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static void test_cr2res_vector_get_int(void);
static void test_cr2res_vector_get_rest(void);
static void test_cr2res_image_cut_rectify(void);
static void test_cr2res_image_insert_rect(void);
static void test_cr2res_polynomial_eval_vector(void);
static void test_cr2res_threshold_spec(void);
static void test_cr2res_get_base_name(void);
static void test_cr2res_get_root_name(void);
static void test_cr2res_extract_filename(void);
static void test_cr2res_extract_frameset(void);
static void test_cr2res_convert_array_to_poly(void);
static void test_cr2res_convert_poly_to_array(void);
static void test_cr2res_detector_shotnoise_model(void);
static void test_cr2res_fit_interorder(void);
static void test_cr2res_slit_pos(void);
static void test_cr2res_slit_pos_img(void);
static void test_cr2res_get_license(void);
static void test_cr2res_slit_curv_compute_order_trace(void);
static void test_cr2res_optimal_filter_2d(void);
static void test_cr2res_polyfit_2d(void);

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_utils-test    Unit test of cr2res_utils
 *
 */
/*----------------------------------------------------------------------------*/

/**@{*/

static void test_cr2res_vector_get_int(void)
{
    int i;
    double d;
    int n = 10;
    cpl_vector *in = cpl_vector_new(n);
    int *res;

    for (i = 0; i < n; i++)
    {
        d = (double)i;
        cpl_vector_set(in, i, d + (d / (n + 1)));
    }

    cpl_test_null(cr2res_vector_get_int(NULL));

    cpl_test(res = cr2res_vector_get_int(in));

    for (i = 0; i < n; i++)
    {
        cpl_test_eq(i, res[i]);
    }

    cpl_vector_delete(in);
    cpl_free(res);

    return;
}
static void test_cr2res_vector_get_rest(void)
{
    int i;
    double d;
    int n = 1000;
    cpl_vector *in = cpl_vector_new(n);
    cpl_vector *out = cpl_vector_new(n);
    double *res;

    for (i = 0; i < n; i++)
    {
        d = (double)i;
        cpl_vector_set(in, i, d + (d / (n + 1)));
        cpl_vector_set(out, i, (d / (n + 1)));
    }

    cpl_test_null(cr2res_vector_get_rest(NULL));

    cpl_test(res = cr2res_vector_get_rest(in));
    cpl_vector_delete(in);
    in = cpl_vector_wrap(n, res);
    cpl_test_vector_abs(in, out, DBL_EPSILON * n);

    cpl_vector_delete(in);
    cpl_vector_delete(out);

    return;
}
static void test_cr2res_image_cut_rectify(void)
{
    cpl_image *res;
    int imdata[] = {1, 2, 3, 2, 1,
                    1, 2, 9, 2, 9,
                    1, 9, 3, 9, 1,
                    9, 2, 3, 2, 1};
    cpl_image *img = cpl_image_wrap_int(5, 4, imdata);
    cpl_image_flip(img, 0); // so that the image looks as formatted above.

    // What result should be
    int cmpdata[] = {9, 9, 9, 9, 9};
    cpl_image *cmp = cpl_image_wrap_int(5, 1, cmpdata);

    double ydata[] = {1.9, 2.1, 3.5, 2.8, 3.99};
    cpl_vector *ycen = cpl_vector_wrap(5, ydata);

    // Run the main function to be tested
    cpl_test_null(cr2res_image_cut_rectify(NULL, ycen, 1));
    cpl_test_null(cr2res_image_cut_rectify(img, NULL, 1));
    cpl_test_null(cr2res_image_cut_rectify(img, ycen, 0));

    cpl_test(res = cr2res_image_cut_rectify(img, ycen, 1));

    // Compare the two
    cpl_test_image_abs(res, cmp, 0);

    cpl_image_unwrap(img);
    cpl_image_unwrap(cmp);
    cpl_vector_unwrap(ycen);
    cpl_image_delete(res);

    return;
}
static void test_cr2res_image_insert_rect(void)
{
    int recdata[] = {1, 2, 3, 2, 1,
                     1, 2, 9, 2, 9,
                     1, 9, 3, 9, 1,
                     9, 2, 3, 2, 1};
    cpl_image *rect_in = cpl_image_wrap_int(5, 4, recdata);
    cpl_image_flip(rect_in, 0);

    double ydata[] = {0.5, 1.1, 6.7, 11.9, 12.1};
    cpl_vector *ycen = cpl_vector_wrap(5, ydata);
    cpl_image *img_out = cpl_image_new(5, 12, CPL_TYPE_INT);
    int cmpdata[] = {0, 0, 0, 2, 9,
                     0, 0, 0, 2, 1,
                     0, 0, 0, 9, 1,
                     0, 0, 0, 2, 0,
                     0, 0, 0, 0, 0,
                     0, 0, 3, 0, 0,
                     0, 0, 9, 0, 0,
                     0, 0, 3, 0, 0,
                     0, 0, 3, 0, 0,
                     0, 0, 0, 0, 0,
                     0, 2, 0, 0, 0,
                     1, 2, 0, 0, 0};
    cpl_image *compare = cpl_image_wrap_int(5, 12, cmpdata);
    cpl_image_flip(compare, 0);

    cpl_test_eq(-1, cr2res_image_insert_rect(NULL, ycen, img_out));
    cpl_test_eq(-1, cr2res_image_insert_rect(rect_in, NULL, img_out));
    cpl_test_eq(-1, cr2res_image_insert_rect(rect_in, ycen, NULL));

    cpl_test_zero(cr2res_image_insert_rect(rect_in, ycen, img_out));

    if (cpl_msg_get_level() == CPL_MSG_DEBUG)
    {
        cpl_image_save(img_out, "TEST_out.fits", CPL_TYPE_INT, NULL,
            CPL_IO_CREATE);
        cpl_image_save(compare, "TEST_cmp.fits", CPL_TYPE_INT, NULL,
            CPL_IO_CREATE);
    }

    // img_out == compare ?
    cpl_test_image_abs(img_out, compare, 0);

    cpl_image_unwrap(rect_in);
    cpl_image_unwrap(compare);
    cpl_vector_unwrap(ycen);
    cpl_image_delete(img_out);

    return;
}
static void test_cr2res_polynomial_eval_vector(void)
{
    int i;
    cpl_size power;
    double p0 = 1.1, p1 = 2.2, p2 = 3.3;
    double d, val;
    int n = 100;
    cpl_vector *in = cpl_vector_new(n);
    cpl_vector *out = cpl_vector_new(n);
    cpl_vector *res;
    cpl_polynomial *poly = cpl_polynomial_new(1);

    power = 0;
    cpl_polynomial_set_coeff(poly, &power, p0);
    power = 1;
    cpl_polynomial_set_coeff(poly, &power, p1);
    power = 2;
    cpl_polynomial_set_coeff(poly, &power, p2);

    for (i = 0; i < n; i++)
    {
        d = (double)i;
        val = d + (d / (n + 1));
        cpl_vector_set(in, i, d);
        val = (p2 * d * d) + (p1 * d) + p0;
        cpl_vector_set(out, i, val);
    }

    cpl_test_null(cr2res_polynomial_eval_vector(NULL, in));
    cpl_test_null(cr2res_polynomial_eval_vector(poly, NULL));

    cpl_test(res = cr2res_polynomial_eval_vector(poly, in));

    cpl_test_vector_abs(res, out, DBL_EPSILON * n * n * 10);

    cpl_vector_delete(in);
    cpl_vector_delete(out);
    cpl_vector_delete(res);
    cpl_polynomial_delete(poly);

    return;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Find the regions with over-average values in a vector
  @param    invector    The vector to be analyzed
  @param    smooth      The size of the boxcar smoothing kernel
  @return   Vector derived as (invector - smoothed_vector - thresh),
            meaning that positive values are at least thresh larger than
            the smoothed vector.
            The returned vector needs to be deallocated by the caller.
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_threshold_spec(void)
{
    //define input
    int n = 10;
    double data[] = {1., 2., 1., 5., 3., 1., 15., 2., 0., 1.};
    cpl_vector *invector = cpl_vector_wrap(n, data);
    // expected data = data - median of boxcar - thresh
    // what is the expected behaviour at the borders?
    // -3, -3, -4, 0, -3, -4, 10, -2, -3, -1
    double outdata[] = {1 - 1 - 3, 2 - 1 - 3, 1 - 2 - 3, 5 - 3 - 3, 3 - 3 - 3,
        1 - 3 - 3, 15 - 2 - 3, 2 - 2 - 3, 0 - 1 - 3, 1 - 1 - 3};
    cpl_vector *outvector = cpl_vector_wrap(n, outdata);

    //boxcar size = smooth + 3, for even values, and smooth + 2 for odd values
    int smooth = 0; //the documentation isn't really right about what smooth is
    double thresh = 3;
    //define output
    cpl_vector *res;

    //run test
    cpl_test_null(cr2res_threshold_spec(NULL, smooth, thresh));
    cpl_test_null(cr2res_threshold_spec(invector, -1, thresh));
    //cpl_test_null(cr2res_threshold_spec(invector, smooth, 4545));

    cpl_test(res = cr2res_threshold_spec(invector, smooth, thresh));

    //cpl_vector_dump(res, "test.log");
    //check output
    cpl_test_vector_abs(outvector, res, DBL_EPSILON * n * n * 10);

    //deallocate memory
    cpl_vector_unwrap(outvector);
    cpl_vector_unwrap(invector);
    cpl_vector_delete(res);

    return;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Find out the base name of a file (i.e. without prefix path)
  @param    filename    Full path name to scan.
  @return   Pointer to char within the input string.
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_get_base_name(void)
{
    //define input
    char *filename = "./../tests/cr2res_trace-test.log";
    char *res;

    //run test
    cpl_test_null(cr2res_get_base_name(NULL));

    cpl_test(res = cr2res_get_base_name(filename));
    //test output
    cpl_test_eq_string(res, "cr2res_trace-test.log");
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Find out the root part of a basename (name without extension).
  @param    filename    File name to scan.
  @return   Pointer to statically allocated string.
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_get_root_name(void)
{
    //define input
    //it only removes the extension for fits, dat, paf, txt, and ascii files
    char *filename = "cr2res_trace-test.fits";
    char *res;

    //run test
    cpl_test_null(cr2res_get_root_name(NULL));

    cpl_test(res = cr2res_get_root_name(filename));
    //test output
    cpl_test_eq_string(res, "cr2res_trace-test");
}

/*----------------------------------------------------------------------------*/
/**
   @brief   Create a frameset two frames, with different tags, and check that 
            only the correct one is recovered
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_extract_filename(void)
{
    //define input
    cpl_frame *frame = cpl_frame_new();
    cpl_frame_set_filename(frame, "bla-test.log");
    cpl_frame_set_tag(frame, "test_correct");

    cpl_frame *other = cpl_frame_new();
    cpl_frame_set_filename(other, "blub-test.log");
    cpl_frame_set_tag(other, "test_wrong");

    cpl_frameset *in = cpl_frameset_new();
    cpl_frameset_insert(in, other);
    cpl_frameset_insert(in, frame);

    char *tag = "test_correct";
    const char *res;

    //run test
    cpl_test_null(cr2res_extract_filename(NULL, tag));
    cpl_test_null(cr2res_extract_filename(in, NULL));

    cpl_test(res = cr2res_extract_filename(in, tag));
    //test output
    cpl_test_eq_string(res, "bla-test.log");

    //deallocate memory
    cpl_frameset_delete(in); //this should also delete the frames
}

/*----------------------------------------------------------------------------*/
/**
   @brief   Extract the frames with the given tag from a frameset
   @param   in      A non-empty frameset
   @param   tag     The tag of the requested frames
   @return  The newly created frameset or NULL on error

   The returned frameset must be de allocated with cpl_frameset_delete
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_extract_frameset(void)
{
    //define input
    cpl_frame *frame = cpl_frame_new();
    cpl_frame_set_filename(frame, "bla-test.log");
    cpl_frame_set_tag(frame, "test_correct");

    cpl_frame *other = cpl_frame_new();
    cpl_frame_set_filename(other, "blub-test.log");
    cpl_frame_set_tag(other, "test_wrong");

    cpl_frameset *in = cpl_frameset_new();
    cpl_frameset_insert(in, frame);
    cpl_frameset_insert(in, other);

    char *tag = "test_correct";
    cpl_frameset *res;

    //run test
    cpl_test_null(cr2res_extract_frameset(NULL, tag));
    cpl_test_null(cr2res_extract_frameset(in, NULL));

    cpl_test(res = cr2res_extract_frameset(in, tag));
    //test output
    //test size
    cpl_test_eq(1, cpl_frameset_get_size(res));
    //check if filenames fit
    const char *fname1 = "bla-test.log";
    const char *fname2 = cpl_frame_get_filename(
        cpl_frameset_get_position(res, 0));
    cpl_test_eq_string(fname1, fname2); //Is that the right comparison?
    //check that the reference was copied as it is supposed to
    cpl_test_noneq_ptr(cpl_frameset_get_position(res, 0), frame);

    //deallocate memory
    //this also deletes the frames
    cpl_frameset_delete(res);
    cpl_frameset_delete(in);
}
/*----------------------------------------------------------------------------*/
/**
   @brief   Check that the coefficients stay the same
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_convert_array_to_poly(void)
{
    //define input
    int n = 10;
    double data[] = {0.9, 1.5, 219.1, 123.8, 18, 123.3, 0.623, 0., 0.9, 1};
    cpl_array *arr = cpl_array_wrap_double(data, n);
    cpl_polynomial *res;
    cpl_size power = 0;
    double poly;

    //run test
    cpl_test_null(cr2res_convert_array_to_poly(NULL));

    cpl_test(res = cr2res_convert_array_to_poly(arr));

    //test output
    for (int i = 0; i < n; i++)
    {
        power = i;
        poly = cpl_polynomial_get_coeff(res, &power);
        cpl_test_eq(data[i], poly);
    }

    //deallocate memory
    cpl_polynomial_delete(res);
    cpl_array_unwrap(arr);
}

/*----------------------------------------------------------------------------*/
/**
   @brief   Check that the coefficients stay the same
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_convert_poly_to_array(void)
{
    //define input
    int n = 10;
    cpl_array *res;
    double data[] = {0.9, 1.5, 219.1, 123.8, 18, 123.3, 0.623, 0., 0.9, 1};
    cpl_polynomial *poly = cpl_polynomial_new(1);
    for (cpl_size i = 0; i < n; i++)
        cpl_polynomial_set_coeff(poly, &i, data[i]);

    //also if size = NULL, no error is raised. Problem?

    //run test
    cpl_test_null(cr2res_convert_poly_to_array(NULL, n));
    cpl_test_null(cr2res_convert_poly_to_array(poly, 0));

    cpl_test(res = cr2res_convert_poly_to_array(poly, n));
    //test output
    for (int j = 0; j < n; j++)
        cpl_test_eq(cpl_array_get(res, j, NULL), data[j]);

    //deallocate memory
    cpl_polynomial_delete(poly);
    cpl_array_delete(res);
}

/*----------------------------------------------------------------------------*/
/**
  @brief   Create a sample image with uniform count and check the result
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_detector_shotnoise_model(void)
{
    //define input
    const double count = 1;
    const double gain = 7;
    const double ron = 3;
    const double err = sqrt(count/gain + ron * ron);
    int width = 5;
    int height = 12;
 
    cpl_image *ima_data = cpl_image_new(width, height, CPL_TYPE_INT);
    cpl_image_add_scalar(ima_data, count);
    // set first pixel to negative to check behaviour
    cpl_image_set(ima_data, 1, 1, -1);

    cpl_image *ima_errs;
    cpl_error_code res;

    cpl_image * compare = cpl_image_new(width, height, CPL_TYPE_INT);
    cpl_image_add_scalar(compare, err);
    cpl_image_set(compare, 1, 1, ron);

    //run test
    cpl_test_eq(CPL_ERROR_NULL_INPUT, 
        cr2res_detector_shotnoise_model(NULL, gain, ron, &ima_errs));
    cpl_test_error(CPL_ERROR_NULL_INPUT);
    cpl_test_eq(CPL_ERROR_ILLEGAL_INPUT, 
        cr2res_detector_shotnoise_model(ima_data, 0, ron, &ima_errs));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_eq(CPL_ERROR_ILLEGAL_INPUT, 
        cr2res_detector_shotnoise_model(ima_data, gain, -1, &ima_errs));
    cpl_test_error(CPL_ERROR_ILLEGAL_INPUT);
    cpl_test_eq(CPL_ERROR_NULL_INPUT, cr2res_detector_shotnoise_model(ima_data,
        gain, ron, NULL));
    cpl_test_error(CPL_ERROR_NULL_INPUT);

    cpl_test_eq(CPL_ERROR_NONE, cr2res_detector_shotnoise_model(ima_data, gain,
        ron, &ima_errs));
    //test output
    cpl_test_image_abs(ima_errs, compare, 0);

    //deallocate memory
    cpl_image_delete(ima_data);
    cpl_image_delete(ima_errs);
    cpl_image_delete(compare);
}

static cpl_table *create_test_table()
{
    int poly_order = 2;
    int n_orders = 2;
    cpl_table * traces = cpl_table_new(n_orders);
    cpl_table_new_column_array(traces, CR2RES_COL_ALL,
        CPL_TYPE_DOUBLE, poly_order);
    cpl_table_new_column_array(traces, CR2RES_COL_UPPER,
        CPL_TYPE_DOUBLE, poly_order);
    cpl_table_new_column_array(traces, CR2RES_COL_LOWER,
        CPL_TYPE_DOUBLE, poly_order);
    cpl_table_new_column(traces, CR2RES_COL_ORDER, CPL_TYPE_INT);
    cpl_table_new_column(traces, CR2RES_COL_TRACENB, CPL_TYPE_INT);

    double all_1[] = {86.6279, 175.5738};
    double all_2[] = {0.01699, 0.07512};
    double upper_1[] = {108.5065, 197.3485};
    double upper_2[] = {0.016601, 0.0184364};
    double lower_1[] = {64.05477, 153.7987};
    double lower_2[] = {0.017355, 0.01659297};

    cpl_array * array = cpl_array_new(poly_order, CPL_TYPE_DOUBLE);
    for (int i = 0; i < n_orders; i++)
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
        cpl_table_set(traces, CR2RES_COL_ORDER, i, 7);
        cpl_table_set(traces, CR2RES_COL_TRACENB, i, i + 1);
    }

    cpl_array_delete(array);
    return traces;
}

static cpl_image *create_test_image()
{
    char *my_path = cpl_sprintf("%s/cr2res_utils_test_image.fits", 
            getenv("srcdir"));
    cpl_image *img = cpl_image_load(my_path, CPL_TYPE_INT, 0, 1);
    cpl_free(my_path) ;
    return img;
}

/*----------------------------------------------------------------------------*/
/**
  @brief Load sample image as input and compare with previous results
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_fit_interorder(void)
{
    // Define all variables
    cpl_image *img = create_test_image();

    // for(int i = 1; i < 100; i++)
    // {
    //     for(int j = 1; j < 100; j++)
    //     {
    //         cpl_image_set(img, i, j, 10000);
    //     }
    // }

    cpl_table *trace_wave = create_test_table();
    cpl_polynomial * res;
    cpl_polynomial *cmp = cpl_polynomial_new(2);
    cpl_size power[] = {0,0};

    // Define comparison polynomial, values from previous run
    // 1.dim.power  2.dim.power  coefficient
    //   0            0      25.1533
    //   1            0      0.245861
    //   0            1      0.519612
    //   1            1      -0.000952474
    cpl_polynomial_set_coeff(cmp, power, 25.1533);
    power[0] = 1;
    cpl_polynomial_set_coeff(cmp, power, 0.245861);
    power[1] = 1;
    cpl_polynomial_set_coeff(cmp, power, -0.000952474);
    power[0] = 0;
    power[1] = 1;
    cpl_polynomial_set_coeff(cmp, power, 0.519612);

    // Run function
    cpl_test(res = cr2res_fit_interorder(img, trace_wave, 1, 1));

    // Compare output
    // thats as precise as it gets, from the polynomial dump
    cpl_test_polynomial_abs(res, cmp, 1e-4);

    // delete cpl objects
    cpl_table_delete(trace_wave);
    cpl_polynomial_delete(res);
    cpl_polynomial_delete(cmp);
    cpl_image_delete(img);
    return;
}

/*----------------------------------------------------------------------------*/
/**
  @brief Load sample data and check if it runs
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_slit_pos()
{
    int chip = 2;
        
    char *my_path = cpl_sprintf("%s/CRIFORS_H24_F_decker1_trace.fits", 
            getenv("srcdir"));
    cpl_table *tw_decker1 = cpl_table_load(my_path, chip, 0);
    cpl_free(my_path) ;

    int nb_orders;
    int *orders = cr2res_trace_get_order_idx_values(tw_decker1, &nb_orders);


    cpl_polynomial **coef_wave = cpl_malloc(nb_orders *
        sizeof(cpl_polynomial *));
    cpl_polynomial **coef_slit = cpl_malloc(nb_orders *
        sizeof(cpl_polynomial *));

    for (cpl_size i = 0; i < nb_orders; i++)
    {
        coef_wave[i] = cpl_polynomial_new(1);
        coef_slit[i] = cpl_polynomial_new(1);
    }
    

    // test NULL input
    //cpl_test_eq(-1, cr2res_slit_pos(NULL, &coef_slit, &coef_wave));
    //cpl_test_eq(-1, cr2res_slit_pos(tw_decker1, NULL, &coef_wave));
    //cpl_test_eq(-1, cr2res_slit_pos(tw_decker1, &coef_slit, NULL));

    // normal run
    //cpl_test_eq(0, cr2res_slit_pos(tw_decker1, &coef_slit, &coef_wave));

    cpl_table_delete(tw_decker1);
    for (int i=0; i < nb_orders; i++){

        // if (i == 3){
        // FILE * file = fopen("slit.txt", "w");
        // char str[3];
        // sprintf(str, "%i\n", orders[i]);
        // fwrite(str, 1, sizeof(str), file);
        // cpl_polynomial_dump(coef_slit[i], file);
        // fclose(file);

        // FILE * file2 = fopen("wave.txt", "w");
        // char str2[3];
        // sprintf(str2, "%i\n", orders[i]);
        // fwrite(str2, 1, sizeof(str2), file2);
        // cpl_polynomial_dump(coef_wave[i], file2);
        // fclose(file2);
        // }

        cpl_polynomial_delete(coef_wave[i]);
        cpl_polynomial_delete(coef_slit[i]);
    }
    cpl_free(orders);
    cpl_free(coef_wave);
    cpl_free(coef_slit);
}

/*----------------------------------------------------------------------------*/
/**
  @brief Load sample data and check if it runs
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_slit_pos_img()
{
    int chip = 2;
    char *my_path = cpl_sprintf("%s/CRIFORS_H24_F_decker1_trace.fits", 
            getenv("srcdir"));
    cpl_table *tw_decker1 = cpl_table_load(my_path, chip, 0);
    cpl_free(my_path) ;
    cpl_image *slitpos = cpl_image_new(CR2RES_DETECTOR_SIZE,
        CR2RES_DETECTOR_SIZE, CPL_TYPE_DOUBLE);
    cpl_image *wavelength = cpl_image_new(CR2RES_DETECTOR_SIZE,
        CR2RES_DETECTOR_SIZE, CPL_TYPE_DOUBLE);

    //cpl_test_eq(-1, cr2res_slit_pos_image(NULL, &slitpos, &wavelength));
    //cpl_test_eq(-1, cr2res_slit_pos_image(tw_decker1, NULL, &wavelength));
    //cpl_test_eq(-1, cr2res_slit_pos_image(tw_decker1, &slitpos, NULL));

    //cpl_test_eq(0, cr2res_slit_pos_image(tw_decker1, &slitpos, &wavelength));

    cpl_image_save(slitpos, "TEST_slit.fits", CPL_TYPE_DOUBLE,
        NULL, CPL_IO_CREATE);
    cpl_image_save(wavelength, "TEST_wave.fits", CPL_TYPE_DOUBLE,
        NULL, CPL_IO_CREATE);

    cpl_table_delete(tw_decker1);
    // cpl_table_delete(tw_decker2);
    cpl_image_delete(slitpos);
    cpl_image_delete(wavelength);
}

static cpl_image * load_etalon_image(){
    char * path = cpl_sprintf("%s/cr2res_slit_curv_test.fits", 
            getenv("srcdir"));
    cpl_image * img = cpl_image_load(path, CPL_TYPE_INT, 0, 1);
    cpl_free(path) ;
    return img;
}

static cpl_table * load_etalon_table(){
    char * path = cpl_sprintf("%s/cr2res_slit_curv_test_tw.fits", 
            getenv("srcdir"));
    cpl_table * trace_wave = cpl_table_load(path, 1, 0);
    cpl_free(path);
    return trace_wave;
}

static void test_cr2res_slit_curv_compute_order_trace(){

    hdrl_image * img_hdrl;
    cpl_image * img_in = load_etalon_image();
    cpl_table * trace_wave = load_etalon_table();

    cpl_polynomial * poly_a = NULL;
    cpl_polynomial * poly_b = NULL;
    cpl_polynomial * poly_c = NULL;

    int order = 1;
    int trace = 1;
    int height = 100; // The height of the order
    int window = 15;  // The spacing between peaks
    int degree = 2;   // That is the default format
    int fit_c = 1;    // Thats what we want most of the time
    int result;

    // cpl_table_save(trace_wave, NULL, NULL, "debug_tw.fits", CPL_IO_CREATE);

    img_hdrl = hdrl_image_create(img_in, NULL);
    
    result = cr2res_slit_curv_compute_order_trace(img_hdrl, trace_wave,
        order, trace, height, window, degree, fit_c,
        &poly_a, &poly_b, &poly_c);
    cpl_test_eq(0, result);

    if (result == 0){
        cpl_polynomial_dump(poly_a, stderr);
        cpl_polynomial_dump(poly_b, stderr);
        cpl_polynomial_dump(poly_c, stderr);
    }

    cpl_image_delete(img_in);
    hdrl_image_delete(img_hdrl);
    cpl_table_delete(trace_wave);

    cpl_polynomial_delete(poly_a);
    cpl_polynomial_delete(poly_b);
    cpl_polynomial_delete(poly_c);
}

static void test_cr2res_optimal_filter_2d(void)
{
    cpl_image * img;
    cpl_image * weight;
    cpl_image * model;
    double lam_x = 1, lam_y = 1;

    int nx = 100;
    int ny = 50;

    img = cpl_image_new(nx, ny, CPL_TYPE_DOUBLE);
    weight = cpl_image_new(nx, ny, CPL_TYPE_DOUBLE);

    for (cpl_size i = 1; i <= nx; i++)
    {
        for (cpl_size j = 1; j <= ny; j++)
        {
            cpl_image_set(img, i, j, 10);
            cpl_image_set(weight, i, j, 1);
        }
    }
    
    cpl_test_nonnull(
        model = cr2res_util_optimal_filter_2d(img, weight, lam_x, lam_y));

    cpl_image_save(model, "TEST_optimal_filter.fits", 
                    CPL_TYPE_DOUBLE, NULL, CPL_IO_CREATE);

    cpl_image_delete(img);
    cpl_image_delete(weight);
    cpl_image_delete(model);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Perform a 2D polynomial fit f(x, y) = z
  @param    x coordinates
  @param    y coordinates
  @param    z values to fit
  @param    degrees to fit with shape (ndegrees, 2) with one entry per 
            xy combination, x degrees in the 1st column, and y degrees in the 2nd,
            Note that you need to include 0, 0 as well for the constant offset
  @return   0 for all good, -1 for error

    Perform a 2D polynomial fit, where the fit degrees are explicitly specified.

 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_polyfit_2d(void)
{
    cpl_vector *x, *y, *z;
    cpl_matrix *degree;
    cpl_polynomial * poly;
    
    cpl_size npoints = 5 * 5;
    cpl_size ngrid = npoints * npoints;
    cpl_size i, j, k;
    cpl_size power[2];

    x = cpl_vector_new(ngrid);
    y = cpl_vector_new(ngrid);
    z = cpl_vector_new(ngrid);

    // f(x, y) = 1 + x + 2 * x y^2
    k = 0;
    for (i = 0; i < npoints ; i++){
        for (j = 0; j  < npoints; j++){
            cpl_vector_set(x, k, i);
            cpl_vector_set(y, k, j);
            cpl_vector_set(z, k, 1 + i + 2 * i * j * j);
            k++;
        }
    }

    // Define the degrees that should be fitted
    degree = cpl_matrix_new(3, 2);
    // Fit constant
    cpl_matrix_set(degree, 0, 0, 0);
    cpl_matrix_set(degree, 0, 1, 0);

    // Fit x linear
    cpl_matrix_set(degree, 1, 0, 1);
    cpl_matrix_set(degree, 1, 1, 0);

    // Fit x y^2
    cpl_matrix_set(degree, 2, 0, 1);
    cpl_matrix_set(degree, 2, 1, 2);

    // Do the fit
    cpl_test(poly = cr2res_polyfit_2d(x, y, z, degree));

    // Check results
    cpl_test_eq(cpl_polynomial_get_dimension(poly), 2);

    power[0] = 0;
    power[1] = 0;
    cpl_test_abs(cpl_polynomial_get_coeff(poly, power), 1, 0.00001);

    power[0] = 1;
    power[1] = 0;
    cpl_test_abs(cpl_polynomial_get_coeff(poly, power), 1, 0.00001);

    power[0] = 1;
    power[1] = 2;
    cpl_test_abs(cpl_polynomial_get_coeff(poly, power), 2, 0.00001);

    // Clean up
    cpl_vector_delete(x);
    cpl_vector_delete(y);
    cpl_vector_delete(z);
    cpl_matrix_delete(degree);
    cpl_polynomial_delete(poly);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the pipeline copyright and license
  @return   The copyright and license string

  The function returns a pointer to the statically allocated license string.
  This string should not be modified using the returned pointer.
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_get_license(void)
{
    const char *license;
    cpl_test(license = cr2res_get_license());
    return;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Run the Unit tests
 */
/*----------------------------------------------------------------------------*/
int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_DEBUG);

    test_cr2res_vector_get_rest();
    test_cr2res_vector_get_int();
    test_cr2res_polynomial_eval_vector();
    test_cr2res_image_cut_rectify();
    test_cr2res_image_insert_rect();
    test_cr2res_threshold_spec();
    test_cr2res_get_base_name();
    test_cr2res_get_root_name();
    test_cr2res_extract_frameset();
    test_cr2res_convert_array_to_poly();
    test_cr2res_convert_poly_to_array();
    test_cr2res_detector_shotnoise_model();
    test_cr2res_get_license();
    test_cr2res_fit_interorder();
    // test_cr2res_slit_pos();
    // test_cr2res_slit_pos_img();
    test_cr2res_slit_curv_compute_order_trace();
    test_cr2res_optimal_filter_2d();
    test_cr2res_polyfit_2d();

    return cpl_test_end(0);
}

/**@}*/

