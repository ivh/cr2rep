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
#include <cr2res_extract.h>
#include <cr2res_trace.h>
#include <cr2res_slit_curv.h>

static void test_cr2res_slitdec_vert_regression(void);
static void test_cr2res_slitdec_vert_edge_cases(void);
static void test_cr2res_slitdec_vert(void);
static void test_cr2res_slitdec_compare_vert_curved(void);

static cpl_table *load_table()
{
    cpl_table * table;
    if (access("./tests/lamp_trace.fits", F_OK) != -1){
        table = cpl_table_load("./tests/lamp_trace.fits", 1, 0);
    } else{
        table = cpl_table_load("./lamp_trace.fits", 1, 0);
    }
    return table;
}

static cpl_image *load_image()
{
    cpl_image * image;
    if (access("./tests/XDGH.fits", F_OK) != -1){
        image = cpl_image_load("./tests/XDGH.fits", CPL_TYPE_INT, 0, 1);
    }
    else{
        image = cpl_image_load("./XDGH.fits", CPL_TYPE_INT, 0, 1);
    }
    return image;
}

static cpl_table *create_test_table()
{
    int poly_order = 2;
    int n_orders = 2;
    cpl_table * traces = cpl_table_new(n_orders);
    cpl_table_new_column_array(traces, CR2RES_COL_ALL, CPL_TYPE_DOUBLE, poly_order);
    cpl_table_new_column_array(traces, CR2RES_COL_UPPER, CPL_TYPE_DOUBLE, poly_order);
    cpl_table_new_column_array(traces, CR2RES_COL_LOWER, CPL_TYPE_DOUBLE, poly_order);
    cpl_table_new_column(traces, CR2RES_COL_ORDER, CPL_TYPE_INT);
    cpl_table_new_column(traces, CR2RES_COL_TRACENB, CPL_TYPE_INT);

    double all_1[] = {86.6279, 175.5738};
    double all_2[] = {0.01699, 0.07512};
    double upper_1[] = {108.5065, 197.3485};
    double upper_2[] = {0.016601, 0.0184364};
    double lower_1[] = {64.05477, 153.7987};
    double lower_2[] = {0.017355, 0.01659297};

    cpl_array * array = cpl_array_new(poly_order, CPL_TYPE_DOUBLE);
    for (int i = 0; i < n_orders; i++) {
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
    cpl_image *img;

    // The debugger runs in a different directory than make check
    // to get the right file in either case check if it exists or not
    if (access("./tests/cr2res_utils_test_image.fits", F_OK) != -1){
        img = cpl_image_load("./tests/cr2res_utils_test_image.fits", CPL_TYPE_INT, 0, 1);
    }
    else{
        img = cpl_image_load("./cr2res_utils_test_image.fits", CPL_TYPE_INT, 0, 1);
    }
    return img;
}

static cpl_image *create_image_linear_increase(int width, int height, double spec_in[width])
{
    int k;
    double min = 5;
    double step = 0.1;
    double spec, slitf;
    double sigma = height/16.;
    double mu = height/2.;
    cpl_image * img = cpl_image_new(width, height, CPL_TYPE_DOUBLE);

    for(int i = 1; i <= width; i++) {
        spec = i * step + min;       
        if (i > width/2 - 5 & i < width/2 + 5)
            spec -= 2 * min * exp(- (i - width/2) * (i - width/2) / (2. * 1 * 1));

        spec_in[i-1] = spec;
        for(int j = 1; j <= height; j++) {
            slitf =  exp(- (j - mu) * (j - mu) / (2. * sigma * sigma));
            cpl_image_set(img, i, j, spec * slitf);
        }
    }

    return img;
}

static cpl_image * create_image_sinusoidal(
        int         width,
        int         height,
        double      spec_in[width])
{
    int k = 0;
    double min = 5.0;
    double sigma = height/8.;
    double mu = height/2.;
    double spec = 0.;
    double slitf = 0.;

    cpl_image * img = cpl_image_new(width, height, CPL_TYPE_DOUBLE);

    for(int i = 1; i <= width; i++) {
        // Sinusoidal spectrum
        spec = 3 * sin(CPL_MATH_2PI * (i-1.)/(0.2 * width)) + min;
        spec_in[i-1] = spec;

        for(int j = 1; j <= height; j++) {
            // Normal distribution as slitfunction
            slitf = exp(- (j - mu) * (j - mu) / (2. * sigma * sigma));
            cpl_image_set(img, i, j, spec * slitf);
        }
    }
    return img;
}

cpl_image * apply_shear(cpl_image * img, int width, int height, double shear){
    int k;
    double mu = height * 0.5;
    double spec;

    cpl_image * img_tilt = cpl_image_new(width, height, CPL_TYPE_DOUBLE);

    for(int i = 1; i <= width; i++) {
        for(int j = 1; j <= height; j++) {
            k = (int)(i - (j - mu) * shear);
            // periodic boundary conditions
            if (k < 1)     k = k + width;
            if (k > width) k = k - width;

            spec = cpl_image_get(img, k, j, &k);
            cpl_image_set(img_tilt, i, j, spec);
        }
    }
    cpl_image_delete(img);

    return img_tilt;
}

static cpl_table * create_table_linear_increase(
        int         width,
        int         height,
        double      shear)
{
    int poly_order = 1;
    cpl_table * traces = cpl_table_new(2);
    cpl_array * array = cpl_array_new(poly_order, CPL_TYPE_DOUBLE);
    cpl_array * wl = cpl_array_new(2, CPL_TYPE_DOUBLE);
    cpl_array_set(wl, 0, 10);
    cpl_array_set(wl, 1, 1);

    cpl_table_new_column_array(traces, CR2RES_COL_ALL, CPL_TYPE_DOUBLE, poly_order);
    cpl_table_new_column_array(traces, CR2RES_COL_UPPER, CPL_TYPE_DOUBLE, poly_order);
    cpl_table_new_column_array(traces, CR2RES_COL_LOWER, CPL_TYPE_DOUBLE, poly_order);
    cpl_table_new_column_array(traces, CR2RES_COL_WAVELENGTH, CPL_TYPE_DOUBLE, 2);
    cpl_table_new_column(traces, CR2RES_COL_ORDER, CPL_TYPE_INT);
    cpl_table_new_column(traces, CR2RES_COL_TRACENB, CPL_TYPE_INT);
    cpl_table_new_column_array(traces, CR2RES_COL_SLIT_CURV_A, CPL_TYPE_DOUBLE, 2);
    cpl_table_new_column_array(traces, CR2RES_COL_SLIT_CURV_B, CPL_TYPE_DOUBLE, 2);
    cpl_table_new_column_array(traces, CR2RES_COL_SLIT_CURV_C, CPL_TYPE_DOUBLE, 2);



    cpl_array_set(array, 0, height * 0.5);
    cpl_table_set_array(traces, CR2RES_COL_ALL, 0, array);
    cpl_array_set(array, 0, height * 0.8);
    cpl_table_set_array(traces, CR2RES_COL_UPPER, 0, array);
    cpl_array_set(array, 0, height * 0.2);
    cpl_table_set_array(traces, CR2RES_COL_LOWER, 0, array);
    cpl_table_set_array(traces, CR2RES_COL_WAVELENGTH, 0, wl);
    cpl_table_set(traces, CR2RES_COL_ORDER, 0, 1);
    cpl_table_set(traces, CR2RES_COL_TRACENB, 0, 1);

    cpl_array_set(array, 0, height * 0.5);
    cpl_table_set_array(traces, CR2RES_COL_ALL, 1, array);
    cpl_array_set(array, 0, height * 0.8);
    cpl_table_set_array(traces, CR2RES_COL_UPPER, 1, array);
    cpl_array_set(array, 0, height * 0.2);
    cpl_table_set_array(traces, CR2RES_COL_LOWER, 1, array);
    cpl_table_set_array(traces, CR2RES_COL_WAVELENGTH, 1, wl);
    cpl_table_set(traces, CR2RES_COL_ORDER, 1, 1);
    cpl_table_set(traces, CR2RES_COL_TRACENB, 1, 2);

    cpl_array_set(wl, 0, 0);
    cpl_array_set(wl, 1, 1);
    cpl_table_set_array(traces, CR2RES_COL_SLIT_CURV_A, 0, wl);
    cpl_table_set_array(traces, CR2RES_COL_SLIT_CURV_A, 1, wl);

    cpl_array_set(wl, 0, shear);
    cpl_array_set(wl, 1, 0);
    cpl_table_set_array(traces, CR2RES_COL_SLIT_CURV_B, 0, wl);
    cpl_table_set_array(traces, CR2RES_COL_SLIT_CURV_B, 1, wl);

    cpl_array_set(wl, 0, 0);
    cpl_array_set(wl, 1, 0);
    cpl_table_set_array(traces, CR2RES_COL_SLIT_CURV_C, 0, wl);
    cpl_table_set_array(traces, CR2RES_COL_SLIT_CURV_C, 1, wl);


    cpl_array_delete(array);
    cpl_array_delete(wl);
    return traces;
}

static void test_cr2res_slitdec_vert_regression(void)
{
    cpl_image * img_in = create_test_image();
    hdrl_image * img_hdrl = hdrl_image_create(img_in, NULL);
    cpl_table * trace_table = create_test_table();
    int order = 7;
    int trace = 1;
    int height = 70;
    int swath = 50;
    int oversample = 2;
    double smooth_slit = 0.1;
    int width = cpl_image_get_size_x(img_in);

    cpl_vector * slit_func;
    cpl_bivector * spec;
    hdrl_image * model;
    int res;

    // Load Expected results
    static double data_slitfunc[143] = {0.000115859, 7.43396e-05, 7.10328e-05, 7.62753e-05, 7.10487e-05, 8.11558e-05, 7.85971e-05, 6.83344e-05, 7.17496e-05, 7.96937e-05, 7.35047e-05, 5.22185e-05, 7.33561e-05, 7.43463e-05, 0.000144974, 0.000102243, 0.000359142, 0.000521881, 0.00104286, 0.00146012, 0.00226715, 0.00301241, 0.00413671, 0.00509677, 0.00643321, 0.00773337, 0.00897916, 0.0104182, 0.0117634, 0.0131721, 0.0144836, 0.0158845, 0.0169486, 0.0181268, 0.0189894, 0.020072, 0.0205192, 0.0213352, 0.021635, 0.0220201, 0.0221268, 0.0223179, 0.0223631, 0.0223449, 0.0224997, 0.0223255, 0.0224773, 0.0224355, 0.0225178, 0.0224838, 0.0224948, 0.0226417, 0.0224643, 0.0228002, 0.0226002, 0.0228257, 0.0225603, 0.0226531, 0.0226004, 0.0224983, 0.0225823, 0.0224062, 0.0225994, 0.0223767, 0.0226528, 0.0224482, 0.0227148, 0.0224825, 0.0226146, 0.0225114, 0.022532, 0.022325, 0.022366, 0.0221547, 0.0222041, 0.0220273, 0.0221068, 0.0220118, 0.022042, 0.0220703, 0.0221329, 0.0221981, 0.0222083, 0.0223475, 0.0221351, 0.0224142, 0.0221839, 0.0224193, 0.0221855, 0.0223566, 0.0221801, 0.0223412, 0.022238, 0.0222933, 0.0222503, 0.0223179, 0.0223925, 0.0223932, 0.0224435, 0.0223841, 0.0225682, 0.0223857, 0.0224443, 0.0224466, 0.0223342, 0.0222627, 0.0218657, 0.0217104, 0.021, 0.0205833, 0.0196169, 0.0188126, 0.017692, 0.0164563, 0.0152199, 0.0138034, 0.012523, 0.0108846, 0.00956393, 0.00814578, 0.00676929, 0.00544948, 0.00429528, 0.00323442, 0.0024136, 0.00157264, 0.00111565, 0.000598061, 0.000371827, 0.000139823, 0.000121443, 5.84376e-05, 8.41633e-05, 7.4026e-05, 5.50628e-05, 7.25114e-05, 6.36848e-05, 6.18246e-05, 7.54592e-05, 7.05881e-05, 6.02867e-05, 6.36684e-05, 6.36684e-05 };
    cpl_vector *cmp_slitfunc = cpl_vector_wrap(143, data_slitfunc);

    static double data_spec[] = {52427.1, 56057, 61339.8, 65980.5, 68602.9, 68408.7, 67261.7, 65948, 64419.7, 63188.3, 61880.6, 60079.3, 58767.2, 56827.8, 55461.9, 54276.3, 53239.5, 52180.1, 51012.6, 49862.9, 49346.3, 49106.3, 49512, 50059.1, 50739.6, 51833, 52867.8, 56163, 58267.9, 59918.2, 61755.8, 63357.5, 64534.8, 66216.8, 67155.8, 68099.7, 69762.6, 71421.4, 72797.3, 74520.5, 75136.2, 75822.4, 75580.1, 74968.4, 73947.4, 72564.9, 71471.5, 70254.3, 69114.7, 67435.2, 66039.3, 64097.1, 62459.1, 61809.6, 60585.5, 59886.6, 56380, 54752.9, 53538.5, 52106.8, 51201.9, 50227.7, 50485.6, 50145.1, 50502.5, 51539.1, 52777.9, 54718.5, 56431.5, 58163.9, 59658.8, 60643.3, 61944.8, 63230.7, 64646.9, 66172.3, 67801.6, 68887.1, 70113.4, 71368.6, 73070.2, 74478.4, 77749.9, 77479.8, 76418.1, 74976.4, 73924.9, 73339.2, 72656.8, 71766, 70424.5, 69250, 68089.5, 67069, 65597.9, 64096.4, 62310.3, 60617.4, 59116.6, 57845.2, 56924, 56161.2, 55271.3, 54702.6, 54136.6, 54096.7, 54553.9, 55668.1, 57046.9, 57986, 59175.6, 58235.6, 59366.1, 61054.1, 62524.4, 64416.5, 66584, 68586.2, 71234.8, 73330.2, 74696.1, 75223.8, 75189.4, 75280.4, 75444.3, 75874.5, 75904, 75838.2, 74891, 74023.9, 72514.4, 70963.7, 69558.8, 67629.4, 66282, 64630.2, 62758.7, 63849.1, 62848.7, 61904.3, 60854.6, 59230, 57784.9, 56456.6, 55437.3, 55184.8, 56000.7, 57246.1, 58924.6, 60388.7, 61565.1, 62219.3, 62859.2, 63542.2, 64846.9, 66776.1, 68786.9, 70464.1, 72447.6, 73847.4, 75108.2, 76561.9, 77626.4, 78978.3, 82490.1, 83193.5, 81041.9, 81479.2, 81559.7, 80826.1, 79666.5, 77906, 76351.9, 75025.8, 73588.5, 72518.2, 71188.3, 69730.8, 68496.9, 66752.2, 64862.7, 63344.7, 61420, 60027.9, 59126.9, 58096, 57417.4, 57173, 57299.7, 58309.4, 59106.8, 60063.1, 61092.4, 60288.4, 61546.3, 63153.4, 65249.9, 67532.7, 70186.9, 72424.3, 74198.5, 75317.4, 75786.2, 75659.8, 76542, 77153.9, 78396.8, 79690.9, 80568.9, 80859.1, 80186.7, 79281.6, 78150.2, 76752.7, 75368, 73920.9, 71622.6, 70165, 68956.5, 67348.3 };
    cpl_vector * cmp_spec = cpl_vector_wrap(width, data_spec);

    static double data_err[] = {523.439, 499.701, 473.877, 445.496, 387.684, 342.322, 329.536, 333.608, 322.94, 304.494, 343.882, 361.009, 302.392, 307.149, 332.682, 323.992, 362.454, 340.802, 359.662, 340.305, 357.704, 383.636, 364.937, 362.076, 365.54, 352.365, 404.616, 386.916, 364.15, 340.36, 326.278, 354.856, 428.537, 436.106, 429.022, 449.069, 413.98, 416.991, 358.546, 288.909, 305.632, 328.581, 416.744, 466.815, 541.072, 511.904, 405.265, 361.332, 385.985, 403.975, 347.848, 490.804, 469.661, 414.215, 378.947, 372.721, 377.953, 395.166, 419.834, 417.746, 435.44, 433.026, 405.965, 355.437, 374.395, 397.416, 363.972, 388.46, 337.569, 314.348, 269.74, 251.195, 262.08, 312.284, 343.845, 339.1, 332.913, 298.026, 248.648, 281.622, 355.802, 352.142, 346.707, 307.605, 344.607, 458.379, 503.564, 575.246, 589.224, 566.922, 444.094, 374.106, 392.472, 429.108, 462.171, 452.507, 402.114, 336.127, 323.267, 274.167, 258.247, 239.089, 255.083, 341.934, 327.816, 356.862, 376.879, 403.153, 391.552, 369.626, 548.004, 534.1, 551.036, 483.023, 407.938, 367.108, 421.151, 424.767, 455.302, 493.192, 477.779, 471.41, 465.075, 474.877, 379.919, 337.061, 418.257, 508.22, 514.854, 504.106, 458.198, 459.895, 387.318, 290.826, 291.78, 269.617, 270.05, 255.638, 300.431, 329.948, 337.552, 365.143, 374.223, 405.786, 402.211, 377.809, 358.27, 345.365, 342.418, 301.744, 340.773, 347.53, 484.977, 568.101, 600.923, 584.527, 555.079, 468.593, 383.9, 329.843, 444.453, 591.095, 623.873, 590.857, 464.361, 399.636, 366.96, 415.829, 422.989, 425.515, 497.843, 520.913, 585.245, 559.485, 455.038, 428.321, 417.609, 423.205, 410.085, 380.559, 458.3, 498.045, 444.246, 343.187, 325.428, 343.338, 374.162, 317.982, 256.119, 262.269, 335.95, 402.572, 393.778, 395.558, 364.628, 270.383, 256.734, 313.647, 396.637, 381.443, 425.537, 384.808, 448.483, 484.637, 479.905, 489.266, 422.141, 335.441, 319.878, 303.223, 327.866, 369.895, 407.826, 433.113, 539.369, 578.454, 574.983, 440.167, 386.225, 382.11};
    cpl_vector * cmp_err = cpl_vector_wrap(width, data_err);

    cpl_image * cmp_model;
    if (access("./tests/debug_model.fits", F_OK) != -1){
        cmp_model = cpl_image_load("./tests/debug_model.fits",  CPL_TYPE_DOUBLE, 0, 0);
    } else {
        cmp_model = cpl_image_load("./debug_model.fits", CPL_TYPE_DOUBLE, 0, 0);
    }

    res = cr2res_extract_slitdec_vert(img_hdrl, trace_table, order, trace, height, swath, oversample, smooth_slit, &slit_func, &spec, &model);

    // Tests compare with expected results
    cpl_test_eq(res, 0);
    cpl_test_vector_abs(slit_func, cmp_slitfunc, FLT_EPSILON);
    cpl_test_vector_abs(cpl_bivector_get_x(spec), cmp_spec, 1e-1);
    cpl_test_vector_abs(cpl_bivector_get_y(spec), cmp_err, 1e-3);
    cpl_test_image_abs(hdrl_image_get_image(model), cmp_model, 1e-4);

    // in case the files are needed again (e.g. when the input is changed)
    // cpl_vector_save(cpl_bivector_get_x(spec), "debug_spec.fits", CPL_TYPE_DOUBLE, NULL,
    //     CPL_IO_CREATE);
    // cpl_vector_save(cpl_bivector_get_y(spec), "debug_err.fits", CPL_TYPE_DOUBLE, NULL,
    //     CPL_IO_CREATE);
    // cpl_vector_save(slit_func, "debug_slitfunc.fits", CPL_TYPE_DOUBLE,
    //     NULL, CPL_IO_CREATE);
    // cpl_image_save(hdrl_image_get_image(model), "debug_model.fits", CPL_TYPE_FLOAT,
    //     NULL, CPL_IO_CREATE);

    // Free memory
    cpl_vector_unwrap(cmp_slitfunc);
    cpl_vector_unwrap(cmp_spec);
    cpl_vector_unwrap(cmp_err);
    cpl_image_delete(cmp_model);

    cpl_image_delete(img_in);
    cpl_table_delete(trace_table);
    cpl_vector_delete(slit_func);
    cpl_bivector_delete(spec);
    hdrl_image_delete(model);
    hdrl_image_delete(img_hdrl);
}

static void test_cr2res_slitdec_vert_edge_cases(void)
{
    int width = 100;
    int height = 20;
    double spec_in[width];
    double shear[width];
    for (int k =0; k < width; k++) shear[k] = 0;
    cpl_image * img_in = create_image_linear_increase(width, height, spec_in);
    hdrl_image * img_hdrl = hdrl_image_create(img_in, NULL);
    cpl_table * trace_table = create_table_linear_increase(width, height, 0);
    int order = 1;
    int trace = 1;
    int swath = 15;
    int oversample = 2;
    double smooth_slit = 0.1;

    cpl_vector * slit_func;
    cpl_bivector * spec;
    hdrl_image * model;
    int res;

    // Test Edge cases
    cpl_test_eq(-1, cr2res_extract_slitdec_vert(NULL, trace_table, order, trace, height, swath,
                oversample, smooth_slit, &slit_func, &spec, &model));
    cpl_test_eq(-1, cr2res_extract_slitdec_vert(img_hdrl, NULL, order, trace, height, swath,
                oversample, smooth_slit, &slit_func, &spec, &model));
    cpl_test_eq(-1, cr2res_extract_slitdec_vert(img_hdrl, trace_table, -10, trace, height, swath,
                oversample, smooth_slit, &slit_func, &spec, &model));
    cpl_test_eq(-1, cr2res_extract_slitdec_vert(img_hdrl, trace_table, order, -3, height, swath,
                oversample, smooth_slit, &slit_func, &spec, &model));

    cpl_test_eq( 0, cr2res_extract_slitdec_vert(img_hdrl, trace_table, order, trace, 0, swath,
                oversample, smooth_slit, &slit_func, &spec, &model));
    cpl_vector_delete(slit_func);
    cpl_bivector_delete(spec);
    hdrl_image_delete(model);

    cpl_test_eq( 0, cr2res_extract_slitdec_vert(img_hdrl, trace_table, order, trace, 10000, swath,
                oversample, smooth_slit, &slit_func, &spec, &model));
    cpl_vector_delete(slit_func);
    cpl_bivector_delete(spec);
    hdrl_image_delete(model);

    cpl_test_eq(-1, cr2res_extract_slitdec_vert(img_hdrl, trace_table, order, trace, height, 0,
                oversample, smooth_slit, &slit_func, &spec, &model));
    cpl_test_eq(-1, cr2res_extract_slitdec_vert(img_hdrl, trace_table, order, trace, height, -1,
                oversample, smooth_slit, &slit_func, &spec, &model));

    cpl_test_eq( 0, cr2res_extract_slitdec_vert(img_hdrl, trace_table, order, trace, height, 10000,
                oversample, smooth_slit, &slit_func, &spec, &model));
    cpl_vector_delete(slit_func);
    cpl_bivector_delete(spec);
    hdrl_image_delete(model);

    cpl_test_eq( 0, cr2res_extract_slitdec_vert(img_hdrl, trace_table, order, trace, height, swath,
                0, smooth_slit, &slit_func, &spec, &model));
    cpl_vector_delete(slit_func);
    cpl_bivector_delete(spec);
    hdrl_image_delete(model);

    cpl_test_eq( 0, cr2res_extract_slitdec_vert(img_hdrl, trace_table, order, trace, height, swath,
                -1, smooth_slit, &slit_func, &spec, &model));
    cpl_vector_delete(slit_func);
    cpl_bivector_delete(spec);
    hdrl_image_delete(model);

    // Free memory
    cpl_image_delete(img_in);
    cpl_table_delete(trace_table);
    hdrl_image_delete(img_hdrl);
}

static void test_cr2res_slitdec_vert(void)
{
    int width = 2000;
    int height = 50;
    double spec_in[width];
    double shear[width];
    for (int k =0; k < width; k++) shear[k] = 0;
    cpl_image * img_in = create_image_linear_increase(width, height, spec_in);
    hdrl_image * img_hdrl = hdrl_image_create(img_in, NULL);
    cpl_table * trace_table = create_table_linear_increase(width, height, 0);
    int order = 1;
    int trace = 1;
    int swath = 400;
    int oversample = 2;
    double smooth_slit = 10;

    cpl_vector * slit_func;
    cpl_bivector * spec;
    hdrl_image * model;
    int res;

    cpl_test_eq( 0, cr2res_extract_slitdec_vert(img_hdrl, trace_table, order, trace, height, swath,
                oversample, smooth_slit, &slit_func, &spec, &model));

    // check results
    double ratio0 = cpl_bivector_get_x_data(spec)[0] / spec_in[0];
    double ratio;
    for(int i = 0; i < width; i++)
    {
        //spectrum, with linear increase
        ratio = cpl_bivector_get_x_data(spec)[i] / spec_in[i];
        cpl_test_abs(ratio, ratio0, 0.2);

        //relative error
        ratio = cpl_bivector_get_y_data(spec)[i] / cpl_bivector_get_x_data(spec)[i];
        cpl_test_lt(ratio,  0.1);
    }

    ratio = cpl_vector_get_sum(slit_func) / oversample;
    cpl_test_abs(ratio, 1, DBL_EPSILON);


    cpl_vector_save(cpl_bivector_get_x(spec), "debug_spec.fits", CPL_TYPE_DOUBLE, NULL,
        CPL_IO_CREATE);
    cpl_vector_save(cpl_bivector_get_y(spec), "debug_err.fits", CPL_TYPE_DOUBLE, NULL,
        CPL_IO_CREATE);
    cpl_vector_save(slit_func, "debug_slitfunc.fits", CPL_TYPE_DOUBLE,
        NULL, CPL_IO_CREATE);
    cpl_image_save(hdrl_image_get_image(model), "debug_model.fits", CPL_TYPE_FLOAT,
        NULL, CPL_IO_CREATE);

    // Free memory
    cpl_vector_delete(slit_func);
    cpl_bivector_delete(spec);
    hdrl_image_delete(model);
    cpl_image_delete(img_in);
    cpl_table_delete(trace_table);
    hdrl_image_delete(img_hdrl);
}

static void test_cr2res_slitdec_curved(void)
{
    int width = 1000;
    int height = 20;
    int order = 1;
    int trace = 1;
    int swath = 400;
    int oversample = 2;
    double smooth_slit = 10;
    double spec_in[width];
    double const_shear = 1;

    cpl_vector * shear = cpl_vector_new(width);
    for (int k = 0; k < width; k++) cpl_vector_set(shear, k, const_shear);

    cpl_image * img_in = create_image_sinusoidal(width, height, spec_in);
    img_in = apply_shear(img_in, width, height, const_shear);
    //cpl_image * img_in = create_image_linear_increase(width, height, cpl_vector_get_data(shear), spec_in);
    hdrl_image * img_hdrl = hdrl_image_create(img_in, NULL);
    cpl_table * trace_table = create_table_linear_increase(width, height, const_shear);

    cpl_vector * slit_func;
    cpl_bivector * spec;
    hdrl_image * model;
    cpl_vector * cmp = cpl_vector_new(width);
    cpl_vector * spec_in_vec = cpl_vector_wrap(width, spec_in);

    int res;
    cpl_test_eq(0,
        cr2res_extract_slitdec_curved(img_hdrl, trace_table, order, trace, height,
                swath, oversample, smooth_slit, &slit_func, &spec, &model));

    // check results
    cpl_vector_copy(cmp, cpl_bivector_get_x(spec));
    cpl_vector_divide(cmp, spec_in_vec);
    double ratio0 = cpl_vector_get_median(cmp);
    double ratio = 0;
    for(int i = 20; i < width-20; i++) {
        // shape of spectrum
        ratio = cpl_bivector_get_x_data(spec)[i] / spec_in[i];
        cpl_test_abs(ratio, ratio0, 2);

        //relative error
        ratio = cpl_bivector_get_y_data(spec)[i] / cpl_bivector_get_x_data(spec)[i];
        cpl_test_lt(ratio,  0.1);
    }

    ratio = cpl_vector_get_sum(slit_func) / oversample;
    cpl_test_abs(ratio, 1, FLT_EPSILON);

    cpl_vector_save(cpl_bivector_get_x(spec), "debug_spec.fits", CPL_TYPE_DOUBLE, NULL,
        CPL_IO_CREATE);
    cpl_vector_save(cpl_bivector_get_y(spec), "debug_err.fits", CPL_TYPE_DOUBLE, NULL,
        CPL_IO_CREATE);
    cpl_vector_save(slit_func, "debug_slitfunc.fits", CPL_TYPE_DOUBLE,
        NULL, CPL_IO_CREATE);
    cpl_image_save(hdrl_image_get_image(model), "debug_model.fits", CPL_TYPE_FLOAT,
        NULL, CPL_IO_CREATE);

    // Free memory
    cpl_vector_delete(slit_func);
    cpl_vector_delete(cmp);
    cpl_vector_unwrap(spec_in_vec);
    cpl_bivector_delete(spec);
    cpl_vector_delete(shear);
    hdrl_image_delete(model);
    cpl_image_delete(img_in);
    cpl_table_delete(trace_table);
    hdrl_image_delete(img_hdrl);
}

static void test_cr2res_slitdec_compare_vert_curved(void)
{
    int width = 1000;
    int height = 21;
    int order = 1;
    int trace = 1;
    int swath = 200;
    int oversample = 5;
    double smooth_slit = 10;
    double spec_in[width];
    double const_shear = 1;


    //cpl_image * img_in = create_image_sinusoidal(width, height * 2, spec_in);
    cpl_image * img_in = create_image_linear_increase(width, height * 2, spec_in);
    //cpl_image * img_in = load_image();
    hdrl_image * img_hdrl = hdrl_image_create(img_in, NULL);
    
    cpl_table * trace_table = create_table_linear_increase(width, height * 2, const_shear);
    //cpl_table * trace_table = load_table();

    cpl_vector * slit_func_vert;
    cpl_bivector * spec_vert;
    hdrl_image * model_vert;
    cpl_vector * slit_func_curved;
    cpl_bivector * spec_curved;
    hdrl_image * model_curved;

    cpl_image_save(img_in, "debug_img_in_vert.fits", CPL_TYPE_DOUBLE, NULL, CPL_IO_CREATE);

    cpl_test_eq(0, cr2res_extract_slitdec_vert(img_hdrl, trace_table, order, trace, height, swath,
                oversample, smooth_slit, &slit_func_vert, &spec_vert, &model_vert));
    

    img_in = apply_shear(img_in, width, height * 2, const_shear);
    hdrl_image_delete(img_hdrl);
    img_hdrl = hdrl_image_create(img_in, NULL);
    cpl_image_save(img_in, "debug_img_in.fits", CPL_TYPE_DOUBLE, NULL, CPL_IO_CREATE);
    cpl_table_save(trace_table, NULL, NULL, "debug_tracetable.fits", CPL_IO_CREATE);

    cpl_test_eq(0, cr2res_extract_slitdec_curved(img_hdrl, trace_table, order, trace, height,
                swath, oversample, smooth_slit * oversample, &slit_func_curved, &spec_curved, &model_curved));



    cpl_vector_save(cpl_bivector_get_x(spec_vert), "debug_spec_vert.fits", CPL_TYPE_DOUBLE, NULL,
        CPL_IO_CREATE);

    cpl_vector_save(slit_func_vert, "debug_slitf_vert.fits", CPL_TYPE_DOUBLE, NULL,
        CPL_IO_CREATE);

    cpl_vector_save(cpl_bivector_get_x(spec_curved), "debug_spec_curved.fits", CPL_TYPE_DOUBLE, NULL,
        CPL_IO_CREATE);


    cpl_vector_save(slit_func_curved, "debug_slitf_curved.fits", CPL_TYPE_DOUBLE, NULL,
        CPL_IO_CREATE);


    // Compare vert and curved
    // cpl_test_vector_abs(cpl_bivector_get_x(spec_vert), cpl_bivector_get_x(spec_curved), FLT_EPSILON);
    // cpl_test_vector_abs(slit_func_vert, slit_func_curved, FLT_EPSILON);

    int delta_x = const_shear * height * 0.5;
    
    cpl_test_abs(cpl_vector_get_sum(slit_func_curved), cpl_vector_get_sum(slit_func_vert), 1e-5);

    for(int i = delta_x; i < width - delta_x - 1; i++)
    {
        cpl_test_abs(cpl_bivector_get_x_data(spec_vert)[i] / cpl_bivector_get_x_data(spec_curved)[i], 1, 0.1);
    }
    

    // Free memory
    cpl_vector_delete(slit_func_vert);
    cpl_bivector_delete(spec_vert);
    hdrl_image_delete(model_vert);

    cpl_vector_delete(slit_func_curved);
    cpl_bivector_delete(spec_curved);
    hdrl_image_delete(model_curved);

    cpl_image_delete(img_in);
    cpl_table_delete(trace_table);
    hdrl_image_delete(img_hdrl);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Run the Unit tests
 */
/*----------------------------------------------------------------------------*/
int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_DEBUG);

    //test_cr2res_slitdec_vert_edge_cases();
    //test_cr2res_slitdec_vert_regression();
    //test_cr2res_slitdec_vert();
    test_cr2res_slitdec_curved();
    //test_cr2res_slitdec_compare_vert_curved();

    return cpl_test_end(0);
}
