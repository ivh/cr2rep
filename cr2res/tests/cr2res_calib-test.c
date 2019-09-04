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
#include <cr2res_calib.h>
#include <cr2res_io.h>
#include "cr2res_pfits.h"

#define MODE_FLAT 0
#define MODE_DARK 1
#define MODE_BPM 2
#define MODE_DETLIN 3

#define pow2(x) (x) * (x)
#define detlin(x, a, b, c) (x) * ((a) + (x) * ((b) + (x) * (c)))
#define deterr(x, a, b, c, sx, sa, sb, sc) sqrt(detlin(pow2(x), pow2(sa), pow2(sb), pow2(sc)) + pow2(sx) * pow2((a) + 2 * (b) * (x) + 3 * (c) * pow2(x)))

static hdrl_image * create_hdrl(int nx, int ny, double value, double error);

static void test_cr2res_calib_image(void);
static void test_cr2res_calib_cosmic(void);
static void test_cr2res_calib_flat(void);
static void test_cr2res_calib_dark(void);
static void test_cr2res_calib_bpm(void);
static void test_cr2res_calib_detlin(void);



/*----------------------------------------------------------------------------*/
/**
  @brief    Save the given hdrl image, in the right way depending on what we want
 */
/*----------------------------------------------------------------------------*/
static void save_hdrl(char * filename, hdrl_image * hdrl, int mode, double dit)
{
    // Empty structures needed for master flat frame, but not for the test
    // Need an empty fits file with header for the framesets
    cpl_frame * empty = cpl_frame_new();
    cpl_frame_set_filename(empty, "./empty.fits");
    cpl_frame_set_tag(empty, "DEBUG");
    cpl_frame_set_group(empty, CPL_FRAME_GROUP_CALIB);

    cpl_parameterlist * parlist = cpl_parameterlist_new();
    cpl_propertylist * ext1 = cpl_propertylist_new();
    cpl_propertylist_append_double(ext1, CR2RES_HEADER_DIT, dit);
    cpl_propertylist * ext[] = {ext1, NULL, NULL};


    cpl_frameset * all = cpl_frameset_new();
    cpl_frameset * in = cpl_frameset_new();
    cpl_frameset_insert(all, empty);
    empty = cpl_frame_duplicate(empty);
    cpl_frameset_insert(in, empty);

    hdrl_image * list[] = {hdrl, NULL, NULL};

    if (mode == MODE_FLAT)
        cr2res_io_save_MASTER_FLAT(filename, all, in, parlist, list, ext1, ext, CR2RES_CAL_FLAT_MASTER_PROCATG, "debug");
    if (mode == MODE_DARK)
        cr2res_io_save_MASTER_DARK(filename, all, in, parlist, list, ext1, ext, CR2RES_CAL_DARK_MASTER_PROCATG, "debug");
    if (mode == MODE_BPM)
    {
        cpl_image * list2[] = {hdrl_image_get_image(hdrl), NULL, NULL};
        cr2res_io_save_BPM(filename, all, in, parlist, list2, ext1, ext, CR2RES_CAL_FLAT_BPM_PROCATG, "debug");    
    }
    if (mode == MODE_DETLIN){
        hdrl_imagelist * list3 = hdrl_imagelist_new();
        hdrl_imagelist_set(list3, hdrl, 0);
        hdrl_imagelist_set(list3, hdrl, 1);
        hdrl_imagelist_set(list3, hdrl, 2);

        hdrl_imagelist * list4[] = {list3, NULL, NULL};

        cr2res_io_save_DETLIN_COEFFS(filename, all, in, parlist, list4, ext1, ext, CR2RES_CAL_DETLIN_COEFFS_PROCATG, "debug");

        // hdrl_imagelist_unset(list3, 0);
        // hdrl_imagelist_unset(list3, 1);
        // hdrl_imagelist_unset(list3, 2);
        hdrl_imagelist_delete(list3);
    }
    cpl_frameset_delete(all);
    cpl_frameset_delete(in);
    cpl_parameterlist_delete(parlist);
    cpl_propertylist_delete(ext1);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Create a Flat Field fits file and cpl_frame
 */
/*----------------------------------------------------------------------------*/
static cpl_frame * create_master_flat(char * filename, int nx, int ny, double value, double error, cpl_mask ** bpm)
{
    cpl_frame * out = cpl_frame_new();
    cpl_frame_set_filename(out, filename);
    cpl_frame_set_tag(out, "FLAT");
    cpl_frame_set_group(out, CPL_FRAME_GROUP_CALIB);

    hdrl_image * hdrl = create_hdrl(nx, ny, value, error);

    if (bpm != NULL){
        hdrl_image_reject(hdrl, 1, 1);
        *bpm = cpl_mask_duplicate(hdrl_image_get_mask(hdrl));
    }

    save_hdrl(filename, hdrl, MODE_FLAT, 0);
    hdrl_image_delete(hdrl);

    return out;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Create a Dark Current fits file and cpl_frame
 */
/*----------------------------------------------------------------------------*/
static cpl_frame * create_master_dark(char * filename, int nx, int ny, double value, double error, double dit, cpl_mask ** bpm)
{

    cpl_frame * out = cpl_frame_new();
    cpl_frame_set_filename(out, filename);
    cpl_frame_set_tag(out, "DARK");
    cpl_frame_set_group(out, CPL_FRAME_GROUP_CALIB);

    hdrl_image * hdrl = create_hdrl(nx, ny, value, error);

    if (bpm != NULL){
        hdrl_image_reject(hdrl, 1, 1);
        *bpm = cpl_mask_duplicate(hdrl_image_get_mask(hdrl));
    }

    save_hdrl(filename, hdrl, MODE_DARK, dit);
    hdrl_image_delete(hdrl);

    return out;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Create a Bad Pixel Map (BPM) fits file and cpl_frame
 */
/*----------------------------------------------------------------------------*/
static cpl_frame * create_bpm(char * filename, int nx, int ny, double value)
{
    cpl_frame * out = cpl_frame_new();
    cpl_frame_set_filename(out, filename);
    cpl_frame_set_tag(out, "BPM");
    cpl_frame_set_group(out, CPL_FRAME_GROUP_CALIB);

    hdrl_image * hdrl = create_hdrl(nx, ny, value, 0);
    save_hdrl(filename, hdrl, MODE_BPM, 0);
    hdrl_image_delete(hdrl);

    return out;   
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Create a Determine Linearity (DetLin) fits file and cpl_frame
 */
/*----------------------------------------------------------------------------*/
static cpl_frame * create_detlin(char * filename, hdrl_image * a, hdrl_image * b, hdrl_image * c)
{
    // Empty structures needed for master flat frame, but not for the test
    // Need an empty fits file with header for the framesets
    cpl_frame * empty = cpl_frame_new();
    cpl_frame_set_filename(empty, "./empty.fits");
    cpl_frame_set_tag(empty, "DEBUG");
    cpl_frame_set_group(empty, CPL_FRAME_GROUP_CALIB);

    cpl_parameterlist * parlist = cpl_parameterlist_new();
    cpl_propertylist * ext1 = cpl_propertylist_new();
    cpl_propertylist * ext[] = {ext1, NULL, NULL};

    cpl_frameset * all = cpl_frameset_new();
    cpl_frameset * in = cpl_frameset_new();
    cpl_frameset_insert(all, empty);
    empty = cpl_frame_duplicate(empty);
    cpl_frameset_insert(in, empty);

    hdrl_imagelist * list3 = hdrl_imagelist_new();
    hdrl_imagelist_set(list3, a, 0);
    hdrl_imagelist_set(list3, b, 1);
    hdrl_imagelist_set(list3, c, 2);

    hdrl_imagelist * list4[] = {list3, NULL, NULL};

    cr2res_io_save_DETLIN_COEFFS(filename, all, in, parlist, list4, ext1, ext, CR2RES_CAL_DETLIN_COEFFS_PROCATG, "debug");

    hdrl_imagelist_unset(list3, 2);
    hdrl_imagelist_unset(list3, 1);
    hdrl_imagelist_unset(list3, 0);

    hdrl_imagelist_delete(list3);
    
    cpl_frameset_delete(all);
    cpl_frameset_delete(in);
    cpl_parameterlist_delete(parlist);
    cpl_propertylist_delete(ext1);

    cpl_frame * out = cpl_frame_new();
    cpl_frame_set_filename(out, filename);
    cpl_frame_set_tag(out, "DETLIN");
    cpl_frame_set_group(out, CPL_FRAME_GROUP_CALIB);

    return out;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Create a hdrl image with constant value and error
 */
/*----------------------------------------------------------------------------*/
static hdrl_image * create_hdrl(int nx, int ny, double value, double error)
{
    hdrl_image * out = hdrl_image_new(nx, ny);
    hdrl_value hv;
    hv.data = value;
    hv.error = error;
    hdrl_image_add_scalar(out, hv);
    return out;
}

static void test_cr2res_calib_image()
{
    int nx = 5;
    int ny = 5;
    int badpix;

    double img_value = 100;
    double img_error = 1;

    hdrl_image * in = create_hdrl(nx, ny, img_value, img_error);
    int chip = 1;
    int cosmics_corr = 1;

    cpl_frame * flat = create_master_flat("master_flat.fits", nx, ny, 1, 0, NULL);
    cpl_frame * dark = create_master_dark("master_dark.fits", nx, ny, 10, 1, 10, NULL);
    cpl_frame * bpm = create_bpm("bpm.fits", nx, ny, 0);
    cpl_frame * detlin = NULL;
    double dit = 10;

    hdrl_image * out;
    hdrl_image * cmp;

    // NULL input / output
    out = cr2res_calib_image(NULL, chip, 0, 0, NULL, NULL, NULL, NULL, dit);
    cpl_test_null(out);

    out = cr2res_calib_image(in, 0, 0, 0, NULL, NULL, NULL, NULL, dit);
    cpl_test_null(out);

    out = cr2res_calib_image(in, CR2RES_NB_DETECTORS + 1, 0, 0, NULL, NULL, NULL, NULL, dit);
    cpl_test_null(out);

    // No correction
    out = cr2res_calib_image(in, chip, 0, 0, NULL, NULL, NULL, NULL, dit);
    cpl_test_nonnull(out);
    cpl_test_image_abs(hdrl_image_get_image(in), hdrl_image_get_image(out), DBL_EPSILON);
    cpl_test_image_abs(hdrl_image_get_error(in), hdrl_image_get_error(out), DBL_EPSILON);
    hdrl_image_delete(out);

    // Test all
    // out = cr2res_calib_image(in, chip, 0, cosmics_corr, flat, dark, bpm, detlin, dit);

    hdrl_image_delete(in);
    cpl_frame_delete(flat);
    cpl_frame_delete(dark);
    cpl_frame_delete(bpm);
    cpl_frame_delete(detlin);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Test Cosmic Correction
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_calib_cosmic()
{
    int nx = 5;
    int ny = 5;
    int badpix;

    double img_value = 100;
    double img_error = 1;
    double out_value, out_error;
    double flat_value, flat_error;

    hdrl_image * in, * out, * cmp;
    cpl_frame * flat;
    int chip = 1;
    int cosmics_corr = 1;
    double dit = 10;

    in = create_hdrl(nx, ny, img_value, img_error);

    // Case 1: No Cosmic Correction, i.e. nothing happens
    out = cr2res_calib_image(in, chip, 0, 0, NULL, NULL, NULL, NULL, dit);
    cpl_test_image_abs(hdrl_image_get_image(in), hdrl_image_get_image(out), DBL_EPSILON);
    cpl_test_image_abs(hdrl_image_get_error(in), hdrl_image_get_error(out), DBL_EPSILON);

    // Case 2: Cosmic Correction
    // TODO when cosmic correction is implemented

    hdrl_image_delete(in);
    hdrl_image_delete(out);

}


/*----------------------------------------------------------------------------*/
/**
  @brief    Only test the flat field calibration
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_calib_flat()
{
    int nx = 5;
    int ny = 5;
    int badpix;

    double img_value = 100;
    double img_error = 1;
    double out_value, out_error;
    double flat_value, flat_error;

    hdrl_image * in, * out, * cmp;
    cpl_mask * bpm;
    cpl_frame * flat;
    int chip = 1;
    int cosmics_corr = 1;
    double dit = 10;
    cpl_error_code error;

    in = create_hdrl(nx, ny, img_value, img_error);

    // Case 1: Flat is just 1 and no error, i.e. no change
    flat_value = 1;
    flat_error = 0;
    flat = create_master_flat("master_flat.fits", nx, ny, flat_value, flat_error, NULL);
    out = cr2res_calib_image(in, chip, 0, 0, flat, NULL, NULL, NULL, dit);
    
    out_value = img_value / flat_value;
    out_error = sqrt( pow(img_error / flat_value, 2) + pow(img_value * flat_error/ (flat_value * flat_value), 2));
    cmp = create_hdrl(nx, ny, out_value, out_error);

    cpl_test_image_abs(hdrl_image_get_image(cmp), hdrl_image_get_image(out), DBL_EPSILON);
    cpl_test_image_abs(hdrl_image_get_error(cmp), hdrl_image_get_error(out), DBL_EPSILON);

    hdrl_image_delete(out);
    hdrl_image_delete(cmp);
    cpl_frame_delete(flat);

    // Case 2: Flat is constant == 2, error 1 
    flat_value = 2;
    flat_error = 1;
    flat = create_master_flat("master_flat.fits", nx, ny, flat_value, flat_error, NULL);
    out = cr2res_calib_image(in, chip, 0, 0, flat, NULL, NULL, NULL, dit);
    
    out_value = img_value / flat_value;
    out_error = sqrt( pow(img_error / flat_value, 2) + pow(img_value * flat_error/ (flat_value * flat_value), 2));
    cmp = create_hdrl(nx, ny, out_value, out_error);

    cpl_test_image_abs(hdrl_image_get_image(cmp), hdrl_image_get_image(out), DBL_EPSILON);
    cpl_test_image_abs(hdrl_image_get_error(cmp), hdrl_image_get_error(out), DBL_EPSILON);

    hdrl_image_delete(out);
    hdrl_image_delete(cmp);
    cpl_frame_delete(flat);

    // Case 3: Flat is 0, i.e. all pixels become bad
    flat_value = 0;
    flat_error = 1;
    flat = create_master_flat("master_flat.fits", nx, ny, flat_value, flat_error, NULL);
    out = cr2res_calib_image(in, chip, 0, 0, flat, NULL, NULL, NULL, dit);
    
    cpl_test_eq(nx * ny, hdrl_image_count_rejected(out));

    hdrl_image_delete(out);
    cpl_frame_delete(flat);

    // Case 4: Flat file does not exist
    flat = cpl_frame_new();
    cpl_frame_set_filename(flat, "tobeornottobe.fits");
    cpl_frame_set_tag(flat, "FLAT");
    cpl_frame_set_group(flat, CPL_FRAME_GROUP_CALIB);

    out = cr2res_calib_image(in, chip, 0, 0, flat, NULL, NULL, NULL, dit);
    cpl_test_error(CPL_ERROR_FILE_IO);
    hdrl_image_delete(out);
    cpl_frame_delete(flat);

    // Case 5: image is in a wrong group
    flat = create_master_dark("master_flat.fits", nx, ny, 1, 0, 10, NULL);
    out = cr2res_calib_image(in, chip, 0, 0, flat, NULL, NULL, NULL, dit);
    cpl_test_null(out);
    cpl_frame_delete(flat);

    // Case 6: No Filename set
    flat = cpl_frame_new();
    out = cr2res_calib_image(in, chip, 0, 0, flat, NULL, NULL, NULL, dit);
    cpl_test_error(CPL_ERROR_DATA_NOT_FOUND);
    cpl_test_null(out);
    cpl_frame_delete(flat);

    // Case 7: DataFile is empty, i.e. only header
    flat = cpl_frame_new();
    cpl_frame_set_filename(flat, "empty.fits");
    out = cr2res_calib_image(in, chip, 0, 0, flat, NULL, NULL, NULL, dit);
    cpl_test_null(out);
    cpl_frame_delete(flat);

    // // Case 8: BPM is set 
    // flat_value = 1;
    // flat_error = 1;
    // flat = create_master_flat("master_flat.fits", nx, ny, flat_value, flat_error, &bpm);
    // out = cr2res_calib_image(in, chip, 0, 0, flat, NULL, NULL, NULL, dit);
    
    // cpl_test_eq_mask(bpm, hdrl_image_get_mask(out));

    // hdrl_image_delete(out);
    // cpl_frame_delete(flat);
    // cpl_mask_delete(bpm);


    hdrl_image_delete(in);
}

static void test_cr2res_calib_dark()
{
    int nx = 5;
    int ny = 5;
    int badpix;

    double img_value = 100;
    double img_error = 1;
    double out_value, out_error;
    double dark_value, dark_error, dark_dit;

    hdrl_image * in, * out, * cmp;
    cpl_mask *bpm;
    cpl_frame * dark;
    int chip = 1;
    int cosmics_corr = 1;
    double dit = 10;

    in = create_hdrl(nx, ny, img_value, img_error);

    // Case 1: Dark is 0, dit is the same as the image, i.e. no correction
    dark_value = 0;
    dark_error = 0;
    dark_dit = dit;
    dark = create_master_dark("master_dark.fits", nx, ny, dark_value, dark_error, dark_dit, NULL);
    out = cr2res_calib_image(in, chip, 0, 0, NULL, dark, NULL, NULL, dit);
    
    out_value = img_value - dit/dark_dit * dark_value;
    out_error = sqrt(pow(img_error, 2) + pow(dit/dark_dit * dark_error, 2));
    cmp = create_hdrl(nx, ny, out_value, out_error);

    cpl_test_image_abs(hdrl_image_get_image(cmp), hdrl_image_get_image(out), DBL_EPSILON);
    cpl_test_image_abs(hdrl_image_get_error(cmp), hdrl_image_get_error(out), DBL_EPSILON);

    hdrl_image_delete(out);
    hdrl_image_delete(cmp);
    cpl_frame_delete(dark);

    // Case 2: Dark is 0, dit is different, i.e. still nothing
    dark_value = 0;
    dark_error = 0;
    dark_dit = dit * 2.3;
    dark = create_master_dark("master_dark.fits", nx, ny, dark_value, dark_error, dark_dit, NULL);
    out = cr2res_calib_image(in, chip, 0, 0, NULL, dark, NULL, NULL, dit);
    
    out_value = img_value - dit/dark_dit * dark_value;
    out_error = sqrt(pow(img_error, 2) + pow(dit/dark_dit * dark_error, 2));
    cmp = create_hdrl(nx, ny, out_value, out_error);

    cpl_test_image_abs(hdrl_image_get_image(cmp), hdrl_image_get_image(out), DBL_EPSILON);
    cpl_test_image_abs(hdrl_image_get_error(cmp), hdrl_image_get_error(out), DBL_EPSILON);

    hdrl_image_delete(out);
    hdrl_image_delete(cmp);
    cpl_frame_delete(dark);

    // Case 3: Dark is 0, dit is same, but dark has error now
    dark_value = 0;
    dark_error = 2;
    dark_dit = dit;
    dark = create_master_dark("master_dark.fits", nx, ny, dark_value, dark_error, dark_dit, NULL);
    out = cr2res_calib_image(in, chip, 0, 0, NULL, dark, NULL, NULL, dit);
    
    out_value = img_value - dit/dark_dit * dark_value;
    out_error = sqrt(pow(img_error, 2) + pow(dit/dark_dit * dark_error, 2));
    cmp = create_hdrl(nx, ny, out_value, out_error);

    cpl_test_image_abs(hdrl_image_get_image(cmp), hdrl_image_get_image(out), DBL_EPSILON);
    cpl_test_image_abs(hdrl_image_get_error(cmp), hdrl_image_get_error(out), DBL_EPSILON);

    hdrl_image_delete(out);
    hdrl_image_delete(cmp);
    cpl_frame_delete(dark);

    // Case 4: Dark is constant 10, dit is same, and has error
    dark_value = 10;
    dark_error = 2;
    dark_dit = dit;
    dark = create_master_dark("master_dark.fits", nx, ny, dark_value, dark_error, dark_dit, NULL);
    out = cr2res_calib_image(in, chip, 0, 0, NULL, dark, NULL, NULL, dit);
    
    out_value = img_value - dit/dark_dit * dark_value;
    out_error = sqrt(pow(img_error, 2) + pow(dit/dark_dit * dark_error, 2));
    cmp = create_hdrl(nx, ny, out_value, out_error);

    cpl_test_image_abs(hdrl_image_get_image(cmp), hdrl_image_get_image(out), DBL_EPSILON);
    cpl_test_image_abs(hdrl_image_get_error(cmp), hdrl_image_get_error(out), DBL_EPSILON);

    hdrl_image_delete(out);
    hdrl_image_delete(cmp);
    cpl_frame_delete(dark);

    // Case 5: Dark is constant 10, dit is different, and has error
    dark_value = 10;
    dark_error = 2;
    dark_dit = dit * 2.3;
    dark = create_master_dark("master_dark.fits", nx, ny, dark_value, dark_error, dark_dit, NULL);
    out = cr2res_calib_image(in, chip, 0, 0, NULL, dark, NULL, NULL, dit);
    
    out_value = img_value - dit/dark_dit * dark_value;
    out_error = sqrt(pow(img_error, 2) + pow(dit/dark_dit * dark_error, 2));
    cmp = create_hdrl(nx, ny, out_value, out_error);

    cpl_test_image_abs(hdrl_image_get_image(cmp), hdrl_image_get_image(out), DBL_EPSILON);
    cpl_test_image_abs(hdrl_image_get_error(cmp), hdrl_image_get_error(out), DBL_EPSILON);

    hdrl_image_delete(out);
    hdrl_image_delete(cmp);
    cpl_frame_delete(dark);

    // Case 6: Dark is negative, dit positive
    dark_value = -10;
    dark_error = 2;
    dark_dit = dit * 2.3;
    dark = create_master_dark("master_dark.fits", nx, ny, dark_value, dark_error, dark_dit, NULL);
    out = cr2res_calib_image(in, chip, 0, 0, NULL, dark, NULL, NULL, dit);
    
    out_value = img_value - dit/dark_dit * dark_value;
    out_error = sqrt(pow(img_error, 2) + pow(dit/dark_dit * dark_error, 2));
    cmp = create_hdrl(nx, ny, out_value, out_error);

    cpl_test_image_abs(hdrl_image_get_image(cmp), hdrl_image_get_image(out), DBL_EPSILON);
    cpl_test_image_abs(hdrl_image_get_error(cmp), hdrl_image_get_error(out), DBL_EPSILON);

    hdrl_image_delete(out);
    hdrl_image_delete(cmp);
    cpl_frame_delete(dark);

    // Case 7: Dark is positive, dit negative
    dark_value = 10;
    dark_error = 2;
    dark_dit = -dit * 2.3;
    dark = create_master_dark("master_dark.fits", nx, ny, dark_value, dark_error, dark_dit, NULL);
    out = cr2res_calib_image(in, chip, 0, 0, NULL, dark, NULL, NULL, dit);
    
    out_value = img_value - dit/dark_dit * dark_value;
    out_error = sqrt(pow(img_error, 2) + pow(dit/dark_dit * dark_error, 2));
    cmp = create_hdrl(nx, ny, out_value, out_error);

    cpl_test_image_abs(hdrl_image_get_image(cmp), hdrl_image_get_image(out), DBL_EPSILON);
    cpl_test_image_abs(hdrl_image_get_error(cmp), hdrl_image_get_error(out), DBL_EPSILON);

    hdrl_image_delete(out);
    hdrl_image_delete(cmp);
    cpl_frame_delete(dark);

    // Case 8: Dark and dit negative
    dark_value = -10;
    dark_error = 2;
    dark_dit = -dit * 2.3;
    dark = create_master_dark("master_dark.fits", nx, ny, dark_value, dark_error, dark_dit, NULL);
    out = cr2res_calib_image(in, chip, 0, 0, NULL, dark, NULL, NULL, dit);
    
    out_value = img_value - dit/dark_dit * dark_value;
    out_error = sqrt(pow(img_error, 2) + pow(dit/dark_dit * dark_error, 2));
    cmp = create_hdrl(nx, ny, out_value, out_error);

    cpl_test_image_abs(hdrl_image_get_image(cmp), hdrl_image_get_image(out), DBL_EPSILON);
    cpl_test_image_abs(hdrl_image_get_error(cmp), hdrl_image_get_error(out), DBL_EPSILON);

    hdrl_image_delete(out);
    hdrl_image_delete(cmp);
    cpl_frame_delete(dark);

    // Case 9: Dark file does not exist
    dark = cpl_frame_new();
    cpl_frame_set_filename(dark, "tobeornottobe.fits");
    cpl_frame_set_tag(dark, "DARK");
    cpl_frame_set_group(dark, CPL_FRAME_GROUP_CALIB);

    out = cr2res_calib_image(in, chip, 0, 0, NULL, dark, NULL, NULL, dit);
    cpl_test_error(CPL_ERROR_FILE_IO);
    cpl_test_null(out);
    cpl_frame_delete(dark);

    // Case 10: image is in a wrong group
    dark = create_master_flat("master_dark.fits", nx, ny, 1, 0, NULL);
    out = cr2res_calib_image(in, chip, 0, 0, NULL, dark, NULL, NULL, dit);
    cpl_test_null(out);
    cpl_frame_delete(dark);

    // Case 11: No Filename set
    dark = cpl_frame_new();
    out = cr2res_calib_image(in, chip, 0, 0, NULL, dark, NULL, NULL, dit);
    cpl_test_error(CPL_ERROR_DATA_NOT_FOUND);
    cpl_test_null(out);
    cpl_frame_delete(dark);

    // Case 12: DataFile is empty, i.e. only header
    dark = cpl_frame_new();
    cpl_frame_set_filename(dark, "empty.fits");
    out = cr2res_calib_image(in, chip, 0, 0, NULL, dark, NULL, NULL, dit);
    cpl_test_null(out);
    cpl_frame_delete(dark);

    // // Case 13: Check that bpm of dark is passed on to out
    // dark_value = 10;
    // dark_error = 0;
    // dark_dit = 10;
    // dark = create_master_dark("master_dark.fits", nx, ny, dark_value, dark_error, dark_dit, &bpm);
    // out = cr2res_calib_image(in, chip, 0, 0, NULL, dark, NULL, NULL, dit);
    
    // cpl_test_eq_mask(bpm, hdrl_image_get_mask(out));

    // cpl_mask_delete(bpm);
    // hdrl_image_delete(out);
    // cpl_frame_delete(dark);

    hdrl_image_delete(in);
}

static void test_cr2res_calib_bpm()
{
    int nx = 5;
    int ny = 5;
    int badpix;

    double img_value = 100;
    double img_error = 1;
    double out_value, out_error;
    double bpm_value;
    hdrl_image * hdrl;
    hdrl_value value;

    hdrl_image * in, * out, * cmp;
    cpl_frame * bpm;
    int chip = 1;
    int cosmics_corr = 1;
    double dit = 10;

    in = create_hdrl(nx, ny, img_value, img_error);
    value.data = 10; value.error = img_error;
    hdrl_image_set_pixel(in, 2, 2, value);

    // Case 1: BPM is all 0, i.e. all good pixels
    bpm = create_bpm("master_bpm.fits", nx, ny, 0);
    out = cr2res_calib_image(in, chip, 0, 0, NULL, NULL, bpm, NULL, dit);
    cpl_test_image_abs(hdrl_image_get_image(out), hdrl_image_get_image(in), DBL_EPSILON);
    cpl_test_image_abs(hdrl_image_get_error(out), hdrl_image_get_error(in), DBL_EPSILON);
    cpl_frame_delete(bpm);
    hdrl_image_delete(out);

    // Case 2: BPM is all 1, i.e all pixels are bad
    bpm = create_bpm("master_bpm.fits", nx, ny, 1);
    out = cr2res_calib_image(in, chip, 0, 0, NULL, NULL, bpm, NULL, dit);
    cpl_test_nonnull(out);
    cpl_test_eq(hdrl_image_count_rejected(out), nx * ny);
    cpl_frame_delete(bpm);
    hdrl_image_delete(out);

    // Case 3: One Bad Pixel
    bpm = cpl_frame_new();
    cpl_frame_set_filename(bpm, "master_bpm.fits");
    cpl_frame_set_tag(bpm, "BPM");
    cpl_frame_set_group(bpm, CPL_FRAME_GROUP_CALIB);

    hdrl = create_hdrl(nx, ny, 0, 0);
    value.data = 1;
    hdrl_image_set_pixel(hdrl, 2, 2, value);
    save_hdrl("master_bpm.fits", hdrl, MODE_BPM, 0);
    hdrl_image_delete(hdrl);

    out = cr2res_calib_image(in, chip, 0, 0, NULL, NULL, bpm, NULL, dit);

    cpl_test_eq(hdrl_image_count_rejected(out), 1);
    cpl_test_image_abs(hdrl_image_get_image(out), hdrl_image_get_image(in), DBL_EPSILON);
    cpl_test_image_abs(hdrl_image_get_error(out), hdrl_image_get_error(in), DBL_EPSILON);

    cpl_frame_delete(bpm);
    hdrl_image_delete(out);

    // Case 4: Only one good pixel
    bpm = cpl_frame_new();
    cpl_frame_set_filename(bpm, "master_bpm.fits");
    cpl_frame_set_tag(bpm, "BPM");
    cpl_frame_set_group(bpm, CPL_FRAME_GROUP_CALIB);

    hdrl = create_hdrl(nx, ny, 1, 0);
    value.data = 0;
    hdrl_image_set_pixel(hdrl, 2, 2, value);
    save_hdrl("master_bpm.fits", hdrl, MODE_BPM, 0);
    hdrl_image_delete(hdrl);

    out = cr2res_calib_image(in, chip, 0, 0, NULL, NULL, bpm, NULL, dit);
    cpl_test_eq(hdrl_image_count_rejected(out), nx*ny - 1);
    cpl_test_image_abs(hdrl_image_get_image(out), hdrl_image_get_image(in), DBL_EPSILON);
    cpl_test_image_abs(hdrl_image_get_error(out), hdrl_image_get_error(in), DBL_EPSILON);

    cpl_frame_delete(bpm);
    hdrl_image_delete(out);

    // Case 5: BPM file does not exist
    bpm = cpl_frame_new();
    cpl_frame_set_filename(bpm, "tobeornottobe.fits");
    cpl_frame_set_tag(bpm, "BPM");
    cpl_frame_set_group(bpm, CPL_FRAME_GROUP_CALIB);

    out = cr2res_calib_image(in, chip, 0, 0, NULL, NULL, bpm, NULL, dit);
    cpl_test_error(CPL_ERROR_FILE_IO);
    cpl_test_null(out);
    cpl_frame_delete(bpm);

    // Case 6: image is in a wrong group
    bpm = create_master_flat("master_bpm.fits", nx, ny, 1, 0, NULL);
    out = cr2res_calib_image(in, chip, 0, 0, NULL, NULL, bpm, NULL, dit);
    cpl_test_null(out);
    cpl_frame_delete(bpm);

    // Case 7: No Filename set
    bpm = cpl_frame_new();
    out = cr2res_calib_image(in, chip, 0, 0, NULL, NULL, bpm, NULL, dit);
    cpl_test_error(CPL_ERROR_DATA_NOT_FOUND);
    cpl_test_null(out);
    cpl_frame_delete(bpm);

    // Case 8: DataFile is empty, i.e. only header
    bpm = cpl_frame_new();
    cpl_frame_set_filename(bpm, "empty.fits");
    out = cr2res_calib_image(in, chip, 0, 0, NULL, NULL, bpm, NULL, dit);
    cpl_test_null(out);
    cpl_frame_delete(bpm);

    hdrl_image_delete(in);
}

static void test_cr2res_calib_detlin()
{
    int nx = 5;
    int ny = 5;
    int badpix;

    double img_value = 100;
    double img_error = 1;
    double out_value, out_error;
    double bpm_value;
    hdrl_image * hdrl;
    hdrl_value value;

    hdrl_image * in, * out, * cmp;
    hdrl_image * ima, * imb, * imc;
    double a, b, c, sa, sb, sc;
    double i, si;
    cpl_frame * detlin;
    int chip = 1;
    int cosmics_corr = 1;
    double dit = 10;

    in = create_hdrl(nx, ny, img_value, img_error);
    i = img_value; si = img_error;

    // DetLin images aren't actually images, but polynomial coefficients for each pixel
    // The new image is the evaluation of that polynomial for each pixel
    // new image = old_image * detlin_poly(old_image)

    // Case 1: DetLin a = c = 0, b = 1, doesn't change anything
    a = 1; c = 0; b = 0;
    sa = sb = sc = 0;
    ima = create_hdrl(nx, ny, a, sa);
    imb = create_hdrl(nx, ny, b, sb);
    imc = create_hdrl(nx, ny, c, sc);

    detlin = create_detlin("master_detlin.fits", ima, imb, imc);
    out = cr2res_calib_image(in, chip, 0, 0, NULL, NULL, NULL, detlin, dit);
    cpl_test_image_abs(hdrl_image_get_image(out), hdrl_image_get_image(in), DBL_EPSILON);
    cpl_test_image_abs(hdrl_image_get_error(out), hdrl_image_get_error(in), DBL_EPSILON);
    hdrl_image_delete(out);
    cpl_frame_delete(detlin);
    hdrl_image_delete(ima);
    hdrl_image_delete(imb);
    hdrl_image_delete(imc);

    // Case 2: DetLin a = 1, c = 0, b = 1, x + x**2
    a = b = 1; c = 0;
    sa = sb = sc = 0;
    ima = create_hdrl(nx, ny, a, sa);
    imb = create_hdrl(nx, ny, b, sb);
    imc = create_hdrl(nx, ny, c, sc);
    detlin = create_detlin("master_detlin.fits", ima, imb, imc);
    out = cr2res_calib_image(in, chip, 0, 0, NULL, NULL, NULL, detlin, dit);
    cmp = create_hdrl(nx, ny, detlin(img_value, a, b, c), deterr(img_value, a, b, c, img_error, sa, sb, sc));
    cpl_test_image_abs(hdrl_image_get_image(out), hdrl_image_get_image(cmp), DBL_EPSILON);
    cpl_test_image_abs(hdrl_image_get_error(out), hdrl_image_get_error(cmp), DBL_EPSILON);
    hdrl_image_delete(out);
    hdrl_image_delete(cmp);
    cpl_frame_delete(detlin);
    hdrl_image_delete(ima);
    hdrl_image_delete(imb);
    hdrl_image_delete(imc);


    // Case 3: DetLin a = b = c = 0, everything is 0
    a = c = b = 0;
    sa = sb = sc = 0;
    ima = create_hdrl(nx, ny, a, sa);
    imb = create_hdrl(nx, ny, b, sb);
    imc = create_hdrl(nx, ny, c, sc);
    detlin = create_detlin("master_detlin.fits", ima, imb, imc);
    out = cr2res_calib_image(in, chip, 0, 0, NULL, NULL, NULL, detlin, dit);
    cmp = create_hdrl(nx, ny, 0, 0);
    cpl_test_image_abs(hdrl_image_get_image(out), hdrl_image_get_image(cmp), DBL_EPSILON);
    cpl_test_image_abs(hdrl_image_get_error(out), hdrl_image_get_error(cmp), DBL_EPSILON);
    hdrl_image_delete(out);
    hdrl_image_delete(cmp);
    cpl_frame_delete(detlin);
    hdrl_image_delete(ima);
    hdrl_image_delete(imb);
    hdrl_image_delete(imc);

    // Case 4: DetLin a = 0.5, c = b = 0, divide all values by 2, including the error
    a = 0.5; b = c = 0;
    sa = sb = sc = 0;
    ima = create_hdrl(nx, ny, a, sa);
    imb = create_hdrl(nx, ny, b, sb);
    imc = create_hdrl(nx, ny, c, sc);
    detlin = create_detlin("master_detlin.fits", ima, imb, imc);
    out = cr2res_calib_image(in, chip, 0, 0, NULL, NULL, NULL, detlin, dit);
    cmp = create_hdrl(nx, ny, img_value * a, img_error * a);
    cpl_test_image_abs(hdrl_image_get_image(out), hdrl_image_get_image(cmp), DBL_EPSILON);
    cpl_test_image_abs(hdrl_image_get_error(out), hdrl_image_get_error(cmp), DBL_EPSILON);
    hdrl_image_delete(out);
    hdrl_image_delete(cmp);
    cpl_frame_delete(detlin);
    hdrl_image_delete(ima);
    hdrl_image_delete(imb);
    hdrl_image_delete(imc);

    // Case 5: DetLin a = 1, c = 0, b = 2
    // x = a + y * b
    a = 2; b = 0; c = 0;
    sa = sb = sc = 0;
    ima = create_hdrl(nx, ny, a, sa);
    imb = create_hdrl(nx, ny, b, sb);
    imc = create_hdrl(nx, ny, c, sc);
    detlin = create_detlin("master_detlin.fits", ima, imb, imc);
    out = cr2res_calib_image(in, chip, 0, 0, NULL, NULL, NULL, detlin, dit);
    cmp = create_hdrl(nx, ny, detlin(i, a, b, c), deterr(i, a, b, c, si, sa, sb, sc));
    cpl_test_image_abs(hdrl_image_get_image(out), hdrl_image_get_image(cmp), DBL_EPSILON);
    cpl_test_image_abs(hdrl_image_get_error(out), hdrl_image_get_error(cmp), DBL_EPSILON);
    hdrl_image_delete(out);
    hdrl_image_delete(cmp);
    cpl_frame_delete(detlin);
    hdrl_image_delete(ima);
    hdrl_image_delete(imb);
    hdrl_image_delete(imc);

    // Case 6: DetLin a = 1, c = 1, b = 2
    // x = a + y * b + y**2 * c
    a = c = 1; b = 2;
    sa = sb = sc = 0;
    ima = create_hdrl(nx, ny, a, sa);
    imb = create_hdrl(nx, ny, b, sb);
    imc = create_hdrl(nx, ny, c, sc);
    detlin = create_detlin("master_detlin.fits", ima, imb, imc);
    out = cr2res_calib_image(in, chip, 0, 0, NULL, NULL, NULL, detlin, dit);
    cmp = create_hdrl(nx, ny, detlin(i, a, b, c), deterr(i, a, b, c, si, sa, sb, sc));
    cpl_test_image_abs(hdrl_image_get_image(out), hdrl_image_get_image(cmp), DBL_EPSILON);
    cpl_test_image_abs(hdrl_image_get_error(out), hdrl_image_get_error(cmp), DBL_EPSILON);
    hdrl_image_delete(out);
    hdrl_image_delete(cmp);
    cpl_frame_delete(detlin);
    hdrl_image_delete(ima);
    hdrl_image_delete(imb);
    hdrl_image_delete(imc);

    // Case 7: DetLin a = 0, c = 0, b = 1, but now with error
    // x = b * y
    a = c = 0; b = 1;
    sa = sc = 0; sb = 1;
    ima = create_hdrl(nx, ny, a, sa);
    imb = create_hdrl(nx, ny, b, sb);
    imc = create_hdrl(nx, ny, c, sc);
    detlin = create_detlin("master_detlin.fits", ima, imb, imc);
    out = cr2res_calib_image(in, chip, 0, 0, NULL, NULL, NULL, detlin, dit);
    cmp = create_hdrl(nx, ny, pow2(img_value), deterr(i, a, b, c, si, sa, sb, sc));
    cpl_test_image_abs(hdrl_image_get_image(out), hdrl_image_get_image(cmp), DBL_EPSILON);
    cpl_test_image_abs(hdrl_image_get_error(out), hdrl_image_get_error(cmp), DBL_EPSILON);
    hdrl_image_delete(out);
    hdrl_image_delete(cmp);
    cpl_frame_delete(detlin);
    hdrl_image_delete(ima);
    hdrl_image_delete(imb);
    hdrl_image_delete(imc);

    // Case 8: DetLin a = 1, c = 0, b = 1, but now with error
    a = 1; c = 0; b = 1;
    sa = 1 ; sc = 0; sb = 1;
    ima = create_hdrl(nx, ny, a, sa);
    imb = create_hdrl(nx, ny, b, sb);
    imc = create_hdrl(nx, ny, c, sc);
    detlin = create_detlin("master_detlin.fits", ima, imb, imc);
    out = cr2res_calib_image(in, chip, 0, 0, NULL, NULL, NULL, detlin, dit);
    cmp = create_hdrl(nx, ny, detlin(i, a, b, c), deterr(i, a, b, c, si, sa, sb, sc));
    cpl_test_image_abs(hdrl_image_get_image(out), hdrl_image_get_image(cmp), DBL_EPSILON);
    cpl_test_image_abs(hdrl_image_get_error(out), hdrl_image_get_error(cmp), DBL_EPSILON);
    hdrl_image_delete(out);
    hdrl_image_delete(cmp);
    cpl_frame_delete(detlin);
    hdrl_image_delete(ima);
    hdrl_image_delete(imb);
    hdrl_image_delete(imc);

    // Case 9: DetLin a = 1, c = 1, b = 0, but now with error
    // x = a + c * y**2
    a = 1; c = 1; b = 0;
    sa = 1 ; sc = 1; sb = 0;
    i = img_value; si = img_error;
    ima = create_hdrl(nx, ny, a, sa);
    imb = create_hdrl(nx, ny, b, sb);
    imc = create_hdrl(nx, ny, c, sc);
    detlin = create_detlin("master_detlin.fits", ima, imb, imc);
    out = cr2res_calib_image(in, chip, 0, 0, NULL, NULL, NULL, detlin, dit);
    cmp = create_hdrl(nx, ny, detlin(i, a, b, c), deterr(i, a, b, c, si, sa, sb, sc));
    cpl_test_image_abs(hdrl_image_get_image(out), hdrl_image_get_image(cmp), DBL_EPSILON);
    cpl_test_image_abs(hdrl_image_get_error(out), hdrl_image_get_error(cmp), DBL_EPSILON);
    hdrl_image_delete(out);
    hdrl_image_delete(cmp);
    cpl_frame_delete(detlin);
    hdrl_image_delete(ima);
    hdrl_image_delete(imb);
    hdrl_image_delete(imc);

    // Case 10: DetLin a = 1, c = 1, b = 1, but now with error
    // x = a + b * y + c * y**2
    a = 1; c = 1; b = 1;
    sa = 1 ; sc = 1; sb = 1;
    i = img_value; si = img_error;
    ima = create_hdrl(nx, ny, a, sa);
    imb = create_hdrl(nx, ny, b, sb);
    imc = create_hdrl(nx, ny, c, sc);
    detlin = create_detlin("master_detlin.fits", ima, imb, imc);
    out = cr2res_calib_image(in, chip, 0, 0, NULL, NULL, NULL, detlin, dit);
    cmp = create_hdrl(nx, ny, detlin(i, a, b, c), deterr(i, a, b, c, si, sa, sb, sc));
    cpl_test_image_abs(hdrl_image_get_image(out), hdrl_image_get_image(cmp), DBL_EPSILON);
    cpl_test_image_abs(hdrl_image_get_error(out), hdrl_image_get_error(cmp), 1e-6);
    hdrl_image_delete(out);
    hdrl_image_delete(cmp);
    cpl_frame_delete(detlin);
    hdrl_image_delete(ima);
    hdrl_image_delete(imb);
    hdrl_image_delete(imc);

    // Case 11: No Filename set
    detlin = cpl_frame_new();
    out = cr2res_calib_image(in, chip, 0, 0, NULL, NULL, NULL, detlin, dit);
    cpl_test_error(CPL_ERROR_DATA_NOT_FOUND);
    cpl_test_null(out);
    cpl_frame_delete(detlin);

    // Case 12: DataFile is empty, i.e. only header
    detlin = cpl_frame_new();
    cpl_frame_set_filename(detlin, "empty.fits");
    out = cr2res_calib_image(in, chip, 0, 0, NULL, NULL, NULL, detlin, dit);
    cpl_test_null(out);
    cpl_frame_delete(detlin);

    hdrl_image_delete(in);

}

/*----------------------------------------------------------------------------*/
/**
  @brief    Run the Unit tests
 */
/*----------------------------------------------------------------------------*/
int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_DEBUG);

    test_cr2res_calib_image();
    test_cr2res_calib_cosmic();
    test_cr2res_calib_flat();
    test_cr2res_calib_dark();
    test_cr2res_calib_bpm();
    test_cr2res_calib_detlin();

    return cpl_test_end(0);
}
