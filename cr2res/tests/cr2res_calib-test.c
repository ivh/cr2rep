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

#define MODE_FLAT 0
#define MODE_DARK 1
#define MODE_BPM 2
#define MODE_DETLIN 3

static hdrl_image * create_hdrl(int nx, int ny, double value, double error);

static void test_cr2res_calib_image(void);
static void test_cr2res_calib_cosmic(void);
static void test_cr2res_calib_flat(void);
static void test_cr2res_calib_dark(void);
static void test_cr2res_calib_bpm(void);



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
    cpl_propertylist_append_double(ext1,"ESO DET SEQ1 DIT", dit);
    cpl_propertylist * ext[] = {ext1, NULL, NULL};


    cpl_frameset * all = cpl_frameset_new();
    cpl_frameset * in = cpl_frameset_new();
    cpl_frameset_insert(all, empty);
    empty = cpl_frame_duplicate(empty);
    cpl_frameset_insert(in, empty);

    hdrl_image * list[] = {hdrl, NULL, NULL};

    if (mode == MODE_FLAT)
        cr2res_io_save_MASTER_FLAT(filename, all, in, parlist, list, ext1, ext, CR2RES_FLAT_MASTER_FLAT_PROCATG, "debug");
    if (mode == MODE_DARK)
        cr2res_io_save_MASTER_DARK(filename, all, in, parlist, list, ext1, ext, CR2RES_MASTER_DARK_PROCATG, "debug");
    if (mode == MODE_BPM)
    {
        cpl_image * list2[] = {hdrl_image_get_image(hdrl), NULL, NULL};
        cr2res_io_save_BPM(filename, all, in, parlist, list2, ext1, ext, CR2RES_FLAT_BPM_PROCATG, "debug");    
    }
    if (mode == MODE_DETLIN){
        hdrl_imagelist * list3 = hdrl_imagelist_new();
        hdrl_imagelist_set(list3, hdrl, 0);
        hdrl_imagelist_set(list3, hdrl, 1);
        hdrl_imagelist_set(list3, hdrl, 2);

        hdrl_imagelist * list4[] = {list3, NULL, NULL};

        cr2res_io_save_DETLIN_COEFFS(filename, all, in, parlist, list4, ext1, ext, CR2RES_DETLIN_COEFFS_PROCATG, "debug");

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
static cpl_frame * create_master_flat(char * filename, int nx, int ny, double value, double error)
{
    cpl_frame * out = cpl_frame_new();
    cpl_frame_set_filename(out, filename);
    cpl_frame_set_tag(out, "FLAT");
    cpl_frame_set_group(out, CPL_FRAME_GROUP_CALIB);

    hdrl_image * hdrl = create_hdrl(nx, ny, value, error);
    save_hdrl(filename, hdrl, MODE_FLAT, 0);
    hdrl_image_delete(hdrl);

    return out;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Create a Dark Current fits file and cpl_frame
 */
/*----------------------------------------------------------------------------*/
static cpl_frame * create_master_dark(char * filename, int nx, int ny, double value, double error, double dit)
{

    cpl_frame * out = cpl_frame_new();
    cpl_frame_set_filename(out, filename);
    cpl_frame_set_tag(out, "DARK");
    cpl_frame_set_group(out, CPL_FRAME_GROUP_CALIB);

    hdrl_image * hdrl = create_hdrl(nx, ny, value, error);
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
static cpl_frame * create_detlin(char * filename, int nx, int ny)
{
    cpl_frame * out = cpl_frame_new();
    cpl_frame_set_filename(out, filename);
    cpl_frame_set_tag(out, "DETLIN");
    cpl_frame_set_group(out, CPL_FRAME_GROUP_CALIB);

    hdrl_image * hdrl = hdrl_image_new(nx, ny);
    hdrl_value one;
    one.data = 10;
    one.error = 10;
    hdrl_image_add_scalar(hdrl, one);
    
    save_hdrl(filename, hdrl, MODE_DETLIN, 0);
    // hdrl_image_delete(hdrl);

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

    cpl_frame * flat = create_master_flat("master_flat.fits", nx, ny, 1, 0);
    cpl_frame * dark = create_master_dark("master_dark.fits", nx, ny, 10, 1, 10);
    cpl_frame * bpm = create_bpm("bpm.fits", nx, ny, 0);
    cpl_frame * detlin = create_detlin("master_detlin.fits", nx, ny);
    // cpl_frame * detlin = NULL;
    double dit = 10;

    hdrl_image * out;
    hdrl_image * cmp;

    // NULL input / output
    out = cr2res_calib_image(NULL, chip, 0, NULL, NULL, NULL, NULL, dit);
    cpl_test_null(out);

    out = cr2res_calib_image(in, 0, 0, NULL, NULL, NULL, NULL, dit);
    cpl_test_null(out);

    out = cr2res_calib_image(in, CR2RES_NB_DETECTORS + 1, 0, NULL, NULL, NULL, NULL, dit);
    cpl_test_null(out);

    // No correction
    out = cr2res_calib_image(in, chip, 0, NULL, NULL, NULL, NULL, dit);
    cpl_test_nonnull(out);
    cpl_test_image_abs(hdrl_image_get_image(in), hdrl_image_get_image(out), DBL_EPSILON);
    cpl_test_image_abs(hdrl_image_get_error(in), hdrl_image_get_error(out), DBL_EPSILON);
    hdrl_image_delete(out);

    // Test all
    // out = cr2res_calib_image(in, chip, cosmics_corr, flat, dark, bpm, detlin, dit);

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
    out = cr2res_calib_image(in, chip, 0, NULL, NULL, NULL, NULL, dit);
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
    cpl_frame * flat;
    int chip = 1;
    int cosmics_corr = 1;
    double dit = 10;
    cpl_error_code error;

    in = create_hdrl(nx, ny, img_value, img_error);

    // Case 1: Flat is just 1 and no error, i.e. no change
    flat_value = 1;
    flat_error = 0;
    flat = create_master_flat("master_flat.fits", nx, ny, flat_value, flat_error);
    out = cr2res_calib_image(in, chip, 0, flat, NULL, NULL, NULL, dit);
    
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
    flat = create_master_flat("master_flat.fits", nx, ny, flat_value, flat_error);
    out = cr2res_calib_image(in, chip, 0, flat, NULL, NULL, NULL, dit);
    
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
    flat = create_master_flat("master_flat.fits", nx, ny, flat_value, flat_error);
    out = cr2res_calib_image(in, chip, 0, flat, NULL, NULL, NULL, dit);
    
    cpl_test_eq(nx * ny, hdrl_image_count_rejected(out));

    hdrl_image_delete(out);
    cpl_frame_delete(flat);

    // Case 4: Flat file does not exist
    // flat = cpl_frame_new();
    // cpl_frame_set_filename(flat, "tobeornottobe.fits");
    // cpl_frame_set_tag(flat, "FLAT");
    // cpl_frame_set_group(flat, CPL_FRAME_GROUP_CALIB);

    // out = cr2res_calib_image(in, chip, 0, flat, NULL, NULL, NULL, dit);
    // cpl_test_error(CPL_ERROR_FILE_NOT_FOUND);

    // Case 5: image is in a wrong group
    flat = create_master_dark("master_flat.fits", nx, ny, 1, 0, 10);
    out = cr2res_calib_image(in, chip, 0, flat, NULL, NULL, NULL, dit);
    cpl_test_null(out);
    cpl_frame_delete(flat);

    // Case 6: No Filename set
    flat = cpl_frame_new();
    out = cr2res_calib_image(in, chip, 0, flat, NULL, NULL, NULL, dit);
    cpl_test_error(CPL_ERROR_DATA_NOT_FOUND);
    cpl_test_null(out);
    cpl_frame_delete(flat);

    // Case 7: DataFile is empty, i.e. only header
    flat = cpl_frame_new();
    cpl_frame_set_filename(flat, "empty.fits");
    out = cr2res_calib_image(in, chip, 0, flat, NULL, NULL, NULL, dit);
    cpl_test_null(out);
    cpl_frame_delete(flat);


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
    cpl_frame * dark;
    int chip = 1;
    int cosmics_corr = 1;
    double dit = 10;

    in = create_hdrl(nx, ny, img_value, img_error);

    // Case 1: Dark is 0, dit is the same as the image, i.e. no correction
    dark_value = 0;
    dark_error = 0;
    dark_dit = dit;
    dark = create_master_dark("master_dark.fits", nx, ny, dark_value, dark_error, dark_dit);
    out = cr2res_calib_image(in, chip, 0, NULL, dark, NULL, NULL, dit);
    
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
    dark = create_master_dark("master_dark.fits", nx, ny, dark_value, dark_error, dark_dit);
    out = cr2res_calib_image(in, chip, 0, NULL, dark, NULL, NULL, dit);
    
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
    dark = create_master_dark("master_dark.fits", nx, ny, dark_value, dark_error, dark_dit);
    out = cr2res_calib_image(in, chip, 0, NULL, dark, NULL, NULL, dit);
    
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
    dark = create_master_dark("master_dark.fits", nx, ny, dark_value, dark_error, dark_dit);
    out = cr2res_calib_image(in, chip, 0, NULL, dark, NULL, NULL, dit);
    
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
    dark = create_master_dark("master_dark.fits", nx, ny, dark_value, dark_error, dark_dit);
    out = cr2res_calib_image(in, chip, 0, NULL, dark, NULL, NULL, dit);
    
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
    dark = create_master_dark("master_dark.fits", nx, ny, dark_value, dark_error, dark_dit);
    out = cr2res_calib_image(in, chip, 0, NULL, dark, NULL, NULL, dit);
    
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
    dark = create_master_dark("master_dark.fits", nx, ny, dark_value, dark_error, dark_dit);
    out = cr2res_calib_image(in, chip, 0, NULL, dark, NULL, NULL, dit);
    
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
    dark = create_master_dark("master_dark.fits", nx, ny, dark_value, dark_error, dark_dit);
    out = cr2res_calib_image(in, chip, 0, NULL, dark, NULL, NULL, dit);
    
    out_value = img_value - dit/dark_dit * dark_value;
    out_error = sqrt(pow(img_error, 2) + pow(dit/dark_dit * dark_error, 2));
    cmp = create_hdrl(nx, ny, out_value, out_error);

    cpl_test_image_abs(hdrl_image_get_image(cmp), hdrl_image_get_image(out), DBL_EPSILON);
    cpl_test_image_abs(hdrl_image_get_error(cmp), hdrl_image_get_error(out), DBL_EPSILON);

    hdrl_image_delete(out);
    hdrl_image_delete(cmp);
    cpl_frame_delete(dark);

    // Case 9: Dark file does not exist
    // dark = cpl_frame_new();
    // cpl_frame_set_filename(dark, "tobeornottobe.fits");
    // cpl_frame_set_tag(dark, "DARK");
    // cpl_frame_set_group(dark, CPL_FRAME_GROUP_CALIB);

    // out = cr2res_calib_image(in, chip, 0, NULL, dark, NULL, NULL, dit);
    // cpl_test_error(CPL_ERROR_FILE_NOT_FOUND);

    // Case 10: image is in a wrong group
    dark = create_master_flat("master_dark.fits", nx, ny, 1, 0);
    out = cr2res_calib_image(in, chip, 0, NULL, dark, NULL, NULL, dit);
    cpl_test_null(out);
    cpl_frame_delete(dark);

    // Case 11: No Filename set
    dark = cpl_frame_new();
    out = cr2res_calib_image(in, chip, 0, NULL, dark, NULL, NULL, dit);
    cpl_test_error(CPL_ERROR_DATA_NOT_FOUND);
    cpl_test_null(out);
    cpl_frame_delete(dark);

    // Case 12: DataFile is empty, i.e. only header
    dark = cpl_frame_new();
    cpl_frame_set_filename(dark, "empty.fits");
    out = cr2res_calib_image(in, chip, 0, NULL, dark, NULL, NULL, dit);
    cpl_test_null(out);
    cpl_frame_delete(dark);

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

    // Case 1: BPM is all 0, i.e. all good pixels
    bpm = create_bpm("master_bpm.fits", nx, ny, 0);
    out = cr2res_calib_image(in, chip, 0, NULL, NULL, bpm, NULL, dit);
    cpl_test_image_abs(hdrl_image_get_image(out), hdrl_image_get_image(in), DBL_EPSILON);
    cpl_test_image_abs(hdrl_image_get_error(out), hdrl_image_get_error(in), DBL_EPSILON);
    cpl_frame_delete(bpm);
    hdrl_image_delete(out);

    // Case 2: BPM is all 1, i.e all pixels are bad
    bpm = create_bpm("master_bpm.fits", nx, ny, 1);
    out = cr2res_calib_image(in, chip, 0, NULL, NULL, bpm, NULL, dit);
    cpl_test_null(out);
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

    out = cr2res_calib_image(in, chip, 0, NULL, NULL, bpm, NULL, dit);

    cmp = hdrl_image_duplicate(in);
    value.data = img_value;
    value.error = 0; // the value is just the mean of the surrounding pixels 
    hdrl_image_set_pixel(cmp, 2, 2, value);

    cpl_test_image_abs(hdrl_image_get_image(out), hdrl_image_get_image(cmp), DBL_EPSILON);
    // cpl_test_image_abs(hdrl_image_get_error(out), hdrl_image_get_error(cmp), DBL_EPSILON);

    cpl_frame_delete(bpm);
    hdrl_image_delete(out);
    hdrl_image_delete(cmp);

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

    out = cr2res_calib_image(in, chip, 0, NULL, NULL, bpm, NULL, dit);

    cmp = hdrl_image_duplicate(in);
    value.data = img_value;
    value.error = 0; // the value is just the mean of the surrounding pixels 
    hdrl_image_set_pixel(cmp, 2, 2, value);

    cpl_test_image_abs(hdrl_image_get_image(out), hdrl_image_get_image(cmp), DBL_EPSILON);
    // cpl_test_image_abs(hdrl_image_get_error(out), hdrl_image_get_error(cmp), DBL_EPSILON);

    cpl_frame_delete(bpm);
    hdrl_image_delete(out);
    hdrl_image_delete(cmp);

    // Case 5: BPM file does not exist
    // bpm = cpl_frame_new();
    // cpl_frame_set_filename(bpm, "tobeornottobe.fits");
    // cpl_frame_set_tag(bpm, "BPM");
    // cpl_frame_set_group(bpm, CPL_FRAME_GROUP_CALIB);

    // out = cr2res_calib_image(in, chip, 0, NULL, NULL, bpm, NULL, dit);
    // cpl_test_error(CPL_ERROR_FILE_NOT_FOUND);

    // Case 6: image is in a wrong group
    bpm = create_master_flat("master_bpm.fits", nx, ny, 1, 0);
    out = cr2res_calib_image(in, chip, 0, NULL, NULL, bpm, NULL, dit);
    cpl_test_null(out);
    cpl_frame_delete(bpm);

    // Case 7: No Filename set
    bpm = cpl_frame_new();
    out = cr2res_calib_image(in, chip, 0, NULL, NULL, bpm, NULL, dit);
    cpl_test_error(CPL_ERROR_DATA_NOT_FOUND);
    cpl_test_null(out);
    cpl_frame_delete(bpm);

    // Case 8: DataFile is empty, i.e. only header
    bpm = cpl_frame_new();
    cpl_frame_set_filename(bpm, "empty.fits");
    out = cr2res_calib_image(in, chip, 0, NULL, NULL, bpm, NULL, dit);
    cpl_test_null(out);
    cpl_frame_delete(bpm);

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

    return cpl_test_end(0);
}
