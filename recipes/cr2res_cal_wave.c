
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

#include <locale.h>
#include <string.h>

#include <cpl.h>
#include <math.h>
#include "hdrl.h"

#include "cr2res_utils.h"
#include "cr2res_calib.h"
#include "cr2res_pfits.h"
#include "cr2res_dfs.h"
#include "cr2res_io.h"
#include "cr2res_extract.h"
#include "cr2res_trace.h"
#include "cr2res_wave.h"

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

#define RECIPE_STRING "cr2res_cal_wave"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_cal_wave_reduce(
        const cpl_frameset  *   rawframes,
        const cpl_frame     *   detlin_frame,
        const cpl_frame     *   master_dark_frame,
        const cpl_frame     *   master_flat_frame,
        const cpl_frame     *   bpm_frame,
        const cpl_frame     *   trace_wave_frame,
        const cpl_frame     *   lines_frame,
        int                     reduce_det,
        int                     reduce_order,
        int                     reduce_trace,
        int                     ext_height,
        int                     ext_swath_width,
        int                     ext_oversample,
        double                  ext_smooth_slit,
        cr2res_wavecal_type     wavecal_type,
        int                     wl_degree,
        double                  wl_start,
        double                  wl_end,
        double                  wl_err_start,
        double                  wl_err_end,
        double                  wl_shift,
        int                     display,
        cpl_table           **  out_trace_wave,
        hdrl_image          **  out_wave_map,
        cpl_propertylist    **  ext_plist) ;

static int cr2res_cal_wave_create(cpl_plugin *);
static int cr2res_cal_wave_exec(cpl_plugin *);
static int cr2res_cal_wave_destroy(cpl_plugin *);
static int cr2res_cal_wave(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_cal_wave_description[] =
"CRIRES+ wavelength calibration recipe\n"
"The files listed in the Set Of Frames (sof-file) must be tagged:\n"
"raw-file.fits " CR2RES_WAVE_RAW "\n"
"trace_wave.fits " CR2RES_FLAT_TRACE_WAVE_PROCATG "\n"
"             or " CR2RES_FLAT_TRACE_WAVE_MERGED_PROCATG "\n"
"             or " CR2RES_UTIL_TRACE_WAVE_PROCATG "\n"
"             or " CR2RES_UTIL_WAVE_TRACE_WAVE_PROCATG "\n"
"             or " CR2RES_CAL_WAVE_TRACE_WAVE_PROCATG "\n"
"             or " CR2RES_UTIL_SLIT_CURV_TRACE_WAVE_PROCATG "\n"
"lines.fits " CR2RES_EMISSION_LINES_PROCATG " (optional) \n"
"detlin_coeffs.fits " CR2RES_DETLIN_COEFFS_PROCATG " (optional) \n"
"master_dark.fits " CR2RES_MASTER_DARK_PROCATG " (optional) \n"
"master_flat.fits " CR2RES_FLAT_MASTER_FLAT_PROCATG " (optional) \n"
"bpm.fits " CR2RES_FLAT_BPM_PROCATG " (optional) \n"
"      or " CR2RES_DETLIN_BPM_PROCATG "\n"
"      or " CR2RES_DARK_BPM_PROCATG "\n"
"      or " CR2RES_UTIL_BPM_SPLIT_PROCATG "\n"
"\n"
"Four different methods are offered by the utility.\n"
"   XCORR:  Cross Correlation with a emission lines catalog (default)\n"
"   LINE1D: Line identification and fitting for each 1D spectra\n"
"   LINE2D: Line identification and fitting for all 1D spectra at once\n"
"   ETALON: Does not require any static calibration file.\n"
"\n"
"The recipe produces the following products:\n"
"cr2res_cal_wave_trace.fits " CR2RES_CAL_WAVE_TRACE_WAVE_PROCATG "\n"
"cr2res_cal_wave_map.fits " CR2RES_CAL_WAVE_MAP_PROCATG "\n"
"\n";

/*-----------------------------------------------------------------------------
                                Function code
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Build the list of available plugins, for this module.
  @param    list    the plugin list
  @return   0 if everything is ok, 1 otherwise
  @note     Only this function is exported

  Create the recipe instance and make it available to the application using the
  interface.
 */
/*----------------------------------------------------------------------------*/
int cpl_plugin_get_info(cpl_pluginlist * list)
{
    cpl_recipe  *   recipe = cpl_calloc(1, sizeof *recipe );
    cpl_plugin  *   plugin = &recipe->interface;

    if (cpl_plugin_init(plugin,
                    CPL_PLUGIN_API,
                    CR2RES_BINARY_VERSION,
                    CPL_PLUGIN_TYPE_RECIPE,
                    RECIPE_STRING,
                    "Wavelength Calibration",
                    cr2res_cal_wave_description,
                    "Ansgar Wehrhahn, Thomas Marquart, Yves Jung",
                    PACKAGE_BUGREPORT,
                    cr2res_get_license(),
                    cr2res_cal_wave_create,
                    cr2res_cal_wave_exec,
                    cr2res_cal_wave_destroy)) {
        cpl_msg_error(cpl_func, "Plugin initialization failed");
        (void)cpl_error_set_where(cpl_func);
        return 1;
    }

    if (cpl_pluginlist_append(list, plugin)) {
        cpl_msg_error(cpl_func, "Error adding plugin to list");
        (void)cpl_error_set_where(cpl_func);
        return 1;
    }

    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Setup the recipe options
  @param    plugin  the plugin
  @return   0 if everything is ok

  Defining the command-line/configuration parameters for the recipe.
 */
/*----------------------------------------------------------------------------*/
static int cr2res_cal_wave_create(cpl_plugin * plugin)
{
    cpl_recipe    * recipe;
    cpl_parameter * p;

    /* Check that the plugin is part of a valid recipe */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else
        return -1;

    /* Create the parameters list in the cpl_recipe object */
    recipe->parameters = cpl_parameterlist_new();

    /* Fill the parameters list */
    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.detector",
            CPL_TYPE_INT, "Only reduce the specified detector",
            "cr2res.cr2res_cal_wave", 0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "detector");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.order",
            CPL_TYPE_INT, "Only reduce the specified order",
            "cr2res.cr2res_cal_wave", -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "order");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.trace_nb",
            CPL_TYPE_INT, "Only reduce the specified trace number",
            "cr2res.cr2res_cal_wave", -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "trace_nb");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.ext_oversample",
            CPL_TYPE_INT, "factor by which to oversample the extraction",
            "cr2res.cr2res_cal_wave", 5);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ext_oversample");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.ext_swath_width",
            CPL_TYPE_INT, "The swath width", "cr2res.cr2res_cal_wave", 32);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ext_swath_width");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.ext_height",
            CPL_TYPE_INT, "Extraction height",
            "cr2res.cr2res_cal_wave", -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ext_height");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.ext_smooth_slit",
            CPL_TYPE_DOUBLE,
            "Smoothing along the slit (1 for high S/N, 5 for low)",
            "cr2res.cr2res_cal_wave", 1.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ext_smooth_slit");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.wl_method",
            CPL_TYPE_STRING, "Data Type (XCORR / LINE1D / LINE2D / ETALON)",
            "cr2res.cr2res_cal_wave", "XCORR");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "wl_method");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.wl_shift",
            CPL_TYPE_DOUBLE, "Wavelength shift (nm) to apply to the guess",
            "cr2res.cr2res_cal_wave", 0.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "wl_shift");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.wl_est",
            CPL_TYPE_STRING, "Estimated wavelength start and end",
            "cr2res.cr2res_cal_wave", "-1.0, -1.0");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "wl_est");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.wl_err",
            CPL_TYPE_STRING,
            "Estimated wavelength error [start_err, end_err] (in nm)",
            "cr2res.cr2res_cal_wave", "-1.0, -1.0");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "wl_err");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.wl_degree",
            CPL_TYPE_INT, "Wavelegth Polynomial degree",
            "cr2res.cr2res_cal_wave", 3);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "wl_degree");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.log",
            CPL_TYPE_BOOL, "Flag for taking the Log() value of the lines",
            "cr2res.cr2res_cal_wave", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "log");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.display",
            CPL_TYPE_BOOL, "Flag for display",
            "cr2res.cr2res_cal_wave", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "display");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);
    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Execute the plugin instance given by the interface
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_cal_wave_exec(cpl_plugin * plugin)
{
    cpl_recipe  *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1;

    return cr2res_cal_wave(recipe->frames, recipe->parameters);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_cal_wave_destroy(cpl_plugin * plugin)
{
    cpl_recipe *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1 ;

    cpl_parameterlist_delete(recipe->parameters);
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Interpret the command line options and execute the data processing
  @param    frameset   the frames list
  @param    parlist    the parameters list
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_cal_wave(
        cpl_frameset            *   frameset,
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter *   param;
    int                     reduce_det, reduce_order, reduce_trace,
                            ext_oversample, ext_swath_width, ext_height,
                            wl_degree, display, log_flag ;
    double                  ext_smooth_slit, wl_start, wl_end, wl_err_start, 
                            wl_err_end, wl_shift ;
    cr2res_wavecal_type     wavecal_type ;
    const char          *   sval ;
    cpl_frameset        *   rawframes ;
    const cpl_frame     *   detlin_frame ;
    const cpl_frame     *   master_dark_frame ;
    const cpl_frame     *   master_flat_frame ;
    const cpl_frame     *   bpm_frame ;
    const cpl_frame     *   trace_wave_frame ;
    const cpl_frame     *   lines_frame ;
    char                *   out_file;
    cpl_table           *   out_trace_wave[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   out_wave_map[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   ext_plist[CR2RES_NB_DETECTORS] ;
    int                     det_nr, order, i ;

    /* Needed for sscanf() */
    setlocale(LC_NUMERIC, "C");

    /* Initialise */
    wavecal_type = CR2RES_UNSPECIFIED ;
    wl_start = wl_end = wl_err_start = wl_err_end = -1.0 ;
    wl_shift = 0.0 ;

    /* RETRIEVE INPUT PARAMETERS */
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_wave.detector");
    reduce_det = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_wave.order");
    reduce_order = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_wave.trace_nb");
    reduce_trace = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
           "cr2res.cr2res_util_extract.ext_oversample");
    ext_oversample = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_extract.ext_swath_width");
    ext_swath_width = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_extract.ext_height");
    ext_height = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_extract.ext_smooth_slit");
    ext_smooth_slit = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_wave.wl_method");
    sval = cpl_parameter_get_string(param) ;
    if (!strcmp(sval, "XCORR"))         wavecal_type = CR2RES_XCORR ;
    else if (!strcmp(sval, "LINE1D"))   wavecal_type = CR2RES_LINE1D ;
    else if (!strcmp(sval, "LINE2D"))   wavecal_type = CR2RES_LINE2D ;
    else if (!strcmp(sval, "ETALON"))   wavecal_type = CR2RES_ETALON ;
    else {
        cpl_msg_error(__func__, "Invalid Data Type specified");
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1;
    }
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_wave.wl_shift");
    wl_shift = cpl_parameter_get_double(param) ;
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_wave.wl_est");
    sval = cpl_parameter_get_string(param) ;
    if (sscanf(sval, "%lg,%lg", &wl_start, &wl_end) != 2) {
        return -1 ;
    }
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_wave.wl_err");
    sval = cpl_parameter_get_string(param) ;
    if (sscanf(sval, "%lg,%lg", &wl_err_start, &wl_err_end) != 2) {
        return -1 ;
    }
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_wave.wl_degree");
    wl_degree = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_wave.log");
    log_flag = cpl_parameter_get_bool(param) ;
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_wave.display");
    display = cpl_parameter_get_bool(param) ;

    /* Check Parameters */
    if (reduce_order > -1 && wavecal_type == CR2RES_LINE2D) {
        cpl_msg_error(__func__, "Limiting to one order with LINE2D impossible");
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }
    if (wl_degree < 0) {
        cpl_msg_error(__func__, "The degree needs to be >= 0");
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Identify the RAW and CALIB frames in the input frameset */
    if (cr2res_dfs_set_groups(frameset) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot identify RAW and CALIB frames") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Retrieve calibration data */
    trace_wave_frame = cr2res_io_find_TRACE_WAVE(frameset) ;
    if (trace_wave_frame == NULL) {
        cpl_msg_error(__func__, "Could not find TRACE_WAVE frame") ;
        return -1 ;
    }
    detlin_frame = cpl_frameset_find_const(frameset,
            CR2RES_DETLIN_COEFFS_PROCATG);
    master_dark_frame = cpl_frameset_find_const(frameset,
            CR2RES_MASTER_DARK_PROCATG) ;
    master_flat_frame = cpl_frameset_find_const(frameset,
            CR2RES_FLAT_MASTER_FLAT_PROCATG) ;
    bpm_frame = cr2res_io_find_BPM(frameset) ;
    lines_frame = cpl_frameset_find_const(frameset,
            CR2RES_EMISSION_LINES_PROCATG) ;

    /* Get the RAW Frames */
    rawframes = cr2res_extract_frameset(frameset, CR2RES_WAVE_RAW) ;
    if (rawframes == NULL) {
        cpl_msg_error(__func__, "Could not find RAW frames") ;
        return -1 ;
    }

    /* Loop over the detectors */
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {

        /* Initialise */
        out_trace_wave[det_nr-1] = NULL ;
        out_wave_map[det_nr-1] = NULL ;
        ext_plist[det_nr-1] = NULL ;

        /* Compute only one detector */
        if (reduce_det != 0 && det_nr != reduce_det) continue ;

        cpl_msg_info(__func__, "Process detector number %d", det_nr) ;
        cpl_msg_indent_more() ;

        /* Call the reduction function */
        if (cr2res_cal_wave_reduce(rawframes, detlin_frame,
                    master_dark_frame, master_flat_frame, bpm_frame,
                    trace_wave_frame, lines_frame, det_nr, reduce_order,
                    reduce_trace,  ext_height, ext_swath_width,
                    ext_oversample, ext_smooth_slit, wavecal_type, wl_degree, 
                    wl_start, wl_end, wl_err_start, wl_err_end, wl_shift, 
                    display,
                    &(out_trace_wave[det_nr-1]),
                    &(out_wave_map[det_nr-1]),
                    &(ext_plist[det_nr-1])) == -1) {
            cpl_msg_warning(__func__, "Failed to reduce detector %d", det_nr);
        }
        cpl_msg_indent_less() ;
    }

    /* Ѕave Products */
    out_file = cpl_sprintf("%s_trace.fits", RECIPE_STRING) ;
    cr2res_io_save_TRACE_WAVE(out_file, frameset, rawframes, parlist, 
            out_trace_wave, NULL, ext_plist, 
            CR2RES_CAL_WAVE_TRACE_WAVE_PROCATG, RECIPE_STRING) ;
    cpl_free(out_file);

    out_file = cpl_sprintf("%s_map.fits", RECIPE_STRING) ;
    cr2res_io_save_WAVE_MAP(out_file, frameset, rawframes, parlist, 
            out_wave_map, NULL, ext_plist, 
            CR2RES_CAL_WAVE_MAP_PROCATG, RECIPE_STRING) ;
    cpl_free(out_file);

    /* Free and return */
    cpl_frameset_delete(rawframes) ;
    for (i=0 ; i<CR2RES_NB_DETECTORS ; i++) {
        if (ext_plist[i] != NULL)
            cpl_propertylist_delete(ext_plist[i]) ;
        if (out_trace_wave[i] != NULL)
            cpl_table_delete(out_trace_wave[i]) ;
        if (out_wave_map[i] != NULL) {
            hdrl_image_delete(out_wave_map[i]) ;
        }
    }
    return (int)cpl_error_get_code();
}
 
/*----------------------------------------------------------------------------*/
/**
  @brief  
  @param 
  @return  
 */
/*----------------------------------------------------------------------------*/
static int cr2res_cal_wave_reduce(
        const cpl_frameset  *   rawframes,
        const cpl_frame     *   detlin_frame,
        const cpl_frame     *   master_dark_frame,
        const cpl_frame     *   master_flat_frame,
        const cpl_frame     *   bpm_frame,
        const cpl_frame     *   trace_wave_frame,
        const cpl_frame     *   lines_frame,
        int                     reduce_det,
        int                     reduce_order,
        int                     reduce_trace,
        int                     ext_height, 
        int                     ext_swath_width,
        int                     ext_oversample, 
        double                  ext_smooth_slit,
        cr2res_wavecal_type     wavecal_type,
        int                     wl_degree,
        double                  wl_start,
        double                  wl_end,
        double                  wl_err_start,
        double                  wl_err_end,
        double                  wl_shift,
        int                     display,
        cpl_table           **  out_trace_wave,
        hdrl_image          **  out_wave_map,
        cpl_propertylist    **  ext_plist)
{
    cpl_vector          *   dits ;
    hdrl_imagelist      *   in ;
    hdrl_imagelist      *   in_calib ;
    hdrl_image          *   collapsed ;
    cpl_image           *   contrib ;
    cpl_table           *   tw_in ;
    cpl_table           *   extracted ;
    cpl_table           *   slit_func ;
    hdrl_image          *   model_master ;

    hdrl_image          *   wl_map ;
    cpl_table           *   tw_out ;
    cpl_propertylist    *   plist ;
    
    /* Check Inputs */
    if (rawframes==NULL || trace_wave_frame==NULL || out_trace_wave==NULL ||
            out_wave_map==NULL || ext_plist==NULL) return -1 ;

    /* Load the DITs if necessary */
    if (master_dark_frame != NULL)  dits = cr2res_io_read_dits(rawframes) ;
    else                            dits = NULL ;
    if (cpl_msg_get_level() == CPL_MSG_DEBUG && dits != NULL)
        cpl_vector_dump(dits, stdout) ;

    /* Load image list */
    if ((in = cr2res_io_load_image_list_from_set(rawframes,
                    reduce_det)) == NULL) {
        cpl_msg_error(__func__, "Cannot load images") ;
        if (dits != NULL) cpl_vector_delete(dits) ;
        return -1 ;
    }
    if (hdrl_imagelist_get_size(in) != cpl_frameset_get_size(rawframes)) {
        cpl_msg_error(__func__, "Inconsistent number of loaded images") ;
        if (dits != NULL) cpl_vector_delete(dits) ;
        hdrl_imagelist_delete(in) ;
        return -1 ;
    }

    /* Calibrate the images */
    if ((in_calib = cr2res_calib_imagelist(in, reduce_det, 0, master_flat_frame,
            master_dark_frame, bpm_frame, detlin_frame, dits)) == NULL) {
        cpl_msg_error(__func__, "Failed to apply the calibrations") ;
        if (dits != NULL) cpl_vector_delete(dits) ;
        hdrl_imagelist_delete(in) ;
        return -1 ;
    }
    hdrl_imagelist_delete(in) ;
    if (dits != NULL) cpl_vector_delete(dits) ;

    /* Collapse */
    cpl_msg_info(__func__, "Collapse the input frames") ;
    cpl_msg_indent_more() ;
    if (hdrl_imagelist_collapse_mean(in_calib, &collapsed, &contrib) !=
            CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Failed to Collapse") ;
        hdrl_imagelist_delete(in_calib) ;
        cpl_msg_indent_less() ;
        return -1 ;
    }
    cpl_image_delete(contrib) ;
    hdrl_imagelist_delete(in_calib) ;
    cpl_msg_indent_less() ;

    /* Load the trace wave */
    cpl_msg_info(__func__, "Load the TRACE WAVE") ;
    if ((tw_in = cr2res_io_load_TRACE_WAVE(cpl_frame_get_filename(
                        trace_wave_frame), reduce_det)) == NULL) {
        cpl_msg_error(__func__, "Failed to Load the traces file") ;
        hdrl_image_delete(collapsed) ;
        return -1 ;
    }

    /* Execute the extraction */
    cpl_msg_info(__func__, "Spectra Extraction") ;
    if (cr2res_extract_traces(collapsed, tw_in, reduce_order, reduce_trace,
                CR2RES_EXTR_OPT_CURV, ext_height, ext_swath_width,
                ext_oversample, ext_smooth_slit,
                &extracted, &slit_func, &model_master) == -1) {
        cpl_msg_error(__func__, "Failed to extract");
        hdrl_image_delete(collapsed) ;
        cpl_table_delete(tw_in) ;
        return -1 ;
    }
    cpl_table_delete(slit_func) ;
    hdrl_image_delete(model_master) ;
    hdrl_image_delete(collapsed);
    
    
    
    cpl_table_delete(tw_in) ;
    cpl_table_delete(extracted) ;



    /* Return */
    *out_wave_map = wl_map;
    *out_trace_wave = tw_out ;
    *ext_plist = plist ;
    return 0 ;
}
