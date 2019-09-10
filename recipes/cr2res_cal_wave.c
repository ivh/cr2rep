
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
        cr2res_collapse         collapse,
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
        int                     log_flag,
        int                     propagate_flag,
        int                     display,
        double                  display_wmin,
        double                  display_wmax,
        cpl_table           **  out_trace_wave,
        cpl_table           **  lines_diagnostics,
        cpl_table           **  out_extracted,
        hdrl_image          **  out_wave_map,
        cpl_propertylist    **  ext_plist) ;
static int cr2res_cal_wave_create(cpl_plugin *);
static int cr2res_cal_wave_exec(cpl_plugin *);
static int cr2res_cal_wave_destroy(cpl_plugin *);
static int cr2res_cal_wave(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_cal_wave_description[] = "\
Spectrum Extraction and Wavelength Calibration                          \n\
  This recipe performs the extraction of the various orders along       \n\
  the provided traces, and the wavelength calibration of these          \n\
  extracted spectra.                                                    \n\
  It can support different methods (--wl_method parameter):             \n\
    XCORR:  Cross Correlation with a emission lines catalog (default)   \n\
    LINE1D: Line identification and fitting for each 1D spectra         \n\
    LINE2D: Line identification and fitting for all 1D spectra at once  \n\
    ETALON: Does not require any static calibration filer               \n\
    AUTO:   Guess the Method from the input file header                 \n\
                                                                        \n\
  Inputs                                                                \n\
    raw.fits " CR2RES_WAVE_RAW " [1 to n]                               \n\
    trace.fits " CR2RES_CAL_FLAT_TW_PROCATG " [1]                       \n\
            or " CR2RES_CAL_FLAT_TW_MERGED_PROCATG "                    \n\
            or " CR2RES_UTIL_TRACE_TW_PROCATG "                         \n\
            or " CR2RES_UTIL_WAVE_TW_PROCATG "                          \n\
            or " CR2RES_CAL_WAVE_TW_PROCATG "                           \n\
            or " CR2RES_UTIL_SLIT_CURV_TW_PROCATG "                     \n\
    detlin.fits " CR2RES_CAL_DETLIN_COEFFS_PROCATG " [0 to 1]           \n\
    bpm.fits " CR2RES_CAL_DARK_BPM_PROCATG " [0 to 1]                   \n\
          or " CR2RES_CAL_FLAT_BPM_PROCATG "                            \n\
          or " CR2RES_CAL_DETLIN_BPM_PROCATG "                          \n\
          or " CR2RES_UTIL_BPM_SPLIT_PROCATG "                          \n\
    master_dark.fits " CR2RES_CAL_DARK_MASTER_PROCATG " [0 to 1]        \n\
    master_flat.fits " CR2RES_CAL_FLAT_MASTER_PROCATG " [0 to 1]        \n\
    lines.fits " CR2RES_EMISSION_LINES_PROCATG " [0 to 1]               \n\
                                                                        \n\
  Outputs                                                               \n\
    cr2res_cal_wave_tw.fits " CR2RES_CAL_WAVE_TW_PROCATG"               \n\
    cr2res_cal_wave_wave_map.fits " CR2RES_CAL_WAVE_MAP_PROCATG"        \n\
    cr2res_cal_wave_lines_diagnostics.fits "
    CR2RES_CAL_WAVE_LINES_DIAGNOSTICS_PROCATG"\n\
    cr2res_cal_wave_extracted.fits " CR2RES_CAL_WAVE_EXTRACT_1D_PROCATG"\n\
                                                                        \n\
  Algorithm                                                             \n\
    loop on detectors d:                                                \n\
      Call cr2res_cal_wave_reduce()                                     \n\
        -> out_trace_wave(d)                                            \n\
        -> lines_diagnostics(d)                                         \n\
        -> out_extracted(d)                                             \n\
        -> out_wave_map(d)                                              \n\
    Save out_trace_wave                                                 \n\
    Save lines_diagnostics                                              \n\
    Save out_extracted                                                  \n\
    Save out_wave_map                                                   \n\
                                                                        \n\
    cr2res_cal_wave_reduce():                                           \n\
        Load the raw image list                                         \n\
        Apply the calibrations to the image list                        \n\
        Collapse the image list                                         \n\
        Extract along the traces from the collapsed image               \n\
        Compute the Wavelength with cr2res_wave_apply()                 \n\
         -> out_trace_wave                                              \n\
         -> lines_diagnostics                                           \n\
         -> out_extracted                                               \n\
        Compute the Wavelength map                                      \n\
         -> out_wave_map                                                \n\
                                                                        \n\
    cr2res_wave_apply()                                                 \n\
      loop on the traces t:                                             \n\
        Get the spectrum                                                \n\
        Get the Initial guess                                           \n\
        Switch on the required method:                                  \n\
          CR2RES_LINE2D: cr2res_wave_2d()                               \n\
          CR2RES_LINE1D: cr2res_wave_line_fitting()                     \n\
          CR2RES_ETALON: cr2res_wave_etalon()                           \n\
          CR2RES_XCORR:  cr2res_wave_xcorr()                            \n\
                                                                        \n\
  Library Functions uѕed                                                \n\
    cr2res_io_find_TRACE_WAVE()                                         \n\
    cr2res_io_find_BPM()                                                \n\
    cr2res_io_read_dits()                                               \n\
    cr2res_io_load_image_list_from_set()                                \n\
    cr2res_calib_imagelist()                                            \n\
    cr2res_io_load_TRACE_WAVE()                                         \n\
    cr2res_extract_traces()                                             \n\
    cr2res_wave_apply()                                                 \n\
    cr2res_extract_EXTRACT1D_get_spectrum()                             \n\
    cr2res_wave_estimate_compute()                                      \n\
    cr2res_wave_2d()                                                    \n\
    cr2res_wave_line_fitting()                                          \n\
    cr2res_wave_etalon()                                                \n\
    cr2res_wave_xcorr()                                                 \n\
    cr2res_wave_gen_wave_map()                                          \n\
    cr2res_io_save_TRACE_WAVE()                                         \n\
    cr2res_io_save_WAVE_MAP()                                           \n\
    cr2res_io_save_LINES_DIAGNOSTICS()                                  \n\
    cr2res_io_save_EXTRACT_1D()                                         \n\
" ;

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

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.collapse_method",
            CPL_TYPE_STRING, "Collapse the input images (MEAN or MEDIAN)",
            "cr2res.cr2res_cal_wave", "MEDIAN");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "collapse_method");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.ext_oversample",
            CPL_TYPE_INT, "factor by which to oversample the extraction",
            "cr2res.cr2res_cal_wave", 3);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ext_oversample");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.ext_swath_width",
            CPL_TYPE_INT, "The swath width", "cr2res.cr2res_cal_wave", 90);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ext_swath_width");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.ext_height",
            CPL_TYPE_INT, "Extraction height",
            "cr2res.cr2res_cal_wave", 25);
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
            CPL_TYPE_STRING, 
            "Wavelength Method (AUTO / XCORR / LINE1D / LINE2D / ETALON)",
            "cr2res.cr2res_cal_wave", "AUTO");
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

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.propagate",
            CPL_TYPE_BOOL, "Flag for using the input WL when no computation",
            "cr2res.cr2res_cal_wave", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "propagate");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.display",
            CPL_TYPE_BOOL, "Flag for display",
            "cr2res.cr2res_cal_wave", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "display");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.display_range",
            CPL_TYPE_STRING,
            "Wavelength range to display [start, end] (in nm)",
            "cr2res.cr2res_cal_wave", "-1.0, -1.0");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "display_range");
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
                            wl_degree, display, log_flag, propagate_flag ;
    double                  ext_smooth_slit, wl_start, wl_end, wl_err_start, 
                            wl_err_end, wl_shift, display_wmin, display_wmax ;
    cr2res_collapse         collapse ;
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
    cpl_table           *   lines_diagnostics[CR2RES_NB_DETECTORS] ;
    cpl_table           *   out_extracted[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   out_wave_map[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   ext_plist[CR2RES_NB_DETECTORS] ;
    int                     det_nr, order, i ;

    /* Needed for sscanf() */
    setlocale(LC_NUMERIC, "C");

    /* Initialise */
    wl_start = wl_end = wl_err_start = wl_err_end = -1.0 ;
    wl_shift = 0.0 ;
    collapse = CR2RES_COLLAPSE_UNSPECIFIED ;
    display_wmin = display_wmax = -1.0 ;

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
            "cr2res.cr2res_cal_wave.collapse_method");
    sval = cpl_parameter_get_string(param);
    if (!strcmp(sval, "MEAN"))          collapse = CR2RES_COLLAPSE_MEAN ;
    else if (!strcmp(sval, "MEDIAN"))   collapse = CR2RES_COLLAPSE_MEDIAN ;
    if (collapse!=CR2RES_COLLAPSE_MEAN && collapse!=CR2RES_COLLAPSE_MEDIAN) {
        cpl_msg_error(__func__, "Unsupported collapse method") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }
    param = cpl_parameterlist_find_const(parlist,
           "cr2res.cr2res_cal_wave.ext_oversample");
    ext_oversample = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_wave.ext_swath_width");
    ext_swath_width = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_wave.ext_height");
    ext_height = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_wave.ext_smooth_slit");
    ext_smooth_slit = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_wave.wl_method");
    sval = cpl_parameter_get_string(param) ;
    if (!strcmp(sval, "XCORR"))         wavecal_type = CR2RES_XCORR ;
    else if (!strcmp(sval, "LINE1D"))   wavecal_type = CR2RES_LINE1D ;
    else if (!strcmp(sval, "LINE2D"))   wavecal_type = CR2RES_LINE2D ;
    else if (!strcmp(sval, "ETALON"))   wavecal_type = CR2RES_ETALON ;
    else if (!strcmp(sval, "AUTO"))     wavecal_type = CR2RES_UNSPECIFIED ;
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
            "cr2res.cr2res_cal_wave.propagate");
    propagate_flag = cpl_parameter_get_bool(param) ;
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_wave.display");
    display = cpl_parameter_get_bool(param) ;
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_wave.display_range");
    sval = cpl_parameter_get_string(param) ;
    if (sscanf(sval, "%lg,%lg", &display_wmin, &display_wmax) != 2) {
        return -1 ;
    }

    /* Check Parameters */
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
            CR2RES_CAL_DETLIN_COEFFS_PROCATG);
    master_dark_frame = cpl_frameset_find_const(frameset,
            CR2RES_CAL_DARK_MASTER_PROCATG) ;
    master_flat_frame = cpl_frameset_find_const(frameset,
            CR2RES_CAL_FLAT_MASTER_PROCATG) ;
    bpm_frame = cr2res_io_find_BPM(frameset) ;
    lines_frame = cpl_frameset_find_const(frameset,
            CR2RES_EMISSION_LINES_PROCATG) ;

    /* Get the RAW Frames */
    rawframes = cr2res_extract_frameset(frameset, CR2RES_WAVE_RAW) ;
    if (rawframes == NULL) {
        cpl_msg_error(__func__, "Could not find RAW frames") ;
        return -1 ;
    }

    /* Guess the method to be used from the RAW frames header */
    if (wavecal_type == CR2RES_UNSPECIFIED) {
        if ((wavecal_type = cr2res_wave_guess_method(
                        cpl_frameset_get_position(rawframes, 0))) == 
                CR2RES_UNSPECIFIED) {
            cpl_frameset_delete(rawframes) ;
            cpl_msg_error(__func__, "Cannot guess the method") ;
            cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
            return -1 ;
        }
        char * method_str = cr2res_wave_method_print(wavecal_type) ;
        cpl_msg_info(__func__, "Method Automatically Guessed : %s",
                method_str) ;
        cpl_free(method_str) ;
    }
    if (reduce_order > -1 && wavecal_type == CR2RES_LINE2D) {
        cpl_msg_error(__func__, "Limiting to one order with LINE2D impossible");
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Loop over the detectors */
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {

        /* Initialise */
        out_trace_wave[det_nr-1] = NULL ;
        lines_diagnostics[det_nr-1] = NULL ;
        out_extracted[det_nr-1] = NULL ;
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
                    reduce_trace, collapse, ext_height, ext_swath_width,
                    ext_oversample, ext_smooth_slit, wavecal_type, wl_degree, 
                    wl_start, wl_end, wl_err_start, wl_err_end, wl_shift, 
                    log_flag, propagate_flag, display, display_wmin,
                    display_wmax, 
                    &(out_trace_wave[det_nr-1]),
                    &(lines_diagnostics[det_nr-1]),
                    &(out_extracted[det_nr-1]),
                    &(out_wave_map[det_nr-1]),
                    &(ext_plist[det_nr-1])) == -1) {
            cpl_msg_warning(__func__, "Failed to reduce detector %d", det_nr);
        }
        cpl_msg_indent_less() ;
    }

    /* Ѕave Products */
    out_file = cpl_sprintf("%s_tw.fits", RECIPE_STRING) ;
    cr2res_io_save_TRACE_WAVE(out_file, frameset, rawframes, parlist, 
            out_trace_wave, NULL, ext_plist, 
            CR2RES_CAL_WAVE_TW_PROCATG, RECIPE_STRING) ;
    cpl_free(out_file);

    out_file = cpl_sprintf("%s_wave_map.fits", RECIPE_STRING) ;
    cr2res_io_save_WAVE_MAP(out_file, frameset, rawframes, parlist, 
            out_wave_map, NULL, ext_plist, 
            CR2RES_CAL_WAVE_MAP_PROCATG, RECIPE_STRING) ;
    cpl_free(out_file);

    out_file = cpl_sprintf("%s_extracted.fits", RECIPE_STRING) ;
    cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
            out_extracted, NULL, ext_plist, 
            CR2RES_CAL_WAVE_EXTRACT_1D_PROCATG, RECIPE_STRING) ;
    cpl_free(out_file);

	if (wavecal_type == CR2RES_LINE2D || wavecal_type == CR2RES_LINE1D) {
		/* Save the Lines Diagnostics */
		out_file = cpl_sprintf("%s_lines_diagnostics.fits", RECIPE_STRING);
		cr2res_io_save_LINES_DIAGNOSTICS(out_file, frameset, rawframes, parlist,
                lines_diagnostics, NULL, ext_plist,
				CR2RES_CAL_WAVE_LINES_DIAGNOSTICS_PROCATG, RECIPE_STRING) ;
		cpl_free(out_file);
	}

    /* Free and return */
    cpl_frameset_delete(rawframes) ;
    for (i=0 ; i<CR2RES_NB_DETECTORS ; i++) {
        if (ext_plist[i] != NULL)
            cpl_propertylist_delete(ext_plist[i]) ;
        if (out_trace_wave[i] != NULL)
            cpl_table_delete(out_trace_wave[i]) ;
        if (lines_diagnostics[i] != NULL)
            cpl_table_delete(lines_diagnostics[i]) ;
        if (out_extracted[i] != NULL)
            cpl_table_delete(out_extracted[i]) ;
        if (out_wave_map[i] != NULL) {
            hdrl_image_delete(out_wave_map[i]) ;
        }
    }
    return (int)cpl_error_get_code();
}

/*----------------------------------------------------------------------------*/
/**
  @brief Compute the Wavelength for a detector
  @param rawframes          Input raw frames 
  @param detlin_frame       Associated detlin coefficients
  @param master_dark_frame  Associated master dark
  @param master_flat_frame  Associated master flat
  @param bpm_frame          Associated BPM
  @param trace_wave_frame   Trace Wave table
  @param lines_frame        Emission lines frame
  @param reduce_det         The detector to compute
  @param reduce_order       The order to compute (-1 for all)
  @param reduce_trace       The trace to compute (-1 for all)
  @param collapse           CR2RES_COLLAPSE_MEAN or CR2RES_COLLAPSE_MEDIAN
  @param ext_height         Extraction related
  @param ext_swath_width    Extraction related
  @param ext_oversample     Extraction related
  @param ext_smooth_slit    Extraction related
  @param wavecal_type       CR2RES_XCORR/LINE1D/LINE2D/ETALON
  @param wl_start           WL estimate of the first pixel
  @param wl_end             WL estimate of the last pixel
  @param wl_err_start       WL error of wl_start
  @param wl_err_end         WL error of wl_end
  @param wl_shift           wavelength shift to apply
  @param log_flag           Flag to apply a log() to the lines intensities
  @param propagate_flag     Flag to copy the input WL to the output when they 
                            are not computed
  @param display            Flag to enable display functionalities
  @param display_wmin       Minimum Wavelength to  display
  @param display_wmax       Maximum Wavelength to  display
  @param out_trace_wave     [out] trace wave table
  @param lines_diagnostics  [out] lines diagnostics table
  @param out_extracted      [out] extracted table with updated WL
  @param out_wave_map       [out] Wave map
  @param ext_plist          [out] the header for saving the products
  @return   0 if ok, -1 otherwise
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
        cr2res_collapse         collapse,
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
        int                     log_flag,
        int                     propagate_flag,
        int                     display,
        double                  display_wmin,
        double                  display_wmax,
        cpl_table           **  out_trace_wave,
        cpl_table           **  lines_diagnostics,
        cpl_table           **  out_extracted,
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
    cpl_table           *   lines_diagnostics_out ;
    cpl_table           *   extracted_out ;
    cpl_propertylist    *   plist ;
    cpl_propertylist    *   qcs_plist ;
    const char          *   first_file ;
    int                     ext_nr ;
    double                  best_xcorr ;
    
    /* Check Inputs */
    if (rawframes==NULL || trace_wave_frame==NULL || out_trace_wave==NULL ||
            lines_diagnostics == NULL || out_extracted == NULL || 
            out_wave_map==NULL || ext_plist==NULL) return -1 ;
    if (collapse!=CR2RES_COLLAPSE_MEAN && collapse!=CR2RES_COLLAPSE_MEDIAN) {
        return -1 ;
    }

    /* Initialise */
    best_xcorr = -1 ;

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
    if ((in_calib = cr2res_calib_imagelist(in, reduce_det, 0, 0,
                    master_flat_frame, master_dark_frame, bpm_frame, 
                    detlin_frame, dits)) == NULL) {
        cpl_msg_error(__func__, "Failed to apply the calibrations") ;
        if (dits != NULL) cpl_vector_delete(dits) ;
        hdrl_imagelist_delete(in) ;
        return -1 ;
    }
    hdrl_imagelist_delete(in) ;
    if (dits != NULL) cpl_vector_delete(dits) ;

    /* Collapse */
    if (collapse == CR2RES_COLLAPSE_MEAN) {
        cpl_msg_info(__func__, "Collapse (Mean) the input frames") ;
        cpl_msg_indent_more() ;
        hdrl_imagelist_collapse_mean(in_calib, &collapsed, &contrib) ;
    } else if (collapse == CR2RES_COLLAPSE_MEDIAN) {
        cpl_msg_info(__func__, "Collapse (Median) the input frames") ;
        cpl_msg_indent_more() ;
        hdrl_imagelist_collapse_median(in_calib, &collapsed, &contrib) ;
    } else {
        /* Should never happen */
        collapsed = NULL ;
        contrib = NULL ;
    }
    hdrl_imagelist_delete(in_calib) ;
    if (contrib != NULL) cpl_image_delete(contrib) ;
    if (cpl_error_get_code() != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Failed to Collapse") ;
        cpl_msg_indent_less() ;
        return -1 ;
    }
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
    
    /* Compute the Wavelength Calibration */
    cpl_msg_info(__func__, "Compute the Wavelength") ;
    if (cr2res_wave_apply(tw_in, extracted, lines_frame, reduce_order, 
                reduce_trace, wavecal_type, wl_degree, wl_start, wl_end, 
                wl_err_start, wl_err_end, wl_shift, log_flag, propagate_flag, 
                display, display_wmin, display_wmax,
                &qcs_plist,
                &lines_diagnostics_out,
                &extracted_out,
                &tw_out)) {
        cpl_msg_error(__func__, "Failed to calibrate");
        cpl_table_delete(tw_in) ;
        cpl_table_delete(extracted) ;
        return -1 ;
    }
    cpl_table_delete(tw_in) ;
    cpl_table_delete(extracted) ;

    /* Generate the Wave Map */
    wl_map = cr2res_wave_gen_wave_map(tw_out) ;

    /* Load the extension header for saving */
    first_file = cpl_frame_get_filename(
            cpl_frameset_get_position_const(rawframes, 0)) ;
    ext_nr = cr2res_io_get_ext_idx(first_file, reduce_det, 1) ;
    plist = cpl_propertylist_load(first_file, ext_nr) ;
    if (plist == NULL) {
        cpl_table_delete(tw_out) ;
        cpl_table_delete(lines_diagnostics_out) ;
        cpl_table_delete(extracted_out) ;
        hdrl_image_delete(wl_map) ;
        cpl_msg_error(__func__, "Failed to load the plist") ;
        return -1 ;
    }

    /* Store the QC parameters in the plist */
    if (qcs_plist != NULL) {
        cpl_propertylist_append(plist, qcs_plist) ;
        cpl_propertylist_delete(qcs_plist) ;
    }

    /* Return */
    *out_trace_wave = tw_out ;
    *lines_diagnostics = lines_diagnostics_out ;
    *out_extracted = extracted_out ;
    *out_wave_map = wl_map;
    *ext_plist = plist ;
    return 0 ;
}
