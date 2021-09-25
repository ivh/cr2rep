
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
        const cpl_frameset  *   rawframes_une,
        const cpl_frameset  *   rawframes_fpet,
        const cpl_frame     *   detlin_frame,
        const cpl_frame     *   master_dark_frame,
        const cpl_frame     *   master_flat_frame,
        const cpl_frame     *   bpm_frame,
        const cpl_frame     *   trace_wave_frame,
        const cpl_frame     *   lines_frame,
        int                     reduce_det,
        int                     reduce_order,
        int                     reduce_trace,
        int                     subtract_nolight_rows,
        cr2res_collapse         collapse,
        int                     ext_height,
        int                     ext_swath_width,
        int                     ext_oversample,
        double                  ext_smooth_slit,
        cr2res_wavecal_type     wavecal_type,
        int                     wl_degree,
        double                  wl_start,
        double                  wl_end,
        double                  wl_err,
        double                  wl_shift,
        int                     log_flag,
        int                     fallback_input_wavecal_flag,
        int                     keep_higher_degrees_flag,
        int                     clean_spectrum,
        int                     display,
        double                  display_wmin,
        double                  display_wmax,
        cpl_table           **  out_trace_wave_une,
        cpl_table           **  lines_diagnostics_une,
        cpl_table           **  out_extracted_une,
        hdrl_image          **  out_wave_map_une,
        cpl_propertylist    **  ext_plist_une,
        cpl_table           **  out_trace_wave_fpet,
        cpl_table           **  lines_diagnostics_fpet,
        cpl_table           **  out_extracted_fpet,
        hdrl_image          **  out_wave_map_fpet,
        cpl_propertylist    **  ext_plist_fpet) ;
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
                                                                        \n\
  Inputs                                                                \n\
    raw.fits " CR2RES_WAVE_UNE_RAW " [1 to n]                           \n\
          or " CR2RES_WAVE_FPET_RAW " [0 to n]                          \n\
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
          or " CR2RES_UTIL_BPM_MERGE_PROCATG "                          \n\
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
        -> out_trace_wave[une|fpet](d)                                  \n\
        -> lines_diagnostics[une|fpet](d)                               \n\
        -> out_extracted[une|fpet](d)                                   \n\
        -> out_wave_map[une|fpet](d)                                    \n\
    Save out_trace_wave[une|fpet]                                       \n\
    Save lines_diagnostics[une|fpet]                                    \n\
    Save out_extracted[une|fpet]                                        \n\
    Save out_wave_map[une|fpet]                                         \n\
                                                                        \n\
    cr2res_cal_wave_reduce():                                           \n\
      Successively for UNE and FPET RAW frames:                         \n\
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
  Library Functions used                                                \n\
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
                    CR2RES_PIPELINE_AUTHORS,
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

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.subtract_nolight_rows",
            CPL_TYPE_BOOL, "Subtract the no-light rows.",
            "cr2res.cr2res_cal_wave", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "keep");
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
            "cr2res.cr2res_cal_wave", 5);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ext_oversample");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.ext_swath_width",
            CPL_TYPE_INT, "The swath width", "cr2res.cr2res_cal_wave", 800);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ext_swath_width");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.ext_height",
            CPL_TYPE_INT, "Extraction height",
            "cr2res.cr2res_cal_wave", 120);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ext_height");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.ext_smooth_slit",
            CPL_TYPE_DOUBLE,
            "Smoothing along the slit",
            "cr2res.cr2res_cal_wave", 3.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ext_smooth_slit");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.wl_method",
            CPL_TYPE_STRING, 
            "Wavelength Method (XCORR / LINE1D / LINE2D)",
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
            CPL_TYPE_DOUBLE, "Estimated wavelength error (in nm)",
            "cr2res.cr2res_cal_wave", -1.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "wl_err");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.wl_degree",
            CPL_TYPE_INT, "Wavelegth Polynomial degree",
            "cr2res.cr2res_cal_wave", 2);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "wl_degree");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.log",
            CPL_TYPE_BOOL, "Flag for taking the Log() value of the lines",
            "cr2res.cr2res_cal_wave", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "log");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.fallback_input_wavecal",
            CPL_TYPE_BOOL, "Flag for using the input WL when no computation",
            "cr2res.cr2res_cal_wave", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "fallback");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.keep_higher_degrees",
            CPL_TYPE_BOOL,
            "Flag for re-using higher degrees of first guess in Cross-Corr.",
            "cr2res.cr2res_cal_wave", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "keep");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.clean_spectrum",
            CPL_TYPE_BOOL, "Flag to automatically clean the missing lines",
            "cr2res.cr2res_cal_wave", TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "clean_spectrum");
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
                            wl_degree, display, log_flag,
                            fallback_input_wavecal_flag,
                            keep_higher_degrees_flag, 
                            clean_spectrum, subtract_nolight_rows ;
    double                  ext_smooth_slit, wl_start, wl_end, wl_err, wl_shift,
                            display_wmin, display_wmax ;
    cr2res_collapse         collapse ;
    cr2res_wavecal_type     wavecal_type ;
    const char          *   sval ;
    cpl_frameset        *   rawframes_une ;
    cpl_frameset        *   rawframes_fpet ;
    const cpl_frame     *   detlin_frame ;
    const cpl_frame     *   master_dark_frame ;
    const cpl_frame     *   master_flat_frame ;
    const cpl_frame     *   bpm_frame ;
    const cpl_frame     *   trace_wave_frame ;
    const cpl_frame     *   lines_frame ;
    char                *   out_file;
    cpl_table           *   out_trace_wave_une[CR2RES_NB_DETECTORS] ;
    cpl_table           *   lines_diagnostics_une[CR2RES_NB_DETECTORS] ;
    cpl_table           *   out_extracted_une[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   out_wave_map_une[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   ext_plist_une[CR2RES_NB_DETECTORS] ;
    cpl_table           *   out_trace_wave_fpet[CR2RES_NB_DETECTORS] ;
    cpl_table           *   lines_diagnostics_fpet[CR2RES_NB_DETECTORS] ;
    cpl_table           *   out_extracted_fpet[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   out_wave_map_fpet[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   ext_plist_fpet[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   plist ;
    int                     det_nr, order, i ;

    /* Needed for sscanf() */
    setlocale(LC_NUMERIC, "C");

    /* Initialise */
    wl_start = wl_end = wl_err = -1.0 ;
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
            "cr2res.cr2res_cal_wave.subtract_nolight_rows");
    subtract_nolight_rows = cpl_parameter_get_bool(param) ;
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
    else {
        cpl_msg_error(__func__, "Invalid Method specified");
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
    wl_err = cpl_parameter_get_double(param) ;
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_wave.wl_degree");
    wl_degree = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_wave.log");
    log_flag = cpl_parameter_get_bool(param) ;
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_wave.fallback_input_wavecal");
    fallback_input_wavecal_flag = cpl_parameter_get_bool(param) ;
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_wave.keep_higher_degrees");
    keep_higher_degrees_flag = cpl_parameter_get_bool(param) ;
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_wave.clean_spectrum");
    clean_spectrum = cpl_parameter_get_bool(param) ;
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
    if (display==TRUE && (reduce_order==-1 || reduce_det==0)){
        cpl_msg_warning(__func__,
                "Option --display can only be used with --order"
                " and --detector");
        display=FALSE;
    }
    if (wl_degree < 0) {
        cpl_msg_error(__func__, "The degree needs to be >= 0");
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }
    if (wl_degree == 0 && !keep_higher_degrees_flag) {
        cpl_msg_error(__func__, "The degree 0 can only be used with --keep");
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
    rawframes_une = cr2res_extract_frameset(frameset, CR2RES_WAVE_UNE_RAW) ;
    rawframes_fpet = cr2res_extract_frameset(frameset, CR2RES_WAVE_FPET_RAW) ;
    if (rawframes_une == NULL && rawframes_fpet == NULL) {
        cpl_msg_error(__func__, "Could not find RAW frames") ;
        return -1 ;
    }
    if (lines_frame == NULL && rawframes_une != NULL) {
        if (rawframes_une !=NULL) cpl_frameset_delete(rawframes_une) ;
        if (rawframes_fpet!=NULL) cpl_frameset_delete(rawframes_fpet) ;
        cpl_msg_error(__func__, "The emission lines file is needed");
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }
    if (reduce_order > -1 && wavecal_type == CR2RES_LINE2D) {
        if (rawframes_une !=NULL) cpl_frameset_delete(rawframes_une) ;
        if (rawframes_fpet!=NULL) cpl_frameset_delete(rawframes_fpet) ;
        cpl_msg_error(__func__, "Limiting to one order with LINE2D impossible");
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Loop over the detectors */
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {

        /* Initialise */
        out_trace_wave_une[det_nr-1] = NULL ;
        lines_diagnostics_une[det_nr-1] = NULL ;
        out_extracted_une[det_nr-1] = NULL ;
        out_wave_map_une[det_nr-1] = NULL ;
        ext_plist_une[det_nr-1] = NULL ;

        out_trace_wave_fpet[det_nr-1] = NULL ;
        lines_diagnostics_fpet[det_nr-1] = NULL ;
        out_extracted_fpet[det_nr-1] = NULL ;
        out_wave_map_fpet[det_nr-1] = NULL ;
        ext_plist_fpet[det_nr-1] = NULL ;

        /* Compute only one detector but not this one */
        if (reduce_det != 0 && det_nr != reduce_det) {
            /* This Detector will not be processed here */
            /* The output trace wave contains the input one ... */
            if (rawframes_une != NULL) 
                out_trace_wave_une[det_nr-1] = cr2res_io_load_TRACE_WAVE(
                        cpl_frame_get_filename(trace_wave_frame), det_nr) ;
            if (rawframes_fpet != NULL) 
                out_trace_wave_fpet[det_nr-1] = cr2res_io_load_TRACE_WAVE(
                        cpl_frame_get_filename(trace_wave_frame), det_nr) ;
            /*    ...  without the WL / WL_ERR */
            /*    ... unless fallback_input_wavecal_flag is set */
            
            if (!fallback_input_wavecal_flag) {
                /* Reset WL / WL_ERR */
                if (rawframes_une != NULL) {
                    cpl_table_erase_column(out_trace_wave_une[det_nr-1],
                            CR2RES_COL_WAVELENGTH) ;
                    cpl_table_new_column_array(out_trace_wave_une[det_nr-1],
                            CR2RES_COL_WAVELENGTH, CPL_TYPE_DOUBLE, 2) ;
                    cpl_table_erase_column(out_trace_wave_une[det_nr-1],
                            CR2RES_COL_WAVELENGTH_ERROR) ;
                    cpl_table_new_column_array(out_trace_wave_une[det_nr-1],  
                            CR2RES_COL_WAVELENGTH_ERROR, CPL_TYPE_DOUBLE, 2) ;
                }
                if (rawframes_fpet != NULL) {
                    cpl_table_erase_column(out_trace_wave_fpet[det_nr-1],
                            CR2RES_COL_WAVELENGTH) ;
                    cpl_table_new_column_array(out_trace_wave_fpet[det_nr-1],
                            CR2RES_COL_WAVELENGTH, CPL_TYPE_DOUBLE, 2) ;
                    cpl_table_erase_column(out_trace_wave_fpet[det_nr-1],
                            CR2RES_COL_WAVELENGTH_ERROR) ;
                    cpl_table_new_column_array(out_trace_wave_fpet[det_nr-1],  
                            CR2RES_COL_WAVELENGTH_ERROR, CPL_TYPE_DOUBLE, 2) ;
                }
            }
            continue ;
        }

        cpl_msg_info(__func__, "Process detector number %d", det_nr) ;
        cpl_msg_indent_more() ;

        /* Call the reduction function */
        if (cr2res_cal_wave_reduce(rawframes_une, rawframes_fpet, detlin_frame,
                    master_dark_frame, master_flat_frame, bpm_frame,
                    trace_wave_frame, lines_frame, det_nr, reduce_order,
                    reduce_trace, subtract_nolight_rows, collapse, ext_height, 
                    ext_swath_width, ext_oversample, ext_smooth_slit, 
                    wavecal_type, wl_degree, wl_start, wl_end, wl_err, 
                    wl_shift, log_flag, fallback_input_wavecal_flag,
                    keep_higher_degrees_flag, clean_spectrum,
                    display, display_wmin, 
                    display_wmax, 
                    &(out_trace_wave_une[det_nr-1]),
                    &(lines_diagnostics_une[det_nr-1]),
                    &(out_extracted_une[det_nr-1]),
                    &(out_wave_map_une[det_nr-1]),
                    &(ext_plist_une[det_nr-1]),
                    &(out_trace_wave_fpet[det_nr-1]),
                    &(lines_diagnostics_fpet[det_nr-1]),
                    &(out_extracted_fpet[det_nr-1]),
                    &(out_wave_map_fpet[det_nr-1]),
                    &(ext_plist_fpet[det_nr-1])) == -1) {
            cpl_msg_warning(__func__, "Failed to reduce detector %d", det_nr);
            cpl_error_reset() ;
        }
        cpl_msg_indent_less() ;
    }

    /* Save Products UNE */
    if (rawframes_une != NULL) {
        out_file = cpl_sprintf("%s_tw_une.fits", RECIPE_STRING) ;
        cr2res_io_save_TRACE_WAVE(out_file, frameset, rawframes_une, parlist, 
                out_trace_wave_une, NULL, ext_plist_une, 
                CR2RES_CAL_WAVE_TW_PROCATG, RECIPE_STRING) ;
        cpl_free(out_file);

        out_file = cpl_sprintf("%s_wave_map_une.fits", RECIPE_STRING) ;
        cr2res_io_save_WAVE_MAP(out_file, frameset, rawframes_une, parlist, 
                out_wave_map_une, NULL, ext_plist_une, 
                CR2RES_CAL_WAVE_MAP_PROCATG, RECIPE_STRING) ;
        cpl_free(out_file);

        out_file = cpl_sprintf("%s_extracted_une.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes_une, parlist, 
                out_extracted_une, NULL, ext_plist_une, 
                CR2RES_CAL_WAVE_EXTRACT_1D_PROCATG, RECIPE_STRING) ;
        cpl_free(out_file);

        out_file = cpl_sprintf("%s_lines_diagnostics_une.fits", 
                RECIPE_STRING);
        cr2res_io_save_LINES_DIAGNOSTICS(out_file, frameset,
                rawframes_une, parlist,
                lines_diagnostics_une, NULL, ext_plist_une,
                CR2RES_CAL_WAVE_LINES_DIAGNOSTICS_PROCATG, RECIPE_STRING) ;
        cpl_free(out_file);
    }

    if (rawframes_fpet != NULL) {
        /* Save Products UNE */
        out_file = cpl_sprintf("%s_tw_fpet.fits", RECIPE_STRING) ;
        cr2res_io_save_TRACE_WAVE(out_file, frameset, rawframes_fpet, parlist, 
                out_trace_wave_fpet, NULL, ext_plist_fpet, 
                CR2RES_CAL_WAVE_TW_PROCATG, RECIPE_STRING) ;
        cpl_free(out_file);

        out_file = cpl_sprintf("%s_wave_map_fpet.fits", RECIPE_STRING) ;
        cr2res_io_save_WAVE_MAP(out_file, frameset, rawframes_fpet, parlist, 
                out_wave_map_fpet, NULL, ext_plist_fpet, 
                CR2RES_CAL_WAVE_MAP_PROCATG, RECIPE_STRING) ;
        cpl_free(out_file);

        out_file = cpl_sprintf("%s_extracted_fpet.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes_fpet, parlist, 
                out_extracted_fpet, NULL, ext_plist_fpet, 
                CR2RES_CAL_WAVE_EXTRACT_1D_PROCATG, RECIPE_STRING) ;
        cpl_free(out_file);

        out_file = cpl_sprintf("%s_lines_diagnostics_fpet.fits", 
                RECIPE_STRING);
        cr2res_io_save_LINES_DIAGNOSTICS(out_file, frameset,
                rawframes_fpet, parlist,
                lines_diagnostics_fpet, NULL, ext_plist_fpet,
                CR2RES_CAL_WAVE_LINES_DIAGNOSTICS_PROCATG, RECIPE_STRING) ;
        cpl_free(out_file);
    }

    /* Free and return */
    cpl_frameset_delete(rawframes_une) ;
    if (rawframes_fpet!=NULL) cpl_frameset_delete(rawframes_fpet) ;
    for (i=0 ; i<CR2RES_NB_DETECTORS ; i++) {
        if (ext_plist_une[i] != NULL)
            cpl_propertylist_delete(ext_plist_une[i]) ;
        if (out_trace_wave_une[i] != NULL)
            cpl_table_delete(out_trace_wave_une[i]) ;
        if (lines_diagnostics_une[i] != NULL)
            cpl_table_delete(lines_diagnostics_une[i]) ;
        if (out_extracted_une[i] != NULL)
            cpl_table_delete(out_extracted_une[i]) ;
        if (out_wave_map_une[i] != NULL)
            hdrl_image_delete(out_wave_map_une[i]) ;
        if (ext_plist_fpet[i] != NULL)
            cpl_propertylist_delete(ext_plist_fpet[i]) ;
        if (out_trace_wave_fpet[i] != NULL)
            cpl_table_delete(out_trace_wave_fpet[i]) ;
        if (lines_diagnostics_fpet[i] != NULL)
            cpl_table_delete(lines_diagnostics_fpet[i]) ;
        if (out_extracted_fpet[i] != NULL)
            cpl_table_delete(out_extracted_fpet[i]) ;
        if (out_wave_map_fpet[i] != NULL) 
            hdrl_image_delete(out_wave_map_fpet[i]) ;
    }
    return (int)cpl_error_get_code();
}

/*----------------------------------------------------------------------------*/
/**
  @brief Compute the Wavelength for a detector
  @param rawframes_une      Input raw frames UNE
  @param rawframes_fpet     Input raw frames FPET (or NULL)
  @param detlin_frame       Associated detlin coefficients
  @param master_dark_frame  Associated master dark
  @param master_flat_frame  Associated master flat
  @param bpm_frame          Associated BPM
  @param trace_wave_frame   Trace Wave table
  @param lines_frame        Emission lines frame
  @param reduce_det         The detector to compute
  @param reduce_order       The order to compute (-1 for all)
  @param reduce_trace       The trace to compute (-1 for all)
  @param subtract_nolight_rows
  @param collapse           CR2RES_COLLAPSE_MEAN or CR2RES_COLLAPSE_MEDIAN
  @param ext_height         Extraction related
  @param ext_swath_width    Extraction related
  @param ext_oversample     Extraction related
  @param ext_smooth_slit    Extraction related
  @param wavecal_type       CR2RES_XCORR/LINE1D/LINE2D
  @param wl_start           WL estimate of the first pixel
  @param wl_end             WL estimate of the last pixel
  @param wl_err             WL error 
  @param wl_shift           wavelength shift to apply
  @param log_flag           Flag to apply a log() to the lines intensities
  @param fallback_input_wavecal_flag Flag to copy the input WL to the output 
                            when they are not computed
  @param keep_higher_degrees_flag  Flag to use higher polynomial degrees
                            from the guess
  @param clean_spectrum     Remove the lines that are not in the catalog (1d)
  @param display            Flag to enable display functionalities
  @param display_wmin       Minimum Wavelength to  display
  @param display_wmax       Maximum Wavelength to  display
  @param out_trace_wave_une      [out] trace wave table
  @param lines_diagnostics_une   [out] lines diagnostics table
  @param out_extracted_une       [out] extracted table with updated WL
  @param out_wave_map_une        [out] Wave map
  @param ext_plist_une           [out] the header for saving the products
  @param out_trace_wave_fpet     [out] trace wave table
  @param lines_diagnostics_fpet  [out] lines diagnostics table
  @param out_extracted_fpet      [out] extracted table with updated WL
  @param out_wave_map_fpet       [out] Wave map
  @param ext_plist_fpet          [out] the header for saving the products
  @return   0 if ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
static int cr2res_cal_wave_reduce(
        const cpl_frameset  *   rawframes_une,
        const cpl_frameset  *   rawframes_fpet,
        const cpl_frame     *   detlin_frame,
        const cpl_frame     *   master_dark_frame,
        const cpl_frame     *   master_flat_frame,
        const cpl_frame     *   bpm_frame,
        const cpl_frame     *   trace_wave_frame,
        const cpl_frame     *   lines_frame,
        int                     reduce_det,
        int                     reduce_order,
        int                     reduce_trace,
        int                     subtract_nolight_rows,
        cr2res_collapse         collapse,
        int                     ext_height, 
        int                     ext_swath_width,
        int                     ext_oversample, 
        double                  ext_smooth_slit,
        cr2res_wavecal_type     wavecal_type,
        int                     wl_degree,
        double                  wl_start,
        double                  wl_end,
        double                  wl_err,
        double                  wl_shift,
        int                     log_flag,
        int                     fallback_input_wavecal_flag,
        int                     keep_higher_degrees_flag,
        int                     clean_spectrum,
        int                     display,
        double                  display_wmin,
        double                  display_wmax,
        cpl_table           **  out_trace_wave_une,
        cpl_table           **  lines_diagnostics_une,
        cpl_table           **  out_extracted_une,
        hdrl_image          **  out_wave_map_une,
        cpl_propertylist    **  ext_plist_une,
        cpl_table           **  out_trace_wave_fpet,
        cpl_table           **  lines_diagnostics_fpet,
        cpl_table           **  out_extracted_fpet,
        hdrl_image          **  out_wave_map_fpet,
        cpl_propertylist    **  ext_plist_fpet)
{
    cpl_vector          *   dits_une ;
    hdrl_imagelist      *   in_une ;
    hdrl_imagelist      *   in_fpet ;
    hdrl_imagelist      *   in_une_calib ;
    hdrl_imagelist      *   in_fpet_calib ;
    hdrl_image          *   collapsed_une ;
    hdrl_image          *   collapsed_fpet ;
    cpl_image           *   contrib ;
    cpl_table           *   tw_in ;
    cpl_table           *   extracted_une ;
    cpl_table           *   slit_func_une ;
    hdrl_image          *   model_master_une ;
    cpl_table           *   extracted_fpet ;
    cpl_table           *   slit_func_fpet ;
    hdrl_image          *   model_master_fpet ;
    hdrl_image          *   wl_map_une_out ;
    cpl_table           *   tw_une_out ;
    cpl_table           *   lines_diagnostics_une_out ;
    cpl_table           *   extracted_une_out ;
    cpl_propertylist    *   qcs_une_out ;
    cpl_propertylist    *   plist_une_out ;
    hdrl_image          *   wl_map_fpet_out ;
    cpl_table           *   tw_fpet_out ;
    cpl_table           *   lines_diagnostics_fpet_out ;
    cpl_table           *   extracted_fpet_out ;
    cpl_propertylist    *   qcs_fpet_out ;
    cpl_propertylist    *   plist_fpet_out ;
    cpl_propertylist    *   plist ;
    const char          *   first_file ;
    int                     ext_nr, zp_order_une, zp_order_fpet,
                            grat1_order_une, grat1_order_fpet ;
    
    /* Check Inputs */
    if (rawframes_fpet==NULL && rawframes_une==NULL) return -1 ;
    if (trace_wave_frame==NULL || out_trace_wave_une==NULL || 
            out_trace_wave_fpet==NULL || lines_diagnostics_une == NULL || 
            lines_diagnostics_fpet == NULL || out_extracted_une == NULL || 
            out_extracted_fpet == NULL || out_wave_map_une == NULL || 
            out_wave_map_fpet == NULL || ext_plist_une == NULL || 
            ext_plist_fpet == NULL) return -1 ;
    if (collapse!=CR2RES_COLLAPSE_MEAN && collapse!=CR2RES_COLLAPSE_MEDIAN) 
        return -1 ;

    /* Reduce the UNE */
    if (rawframes_une != NULL) {
        cpl_msg_info(__func__, "Reduce %"CPL_SIZE_FORMAT" UNE Frames",
                cpl_frameset_get_size(rawframes_une)) ;
        cpl_msg_indent_more() ;

        /* Load the UNE DITs if necessary */
        if (master_dark_frame != NULL)  
            dits_une = cr2res_io_read_dits(rawframes_une) ;
        else                            
            dits_une = NULL ;
        if (cpl_msg_get_level() == CPL_MSG_DEBUG && dits_une != NULL)
            cpl_vector_dump(dits_une, stdout) ;

        /* Load UNE image list */
        if ((in_une = cr2res_io_load_image_list_from_set(rawframes_une,
                        reduce_det)) == NULL) {
            cpl_msg_error(__func__, "Cannot load images") ;
            if (dits_une != NULL) cpl_vector_delete(dits_une) ;
            cpl_msg_indent_less() ;
            return -1 ;
        }
        if (hdrl_imagelist_get_size(in_une) !=
                cpl_frameset_get_size(rawframes_une)) {
            cpl_msg_error(__func__, "Inconsistent number of loaded images") ;
            if (dits_une != NULL) cpl_vector_delete(dits_une) ;
            hdrl_imagelist_delete(in_une) ;
            cpl_msg_indent_less() ;
            return -1 ;
        }
        /* Calibrate the UNE images */
        if ((in_une_calib = cr2res_calib_imagelist(in_une, reduce_det, 0,
                        subtract_nolight_rows, 0, master_flat_frame, 
                        master_dark_frame, bpm_frame, detlin_frame, 
                        dits_une)) == NULL) {
            cpl_msg_error(__func__, "Failed to apply the calibrations") ;
            if (dits_une != NULL) cpl_vector_delete(dits_une) ;
            hdrl_imagelist_delete(in_une) ;
            cpl_msg_indent_less() ;
            return -1 ;
        }
        hdrl_imagelist_delete(in_une) ;
        if (dits_une != NULL) cpl_vector_delete(dits_une) ;

        /* Collapse UNE */
        contrib = NULL ;
        if (collapse == CR2RES_COLLAPSE_MEAN) {
            cpl_msg_info(__func__, "Collapse (Mean) the input UNE frames") ;
            cpl_msg_indent_more() ;
            hdrl_imagelist_collapse_mean(in_une_calib, &collapsed_une,
                    &contrib) ;
            cpl_msg_indent_less() ;
        } else if (collapse == CR2RES_COLLAPSE_MEDIAN) {
            cpl_msg_info(__func__, "Collapse (Median) the input UNE frames") ;
            cpl_msg_indent_more() ;
            hdrl_imagelist_collapse_median(in_une_calib, &collapsed_une,
                    &contrib) ;
            cpl_msg_indent_less() ;
        } else {
            /* Should never happen */
            collapsed_une = NULL ;
            contrib = NULL ;
        }
        hdrl_imagelist_delete(in_une_calib) ;
        if (contrib != NULL) cpl_image_delete(contrib) ;
        if (cpl_error_get_code() != CPL_ERROR_NONE) {
            cpl_msg_error(__func__, "Failed Collapse: %d",cpl_error_get_code());
            cpl_msg_indent_less() ;
            return -1 ;
        }

        /* Load the trace wave */
        cpl_msg_info(__func__, "Load the TRACE WAVE") ;
        if ((tw_in = cr2res_io_load_TRACE_WAVE(cpl_frame_get_filename(
                            trace_wave_frame), reduce_det)) == NULL) {
            cpl_msg_error(__func__, "Failed to Load the traces file") ;
            hdrl_image_delete(collapsed_une) ;
            cpl_msg_indent_less() ;
            return -1 ;
        }

        /* Execute the extraction for UNE */
        cpl_msg_info(__func__, "Spectra Extraction UNE") ;
        cpl_msg_indent_more() ;
        if (cr2res_extract_traces(collapsed_une, tw_in, NULL, reduce_order, 
                    reduce_trace, CR2RES_EXTR_OPT_CURV, ext_height, 
                    ext_swath_width, ext_oversample, ext_smooth_slit, 0.0, 
                    0, 0, 0, // display flags
                    &extracted_une, &slit_func_une, &model_master_une) == -1) {
            cpl_msg_error(__func__, "Failed to extract");
            hdrl_image_delete(collapsed_une) ;
            cpl_table_delete(tw_in) ;
            cpl_msg_indent_less() ;
            cpl_msg_indent_less() ;
            return -1 ;
        }
        cpl_msg_indent_less() ;
        cpl_table_delete(slit_func_une) ;
        hdrl_image_delete(model_master_une) ;
        hdrl_image_delete(collapsed_une);

        /* Compute the Wavelength Calibration for UNE */
        first_file = cpl_frame_get_filename(
                cpl_frameset_get_position_const(rawframes_une, 0)) ;
        plist = cpl_propertylist_load(first_file, 0) ;
        zp_order_une = cr2res_pfits_get_order_zp(plist) ;
        grat1_order_une = cr2res_pfits_get_order(plist) ;
        cpl_propertylist_delete(plist);
       
        cpl_msg_info(__func__, "Compute the Wavelength for UNE") ;
        cpl_msg_indent_more() ;
        if (cr2res_wave_apply(tw_in, extracted_une, lines_frame, reduce_order, 
                    reduce_trace, wavecal_type, wl_degree, wl_start, wl_end, 
                    wl_err, wl_shift, log_flag, fallback_input_wavecal_flag, 
                    keep_higher_degrees_flag, clean_spectrum, 
                    display, display_wmin, display_wmax, zp_order_une,
                    grat1_order_une,
                    &qcs_une_out,
                    &lines_diagnostics_une_out,
                    &extracted_une_out,
                    &tw_une_out) || cpl_error_get_code()) {
            cpl_msg_error(__func__, "Failed to calibrate");
            cpl_table_delete(tw_in) ;
            cpl_table_delete(extracted_une) ;
            cpl_msg_indent_less() ;
            cpl_msg_indent_less() ;
            return -1 ;
        }
        cpl_msg_indent_less() ;
        cpl_table_delete(tw_in) ;
        cpl_table_delete(extracted_une) ;

        /* Generate the Wave Map */
        wl_map_une_out = cr2res_wave_gen_wave_map(tw_une_out) ;

        /* Load the extension header for saving UNE*/
        first_file = cpl_frame_get_filename(
                cpl_frameset_get_position_const(rawframes_une, 0)) ;
        ext_nr = cr2res_io_get_ext_idx(first_file, reduce_det, 1) ;
        plist_une_out = cpl_propertylist_load(first_file, ext_nr) ;
        if (plist_une_out == NULL) {
            cpl_propertylist_delete(qcs_une_out) ;
            cpl_table_delete(tw_une_out) ;
            cpl_table_delete(lines_diagnostics_une_out) ;
            cpl_table_delete(extracted_une_out) ;
            hdrl_image_delete(wl_map_une_out) ;
            cpl_msg_error(__func__, "Failed to load the plist") ;
            cpl_msg_indent_less() ;
            return -1 ;
        }
        if (qcs_une_out != NULL) {
            cpl_propertylist_append(plist_une_out, qcs_une_out) ;
            cpl_propertylist_delete(qcs_une_out) ;
        }
        cpl_msg_indent_less() ;
    } else {
        tw_une_out = NULL ;
        lines_diagnostics_une_out = NULL ;
        extracted_une_out = NULL ;
        wl_map_une_out = NULL ;
        plist_une_out = NULL ;
    }

    /* Reduce the FPET */
    if (rawframes_fpet != NULL) {
        cpl_msg_info(__func__, "Reduce %"CPL_SIZE_FORMAT" FPET Frames",
                cpl_frameset_get_size(rawframes_fpet)) ;
        cpl_msg_indent_more() ;

        /* Load FPET image list */
        if ((in_fpet = cr2res_io_load_image_list_from_set(rawframes_fpet,
                        reduce_det)) == NULL) {
            cpl_msg_error(__func__, "Cannot load images") ;
            if (plist_une_out != NULL) cpl_propertylist_delete(plist_une_out) ;
            if (tw_une_out != NULL) cpl_table_delete(tw_une_out) ;
            if (lines_diagnostics_une_out != NULL) 
                cpl_table_delete(lines_diagnostics_une_out) ;
            if (extracted_une_out != NULL) cpl_table_delete(extracted_une_out) ;
            if (wl_map_une_out != NULL) hdrl_image_delete(wl_map_une_out) ;
            return -1 ;
        }
        if (hdrl_imagelist_get_size(in_fpet) !=
                cpl_frameset_get_size(rawframes_fpet)) {
            cpl_msg_error(__func__, "Inconsistent number of loaded images") ;
            hdrl_imagelist_delete(in_fpet) ;
            if (plist_une_out != NULL) cpl_propertylist_delete(plist_une_out) ;
            if (tw_une_out != NULL) cpl_table_delete(tw_une_out) ;
            if (lines_diagnostics_une_out != NULL) 
                cpl_table_delete(lines_diagnostics_une_out) ;
            if (extracted_une_out != NULL) cpl_table_delete(extracted_une_out) ;
            if (wl_map_une_out != NULL) hdrl_image_delete(wl_map_une_out) ;
            return -1 ;
        }

        /* Calibrate the FPET images */
        if ((in_fpet_calib = cr2res_calib_imagelist(in_fpet, reduce_det, 0, 
                        subtract_nolight_rows, 0, master_flat_frame, NULL, 
                        bpm_frame, detlin_frame, NULL)) == NULL) {
            cpl_msg_error(__func__, "Failed to apply the calibrations") ;
            if (in_fpet != NULL) hdrl_imagelist_delete(in_fpet) ;
            if (plist_une_out != NULL) cpl_propertylist_delete(plist_une_out) ;
            if (tw_une_out != NULL) cpl_table_delete(tw_une_out) ;
            if (lines_diagnostics_une_out != NULL) 
                cpl_table_delete(lines_diagnostics_une_out) ;
            if (extracted_une_out != NULL) cpl_table_delete(extracted_une_out) ;
            if (wl_map_une_out != NULL) hdrl_image_delete(wl_map_une_out) ;
            return -1 ;
        }
        hdrl_imagelist_delete(in_fpet) ;

        /* Collapse FPET */
        contrib = NULL ;
        if (collapse == CR2RES_COLLAPSE_MEAN) {
            cpl_msg_info(__func__, "Collapse (Mean) the input FPET frames") ;
            cpl_msg_indent_more() ;
            hdrl_imagelist_collapse_mean(in_fpet_calib, &collapsed_fpet,
                    &contrib) ;
            cpl_msg_indent_less() ;
        } else if (collapse == CR2RES_COLLAPSE_MEDIAN) {
            cpl_msg_info(__func__, "Collapse (Median) the input FPET frames") ;
            cpl_msg_indent_more() ;
            hdrl_imagelist_collapse_median(in_fpet_calib, &collapsed_fpet,
                    &contrib) ;
            cpl_msg_indent_less() ;
        } else {
            /* Should never happen */
            collapsed_fpet = NULL ;
            contrib = NULL ;
        }
        hdrl_imagelist_delete(in_fpet_calib) ;
        if (contrib != NULL) cpl_image_delete(contrib) ;
        if (cpl_error_get_code() != CPL_ERROR_NONE) {
            cpl_msg_error(__func__, 
                    "Failed Collapse: %d", cpl_error_get_code()) ;
            if (collapsed_fpet != NULL) hdrl_image_delete(collapsed_fpet) ;
            if (plist_une_out != NULL) cpl_propertylist_delete(plist_une_out) ;
            if (tw_une_out != NULL) cpl_table_delete(tw_une_out) ;
            if (lines_diagnostics_une_out != NULL) 
                cpl_table_delete(lines_diagnostics_une_out) ;
            if (extracted_une_out != NULL) cpl_table_delete(extracted_une_out) ;
            if (wl_map_une_out != NULL) hdrl_image_delete(wl_map_une_out) ;
            cpl_msg_indent_less() ;
            return -1 ;
        }

        /* Get the trace wave */
        if (tw_une_out != NULL) {
            cpl_msg_info(__func__, "Use the UNE output TRACE WAVE") ;
            tw_in = cpl_table_duplicate(tw_une_out) ;
        } else {
            cpl_msg_info(__func__, "Load the TRACE WAVE") ;
            if ((tw_in = cr2res_io_load_TRACE_WAVE(cpl_frame_get_filename(
                                trace_wave_frame), reduce_det)) == NULL) {
                cpl_msg_error(__func__, "Failed to Load the traces file") ;
                if (collapsed_fpet != NULL) hdrl_image_delete(collapsed_fpet) ;
                if (plist_une_out != NULL) 
                    cpl_propertylist_delete(plist_une_out) ;
                if (tw_une_out != NULL) cpl_table_delete(tw_une_out) ;
                if (lines_diagnostics_une_out != NULL) 
                    cpl_table_delete(lines_diagnostics_une_out) ;
                if (extracted_une_out != NULL) 
                    cpl_table_delete(extracted_une_out) ;
                if (wl_map_une_out != NULL) hdrl_image_delete(wl_map_une_out) ;
                cpl_msg_indent_less() ;
                return -1 ;
            }
        }

        /* Execute the extraction for FPET */
        cpl_msg_info(__func__, "Spectra Extraction FPET") ;
        cpl_msg_indent_more() ;
        if (cr2res_extract_traces(collapsed_fpet, tw_in, NULL, reduce_order, 
                    reduce_trace, CR2RES_EXTR_OPT_CURV, ext_height, 
                    ext_swath_width, ext_oversample, ext_smooth_slit, 0.0, 
                    0, 0, 0, // display flags
                    &extracted_fpet, &slit_func_fpet, &model_master_fpet)==-1) {
            cpl_msg_error(__func__, "Failed to extract");
            cpl_table_delete(tw_in) ;
            if (collapsed_fpet != NULL) hdrl_image_delete(collapsed_fpet) ;
            if (plist_une_out != NULL) cpl_propertylist_delete(plist_une_out) ;
            if (tw_une_out != NULL) cpl_table_delete(tw_une_out) ;
            if (lines_diagnostics_une_out != NULL) 
                cpl_table_delete(lines_diagnostics_une_out) ;
            if (extracted_une_out != NULL) cpl_table_delete(extracted_une_out) ;
            if (wl_map_une_out != NULL) hdrl_image_delete(wl_map_une_out) ;
            cpl_msg_indent_less() ;
            cpl_msg_indent_less() ;
            return -1 ;
        }
        cpl_msg_indent_less() ;
        cpl_table_delete(slit_func_fpet) ;
        hdrl_image_delete(model_master_fpet) ;
        hdrl_image_delete(collapsed_fpet);

        /* Compute the Wavelength Calibration for FPET */
        first_file = cpl_frame_get_filename(
                cpl_frameset_get_position_const(rawframes_fpet, 0)) ;
        plist = cpl_propertylist_load(first_file, 0) ;
        zp_order_fpet = cr2res_pfits_get_order_zp(plist) ;
        grat1_order_fpet = cr2res_pfits_get_order(plist) ;
        cpl_propertylist_delete(plist);
       
        cpl_msg_info(__func__, "Compute the Wavelength for FPET") ;
        cpl_msg_indent_more() ;
        if (cr2res_wave_apply(tw_une_out, extracted_fpet, NULL, reduce_order, 
                    reduce_trace, CR2RES_ETALON, wl_degree, wl_start, wl_end, 
                    wl_err, wl_shift, log_flag, fallback_input_wavecal_flag, 
                    keep_higher_degrees_flag, clean_spectrum, 
                    display, display_wmin, display_wmax, zp_order_fpet,
                    grat1_order_fpet,
                    &qcs_fpet_out,
                    &lines_diagnostics_fpet_out,
                    &extracted_fpet_out,
                    &tw_fpet_out) || cpl_error_get_code()) {
            cpl_msg_error(__func__, "Failed to calibrate");
            cpl_table_delete(extracted_fpet) ;
            if (plist_une_out != NULL) cpl_propertylist_delete(plist_une_out) ;
            if (tw_une_out != NULL) cpl_table_delete(tw_une_out) ;
            if (lines_diagnostics_une_out != NULL) 
                cpl_table_delete(lines_diagnostics_une_out) ;
            if (extracted_une_out != NULL) cpl_table_delete(extracted_une_out) ;
            if (wl_map_une_out != NULL) hdrl_image_delete(wl_map_une_out) ;
            cpl_msg_indent_less() ;
            cpl_msg_indent_less() ;
            return -1 ;
        }
        cpl_msg_indent_less() ;
        cpl_table_delete(tw_in) ;
        cpl_table_delete(extracted_fpet) ;

        /* Generate the Wave Map */
        wl_map_fpet_out = cr2res_wave_gen_wave_map(tw_fpet_out) ;

        /* Load the extension header for saving FPET */
        first_file = cpl_frame_get_filename(
                cpl_frameset_get_position_const(rawframes_fpet, 0)) ;
        ext_nr = cr2res_io_get_ext_idx(first_file, reduce_det, 1) ;
        plist_fpet_out = cpl_propertylist_load(first_file, ext_nr) ;
        if (plist_fpet_out == NULL) {
            cpl_propertylist_delete(qcs_fpet_out) ;
            cpl_table_delete(tw_fpet_out) ;
            cpl_table_delete(lines_diagnostics_fpet_out) ;
            cpl_table_delete(extracted_fpet_out) ;
            hdrl_image_delete(wl_map_fpet_out) ;
            if (plist_une_out != NULL) cpl_propertylist_delete(plist_une_out) ;
            if (tw_une_out != NULL) cpl_table_delete(tw_une_out) ;
            if (lines_diagnostics_une_out != NULL) 
                cpl_table_delete(lines_diagnostics_une_out) ;
            if (extracted_une_out != NULL) cpl_table_delete(extracted_une_out) ;
            if (wl_map_une_out != NULL) hdrl_image_delete(wl_map_une_out) ;
            cpl_msg_indent_less() ;
            cpl_msg_error(__func__, "Failed to load the plist") ;
            return -1 ;
        }
        if (qcs_fpet_out != NULL) {
            cpl_propertylist_append(plist_fpet_out, qcs_fpet_out) ;
            cpl_propertylist_delete(qcs_fpet_out) ;
        }
    } else {
        tw_fpet_out = NULL ;
        lines_diagnostics_fpet_out = NULL ;
        extracted_fpet_out = NULL ;
        wl_map_fpet_out = NULL ;
        plist_fpet_out = NULL ;
    }

    /* Return */
    *out_trace_wave_une = tw_une_out ;
    *lines_diagnostics_une = lines_diagnostics_une_out ;
    *out_extracted_une = extracted_une_out ;
    *out_wave_map_une = wl_map_une_out ;
    *ext_plist_une = plist_une_out ;

    *out_trace_wave_fpet = tw_fpet_out ;
    *lines_diagnostics_fpet = lines_diagnostics_fpet_out ;
    *out_extracted_fpet = extracted_fpet_out ;
    *out_wave_map_fpet = wl_map_fpet_out ;
    *ext_plist_fpet = plist_fpet_out ;
    return 0 ;
}
