/*

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
#include "cr2res_qc.h"
#include "cr2res_slit_curv.h"
#include "cr2res_extract.h"
#include "cr2res_trace.h"
#include "cr2res_wave.h"

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

#define RECIPE_STRING "cr2res_cal_wave"
#define CRIRES_PI       3.1415926535897932384626433832795029L

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_cal_wave_qc_une_flux(
        const cpl_table         **      extracted,
        const cpl_propertylist  **      ext_plist,
        double                          cwlen,
        cpl_propertylist        *       plist) ;
static int cr2res_cal_wave_qc_fpi(
        const cpl_table         **      extracted,
        const cpl_table         **      tw,
        const cpl_propertylist  **      ext_plist,
        double                          cwlen,
        cpl_propertylist        *       plist) ;
static int cr2res_cal_wave_qc_tilt(
        const cpl_table     **      tws,
        double                      cwlen,
        cpl_propertylist    *       plist) ;
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
        int                     slit_degree,
        int                     subtract_nolight_rows,
        int                     cosmics,
        cr2res_collapse         collapse,
        int                     ext_height,
        int                     ext_swath_width,
        int                     ext_oversample,
        double                  ext_smooth_slit,
        double                  ext_smooth_spec,
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
      If FPET RAW frame is present, compute the slit_curvature from it, \n\
      use it in the following, and store it in the out_trace_wave_fpet. \n\
                                                                        \n\
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

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.slit_degree",
            CPL_TYPE_INT, "Slit fitting Polynomial degree (1 or 2)",
            "cr2res.cr2res_cal_wave", 2);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "slit_degree");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.subtract_nolight_rows",
            CPL_TYPE_BOOL, "Subtract the no-light rows.",
            "cr2res.cr2res_cal_wave", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "subtract_nolight_rows");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.cosmics",
            CPL_TYPE_BOOL, "Find and mark cosmic rays hits as bad",
            "cr2res.cr2res_cal_wave", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "cosmics");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.collapse_method",
            CPL_TYPE_STRING, "Collapse the input images (MEAN or MEDIAN)",
            "cr2res.cr2res_cal_wave", "MEDIAN");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "collapse_method");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.extract_oversample",
            CPL_TYPE_INT, "factor by which to oversample the extraction",
            "cr2res.cr2res_cal_wave", 7);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_oversample");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.extract_swath_width",
            CPL_TYPE_INT, "The swath width", "cr2res.cr2res_cal_wave", 2048);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_swath_width");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.extract_height",
            CPL_TYPE_INT, "Extraction height",
            "cr2res.cr2res_cal_wave", 160);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_height");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.extract_smooth_slit",
            CPL_TYPE_DOUBLE, "Smoothing along the slit",
            "cr2res.cr2res_cal_wave", 15.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_smooth_slit");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.extract_smooth_spec",
            CPL_TYPE_DOUBLE, "Smoothing along the spectrum",
            "cr2res.cr2res_cal_wave", 2.0E-7);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_smooth_spec");
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
            "cr2res.cr2res_cal_wave", 0.04);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "wl_err");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.wl_degree",
            CPL_TYPE_INT, "Wavelength Polynomial degree",
            "cr2res.cr2res_cal_wave", 0);
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
            "cr2res.cr2res_cal_wave", TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "keep");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_wave.save_intermediate",
            CPL_TYPE_BOOL,
            "Flag to save UNE results (if UNE and FPET are used as inputs)",
            "cr2res.cr2res_cal_wave", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "save");
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
                            slit_degree, ext_oversample, ext_swath_width, 
                            ext_height, wl_degree, display, log_flag,
                            fallback_input_wavecal_flag, 
                            save_intermediate_flag, keep_higher_degrees_flag, 
                            clean_spectrum, subtract_nolight_rows,
                            cosmics;
    double                  ext_smooth_slit, ext_smooth_spec, wl_start, wl_end,
                            wl_err, wl_shift, display_wmin,
                            display_wmax, central_wlen;
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
    cpl_propertylist    *   plist ;
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
    cpl_propertylist    *   qc_main ;
    int                     det_nr, i;

    /* Needed for sscanf() */
    setlocale(LC_NUMERIC, "C");

    /* Initialise */
    wl_start = wl_end = -1.0 ;

    collapse = CR2RES_COLLAPSE_UNSPECIFIED ;
    display_wmin = display_wmax = -1.0 ;

    int qc_sizes[] = {CR2RES_NB_DETECTORS,CR2RES_NB_DETECTORS};
    cpl_propertylist** qc_plists[] = {ext_plist_une, ext_plist_fpet} ;
    cpl_propertylist** qc_plists_une[] = {ext_plist_une} ;

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
            "cr2res.cr2res_cal_wave.slit_degree");
    slit_degree = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_wave.subtract_nolight_rows");
    subtract_nolight_rows = cpl_parameter_get_bool(param) ;
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_wave.cosmics");
    cosmics = cpl_parameter_get_bool(param);
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
           "cr2res.cr2res_cal_wave.extract_oversample");
    ext_oversample = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_wave.extract_swath_width");
    ext_swath_width = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_wave.extract_height");
    ext_height = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_wave.extract_smooth_slit");
    ext_smooth_slit = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_wave.extract_smooth_spec");
    ext_smooth_spec = cpl_parameter_get_double(param);
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
            "cr2res.cr2res_cal_wave.save_intermediate");
    save_intermediate_flag = cpl_parameter_get_bool(param) ;
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
    if (slit_degree != 1 && slit_degree != 2) {
        cpl_msg_error(__func__, "The slit fit degree must be 1 or 2");
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
        cpl_frameset_delete(rawframes_une) ;
        if (rawframes_fpet!=NULL) cpl_frameset_delete(rawframes_fpet) ;
        cpl_msg_error(__func__, "The emission lines catalog is needed");
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
    if (reduce_order > -1 && wavecal_type == CR2RES_ETALON) {
        if (rawframes_une !=NULL) cpl_frameset_delete(rawframes_une) ;
        if (rawframes_fpet!=NULL) cpl_frameset_delete(rawframes_fpet) ;
        cpl_msg_error(__func__, "Limiting to one order with ETALON impossible");
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

        /* For Testing the QCs without re-running all the time */
        if (0) {
            out_trace_wave_une[det_nr-1] = 
                cr2res_io_load_TRACE_WAVE(
"/home/yjung/P_cr2res/test_cr2res_cal_wave/res/cr2res_cal_wave_tw_une.fits", 
                det_nr) ;
            out_extracted_une[det_nr-1] = cr2res_io_load_EXTRACT_1D(
"/home/yjung/P_cr2res/test_cr2res_cal_wave/res/cr2res_cal_wave_extracted_une.fits", 
                det_nr) ;
            ext_plist_une[det_nr-1] = cpl_propertylist_load(
"/home/yjung/P_cr2res/test_cr2res_cal_wave/res/cr2res_cal_wave_extracted_une.fits", 
                det_nr) ;
            out_trace_wave_fpet[det_nr-1] = 
                cr2res_io_load_TRACE_WAVE(
"/home/yjung/P_cr2res/test_cr2res_cal_wave/res/cr2res_cal_wave_tw_fpet.fits", 
                det_nr) ;
            out_extracted_fpet[det_nr-1] = 
                cr2res_io_load_EXTRACT_1D(
"/home/yjung/P_cr2res/test_cr2res_cal_wave/res/cr2res_cal_wave_extracted_fpet.fits", 
                det_nr) ;
            ext_plist_fpet[det_nr-1] = 
                cpl_propertylist_load(
"/home/yjung/P_cr2res/test_cr2res_cal_wave/res/cr2res_cal_wave_extracted_fpet.fits", 
                det_nr) ;
            lines_diagnostics_une[det_nr-1] = NULL ;
            out_wave_map_une[det_nr-1] = NULL ;
            out_wave_map_fpet[det_nr-1] = NULL ;
            lines_diagnostics_fpet[det_nr-1] = NULL ;
        } else {

        /* Call the reduction function */
        if (cr2res_cal_wave_reduce(rawframes_une, rawframes_fpet, detlin_frame,
                    master_dark_frame, master_flat_frame, bpm_frame,
                    trace_wave_frame, lines_frame, det_nr, reduce_order,
                    reduce_trace, slit_degree, subtract_nolight_rows, cosmics,
                    collapse, ext_height, ext_swath_width, ext_oversample, 
                    ext_smooth_slit, ext_smooth_spec, wavecal_type, wl_degree, 
                    wl_start, wl_end, wl_err, wl_shift, log_flag, 
                    fallback_input_wavecal_flag, keep_higher_degrees_flag, 
                    clean_spectrum, display, display_wmin, 
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
        }

        cpl_msg_indent_less() ;
    }

    /* QC Computation */
    qc_main = cpl_propertylist_new() ;

    cr2res_qc_calculate_mean_and_rmsd(qc_plists, 2, qc_sizes, 
                             CR2RES_HEADER_QC_WAVE_CENTWL, qc_main,
                             CR2RES_HEADER_QC_WAVE_CENTWL_AVG, CR2RES_HEADER_QC_WAVE_CENTWL_RMS);

    cr2res_qc_calculate_mean_and_rmsd(qc_plists, 2, qc_sizes, 
                             CR2RES_HEADER_QC_WAVE_DISPWL, qc_main,
                             CR2RES_HEADER_QC_WAVE_DISPWL_AVG, CR2RES_HEADER_QC_WAVE_DISPWL_RMS);

    char* qc_une_keywords[] = {
    CR2RES_HEADER_QC_WAVE_LAMP_EFFIC,
    CR2RES_HEADER_QC_WAVE_RESOL,
    CR2RES_HEADER_QC_WAVE_RESOL_FWHM,
    CR2RES_HEADER_QC_WAVE_POS,
    CR2RES_HEADER_QC_WAVE_BESTXCORR,
    CR2RES_HEADER_QC_OVEREXPOSED
};


    int qc_une_keywords_size = sizeof(qc_une_keywords) / sizeof(qc_une_keywords[0]);

    if (rawframes_une != NULL) {
        plist = cpl_propertylist_load(cpl_frame_get_filename(
                    cpl_frameset_get_position_const(rawframes_une, 0)), 0) ;
        central_wlen = cr2res_pfits_get_cwlen(plist) ;
        cpl_propertylist_delete(plist);

        /* QC.UNE.FLUX  Computation */
        if (cr2res_cal_wave_qc_une_flux(
                    (const cpl_table **)out_extracted_une,
                    (const cpl_propertylist **)ext_plist_une, 
                    central_wlen, qc_main)) {
            cpl_msg_warning(__func__, "QC.UNE.FLUX Computation failed") ;
        }

        /* Should get min and max order_id's from column of all detectors here*/
        /*int nb_traces = cpl_table_get_nrow(out_trace_wave_une[0]);*/

        int min_order = INT_MAX;
        int max_order = INT_MIN;
        int current_min, current_max;

        for (int qc_i = 0; qc_i < CR2RES_NB_DETECTORS; qc_i++) {
            if (out_trace_wave_une[qc_i] == NULL) {
                continue;
            }
            current_min = cpl_table_get_column_min(out_trace_wave_une[qc_i], "Order");
            if (cpl_error_get_code() != CPL_ERROR_NONE) {
                continue;
            }
            if (current_min < min_order) {
                min_order = current_min;
            }

           current_max = cpl_table_get_column_max(out_trace_wave_une[qc_i], "Order");
            if (cpl_error_get_code() != CPL_ERROR_NONE) {
               continue;
            }
            if (current_max > max_order) {
                max_order = current_max;
            }
        }

        int qc_j;
        for (qc_j=0; qc_j<qc_une_keywords_size ; qc_j++) {
            for (int order_id=min_order ; order_id<=max_order ; order_id++) {

                int qc_sizes_une[] = {CR2RES_NB_DETECTORS};
                char * qc_ref, * qc_res_avg, * qc_res_rmsd ;

                if(strcmp(qc_une_keywords[qc_j], CR2RES_HEADER_QC_OVEREXPOSED)==0) {
                    qc_ref = cpl_sprintf("%s%02d",
                            qc_une_keywords[qc_j], order_id);
                    qc_res_avg = cpl_sprintf("%s%02d %s",
                            qc_une_keywords[qc_j], order_id, "AVG");
                    qc_res_rmsd = cpl_sprintf("%s%02d %s",
                            qc_une_keywords[qc_j], order_id, "RMS");;
                }
                else {
                    qc_ref = cpl_sprintf("%s-%02d-%02d",
                            qc_une_keywords[qc_j], order_id, 1);
                    qc_res_avg = cpl_sprintf("%s-%02d-%02d %s",
                            qc_une_keywords[qc_j], order_id, 1, "AVG");
                    qc_res_rmsd = cpl_sprintf("%s-%02d-%02d %s",
                            qc_une_keywords[qc_j], order_id, 1, "RMS");
                }

                cr2res_qc_calculate_mean_and_rmsd(qc_plists_une, 1, qc_sizes_une, 
                                 qc_ref, qc_main,
                                 qc_res_avg, qc_res_rmsd);

                cpl_free(qc_ref);
                cpl_free(qc_res_avg);
                cpl_free(qc_res_rmsd);
            }
        }

    if (cpl_error_get_code() != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Failed to set QC params for UNE: %s",
                cpl_error_get_message_default ( cpl_error_get_code())) ;
        return -1 ;
    }

    }

    if (rawframes_fpet != NULL) {
        plist = cpl_propertylist_load(cpl_frame_get_filename(
                    cpl_frameset_get_position_const(rawframes_fpet, 0)), 0) ;
        central_wlen = cr2res_pfits_get_cwlen(plist) ;
        cpl_propertylist_delete(plist);

        /* QC.TILT Computation */
        if (cr2res_cal_wave_qc_tilt((const cpl_table **)out_trace_wave_fpet,
                    central_wlen, qc_main)) {
            cpl_msg_warning(__func__, "QC.TILT Computation failed") ;
        }
        /* QC.FPI Computation */
        if (cr2res_cal_wave_qc_fpi(
                    (const cpl_table **)out_extracted_fpet,
                    (const cpl_table **)out_trace_wave_fpet,
                    (const cpl_propertylist **)ext_plist_fpet, 
                    central_wlen, qc_main)) {
            cpl_msg_warning(__func__, "QC.FPI Computation failed") ;
        }
    }

    /* Save only the used RAW ? : rawframes_xxx instead of 2nd frameset */
    /* Beware that the calibration PRO RECi CAL will be missing */ 

    /* Save Products UNE */
    if (rawframes_une != NULL &&
            (rawframes_fpet == NULL ||   // UNE is the main product
             save_intermediate_flag)) { 
        out_file = cpl_sprintf("%s_tw_une.fits", RECIPE_STRING) ;
        cr2res_io_save_TRACE_WAVE(out_file, frameset, frameset, parlist, 
                out_trace_wave_une, qc_main, ext_plist_une, 
                CR2RES_CAL_WAVE_TW_PROCATG, RECIPE_STRING) ;
        cpl_free(out_file);

        out_file = cpl_sprintf("%s_wave_map_une.fits", RECIPE_STRING) ;
        cr2res_io_save_WAVE_MAP(out_file, frameset, frameset, parlist, 
                out_wave_map_une, qc_main, ext_plist_une, 
                CR2RES_CAL_WAVE_MAP_PROCATG, RECIPE_STRING) ;
        cpl_free(out_file);

        out_file = cpl_sprintf("%s_extracted_une.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, frameset, parlist, 
                out_extracted_une, qc_main, ext_plist_une, 
                CR2RES_CAL_WAVE_EXTRACT_1D_PROCATG, RECIPE_STRING) ;
        cpl_free(out_file);

        out_file = cpl_sprintf("%s_lines_diagnostics_une.fits", 
                RECIPE_STRING);
        cr2res_io_save_LINES_DIAGNOSTICS(out_file, frameset, frameset, 
                parlist, lines_diagnostics_une, qc_main, ext_plist_une,
                CR2RES_CAL_WAVE_LINES_DIAGNOSTICS_PROCATG, RECIPE_STRING) ;
        cpl_free(out_file);
    }

    if (rawframes_fpet != NULL) {
        /* Save Products FPET */
        out_file = cpl_sprintf("%s_tw_fpet.fits", RECIPE_STRING) ;
        cr2res_io_save_TRACE_WAVE(out_file, frameset, frameset, parlist, 
                out_trace_wave_fpet, qc_main, ext_plist_fpet, 
                CR2RES_CAL_WAVE_TW_PROCATG, RECIPE_STRING) ;
        cpl_free(out_file);

        out_file = cpl_sprintf("%s_wave_map_fpet.fits", RECIPE_STRING) ;
        cr2res_io_save_WAVE_MAP(out_file, frameset, frameset, parlist, 
                out_wave_map_fpet, qc_main, ext_plist_fpet, 
                CR2RES_CAL_WAVE_MAP_PROCATG, RECIPE_STRING) ;
        cpl_free(out_file);

        out_file = cpl_sprintf("%s_extracted_fpet.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, frameset, parlist, 
                out_extracted_fpet, qc_main, ext_plist_fpet, 
                CR2RES_CAL_WAVE_EXTRACT_1D_PROCATG, RECIPE_STRING) ;
        cpl_free(out_file);

        out_file = cpl_sprintf("%s_lines_diagnostics_fpet.fits", 
                RECIPE_STRING);
        cr2res_io_save_LINES_DIAGNOSTICS(out_file, frameset, frameset, parlist,
                lines_diagnostics_fpet, qc_main, ext_plist_fpet,
                CR2RES_CAL_WAVE_LINES_DIAGNOSTICS_PROCATG, RECIPE_STRING) ;
        cpl_free(out_file);
    }

    /* Free and return */
    cpl_propertylist_delete(qc_main) ;
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
  @param slit_degree        The polynomial degree for the slit fit 
  @param subtract_nolight_rows
  @param cosmics            Flag to correct for cosmics
  @param collapse           CR2RES_COLLAPSE_MEAN or CR2RES_COLLAPSE_MEDIAN
  @param ext_height         Extraction related
  @param ext_swath_width    Extraction related
  @param ext_oversample     Extraction related
  @param ext_smooth_slit    Extraction related
  @param ext_smooth_spec    Extraction related
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
        int                     slit_degree,
        int                     subtract_nolight_rows,
        int                     cosmics,
        cr2res_collapse         collapse,
        int                     ext_height, 
        int                     ext_swath_width,
        int                     ext_oversample, 
        double                  ext_smooth_slit,
        double                  ext_smooth_spec,
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
    cpl_vector          *   ndits_fpet ;
    hdrl_imagelist      *   in_fpet ;
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
    cpl_polynomial      *   slit_polya ;
    cpl_polynomial      *   slit_polyb ;
    cpl_polynomial      *   slit_polyc ;
    double                  gain, error_factor ;
    int                     i, order, trace_id, ext_nr, 
                            zp_order_fpet, nb_traces,
                            grat1_order_fpet ;
    
    /* TODO, make parameters */
    int extract_niter = 10;
    double extract_kappa = 10;

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

    /* Get the Gain */
    if (reduce_det == 1) gain = CR2RES_GAIN_CHIP1 ;
    else if (reduce_det == 2) gain = CR2RES_GAIN_CHIP2 ;
    else if (reduce_det == 3) gain = CR2RES_GAIN_CHIP3 ;
    else {
        cpl_msg_error(__func__, "Failed to get the Gain value") ;
        return -1 ;
    }


    /* Compute the Slit Curvature using the first passed FPET frame */
    tw_in = NULL ;
    if (rawframes_fpet != NULL) {

        const cpl_frame     *   fpet_frame ;
        const char          *   fpet_fname ;
        hdrl_image          *   fpet_image ;

        fpet_frame = cpl_frameset_get_position_const(rawframes_fpet, 0) ;
        fpet_fname = cpl_frame_get_filename(fpet_frame) ;
        cpl_msg_info(__func__, 
                "Compute the Slit Curvature from the first FPET frame %s",
                fpet_fname) ;
        cpl_msg_indent_more() ;

        /* Load the trace wave */
        cpl_msg_info(__func__, "Load the TRACE WAVE") ;
        if ((tw_in = cr2res_io_load_TRACE_WAVE(cpl_frame_get_filename(
                            trace_wave_frame), reduce_det)) == NULL) {
            cpl_msg_error(__func__, "Failed to Load the traces file") ;
            cpl_msg_indent_less() ;
            return -1 ;
        }

        /* Load the FPET image of this detector */
        cpl_msg_info(__func__, "Load the FPET image") ;
        if ((fpet_image = cr2res_io_load_image(fpet_fname, reduce_det))==NULL) {
            cpl_msg_error(__func__, "Failed to load image") ;
            cpl_table_delete(tw_in) ;
            cpl_msg_indent_less() ;
            return -1 ;
        }

        /* Loop over the traces and get the slit curvature */
        nb_traces = cpl_table_get_nrow(tw_in) ;
        for (i=0 ; i<nb_traces ; i++) {
			/* Get Order and trace id */
			order = cpl_table_get(tw_in, CR2RES_COL_ORDER, i, NULL) ;
			trace_id = cpl_table_get(tw_in, CR2RES_COL_TRACENB, i, NULL) ;

            /* Check if this order needs to be skipped */
            if (reduce_order > -1 && order != reduce_order) continue ;

            /* Check if this trace needs to be skipped */
            if (reduce_trace > -1 && trace_id != reduce_trace) continue ;

            cpl_msg_info(__func__, "Process Order %d/Trace %d",order,trace_id) ;
            cpl_msg_indent_more() ;
    
            /* Call the Slit Curvature Computation */
            /* TODO : Should those become parameters ? */
            int height = 100 ;
            int window = 15 ;
            int change_degree = 1;
            if (cr2res_slit_curv_compute_order_trace(fpet_image, tw_in, order, 
                        trace_id, height, window, change_degree, slit_degree, 
                        &slit_polya, &slit_polyb, &slit_polyc)) {
                cpl_msg_error(__func__, "Failed to compute slit curvature") ;
                cpl_table_delete(tw_in) ;
                hdrl_image_delete(fpet_image) ;
                cpl_msg_indent_less() ;
                return -1 ;
            }

            cpl_array           *   slit_array ;
            /* Fill the SLIT_CURVE_A/B/C for the current trace */
            slit_array = cr2res_convert_poly_to_array(slit_polya, 3) ;
            cpl_polynomial_delete(slit_polya) ;
            cpl_table_set_array(tw_in, CR2RES_COL_SLIT_CURV_A, i, slit_array) ;
            cpl_array_delete(slit_array) ;
            slit_array = cr2res_convert_poly_to_array(slit_polyb, 3) ;
            cpl_polynomial_delete(slit_polyb) ;
            cpl_table_set_array(tw_in, CR2RES_COL_SLIT_CURV_B, i, slit_array) ;
            cpl_array_delete(slit_array) ;
            slit_array = cr2res_convert_poly_to_array(slit_polyc, 3) ;
            cpl_polynomial_delete(slit_polyc) ;
            cpl_table_set_array(tw_in, CR2RES_COL_SLIT_CURV_C, i, slit_array) ;
            cpl_array_delete(slit_array) ;

            cpl_msg_indent_less() ;
        }
        hdrl_image_delete(fpet_image) ;
        cpl_msg_indent_less() ;
    }

    /*
        At this point : tw_in is NULL or it is the input TW with the
        newly computed slit curvature
    */

    /* Reduce the UNE */
    if (rawframes_une != NULL) {

        cpl_vector          *   dits_une ;
        cpl_vector          *   ndits_une ;
        hdrl_imagelist      *   in_une ;
        hdrl_imagelist      *   in_une_calib ;
        hdrl_image          *   first_image ;
        int      zp_order_une, grat1_order_une ;

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
        
        /*Load the NDITs */
        ndits_une = cr2res_io_read_ndits(rawframes_une);

        /* Set the error factor. */
        error_factor = gain * cpl_vector_get(ndits_une, 0) *
                                    cpl_frameset_get_size(rawframes_une) ;

        for (i=0; i<cpl_vector_get_size(ndits_une); i++){
            if (cpl_vector_get(ndits_une,i) != cpl_vector_get(ndits_une, 0))
                cpl_msg_warning(__func__, "UNE raw frames have different NDIT! "
                    "Error spectrum will likely be scaled incorrectly.");
        }

        /* Load UNE image list */
        if ((in_une = cr2res_io_load_image_list_from_set(rawframes_une,
                        reduce_det)) == NULL) {
            cpl_msg_error(__func__, "Cannot load images") ;
            if (dits_une != NULL) cpl_vector_delete(dits_une) ;
            if (tw_in != NULL) cpl_table_delete(tw_in); 
            cpl_msg_indent_less() ;
            return -1 ;
        }
        if (hdrl_imagelist_get_size(in_une) !=
                cpl_frameset_get_size(rawframes_une)) {
            cpl_msg_error(__func__, "Inconsistent number of loaded images") ;
            if (dits_une != NULL) cpl_vector_delete(dits_une) ;
            hdrl_imagelist_delete(in_une) ;
            if (tw_in != NULL) cpl_table_delete(tw_in); 
            cpl_msg_indent_less() ;
            return -1 ;
        }
        /* Calibrate the UNE images */
        if ((in_une_calib = cr2res_calib_imagelist(in_une, reduce_det, 0,
                        subtract_nolight_rows, 1, cosmics, master_flat_frame, 
                        master_dark_frame, bpm_frame, detlin_frame, 
                        dits_une, ndits_une)) == NULL) {
            cpl_msg_error(__func__, "Failed to apply the calibrations") ;
            if (dits_une != NULL) cpl_vector_delete(dits_une) ;
            if (ndits_une != NULL) cpl_vector_delete(ndits_une) ;
            hdrl_imagelist_delete(in_une) ;
            if (tw_in != NULL) cpl_table_delete(tw_in); 
            cpl_msg_indent_less() ;
            return -1 ;
        }
        hdrl_imagelist_delete(in_une) ;
        if (dits_une != NULL) cpl_vector_delete(dits_une) ;
        if (ndits_une != NULL) cpl_vector_delete(ndits_une) ;

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
            cpl_msg_error(__func__, "Failed Collapse: %d",
                    cpl_error_get_code());
            cpl_msg_indent_less() ;
            if (tw_in != NULL) cpl_table_delete(tw_in); 
            return -1 ;
        }

        /* Load the trace wave */
        cpl_msg_info(__func__, "Load the TRACE WAVE") ;
        if (tw_in == NULL) 
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
        if (cr2res_extract_traces(collapsed_une, tw_in, NULL, NULL, 0,
                    reduce_order, reduce_trace, CR2RES_EXTR_OPT_CURV, 
                    ext_height, ext_swath_width, ext_oversample, 
                    ext_smooth_slit, ext_smooth_spec,
                    extract_niter, extract_kappa, error_factor, 0, 0, 0,
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
                    reduce_trace, wavecal_type, wl_degree, -1,  // xdegree
                    wl_start, wl_end, 
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
			cpl_table_delete(tw_in) ;
            cpl_table_delete(tw_une_out) ;
            cpl_table_delete(lines_diagnostics_une_out) ;
            cpl_table_delete(extracted_une_out) ;
            hdrl_image_delete(wl_map_une_out) ;
            cpl_msg_error(__func__, "Failed to load the plist") ;
            cpl_msg_indent_less() ;
            return -1 ;
        }

        /* Load the first RAW image and compute QC OVEREXPOSED */
        first_image = cr2res_io_load_image(first_file,reduce_det) ;

		/* QC.OVEREXPOSED */
		/* Loop on the traces */
        nb_traces = cpl_table_get_nrow(tw_in) ;
		for (i=0 ; i<nb_traces ; i++) {
			/* Get Order and trace id */
			order = cpl_table_get(tw_in, CR2RES_COL_ORDER, i, NULL) ;
			trace_id = cpl_table_get(tw_in, CR2RES_COL_TRACENB, i, NULL) ;

			/* Check if this order needs to be skipped */
			if (reduce_order > -1 && order != reduce_order) continue ;

			/* Check if this trace needs to be skipped */
			if (reduce_trace > -1 && trace_id != reduce_trace) continue ;

            double qc_overexposed;
			qc_overexposed = cr2res_qc_overexposed(
					hdrl_image_get_image(first_image), tw_in, order) ;
			char * qc_name = cpl_sprintf("%s%02d",
					CR2RES_HEADER_QC_OVEREXPOSED, order) ;
			if (qcs_une_out == NULL) qcs_une_out = cpl_propertylist_new() ;
			cpl_propertylist_append_double(qcs_une_out, qc_name,qc_overexposed);
			cpl_free(qc_name) ;
        }
        hdrl_image_delete(first_image) ;

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

    /*
        At this point : tw_in is NULL or it is the input TW with the
        newly computed slit curvature
    */

    /* Reduce the FPET */
    if (rawframes_fpet != NULL && reduce_order != -1){
        cpl_msg_warning(__func__, 
            "Etalon wavecal requires all orders to be used, skipping it.");
        if (tw_in != NULL) cpl_table_delete(tw_in); 
        tw_fpet_out = NULL ;
        lines_diagnostics_fpet_out = NULL ;
        extracted_fpet_out = NULL ;
        wl_map_fpet_out = NULL ;
        plist_fpet_out = NULL ;
    } else if (rawframes_fpet != NULL) {
        cpl_msg_info(__func__, "Reduce %"CPL_SIZE_FORMAT" FPET Frames",
                cpl_frameset_get_size(rawframes_fpet)) ;
        cpl_msg_indent_more() ;

        /* Load FPET image list */
        if ((in_fpet = cr2res_io_load_image_list_from_set(rawframes_fpet,
                        reduce_det)) == NULL) {
            cpl_msg_error(__func__, "Cannot load images") ;
            cpl_msg_indent_less() ;
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
            cpl_msg_indent_less() ;
            hdrl_imagelist_delete(in_fpet) ;
            if (plist_une_out != NULL) cpl_propertylist_delete(plist_une_out) ;
            if (tw_une_out != NULL) cpl_table_delete(tw_une_out) ;
            if (lines_diagnostics_une_out != NULL) 
                cpl_table_delete(lines_diagnostics_une_out) ;
            if (extracted_une_out != NULL) cpl_table_delete(extracted_une_out) ;
            if (wl_map_une_out != NULL) hdrl_image_delete(wl_map_une_out) ;
            return -1 ;
        }

        /* Load FPET NDITs*/
        ndits_fpet = cr2res_io_read_ndits(rawframes_fpet);

        /* Set the error factor. */
        error_factor = gain * cpl_vector_get(ndits_fpet, 0) *
                                    cpl_frameset_get_size(rawframes_une) ;

        for (i=0; i<cpl_vector_get_size(ndits_fpet); i++){
            if (cpl_vector_get(ndits_fpet,i) != cpl_vector_get(ndits_fpet, 0))
                cpl_msg_warning(__func__, "FPET raw frames have different NDIT! "
                    "Error spectrum will likely be scaled incorrectly.");
        }

        /* Calibrate the FPET images */
        if ((in_fpet_calib = cr2res_calib_imagelist(in_fpet, reduce_det, 0, 
                        subtract_nolight_rows, 0, cosmics, master_flat_frame, 
                        NULL, bpm_frame, detlin_frame, NULL, 
                        ndits_fpet)) == NULL) {
            cpl_msg_error(__func__, "Failed to apply the calibrations") ;
            cpl_msg_indent_less() ;
            if (in_fpet != NULL) hdrl_imagelist_delete(in_fpet) ;
            if (plist_une_out != NULL) cpl_propertylist_delete(plist_une_out) ;
            if (tw_une_out != NULL) cpl_table_delete(tw_une_out) ;
            if (lines_diagnostics_une_out != NULL) 
                cpl_table_delete(lines_diagnostics_une_out) ;
            if (extracted_une_out != NULL) cpl_table_delete(extracted_une_out) ;
            if (wl_map_une_out != NULL) hdrl_image_delete(wl_map_une_out) ;
            if (ndits_fpet != NULL) cpl_vector_delete(ndits_fpet) ;
            return -1 ;
        }
        hdrl_imagelist_delete(in_fpet) ;
        if (ndits_fpet != NULL) cpl_vector_delete(ndits_fpet) ;

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
            if (tw_in != NULL) cpl_table_delete(tw_in); 
            tw_in = cpl_table_duplicate(tw_une_out) ;
        } else {
            cpl_msg_info(__func__, "Load the TRACE WAVE") ;
            if (tw_in == NULL)
                if ((tw_in = cr2res_io_load_TRACE_WAVE(cpl_frame_get_filename(
                                    trace_wave_frame), reduce_det)) == NULL) {
                    cpl_msg_error(__func__, "Failed to Load the traces file") ;
                    if (collapsed_fpet != NULL) 
                        hdrl_image_delete(collapsed_fpet) ;
                    if (plist_une_out != NULL) 
                        cpl_propertylist_delete(plist_une_out) ;
                    if (tw_une_out != NULL) cpl_table_delete(tw_une_out) ;
                    if (lines_diagnostics_une_out != NULL) 
                        cpl_table_delete(lines_diagnostics_une_out) ;
                    if (extracted_une_out != NULL) 
                        cpl_table_delete(extracted_une_out) ;
                    if (wl_map_une_out != NULL) 
                        hdrl_image_delete(wl_map_une_out) ;
                    cpl_msg_indent_less() ;
                    return -1 ;
                }
        }

        /* Execute the extraction for FPET */
        cpl_msg_info(__func__, "Spectra Extraction FPET") ;
        cpl_msg_indent_more() ;
        if (cr2res_extract_traces(collapsed_fpet, tw_in, NULL, NULL, 0,
                    reduce_order, reduce_trace, CR2RES_EXTR_OPT_CURV, 
                    ext_height, ext_swath_width, ext_oversample, 
                    ext_smooth_slit, ext_smooth_spec, 
                    extract_niter, extract_kappa, error_factor, 0, 0, 0,
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
        cpl_msg_info(__func__, "Running FPET with degrees=4,6") ;

        if (cr2res_wave_apply(tw_in, extracted_fpet, NULL, reduce_order, 
                    reduce_trace, CR2RES_ETALON, 4, 6, wl_start, wl_end, 
                    wl_err, wl_shift, log_flag, fallback_input_wavecal_flag, 
                    FALSE, // --keep=FALSE 
                    clean_spectrum, 
                    display, display_wmin, display_wmax, zp_order_fpet,
                    grat1_order_fpet,
                    &qcs_fpet_out,
                    &lines_diagnostics_fpet_out,
                    &extracted_fpet_out,
                    &tw_fpet_out) || cpl_error_get_code()) {
            cpl_msg_error(__func__, "Failed to calibrate");
            cpl_table_delete(extracted_fpet) ;
            cpl_table_delete(tw_in) ;
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
        cpl_msg_indent_less() ;
    } else {
        if (tw_in != NULL) cpl_table_delete(tw_in); 
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

/*----------------------------------------------------------------------------*/
/**
  @brief Compute QC.UNE.FLUX
  @param    extracted   3 extracted tables (1 per detector)
  @param    ext_plist   3 Extension headers (1 per detector)
  @param    cwlen
  @param    plist       Header for holding the QC values
  @return   0 if ok
  
  Use Detector 1 if CWLEN (from input file headers) is < 1100.0,
  detector 2 otherwise.
  The order that is the closest to the center of the detector is used.

  QC.UNE.FLUX is the average of the lines between 1000 and 37000 intensity

 */
/*----------------------------------------------------------------------------*/
static int cr2res_cal_wave_qc_une_flux(
        const cpl_table         **      extracted,
        const cpl_propertylist  **      ext_plist,
        double                          cwlen,
        cpl_propertylist        *       plist)
{
    const cpl_table         *   ref_extr ;
    const cpl_propertylist  *   ref_plist ;
    cpl_bivector            *   spec ;
    cpl_bivector            *   spec_err ;
    cpl_bivector            *   intensities ;
    double                  *   pintens;
    double                      qc_flux, min_intens, max_intens;
    int                         ref_det, corder_idx ;
	cpl_size					k;

    /* Check Inputs */
    if (extracted==NULL || ext_plist == NULL) return -1 ;
    for (k=0 ; k< CR2RES_NB_DETECTORS ; k++) {
        if (extracted[k] == NULL) return -1 ;
        if (ext_plist[k] == NULL) return -1 ;
    }

    /* Initialize */
    qc_flux = -1.0 ;
    min_intens = 1000.0 ;
    max_intens = 37000.0 ;

    /* Detector used depends on the cwlen */
    if (cwlen <= 1100.0)    ref_det = 0 ;
    else                    ref_det = 1 ;
    ref_extr = extracted[ref_det] ;

    ref_plist = ext_plist[ref_det] ;

    /* Get central order index */
    corder_idx = cr2res_pfits_get_order_idx(ref_plist,
            CR2RES_DETECTOR_SIZE/2.0) ;

    cpl_msg_debug(__func__, "Central wl: %g  / Cent. Order IDX:  %d", 
            cwlen, corder_idx) ;

    /* Load the spectrum */
    cr2res_extract_EXTRACT1D_get_spectrum(ref_extr, corder_idx, 1, &spec,
            &spec_err) ;
    cpl_bivector_delete(spec_err) ;

    /* Plot the spectrum */
    if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
        cpl_errorstate status;
        status = cpl_errorstate_get();
        cpl_plot_bivector(
                "set grid;set xlabel 'Wavelength (nm)';set ylabel 'Intensity';",
                "t 'UNE Spectrum' w lines", "",
                spec);
        cpl_error_reset();
        cpl_errorstate_set(status);
    }

    /* Compute the lines statistics */
    if ((intensities = cr2res_qc_lines_intens_bgd(spec)) == NULL) {
        cpl_msg_warning(__func__, "Cannot get the lines statistics") ;
    }
    else {
        cpl_size nintensities, nval;
        double * pbgs ;
        nintensities = cpl_bivector_get_size(intensities) ;
        pintens = cpl_bivector_get_x_data(intensities) ;
        pbgs = cpl_bivector_get_y_data(intensities) ;
        /* Compute qc_flux of lines in [min_intens,max_intens] */
        nval = 0 ;
        for (k=0 ; k<nintensities ; k++) {
            double val;
            val = pintens[k] - pbgs[k] ;
            if (val < max_intens && val > min_intens) {
                qc_flux += val ;
                nval++;
            }
        }
        if (nval > 0) qc_flux /= nval ;
        cpl_bivector_delete(intensities) ;
    }
    cpl_bivector_delete(spec) ;

    /* Store the QCs */
    cpl_propertylist_append_double(plist, CR2RES_HEADER_QC_UNE_FLUX, 
            qc_flux) ;

    return 0 ;
}
/*----------------------------------------------------------------------------*/
/**
  @brief Compute QC.FPI CONTRAST and SEPARATION
  @param    extracted   3 extracted tables (1 per detector)
  @param    tw          3 TW tables (1 per detector)
  @param    ext_plist   3 Extension headers (1 per detector)
  @param    cwlen
  @param    plist       Header for holding the QC values
  @return   0 if ok
  
  Use Detector 1 if CWLEN (from input file headers) is < 1100.0,
  detector 2 otherwise.
  The order that is the closest to the center of the detector is used.

  QC.FPI.CONTRAST
    Detect all fringes, the peak/valley values, their positions in pix and wl
    -> average of this value (computed for each fringe) :
        (peak_wl-valley_wl) / (peak_wl+valley_wl) 

  QC.FPI.SEPARATION
    First and Last complete fringes:    
        -> (last_wl - first_wl)/NumberOfFringesInBetween
 */
/*----------------------------------------------------------------------------*/
static int cr2res_cal_wave_qc_fpi(
        const cpl_table         **      extracted,
        const cpl_table         **      tw,
        const cpl_propertylist  **      ext_plist,
        double                          cwlen,
        cpl_propertylist        *       plist)
{
    const cpl_table         *   ref_extr ;
    const cpl_table         *   ref_tw ;
    const cpl_propertylist  *   ref_plist ;
    double                      qc_contrast, qc_flux, qc_separation;
    int                         corder_idx, ref_det, first_fringe_idx,
                                last_fringe_idx, nb_fringes_between ;
    cpl_bivector            *   spec ;
    cpl_bivector            *   spec_err ;
    cpl_vector              *   xpeaks ;
    cpl_vector              *   lpeaks ;
    cpl_vector              *   fpeaks ;
    cpl_vector              *   xmins ;
    cpl_vector              *   fmins ;
    cpl_vector              *   lmins ;
    cpl_vector              *   tmp_vec ;
    cpl_vector              *   flux_vec ;
    cpl_polynomial          *   wlpoly ;
    cpl_size                    k, npeaks, end_val_search;

    

    /* Check Inputs */
    if (extracted==NULL || ext_plist == NULL) return -1 ;
    for (k=0 ; k< CR2RES_NB_DETECTORS ; k++) {
        if (extracted[k] == NULL) return -1 ;
        if (ext_plist[k] == NULL) return -1 ;
    }

    /* Detector used depends on the cwlen */
    if (cwlen <= 1100.0)    ref_det = 0 ;
    else                    ref_det = 1 ;
    ref_extr = extracted[ref_det] ;
    ref_tw = tw[ref_det] ;
    ref_plist = ext_plist[ref_det] ;

    /* Get central order index */
    corder_idx = cr2res_pfits_get_order_idx(ref_plist,
            CR2RES_DETECTOR_SIZE/2.0) ;

    cpl_msg_debug(__func__, "Central wl: %g  / Cent. Order IDX:  %d", 
            cwlen, corder_idx) ;

    /* Load the spectrum */
    cr2res_extract_EXTRACT1D_get_spectrum(ref_extr, corder_idx, 1, &spec,
            &spec_err) ;
    cpl_bivector_delete(spec_err) ;

    /* Plot the spectrum */
    if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
        cpl_errorstate status;

        status = cpl_errorstate_get();
        cpl_plot_bivector(
                "set grid;set xlabel 'Wavelength (nm)';set ylabel 'Intensity';",
                "t 'FPI Spectrum' w lines", "",
                spec);
        cpl_error_reset();
        cpl_errorstate_set(status);
    }

    /* Identify the peaks */
    if ((xpeaks=cr2res_wave_etalon_measure_fringes(cpl_bivector_get_y(spec))) 
            == NULL) {
        cpl_msg_warning(__func__, "Cannot identify the peaks") ;
        cpl_bivector_delete(spec) ;
        return -1 ;
    }
    npeaks = cpl_vector_get_size(xpeaks) ;

    if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
        cpl_vector_dump(xpeaks, stdout) ;
    }

    /* Get the wavelength solution for this order */
    wlpoly = cr2res_get_trace_wave_poly(ref_tw, CR2RES_COL_WAVELENGTH,
            corder_idx, 1) ;
    if (wlpoly == NULL) {
        cpl_msg_warning(__func__, "Cannot get the wavelength") ;
        cpl_bivector_delete(spec) ;
        cpl_vector_delete(xpeaks) ;
        return -1 ;
    }

    /* Get the wavelengths at the peak positions */
    lpeaks = cr2res_polynomial_eval_vector(wlpoly, xpeaks);

    /* Get the Min after each peak */
    xmins = cpl_vector_new(npeaks) ;
    fpeaks = cpl_vector_new(npeaks);
    fmins = cpl_vector_new(npeaks);
    for (k=0 ; k<npeaks ; k++) {
        cpl_size start_val_search, xminpos;

        /* Min of the peak */
        start_val_search = (int)cpl_vector_get(xpeaks, k) ;
        if (k == npeaks-1) {
            end_val_search = start_val_search ;// Not used
        } else {
            end_val_search = (int)cpl_vector_get(xpeaks, k+1) ;
        }
        tmp_vec = cpl_vector_extract(cpl_bivector_get_y(spec), 
                start_val_search-1, end_val_search-1, 1) ;
        xminpos = cpl_vector_get_minpos(tmp_vec) ;
        cpl_vector_set(fpeaks,k,cpl_vector_get(tmp_vec,0));
        cpl_vector_set(fmins,k,cpl_vector_get(tmp_vec,xminpos));
        cpl_vector_delete(tmp_vec) ;

        cpl_vector_set(xmins, k, start_val_search + xminpos) ;
    }
    
    /* Get the wavelengths at the min positions */
    lmins = cr2res_polynomial_eval_vector(wlpoly, xmins);
    cpl_polynomial_delete(wlpoly) ;

    for (k=0 ; k<npeaks ; k++) {
        cpl_msg_debug(__func__, 
                "Peak #%"CPL_SIZE_FORMAT" Pos: %g / WL: %g / Val :%g",
                k+1, 
                cpl_vector_get(xpeaks, k), 
                cpl_vector_get(lpeaks, k), 
                cpl_vector_get(cpl_bivector_get_y(spec), 
                    (int)cpl_vector_get(xpeaks, k))) ;
        cpl_msg_debug(__func__, 
                "        Min Pos: %g / WL: %g / Val :%g",
                cpl_vector_get(xmins, k), 
                cpl_vector_get(lmins, k), 
                cpl_vector_get(cpl_bivector_get_y(spec), 
                    (int)cpl_vector_get(xmins, k))) ;
    }
    cpl_bivector_delete(spec) ;
    cpl_vector_delete(xpeaks) ;
    cpl_vector_delete(xmins) ;

    /* CONTRAST */
    /* Ignore first and last peaks */
    tmp_vec = cpl_vector_new(npeaks-2) ;
    flux_vec = cpl_vector_new(npeaks-2) ;
    for (k=1 ; k<npeaks-1 ; k++) {
        double fpeak, fmin, contrast_val ;
        fpeak = cpl_vector_get(fpeaks, k) ;
        fmin = cpl_vector_get(fmins, k) ;
        contrast_val = (fpeak-fmin) / (fmin+fpeak) ;
        cpl_vector_set(tmp_vec, k-1, contrast_val) ;
        cpl_vector_set(flux_vec, k-1, fpeak) ;
    }
    qc_contrast = cpl_vector_get_mean(tmp_vec) ;
    qc_flux = cpl_vector_get_mean(flux_vec) ;
    cpl_vector_delete(tmp_vec) ;
    cpl_vector_delete(flux_vec) ;
    cpl_vector_delete(lmins) ;

    /* SEPARATION */
    first_fringe_idx = 1 ;          // Skip first
    last_fringe_idx = npeaks-2 ;    // Skip Last
    nb_fringes_between = last_fringe_idx - first_fringe_idx + 1 - 2 ;

    qc_separation = (cpl_vector_get(lpeaks, last_fringe_idx) -
            cpl_vector_get(lpeaks, first_fringe_idx)) /
        nb_fringes_between ;

    cpl_vector_delete(lpeaks) ;
    cpl_vector_delete(fpeaks) ;
    cpl_vector_delete(fmins) ;

    /* Store the QCs */
    cpl_propertylist_append_double(plist, CR2RES_HEADER_QC_FPI_CONTRAST, 
            qc_contrast) ;
    cpl_propertylist_append_double(plist, CR2RES_HEADER_QC_FPI_FLUX, 
            qc_flux) ;
    cpl_propertylist_append_double(plist, CR2RES_HEADER_QC_FPI_SEPARATION, 
            qc_separation) ;

    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief Compute QC.TILT
  @param    tws     3 tw tables (1 per detector)
  @param    cwlen
  @param    plist   Header for holding the QC.TILTn values
  @return   0 if ok
  
  Use Detector 1 if CWLEN (from input file headers) is < 1100.0,
  detector 2 otherwise.

  For each order n of the detector, evaluate the slit curvature at
  x=1024 pixels (middle of the detector), evaluate the slit position at
  the bottom (Xbot, Ybot) and the top (Xtop, Ytop) of the trace, and compute 
  QC.TILTn = atan [(Xtop-Xbot)/(Ytop-Ybot)] x 180 / pi 
 
 */
/*----------------------------------------------------------------------------*/
static int cr2res_cal_wave_qc_tilt(
        const cpl_table     **      tws,
        double                      cwlen,
        cpl_propertylist    *       plist)
{
    const cpl_table *   ref_tw ;
    double              tilt_avg;
    int                 ref_x, ntilt ;
    cpl_size            nrows, k ;

    /* Check Inputs */
    if (tws==NULL) return -1 ;
    for (k=0 ; k< CR2RES_NB_DETECTORS ; k++) 
        if (tws[k] == NULL) return -1 ;

    /* Initialize */
    ref_x = 1024 ;
    tilt_avg = 0.0 ;
    ntilt = 0 ;

    /* Detector used depends on the cwlen */
    if (cwlen <= 1100.0)    ref_tw = tws[0] ;
    else                    ref_tw = tws[1] ;

    nrows = cpl_table_get_nrow(ref_tw) ;

    /* Loop on rows */
    for (k=0 ; k<nrows ; k++) {

        const cpl_array *   tmp_array ;
        cpl_polynomial  *   slit_poly_a ;
        cpl_polynomial  *   slit_poly_b ;
        cpl_polynomial  *   slit_poly_c ;
        cpl_polynomial  *   upper_poly ;
        cpl_polynomial  *   lower_poly ;
        cpl_polynomial  *   slit_curv_poly ;
        char            *   qc_name;
        double              top_x, top_y, bot_x, bot_y;
        double              tilt;
        int                 cur_order ;

        cur_order = cpl_table_get(ref_tw, CR2RES_COL_ORDER, k, NULL) ;

        /* Get the Slit curvature for the trace */
        tmp_array = cpl_table_get_array(ref_tw, CR2RES_COL_SLIT_CURV_A, k) ;
        slit_poly_a = cr2res_convert_array_to_poly(tmp_array) ;
        tmp_array = cpl_table_get_array(ref_tw, CR2RES_COL_SLIT_CURV_B, k) ;
        slit_poly_b = cr2res_convert_array_to_poly(tmp_array) ;
        tmp_array = cpl_table_get_array(ref_tw, CR2RES_COL_SLIT_CURV_C, k) ;
        slit_poly_c = cr2res_convert_array_to_poly(tmp_array) ;

		/* Get the Upper Polynomial */
		tmp_array = cpl_table_get_array(ref_tw, CR2RES_COL_UPPER, k) ;
		upper_poly = cr2res_convert_array_to_poly(tmp_array) ;

		/* Get the Lower Polynomial */
		tmp_array = cpl_table_get_array(ref_tw, CR2RES_COL_LOWER, k) ;
		lower_poly = cr2res_convert_array_to_poly(tmp_array) ;

        /* Create the slit curvature polynomial at position ref_x */
        slit_curv_poly = cr2res_slit_curv_build_poly(slit_poly_a,
                slit_poly_b, slit_poly_c, ref_x) ;
        cpl_polynomial_delete(slit_poly_a) ;
        cpl_polynomial_delete(slit_poly_b) ;
        cpl_polynomial_delete(slit_poly_c) ;
                
        top_y = cpl_polynomial_eval_1d(upper_poly, ref_x, NULL) ;
        bot_y = cpl_polynomial_eval_1d(lower_poly, ref_x, NULL) ;
        cpl_polynomial_delete(upper_poly) ;
        cpl_polynomial_delete(lower_poly) ;

        top_x = cpl_polynomial_eval_1d(slit_curv_poly, top_y, NULL) ;
        bot_x = cpl_polynomial_eval_1d(slit_curv_poly, bot_y, NULL) ;
        cpl_polynomial_delete(slit_curv_poly) ;

		/* Compute the tilt */
        tilt = atan((top_x-bot_x)/(top_y-bot_y))*180.0/CRIRES_PI;
        cpl_msg_debug(__func__, 
                "TILT : %g / Upper (x,y)=(%g,%g) Lower (x,y)=(%g,%g)",
                tilt, top_x, top_y, bot_x, bot_y) ;

        /* Store the tilt */
        qc_name = cpl_sprintf(CR2RES_HEADER_QC_TILT, cur_order) ;
		cpl_propertylist_append_double(plist, qc_name, tilt) ;
		cpl_free(qc_name) ;

		/* update the tilt average */
        tilt_avg += tilt ;
        ntilt++ ;
    }

    /* Store global tilt */
    if (ntilt > 0) {
        cpl_msg_debug(__func__, "TILT GLOBAL: %g", tilt_avg/ntilt) ;
		cpl_propertylist_append_double(plist, CR2RES_HEADER_QC_TILT_GLOBAL, 
                tilt_avg/ntilt) ;
    }

    return 0 ;
}
