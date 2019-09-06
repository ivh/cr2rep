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
#include "cr2res_pfits.h"
#include "cr2res_dfs.h"
#include "cr2res_io.h"
#include "cr2res_extract.h"
#include "cr2res_trace.h"
#include "cr2res_wave.h"

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

#define RECIPE_STRING "cr2res_util_wave"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_util_wave_create(cpl_plugin *);
static int cr2res_util_wave_exec(cpl_plugin *);
static int cr2res_util_wave_destroy(cpl_plugin *);
static int cr2res_util_wave(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_util_wave_description[] = "\
Wavelength Calibration                                                  \n\
  This utility performs the wavelength calibration on already extracted \n\
  spectra. It can support different methods (--wl_method parameter):    \n\
    XCORR:  Cross Correlation with a emission lines catalog (default)   \n\
    LINE1D: Line identification and fitting for each 1D spectra         \n\
    LINE2D: Line identification and fitting for all 1D spectra at once  \n\
    ETALON: Does not require any static calibration file                \n\
    AUTO:   Guess the Method from the input file header                 \n\
                                                                        \n\
  Inputs                                                                \n\
    raw.fits " CR2RES_EXTRACT_1D_PROTYPE " [1 to n]                     \n\
    trace.fits " CR2RES_CAL_FLAT_TW_PROCATG " [1]                       \n\
            or " CR2RES_CAL_FLAT_TW_MERGED_PROCATG "                    \n\
            or " CR2RES_UTIL_TRACE_TW_PROCATG "                         \n\
            or " CR2RES_UTIL_WAVE_TW_PROCATG "                          \n\
            or " CR2RES_CAL_WAVE_TW_PROCATG "                           \n\
            or " CR2RES_UTIL_SLIT_CURV_TW_PROCATG "                     \n\
    lines.fits " CR2RES_EMISSION_LINES_PROCATG " [0 to 1]               \n\
                                                                        \n\
  Outputs                                                               \n\
    <input_name>_tw.fits " 
    CR2RES_UTIL_WAVE_TW_PROCATG"\n\
    <input_name>_wave_map.fits " 
    CR2RES_UTIL_WAVE_MAP_PROCATG "\n\
    <input_name>_lines_diagnostics.fits " 
    CR2RES_UTIL_WAVE_LINES_DIAGNOSTICS_PROCATG "\n\
    <input_name>_extracted.fits " 
    CR2RES_UTIL_WAVE_EXTRACT_1D_PROCATG "\n\
                                                                        \n\
  Algorithm                                                             \n\
    loop on raw frames f:                                               \n\
      loop on detectors d:                                              \n\
        Load the trace wave tw(f,d)                                     \n\
        Load the extracted spectra table ext(f,d)                       \n\
        Call cr2res_wave_apply(tw(f,d),ext(f,d),emission_lines)         \n\
            -> lines diagnostics(f,d)                                   \n\
            -> updated_extracted(f,d)                                   \n\
            -> trace_wave_out(f,d)                                      \n\
        Create the wavelength map wave_map(f,d)                         \n\
      Save lines diagnostics(f)                                         \n\
      Save updated_extracted(f)                                         \n\
      Save wave_map(f)                                                  \n\
      Save trace_wave_out(f)                                            \n\
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
    cr2res_io_load_TRACE_WAVE()                                         \n\
    cr2res_io_load_EXTRACT_1D()                                         \n\
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
                    cr2res_util_wave_description,
                    "Thomas Marquart, Yves Jung",
                    PACKAGE_BUGREPORT,
                    cr2res_get_license(),
                    cr2res_util_wave_create,
                    cr2res_util_wave_exec,
                    cr2res_util_wave_destroy)) {
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
static int cr2res_util_wave_create(cpl_plugin * plugin)
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
    p = cpl_parameter_new_value("cr2res.cr2res_util_wave.detector",
            CPL_TYPE_INT, "Only reduce the specified detector",
            "cr2res.cr2res_util_wave", 0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "detector");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_wave.order",
            CPL_TYPE_INT, "Only reduce the specified order",
            "cr2res.cr2res_util_wave", -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "order");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_wave.trace_nb",
            CPL_TYPE_INT, "Only reduce the specified trace number",
            "cr2res.cr2res_util_wave", -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "trace_nb");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_wave.wl_method",
            CPL_TYPE_STRING, 
            "Wavelength Method (AUTO / XCORR / LINE1D / LINE2D / ETALON)",
            "cr2res.cr2res_util_wave", "AUTO");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "wl_method");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_wave.wl_shift",
            CPL_TYPE_DOUBLE, "Wavelength shift (nm) to apply to the guess",
            "cr2res.cr2res_util_wave", 0.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "wl_shift");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_wave.wl_est",
            CPL_TYPE_STRING,
            "Estimated wavelength [start, end] (in nm)",
            "cr2res.cr2res_util_wave", "-1.0, -1.0");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "wl_est");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_wave.wl_err",
            CPL_TYPE_STRING,
            "Estimated wavelength error [start_err, end_err] (in nm)",
            "cr2res.cr2res_util_wave", "-1.0, -1.0");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "wl_err");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_wave.wl_degree",
            CPL_TYPE_INT, "Wavelength Polynomial degree",
            "cr2res.cr2res_util_wave", 3);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "wl_degree");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_wave.log",
            CPL_TYPE_BOOL, "Flag for taking the Log() value of the lines",
            "cr2res.cr2res_util_wave", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "log");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_wave.propagate",
            CPL_TYPE_BOOL, "Flag for using the input WL when no computation",
            "cr2res.cr2res_util_wave", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "propagate");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_wave.display",
            CPL_TYPE_BOOL, "Flag for display",
            "cr2res.cr2res_util_wave", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "display");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_wave.display_range",
            CPL_TYPE_STRING,
            "Wavelength range to display [start, end] (in nm)",
            "cr2res.cr2res_util_wave", "-1.0, -1.0");
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
static int cr2res_util_wave_exec(cpl_plugin * plugin)
{
    cpl_recipe  *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1;

    return cr2res_util_wave(recipe->frames, recipe->parameters);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_util_wave_destroy(cpl_plugin * plugin)
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
static int cr2res_util_wave(
        cpl_frameset            *   frameset,
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter *   param;
    int                     reduce_det, reduce_order, reduce_trace,
                            wl_degree, display, log_flag, propagate_flag ;
    double                  wl_start, wl_end, wl_err_start, wl_err_end, 
                            wl_shift, display_wmin, display_wmax ;
    cr2res_wavecal_type     wavecal_type ;
    const char          *   sval ;
    cpl_frameset        *   rawframes ;
    const cpl_frame     *   cur_frame ;
    const char          *   cur_fname ;
    cpl_frameset        *   cur_fset ;
    const cpl_frame     *   trace_wave_frame ;
    const cpl_frame     *   lines_frame ;
    cpl_table           *   trace_wave ;
    cpl_table           *   extracted_table ;
    char                *   out_file;
    cpl_table           *   out_trace_wave[CR2RES_NB_DETECTORS] ;
    cpl_table           *   lines_diagnostics[CR2RES_NB_DETECTORS] ;
    cpl_table           *   updated_extracted_table[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   out_wave_map[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   ext_plist[CR2RES_NB_DETECTORS] ;
    int                     det_nr, order, i, j ;

    /* Needed for sscanf() */
    setlocale(LC_NUMERIC, "C");

    /* Initialise */
    wavecal_type = CR2RES_UNSPECIFIED ;
    wl_start = wl_end = wl_err_start = wl_err_end = -1.0 ;
    wl_shift = 0.0 ;
    display_wmin = display_wmax = -1.0 ;

    /* RETRIEVE INPUT PARAMETERS */
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_wave.detector");
    reduce_det = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_wave.order");
    reduce_order = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_wave.trace_nb");
    reduce_trace = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_wave.wl_method");
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
            "cr2res.cr2res_util_wave.wl_shift");
    wl_shift = cpl_parameter_get_double(param) ;
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_wave.wl_est");
    sval = cpl_parameter_get_string(param) ;
    if (sscanf(sval, "%lg,%lg", &wl_start, &wl_end) != 2) {
        return -1 ;
    }
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_wave.wl_err");
    sval = cpl_parameter_get_string(param) ;
    if (sscanf(sval, "%lg,%lg", &wl_err_start, &wl_err_end) != 2) {
        return -1 ;
    }
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_wave.wl_degree");
    wl_degree = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_wave.log");
    log_flag = cpl_parameter_get_bool(param) ;
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_wave.propagate");
    propagate_flag = cpl_parameter_get_bool(param) ;
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_wave.display");
    display = cpl_parameter_get_bool(param) ;
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_wave.display_range");
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
    lines_frame = cpl_frameset_find_const(frameset,
            CR2RES_EMISSION_LINES_PROCATG) ;
    if ((wavecal_type == CR2RES_XCORR || wavecal_type == CR2RES_LINE1D ||
                wavecal_type == CR2RES_LINE2D) && lines_frame == NULL) {
        cpl_msg_error(__func__,
                "The catalog file is needed for XCORR/LINE1D/LINE2D");
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Get the RAW Frames */
    rawframes = cr2res_extract_frameset(frameset, CR2RES_EXTRACT_1D_PROTYPE) ;
    if (rawframes == NULL) {
        cpl_msg_error(__func__, "Could not find RAW frames") ;
        return -1 ;
    }

    /* Loop on the RAW frames */
    for (i=0 ; i<cpl_frameset_get_size(rawframes) ; i++) {
        /* Get the Current Frame */
        cur_frame = cpl_frameset_get_position(rawframes, i) ;
        cur_fname = cpl_frame_get_filename(cur_frame) ;
        cpl_msg_info(__func__, "Reduce Frame %s", cur_fname) ;
        cpl_msg_indent_more() ;

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
			cpl_msg_error(__func__, 
                    "Limiting to one order with LINE2D impossible");
			cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
			return -1 ;
		}

        /* Loop over the detectors */
        for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {

            /* Initialise */
            out_trace_wave[det_nr-1] = NULL ;
            lines_diagnostics[det_nr-1] = NULL ;
            updated_extracted_table[det_nr-1] = NULL ;
            out_wave_map[det_nr-1] = NULL ;
            ext_plist[det_nr-1] = NULL ;

            /* Store the extenѕion header for product saving */
            ext_plist[det_nr-1] = cpl_propertylist_load(cur_fname,
                    cr2res_io_get_ext_idx(cur_fname, det_nr, 1)) ;

            /* Compute only one detector */
            if (reduce_det != 0 && det_nr != reduce_det) continue ;

            cpl_msg_info(__func__, "Process detector number %d", det_nr) ;
            cpl_msg_indent_more() ;

            /* Load the TRACE_WAVE table of this detector */
            cpl_msg_info(__func__, "Load the TRACE_WAVE table") ;
            if ((trace_wave = cr2res_io_load_TRACE_WAVE(
                            cpl_frame_get_filename(trace_wave_frame), 
                            det_nr)) == NULL) {
                cpl_msg_error(__func__,"Failed to load table - skip detector");
                cpl_error_reset() ;
                cpl_msg_indent_less() ;
                continue ;
            }

            /* Load the EXTRACT1D table of this detector */
            cpl_msg_info(__func__, "Load the EXTRACT1D table") ;
            if ((extracted_table = cr2res_io_load_EXTRACT_1D(cur_fname,
                            det_nr)) == NULL) {
                cpl_msg_error(__func__,"Failed to load table - skip detector");
                cpl_table_delete(trace_wave) ;
                cpl_error_reset() ;
                cpl_msg_indent_less() ;
                continue ;
            }

            /* Compute the Wavelength Calibration */
            cpl_msg_info(__func__, "Compute the Wavelength") ;
            if (cr2res_wave_apply(trace_wave, extracted_table,
                        lines_frame, reduce_order, reduce_trace, wavecal_type,
                        wl_degree, wl_start, wl_end, wl_err_start, wl_err_end, 
                        wl_shift, log_flag, propagate_flag, display,
                        display_wmin, display_wmax,
                        NULL,
                        &(lines_diagnostics[det_nr-1]),
                        &(updated_extracted_table[det_nr-1]),
                        &(out_trace_wave[det_nr-1]))) {
                cpl_msg_error(__func__, "Failed to calibrate - skip detector");
                cpl_table_delete(trace_wave) ;
                cpl_table_delete(extracted_table) ;
                cpl_error_reset() ;
                cpl_msg_indent_less() ;
                continue ;
            }
            cpl_table_delete(trace_wave) ;
            cpl_table_delete(extracted_table) ;

            /* Generate the Wave Map */
            out_wave_map[det_nr-1] =
                cr2res_wave_gen_wave_map(out_trace_wave[det_nr-1]) ;
            cpl_msg_indent_less() ;
        }

        /* Generate the currently used frameset */
        /* TODO : add calibrations */
        cur_fset = cpl_frameset_new() ;
        cpl_frameset_insert(cur_fset, cpl_frame_duplicate(cur_frame)) ;

        /* Save the new trace_wave table */
        out_file = cpl_sprintf("%s_tw.fits",
                cr2res_get_base_name(cr2res_get_root_name(cur_fname)));
        cr2res_io_save_TRACE_WAVE(out_file, frameset, cur_fset, parlist,
                out_trace_wave, NULL, ext_plist, CR2RES_UTIL_WAVE_TW_PROCATG,
                RECIPE_STRING) ;
        cpl_free(out_file);

        /* Save the Wave Map */
        out_file = cpl_sprintf("%s_wave_map.fits",
                cr2res_get_base_name(cr2res_get_root_name(cur_fname)));
        cr2res_io_save_WAVE_MAP(out_file, frameset, cur_fset, parlist, 
                out_wave_map, NULL, ext_plist, 
                CR2RES_UTIL_WAVE_MAP_PROCATG, RECIPE_STRING) ;
        cpl_free(out_file);

        /* Save the Wave Map */
        out_file = cpl_sprintf("%s_extracted.fits",
                cr2res_get_base_name(cr2res_get_root_name(cur_fname)));
        cr2res_io_save_EXTRACT_1D(out_file, frameset, cur_fset, parlist, 
                updated_extracted_table, NULL, ext_plist, 
                CR2RES_UTIL_WAVE_EXTRACT_1D_PROCATG, RECIPE_STRING) ;
        cpl_free(out_file);

        if (wavecal_type == CR2RES_LINE2D || wavecal_type == CR2RES_LINE1D) {
            /* Save the Lines Diagnostics */
            out_file = cpl_sprintf("%s_lines_diagnostics.fits",
                    cr2res_get_base_name(cr2res_get_root_name(cur_fname)));
            cr2res_io_save_LINES_DIAGNOSTICS(out_file, frameset, cur_fset, 
                    parlist, lines_diagnostics, NULL, ext_plist,
                    CR2RES_UTIL_WAVE_LINES_DIAGNOSTICS_PROCATG, RECIPE_STRING) ;
            cpl_free(out_file);
        }
        cpl_frameset_delete(cur_fset) ;

        /* Free and return */
        for (j=0 ; j<CR2RES_NB_DETECTORS ; j++) {
            if (ext_plist[j] != NULL)
                cpl_propertylist_delete(ext_plist[j]) ;
            if (out_trace_wave[j] != NULL)
                cpl_table_delete(out_trace_wave[j]) ;
            if (lines_diagnostics[j] != NULL)
                cpl_table_delete(lines_diagnostics[j]) ;
            if (updated_extracted_table[j] != NULL)
                cpl_table_delete(updated_extracted_table[j]) ;
            if (out_wave_map[j] != NULL) 
                hdrl_image_delete(out_wave_map[j]) ;
        }
        cpl_msg_indent_less() ;
    }
    cpl_frameset_delete(rawframes) ;
    return (int)cpl_error_get_code();
}
