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

#include <cpl.h>
#include "hdrl.h"

#include "cr2res_utils.h"
#include "cr2res_pfits.h"
#include "cr2res_dfs.h"
#include "cr2res_io.h"
#include "cr2res_trace.h"
#include "cr2res_wave.h"

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

#define RECIPE_STRING "cr2res_util_trace"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static cpl_frameset * cr2res_util_trace_find_RAW(const cpl_frameset * in) ;
static int cr2res_util_trace_create(cpl_plugin *);
static int cr2res_util_trace_exec(cpl_plugin *);
static int cr2res_util_trace_destroy(cpl_plugin *);
static int cr2res_util_trace(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_util_trace_description[] = "\
Traces detection                                                        \n\
  This utility detects the traces, fits polynomials on their edges      \n\
  (Upper and Lower) and in their centers (All), and stores these        \n\
  informations in the TRACE_WAVE file.                                  \n\
  Each trace is uniquely identified by its Order/TraceNb values.        \n\
  The Order values refer to the keywords indices (e.g. HIERARCH ESO INS \n\
  WLEN CENY4) in the product headers.                                   \n\
  The TraceNb starts with 1, identifies traces within the same order.   \n\
  The additional columns :                                              \n\
    "CR2RES_COL_WAVELENGTH"                                             \n\
    "CR2RES_COL_WAVELENGTH_ERROR"                                       \n\
    "CR2RES_COL_SLIT_CURV_A"                                            \n\
    "CR2RES_COL_SLIT_CURV_B"                                            \n\
    "CR2RES_COL_SLIT_CURV_C"                                            \n\
    "CR2RES_COL_SLIT_FRACTION"                                          \n\
  are filled with default values.                                       \n\
                                                                        \n\
  Inputs                                                                \n\
    raw.fits " CR2RES_FLAT_RAW " [1 to n]                               \n\
          or " CR2RES_CALIBRATED_PROTYPE "                              \n\
                                                                        \n\
  Outputs                                                               \n\
    <input_name>_tw.fits " CR2RES_UTIL_TRACE_TW_PROCATG"                \n\
                                                                        \n\
  Algorithm                                                             \n\
    loop on input raw frames f:                                         \n\
      loop on detectors d:                                              \n\
        Use cr2res_trace(--degree, --min_cluster, --smooth, --opening)  \n\
                to measure the traces                                   \n\
        if --split_traces, call cr2res_trace_split_traces() to split    \n\
               the traces                                               \n\
        Use cr2res_trace_add_extra_columns() to add the additional      \n\
                columns (slit fraction, wl, slit curvature)             \n\
      Save the trace wave table                                         \n\
                                                                        \n\
  Library functions uѕed                                                \n\
    cr2res_io_load_image()                                              \n\
    cr2res_trace()                                                      \n\
    cr2res_trace_add_extra_columns()                                    \n\
    cr2res_trace_split_traces()                                         \n\
    cr2res_io_save_TRACE_WAVE()                                         \n\
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
                    "Trace utility",
                    cr2res_util_trace_description,
                    "Thomas Marquart, Yves Jung",
                    PACKAGE_BUGREPORT,
                    cr2res_get_license(),
                    cr2res_util_trace_create,
                    cr2res_util_trace_exec,
                    cr2res_util_trace_destroy)) {
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
static int cr2res_util_trace_create(cpl_plugin * plugin)
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
    p = cpl_parameter_new_value("cr2res.cr2res_util_trace.degree",
            CPL_TYPE_INT, "polynomial degree for the fit to the orders",
            "cr2res.cr2res_util_trace", 5);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "degree");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_trace.min_cluster",
            CPL_TYPE_INT, "size in pixels of the smallest allowed cluster",
            "cr2res.cr2res_util_trace", 40000);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "min_cluster");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_trace.smooth_x",
            CPL_TYPE_INT, "Length of the smoothing kernel in x",
            "cr2res.cr2res_util_trace", 111);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "smooth_x");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_trace.smooth_y",
            CPL_TYPE_INT, "Length of the smoothing kernel in y",
            "cr2res.cr2res_util_trace", 401);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "smooth_y");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_trace.threshold",
            CPL_TYPE_DOUBLE, "Detection Threshold",
            "cr2res.cr2res_util_trace", 300.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "threshold");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_trace.opening",
            CPL_TYPE_BOOL, "Use a morphological opening to rejoin clusters",
            "cr2res.cr2res_util_trace", TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "opening");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_trace.split_traces",
            CPL_TYPE_INT, "Split the full slit traces",
            "cr2res.cr2res_util_trace", 0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "split_traces");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_trace.detector",
            CPL_TYPE_INT, "Only reduce the specified detector",
            "cr2res.cr2res_util_trace", 0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "detector");
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
static int cr2res_util_trace_exec(cpl_plugin * plugin)
{
    cpl_recipe  *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1;

    return cr2res_util_trace(recipe->frames, recipe->parameters);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_util_trace_destroy(cpl_plugin * plugin)
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
static int cr2res_util_trace(
        cpl_frameset            *   frameset,
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter *   param;
    int                     min_cluster, degree, opening, reduce_det,
                            split_traces, smooth_x, smooth_y ;
    double                  threshold ;
    cpl_frameset        *   rawframes ;
    const cpl_frame     *   cur_frame ;
    const char          *   cur_fname ;
    cpl_frameset        *   cur_fset ;
    char                *   out_file;
    hdrl_image          *   flat_ima ;
    cpl_image           *   debug_ima ;
    int                     det_nr ;
    cpl_table           *   traces_tmp ;
    cpl_table           *   traces[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   ext_plist[CR2RES_NB_DETECTORS] ;
    int                     i ;

    /* RETRIEVE INPUT PARAMETERS */
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_trace.min_cluster");
    min_cluster = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_trace.degree");
    degree = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_trace.smooth_x");
    smooth_x = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_trace.smooth_y");
    smooth_y = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_trace.threshold");
    threshold = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_trace.opening");
    opening = cpl_parameter_get_bool(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_trace.split_traces");
    split_traces = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_trace.detector");
    reduce_det = cpl_parameter_get_int(param);

    /* Check Parameters */

    /* Identify the RAW and CALIB frames in the input frameset */
    if (cr2res_dfs_set_groups(frameset) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot identify RAW and CALIB frames") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Get Calibration frames */

    /* Get the rawframes */
    rawframes = cr2res_util_trace_find_RAW(frameset) ;
    if (rawframes==NULL || cpl_frameset_get_size(rawframes) <= 0) {
        cpl_msg_error(__func__, "Cannot find any RAW file") ;
        cpl_error_set(__func__, CPL_ERROR_DATA_NOT_FOUND) ;
        return -1 ;
    }

    /* Loop on the RAW frames */
    for (i=0 ; i<cpl_frameset_get_size(rawframes) ; i++) {
        /* Get the Current Frame */
        cur_frame = cpl_frameset_get_position(rawframes, i) ;
        cur_fname = cpl_frame_get_filename(cur_frame) ;
        cpl_msg_info(__func__, "Reduce Frame %s", cur_fname) ;
        cpl_msg_indent_more() ;

        /* Loop over the detectors */
        for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {

            /* Initialise */
            traces[det_nr-1] = NULL ;
            ext_plist[det_nr-1] = NULL ;

            /* Store the extenѕion header for product saving */
            ext_plist[det_nr-1] = cpl_propertylist_load(cur_fname,
                    cr2res_io_get_ext_idx(cur_fname, det_nr, 1)) ;

            /* Compute only one detector */
            if (reduce_det != 0 && det_nr != reduce_det) continue ;

            cpl_msg_info(__func__, "Process detector number %d", det_nr) ;
            cpl_msg_indent_more() ;

            /* Load the image in which the orders are to extract*/
            cpl_msg_info(__func__, "Load the Image") ;
            if ((flat_ima = cr2res_io_load_image(cur_fname, det_nr)) == NULL) {
                cpl_msg_warning(__func__, 
                        "Cannot load the image - skip detector");
                cpl_error_reset() ;
                cpl_msg_indent_less() ;
                continue ;
            }

            /* Get the traces */
            cpl_msg_info(__func__, "Compute the traces") ;
            cpl_msg_indent_more() ;
            if ((traces[det_nr-1] = cr2res_trace(hdrl_image_get_image(flat_ima),
                            smooth_x, smooth_y, threshold, opening, degree, 
                            min_cluster)) == NULL) {
                cpl_msg_warning(__func__, 
                        "Cannot compute trace - skip detector");
                cpl_error_reset() ;
                hdrl_image_delete(flat_ima) ;
                cpl_msg_indent_less() ;
                cpl_msg_indent_less() ;
                continue ;
            }
            cpl_msg_indent_less() ;

            /* Add The remaining Columns to the trace table */
            if (cr2res_trace_add_extra_columns(traces[det_nr-1],
                    cur_fname, det_nr) != 0) {
                cpl_msg_warning(__func__, 
                        "Cannot complete the trace table - skip detector");
                cpl_error_reset() ;
                hdrl_image_delete(flat_ima) ;
                cpl_table_delete(traces[det_nr-1]) ;
                traces[det_nr-1] = NULL ;
                cpl_msg_indent_less() ;
                continue ;
            }

            /* Split the traces when required */
            if (split_traces) {
                cpl_msg_info(__func__, 
                        "Split the full slit traces in %d traces", 
                        split_traces);
                /* Split the full slit traces */
                if ((traces_tmp = cr2res_trace_split(traces[det_nr-1], -1, 
                                split_traces)) == NULL) {
                    cpl_msg_warning(__func__, 
                            "Failed splitting the traces - skip detector") ;
                    cpl_error_reset() ;
                    hdrl_image_delete(flat_ima) ;
                    cpl_table_delete(traces[det_nr-1]) ;
                    traces[det_nr-1] = NULL ;
                    cpl_msg_indent_less() ;
                    continue ;
                } else {
                    cpl_table_delete(traces[det_nr-1]) ;
                    traces[det_nr-1] = traces_tmp ;
                    traces_tmp = NULL ;
                } 
            }

            /* Debug Image */
            if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
                debug_ima = cr2res_trace_gen_image(traces[det_nr-1],
                        hdrl_image_get_size_x(flat_ima),
                        hdrl_image_get_size_y(flat_ima)) ;
                out_file = cpl_sprintf("debug_%s_trace_map_%d.fits", 
                        cr2res_get_base_name(cr2res_get_root_name(cur_fname)),
                        det_nr);
                cpl_image_save(debug_ima, out_file, CPL_BPP_IEEE_DOUBLE, NULL, 
                        CPL_IO_CREATE) ;
                cpl_free(out_file);
                cpl_image_delete(debug_ima) ;
            }
            hdrl_image_delete(flat_ima) ;
            cpl_msg_indent_less() ;
        }

        /* Save the Products */
        out_file = cpl_sprintf("%s_tw.fits", 
                cr2res_get_base_name(cr2res_get_root_name(cur_fname)));
        cur_fset = cpl_frameset_new() ;
        cpl_frameset_insert(cur_fset, cpl_frame_duplicate(cur_frame)) ;
        cr2res_io_save_TRACE_WAVE(out_file, frameset, cur_fset, parlist, traces,
                NULL, ext_plist, CR2RES_UTIL_TRACE_TW_PROCATG, RECIPE_STRING);
        cpl_frameset_delete(cur_fset) ;
        cpl_free(out_file);

        /* Free and return */
        for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
            if (ext_plist[det_nr-1] != NULL)
                cpl_propertylist_delete(ext_plist[det_nr-1]) ;
            if (traces[det_nr-1] != NULL)
                cpl_table_delete(traces[det_nr-1]) ;
        }
    }
    cpl_frameset_delete(rawframes) ;
    return (int)cpl_error_get_code();
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the RAW frames from a frameset
  @param    set     Input frame set
  @return   the RAW frameset or NULL in error case or if it is missing
    Allowed RAW types : CR2RES_FLAT_RAW
                        CR2RES_CALIBRATED_PROTYPE
 */
/*----------------------------------------------------------------------------*/
static cpl_frameset * cr2res_util_trace_find_RAW(const cpl_frameset * in)
{
    cpl_frameset    *   out ;

    /* Check entries */
    if (in == NULL) return NULL ;

    out = cr2res_extract_frameset(in, CR2RES_FLAT_RAW) ;
    if (out == NULL)
        out = cr2res_extract_frameset(in, CR2RES_CALIBRATED_PROTYPE) ;
    return out ;
}



