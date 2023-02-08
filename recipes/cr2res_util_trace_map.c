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
#include "hdrl.h"

#include "cr2res_utils.h"
#include "cr2res_pfits.h"
#include "cr2res_dfs.h"
#include "cr2res_io.h"
#include "cr2res_trace.h"
#include "cr2res_wave.h"
#include "cr2res_slit_curv.h"

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

#define RECIPE_STRING "cr2res_util_trace_map"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_util_trace_map_create(cpl_plugin *);
static int cr2res_util_trace_map_exec(cpl_plugin *);
static int cr2res_util_trace_map_destroy(cpl_plugin *);
static int cr2res_util_trace_map(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_util_trace_map_description[] = "\
Maps creation                                                           \n\
  Each input TRACE_WAVE file is converted into maps to visualize the    \n\
  traces, wavelengths and the slit curvature                            \n\
                                                                        \n\
  Inputs                                                                \n\
    raw.fits " CR2RES_TW_DRSTYPE " [1 to n]                             \n\
                                                                        \n\
  Outputs                                                               \n\
    <input_name>_slit_curve.fits " 
    CR2RES_UTIL_TRACE_MAP_SLIT_CURVE_PROCATG "\n\
    <input_name>_wave.fits " 
    CR2RES_UTIL_TRACE_MAP_WL_PROCATG "\n\
    <input_name>_trace.fits " 
    CR2RES_UTIL_TRACE_MAP_TRACE_PROCATG "\n\
                                                                        \n\
  Algorithm                                                             \n\
    loop on input raw frames f:                                         \n\
      loop on detectors d:                                              \n\
        Load the trace_wave extension                                   \n\
        Call cr2res_wave_gen_wave_map() to generate wave_map(d)         \n\
        Call cr2res_trace_gen_image() to generate traces_map(d)         \n\
        Call cr2res_slit_curv_gen_map() to generate slit_curv_map(d)    \n\
    Save wave_map, traces_map and slit_curv_map                         \n\
                                                                        \n\
  Library Functions used                                                \n\
    cr2res_io_load_TRACE_WAVE()                                         \n\
    cr2res_trace_gen_image()                                            \n\
    cr2res_slit_curv_gen_map()                                          \n\
    cr2res_io_save_SLIT_CURV_MAP()                                      \n\
    cr2res_io_save_WAVE_MAP()                                           \n\
    cr2res_io_save_TRACE_MAP()                                          \n\
";

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
                    "TRACE_WAVE maps creation",
                    cr2res_util_trace_map_description,
                    CR2RES_PIPELINE_AUTHORS,
                    PACKAGE_BUGREPORT,
                    cr2res_get_license(),
                    cr2res_util_trace_map_create,
                    cr2res_util_trace_map_exec,
                    cr2res_util_trace_map_destroy)) {
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
static int cr2res_util_trace_map_create(cpl_plugin * plugin)
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
    p = cpl_parameter_new_value("cr2res.cr2res_util_trace_map.detector",
            CPL_TYPE_INT, "Only reduce the specified detector",
            "cr2res.cr2res_util_trace_map", 0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "detector");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_trace_map.order",
            CPL_TYPE_INT, "Only reduce the specified order",
            "cr2res.cr2res_util_trace_map", -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "order");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_trace_map.trace_nb",
            CPL_TYPE_INT, "Only reduce the specified trace number",
            "cr2res.cr2res_util_trace_map", -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "trace_nb");
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
static int cr2res_util_trace_map_exec(cpl_plugin * plugin)
{
    cpl_recipe  *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1;

    return cr2res_util_trace_map(recipe->frames, recipe->parameters);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_util_trace_map_destroy(cpl_plugin * plugin)
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
static int cr2res_util_trace_map(
        cpl_frameset            *   frameset,
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter *   param;
    int                     reduce_det, reduce_order, reduce_trace ;
    cpl_frameset        *   rawframes ;
    const cpl_frame     *   cur_frame ;
    const char          *   cur_fname ;
    cpl_frameset        *   cur_fset ;
    char                *   out_file;
    cpl_image           *   img_tmp ;
    hdrl_image          *   wl_maps[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   trace_maps[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   slit_curve_maps[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   ext_plist[CR2RES_NB_DETECTORS] ;
    cpl_table           *   trace_table ;
    int                     det_nr, ext_nr, order, i;

    /* Initialise */

    /* RETRIEVE INPUT PARAMETERS */
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_trace_map.detector");
    reduce_det = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_trace_map.order");
    reduce_order = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_trace_map.trace_nb");
    reduce_trace = cpl_parameter_get_int(param);

    /* Check Parameters */

    /* Identify the RAW and CALIB frames in the input frameset */
    if (cr2res_dfs_set_groups(frameset) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot identify RAW and CALIB frames") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Get Calibration frames */

    /* Get the Rawframes */
    rawframes = cr2res_extract_frameset(frameset, CR2RES_TW_DRSTYPE) ;
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
            wl_maps[det_nr-1] = NULL ;
            trace_maps[det_nr-1] = NULL ;
            slit_curve_maps[det_nr-1] = NULL ;

            /* Get Extension Numbers */
            ext_nr = cr2res_io_get_ext_idx(cur_fname, det_nr, 1) ;
            if (ext_nr < 0) continue ;

            /* Store the extension header for product saving */
            ext_plist[det_nr-1] = cpl_propertylist_load(cur_fname, ext_nr) ;

            /* Compute only one detector */
            if (reduce_det != 0 && det_nr != reduce_det) continue ;

            cpl_msg_info(__func__, "Process detector number %d", det_nr) ;
            cpl_msg_indent_more() ;

            /* Load the trace table of this detector */
            cpl_msg_info(__func__, "Load the trace table") ;
            if ((trace_table = cr2res_io_load_TRACE_WAVE(cur_fname,
                            det_nr)) == NULL) {
                cpl_msg_error(__func__,
                        "Failed to get trace table - skip detector");
                cpl_error_reset() ;
                cpl_msg_indent_less() ;
                continue ;
            }

            /* Create WAVE_MAP */
            wl_maps[det_nr-1] =
                cr2res_wave_gen_wave_map(trace_table) ;

            /* Create TRACE MAP */
            img_tmp = cr2res_trace_gen_image(trace_table,
                    CR2RES_DETECTOR_SIZE, CR2RES_DETECTOR_SIZE) ;
            trace_maps[det_nr-1] = hdrl_image_create(img_tmp, NULL);
            cpl_image_delete(img_tmp);
           
            /* Create SLIT_CURVE MAP */
            slit_curve_maps[det_nr-1] =
                cr2res_slit_curv_gen_map(trace_table, reduce_order,
                        reduce_trace, 50, 0) ;

            cpl_table_delete(trace_table) ;
            cpl_msg_indent_less() ;
        }

        /* Save the Products */
        cur_fset = cpl_frameset_new() ;
        cpl_frameset_insert(cur_fset, cpl_frame_duplicate(cur_frame)) ;

        out_file = cpl_sprintf("%s_slit_curve.fits",
                        cr2res_get_base_name(cr2res_get_root_name(cur_fname)));
        cr2res_io_save_SLIT_CURV_MAP(out_file, frameset, cur_fset, parlist, 
                slit_curve_maps, NULL, ext_plist, 
                CR2RES_UTIL_TRACE_MAP_SLIT_CURVE_PROCATG, RECIPE_STRING);
        cpl_free(out_file);

        out_file = cpl_sprintf("%s_wave.fits",
                        cr2res_get_base_name(cr2res_get_root_name(cur_fname)));
        cr2res_io_save_WAVE_MAP(out_file, frameset, cur_fset, parlist, wl_maps, 
                NULL, ext_plist, CR2RES_UTIL_TRACE_MAP_WL_PROCATG, 
                RECIPE_STRING);
        cpl_free(out_file);

        out_file = cpl_sprintf("%s_trace.fits",
                        cr2res_get_base_name(cr2res_get_root_name(cur_fname)));
        cr2res_io_save_TRACE_MAP(out_file, frameset, cur_fset, parlist, 
                trace_maps, NULL, ext_plist, 
                CR2RES_UTIL_TRACE_MAP_TRACE_PROCATG, RECIPE_STRING);
        cpl_free(out_file);
        cpl_frameset_delete(cur_fset) ;

        /* Free and return */
        for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
            if (ext_plist[det_nr-1] != NULL)
                cpl_propertylist_delete(ext_plist[det_nr-1]) ;
            if (slit_curve_maps[det_nr-1] != NULL) 
                hdrl_image_delete(slit_curve_maps[det_nr-1]) ;
            if (trace_maps[det_nr-1] != NULL)
                hdrl_image_delete(trace_maps[det_nr-1]) ;
            if (wl_maps[det_nr-1] != NULL)
                hdrl_image_delete(wl_maps[det_nr-1]) ;
        }
        cpl_msg_indent_less() ;
    }
    cpl_frameset_delete(rawframes) ;
    return (int)cpl_error_get_code();
}

