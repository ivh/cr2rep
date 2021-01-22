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

#include <string.h>

#include <cpl.h>
#include "hdrl.h"

#include "cr2res_slit_curv.h"
#include "cr2res_utils.h"
#include "cr2res_pfits.h"
#include "cr2res_dfs.h"
#include "cr2res_io.h"

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

#define RECIPE_STRING "cr2res_util_slit_curv"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_util_slit_curv_create(cpl_plugin *);
static int cr2res_util_slit_curv_exec(cpl_plugin *);
static int cr2res_util_slit_curv_destroy(cpl_plugin *);
static int cr2res_util_slit_curv(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_util_slit_curv_description[] = "\
Slit curvature computation                                              \n\
  For each input trace_wave file, the slit curvature is derived from    \n\
  each order that has at least 2 traces.                                \n\
                                                                        \n\
  Inputs                                                                \n\
    trace.fits " CR2RES_CAL_FLAT_TW_PROCATG " [1 to n]                  \n\
            or " CR2RES_CAL_FLAT_TW_MERGED_PROCATG "                    \n\
            or " CR2RES_UTIL_TRACE_TW_PROCATG "                         \n\
            or " CR2RES_UTIL_WAVE_TW_PROCATG "                          \n\
            or " CR2RES_CAL_WAVE_TW_PROCATG "                           \n\
            or " CR2RES_UTIL_SLIT_CURV_TW_PROCATG "                     \n\
    lamp.fits " CR2RES_WAVE_RAW " [1 to n]                              \n\
                                                                        \n\
  Outputs                                                               \n\
    <input_lamp_name>_map.fits " CR2RES_UTIL_SLIT_CURV_MAP_PROCATG"     \n\
    <input_lamp_name>_slit_curv.fits " CR2RES_UTIL_SLIT_CURV_PROCATG "  \n\
    <input_lamp_name>_tw.fits " CR2RES_UTIL_SLIT_CURV_TW_PROCATG"       \n\
                                                                        \n\
                     TODO\n\
  Algorithm                                                             \n\
    loop on input raw files pairs (t,l):                                \n\
      loop on detectors d:                                              \n\
        Load the TRACE_WAVE table for the current detector              \n\
        Loop over the traces t:                                         \n\
          Call cr2res_slit_curv_compute_order_trace() to get the        \n\
                current trace curvatures                                \n\
          Update the slit curvature in the TRACE_WAVE table             \n\
        Generate a slit curve map                                       \n\
      Save the TRACE_WAVE                                               \n\
      Save the SLIT_CURVE_MAP                                           \n\
                                                                        \n\
  Library functions uѕed                                                \n\
    cr2res_io_load_TRACE_WAVE()                                         \n\
    cr2res_slit_curv_compute_order_trace()                              \n\
    cr2res_slit_curv_gen_map()                                          \n\
    cr2res_io_save_SLIT_CURV_MAP()                                      \n\
    cr2res_io_save_TRACE_WAVE()                                         \n\
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
                    "Slit Curvature utility",
                    cr2res_util_slit_curv_description,
                    CR2RES_PIPELINE_AUTHORS,
                    PACKAGE_BUGREPORT,
                    cr2res_get_license(),
                    cr2res_util_slit_curv_create,
                    cr2res_util_slit_curv_exec,
                    cr2res_util_slit_curv_destroy)) {
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
static int cr2res_util_slit_curv_create(cpl_plugin * plugin)
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
    p = cpl_parameter_new_value("cr2res.cr2res_util_slit_curv.detector",
            CPL_TYPE_INT, "Only reduce the specified detector",
            "cr2res.cr2res_util_slit_curv", 0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "detector");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_slit_curv.order",
            CPL_TYPE_INT, "Only reduce the specified order",
            "cr2res.cr2res_util_slit_curv", -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "order");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_slit_curv.trace_nb",
            CPL_TYPE_INT, "Only reduce the specified trace number",
            "cr2res.cr2res_util_slit_curv", -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "trace_nb");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_slit_curv.display",
            CPL_TYPE_INT, "X value to display (1->2048)",
            "cr2res.cr2res_util_slit_curv", 0);
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
static int cr2res_util_slit_curv_exec(cpl_plugin * plugin)
{
    cpl_recipe  *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1;

    return cr2res_util_slit_curv(recipe->frames, recipe->parameters);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_util_slit_curv_destroy(cpl_plugin * plugin)
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
static int cr2res_util_slit_curv(
        cpl_frameset            *   frameset,
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter *   param;
    int                     reduce_det, reduce_order, reduce_trace, display ;
    cpl_frameset        *   rawframes_tw ;
    cpl_frameset        *   rawframes_lamp ;
    const cpl_frame     *   cur_frame_tw ;
    const char          *   cur_fname_tw ;
    const cpl_frame     *   cur_frame_lamp ;
    const char          *   cur_fname_lamp ;
    cpl_table           *   trace_wave[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   lamp_image[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   slit_curv_map[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   ext_plist[CR2RES_NB_DETECTORS] ;
    cpl_polynomial      *   slit_polya ;
    cpl_polynomial      *   slit_polyb ;
    cpl_polynomial      *   slit_polyc ;
    cpl_array           *   slit_array ;
    int                     det_nr, order, trace_id, nb_traces, curv_degree ;
    char                *   col_name ;
    char                *   out_file;
    cpl_array           *   curv_array ;
    int                     i, j, k ;

    /* Initialise */
    curv_degree = 2 ;

    /* RETRIEVE INPUT PARAMETERS */
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_slit_curv.detector");
    reduce_det = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_slit_curv.order");
    reduce_order = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_slit_curv.trace_nb");
    reduce_trace = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_slit_curv.display");
    display = cpl_parameter_get_int(param) ;
 
    /* Identify the RAW and CALIB frames in the input frameset */
    if (cr2res_dfs_set_groups(frameset) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot identify RAW and CALIB frames") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Get Calibration frames */

    /* Get the rawframes */
    rawframes_tw = cr2res_io_find_TRACE_WAVE_all(frameset);
    rawframes_lamp = cr2res_extract_frameset(frameset, CR2RES_WAVE_RAW) ;
    if (rawframes_tw == NULL || rawframes_lamp == NULL || 
            cpl_frameset_get_size(rawframes_lamp) <= 0 ||
            cpl_frameset_get_size(rawframes_tw) <= 0 ||
            cpl_frameset_get_size(rawframes_lamp) != 
            cpl_frameset_get_size(rawframes_tw)) {
        cpl_msg_error(__func__, "Input files not consistent") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        if (rawframes_tw != NULL) cpl_frameset_delete(rawframes_tw) ; 
        if (rawframes_lamp != NULL) cpl_frameset_delete(rawframes_lamp) ; 
        return -1 ;
    }

    /* Loop on the RAW frames */
    for (i=0 ; i<cpl_frameset_get_size(rawframes_tw) ; i++) {
        /* Get the Current Frame */
        cur_frame_tw = cpl_frameset_get_position(rawframes_tw, i) ;
        cur_fname_tw = cpl_frame_get_filename(cur_frame_tw) ;
        cur_frame_lamp = cpl_frameset_get_position(rawframes_lamp, i) ;
        cur_fname_lamp = cpl_frame_get_filename(cur_frame_lamp) ;
        cpl_msg_info(__func__, "Reduce Frames %s (%s) ", 
                cr2res_get_base_name(cur_fname_lamp),
                cr2res_get_base_name(cur_fname_tw)) ;
        cpl_msg_indent_more() ;

        /* Loop over the detectors */
        for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {

            /* Initialise */
            slit_curv_map[det_nr-1] = NULL ;
            trace_wave[det_nr-1] = NULL ;
            lamp_image[det_nr-1] = NULL ;

            /* Store the extenѕion header for product saving */
            ext_plist[det_nr-1] = cpl_propertylist_load(cur_fname_tw,
                    cr2res_io_get_ext_idx(cur_fname_tw, det_nr, 1)) ;

            /* Compute only one detector */
            if (reduce_det != 0 && det_nr != reduce_det) continue ;

            cpl_msg_info(__func__, "Process detector number %d", det_nr) ;
            cpl_msg_indent_more() ;

            /* Load the TRACE_WAVE table of this detector */
            cpl_msg_info(__func__, "Load the TRACE_WAVE table") ;
            if ((trace_wave[det_nr-1] = cr2res_io_load_TRACE_WAVE(cur_fname_tw,
                            det_nr)) == NULL) {
                cpl_msg_error(__func__,"Failed to load table - skip detector");
                cpl_error_reset() ;
                cpl_msg_indent_less() ;
                continue ;
            }
            nb_traces = cpl_table_get_nrow(trace_wave[det_nr-1]) ;

            /* Load the lamp image of this detector */
            cpl_msg_info(__func__, "Load the LAMP image") ;
            if ((lamp_image[det_nr-1] = cr2res_io_load_image(cur_fname_lamp,
                            det_nr)) == NULL) {
                cpl_msg_error(__func__,"Failed to load image - skip detector");
                cpl_error_reset() ;
                cpl_msg_indent_less() ;
                continue ;
            }

            /* Loop over the traces and get the slit curvature */
            for (j=0 ; j<nb_traces ; j++) {
                /* Get Order and trace id */
                order = cpl_table_get(trace_wave[det_nr-1], 
                        CR2RES_COL_ORDER, j, NULL) ;
                trace_id = cpl_table_get(trace_wave[det_nr-1], 
                        CR2RES_COL_TRACENB, j, NULL) ;

                /* Check if this order needs to be skipped */
                if (reduce_order > -1 && order != reduce_order) continue ;

                /* Check if this trace needs to be skipped */
                if (reduce_trace > -1 && trace_id != reduce_trace) continue ;

                cpl_msg_info(__func__, "Process Order %d/Trace %d",
                        order, trace_id) ;
                cpl_msg_indent_more() ;
        
                /* Call the Slit Curvature Computation */
                /* TODO : Should those become parameters ? */
                int height = 100 ;
                int window = 15 ;
                int fit_c = 0 ;
                if (cr2res_slit_curv_compute_order_trace(lamp_image[det_nr-1], 
                            trace_wave[det_nr-1], order, trace_id,
                            height, window, curv_degree, fit_c, 
                            &slit_polya, &slit_polyb, &slit_polyc)) {
                    cpl_msg_warning(__func__, 
                            "Cannot Compute Slit curvature for Order %d",
                            order);
                    cpl_error_reset() ;
                    cpl_msg_indent_less() ;
                    continue ;
                }

                /* Fill the SLIT_CURVE_A/B/C for the current trace */
                slit_array = cr2res_convert_poly_to_array(slit_polya, 3) ;
                cpl_polynomial_delete(slit_polya) ;
                cpl_table_set_array(trace_wave[det_nr-1],
                        CR2RES_COL_SLIT_CURV_A, j, slit_array) ;
                cpl_array_delete(slit_array) ;
                slit_array = cr2res_convert_poly_to_array(slit_polyb, 3) ;
                cpl_polynomial_delete(slit_polyb) ;
                cpl_table_set_array(trace_wave[det_nr-1],
                        CR2RES_COL_SLIT_CURV_B, j, slit_array) ;
                cpl_array_delete(slit_array) ;
                slit_array = cr2res_convert_poly_to_array(slit_polyc, 3) ;
                cpl_polynomial_delete(slit_polyc) ;
                cpl_table_set_array(trace_wave[det_nr-1],
                        CR2RES_COL_SLIT_CURV_C, j, slit_array) ;
                cpl_array_delete(slit_array) ;

                cpl_msg_indent_less() ;
            }

            /* Generate the SLIT CURV Map */
            slit_curv_map[det_nr-1] =
                cr2res_slit_curv_gen_map(trace_wave[det_nr-1], reduce_order,
                        reduce_trace, 50, 0) ;
            cpl_msg_indent_less() ;
        }

        /* Save the SLIT_CURV_MAP */
        out_file = cpl_sprintf("%s_map.fits",
                cr2res_get_base_name(cr2res_get_root_name(cur_fname_tw)));
        cr2res_io_save_SLIT_CURV_MAP(out_file, frameset, frameset, parlist, 
                slit_curv_map, NULL, ext_plist, 
                CR2RES_UTIL_SLIT_CURV_MAP_PROCATG, RECIPE_STRING) ;
        cpl_free(out_file);

        /* Save the new TRACE_WAVE table */
        out_file=cpl_sprintf("%s_tw.fits", 
                cr2res_get_base_name(cr2res_get_root_name(cur_fname_tw)));
        cr2res_io_save_TRACE_WAVE(out_file, frameset, frameset, parlist, 
                trace_wave, NULL, ext_plist, CR2RES_UTIL_SLIT_CURV_TW_PROCATG, 
                RECIPE_STRING) ;
        cpl_free(out_file);

        /* Free and return */
        for (k=0 ; k<CR2RES_NB_DETECTORS ; k++) {
            if (trace_wave[k] != NULL) cpl_table_delete(trace_wave[k]) ;
            if (lamp_image[k] != NULL) hdrl_image_delete(lamp_image[k]) ;
            if (slit_curv_map[k] != NULL) 
                hdrl_image_delete(slit_curv_map[k]) ;
            if (ext_plist[k] != NULL)
                cpl_propertylist_delete(ext_plist[k]) ;
        }
    }
    cpl_frameset_delete(rawframes_lamp) ;
    cpl_frameset_delete(rawframes_tw) ;
    return (int)cpl_error_get_code();
}

