
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

static char cr2res_util_slit_curv_description[] =
"The utility expects 1 file as input:\n"
"   * trace_wave.fits " CR2RES_COMMAND_LINE "\n"
"The slit curvature is derived from each order with more than 1 trace.\n"
"The recipe produces the following products:\n"
"   * SLIT_CURV\n"
"\n";

/*-----------------------------------------------------------------------------
                                Function code
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @defgroup cr2res_util_slit_curv    Slit Curvature
 */
/*----------------------------------------------------------------------------*/

/**@{*/

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
                    "cr2res_util_slit_curv",
                    "Slit Curvature utility",
                    cr2res_util_slit_curv_description,
                    "Thomas Marquart, Yves Jung",
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
            CPL_TYPE_BOOL, "Flag for display",
            "cr2res.cr2res_util_slit_curv", FALSE);
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
    cpl_frame           *   fr ;
    cpl_polynomial      **  curvatures ;
    const char          *   trace_wave_file ;
    cpl_table           *   trace_wave[CR2RES_NB_DETECTORS] ;
    cpl_table           *   slit_curv[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   ext_plist[CR2RES_NB_DETECTORS] ;
    cpl_polynomial      *   slit_polya ;
    cpl_polynomial      *   slit_polyb ;
    cpl_polynomial      *   slit_polyc ;
    cpl_array           *   slit_array ;
    int                     det_nr, order, trace_id, nb_traces, curv_degree ;
    char                *   col_name ;
    char                *   out_file;
    cpl_array           *   curv_array ;
    int                     i, j ;

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
    display = cpl_parameter_get_bool(param) ;
 
    /* Check Parameters */
    /* TODO */

    /* Identify the RAW and CALIB frames in the input frameset */
    if (cr2res_dfs_set_groups(frameset) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot identify RAW and CALIB frames") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Get Inputs */
    fr = cpl_frameset_get_position(frameset, 0);
    trace_wave_file = cpl_frame_get_filename(fr) ;
    if (trace_wave_file==NULL) {
        cpl_msg_error(__func__, "The utility needs at least 1 file as input");
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Loop over the detectors */
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {

        /* Initialise */
        slit_curv[det_nr-1] = cpl_table_new(CR2RES_DETECTOR_SIZE) ;
        ext_plist[det_nr-1] = NULL ;
        trace_wave[det_nr-1] = NULL ;

        /* Store the extenѕion header for product saving */
        ext_plist[det_nr-1] = cpl_propertylist_load(trace_wave_file,
                cr2res_io_get_ext_idx(trace_wave_file, det_nr, 1)) ;

        /* Compute only one detector */
        if (reduce_det != 0 && det_nr != reduce_det) continue ;

        cpl_msg_info(__func__, "Process detector number %d", det_nr) ;
        cpl_msg_indent_more() ;

        /* Load the TRACE_WAVE table of this detector */
        cpl_msg_info(__func__, "Load the TRACE_WAVE table") ;
        if ((trace_wave[det_nr-1] = cr2res_io_load_TRACE_WAVE(trace_wave_file,
                        det_nr)) == NULL) {
            cpl_msg_error(__func__,"Failed to load table - skip detector");
            cpl_error_reset() ;
            cpl_msg_indent_less() ;
            continue ;
        }
        nb_traces = cpl_table_get_nrow(trace_wave[det_nr-1]) ;

        /* Allocate Data containers */

        /* Loop over the traces and get the slit curvature */
        for (i=0 ; i<nb_traces ; i++) {
            /* Get Order and trace id */
            order = cpl_table_get(trace_wave[det_nr-1], CR2RES_COL_ORDER, i, 
                    NULL) ;
            trace_id = cpl_table_get(trace_wave[det_nr-1], CR2RES_COL_TRACENB, 
                    i, NULL) ;

            /* Check if this order needs to be skipped */
            if (reduce_order > -1 && order != reduce_order) {
                cpl_msg_indent_less() ;
                continue ;
            }

            /* Check if this trace needs to be skipped */
            if (reduce_trace > -1 && trace_id != reduce_trace) {
                cpl_msg_indent_less() ;
                continue ;
            }

            cpl_msg_info(__func__, "Process Order %d/Trace %d",order,trace_id) ;
            cpl_msg_indent_more() ;
	
            /* Call the Slit Curvature Computation */
            if ((curvatures = cr2res_slit_curv_compute_order_trace(
                            trace_wave[det_nr-1], order, trace_id, 
                            display, curv_degree)) == NULL) {
                cpl_msg_warning(__func__, 
                        "Cannot Compute Slit curvature for Order %d",
                        order);
                cpl_error_reset() ;
                cpl_msg_indent_less() ;
                continue ;
            }

            /* Fill the SLIT_CURVE_A/B/C for the current trace */
            if (cr2res_slit_curv_fit_coefficients(curvatures,
                        CR2RES_DETECTOR_SIZE,
                        &slit_polya, &slit_polyb, &slit_polyc) == 0) {
                slit_array = cr2res_convert_poly_to_array(slit_polya, 3) ;
                cpl_polynomial_delete(slit_polya) ;
                cpl_table_set_array(trace_wave[det_nr-1],
                        CR2RES_COL_SLIT_CURV_A, i, slit_array) ;
                cpl_array_delete(slit_array) ;
                slit_array = cr2res_convert_poly_to_array(slit_polyb, 3) ;
                cpl_polynomial_delete(slit_polyb) ;
                cpl_table_set_array(trace_wave[det_nr-1],
                        CR2RES_COL_SLIT_CURV_B, i, slit_array) ;
                cpl_array_delete(slit_array) ;
                slit_array = cr2res_convert_poly_to_array(slit_polyc, 3) ;
                cpl_polynomial_delete(slit_polyc) ;
                cpl_table_set_array(trace_wave[det_nr-1],
                        CR2RES_COL_SLIT_CURV_C, i, slit_array) ;
                cpl_array_delete(slit_array) ;
            } 

            /* Store the Solution in the table SLIT_CURV */
            col_name = cr2res_dfs_SLIT_CURV_colname(order, trace_id) ;
            cpl_table_new_column_array(slit_curv[det_nr-1], col_name, 
                    CPL_TYPE_DOUBLE, curv_degree+1) ; 

            /* Loop on the Column rows */
            for (j=0 ; j<CR2RES_DETECTOR_SIZE ; j++) {
                if (curvatures[j] != NULL) {
                    /* Ѕtore the polynomial in the table */
                    curv_array=cr2res_convert_poly_to_array(curvatures[j],
                            curv_degree+1) ;
                    cpl_polynomial_delete(curvatures[j]) ;
                    if (curv_array != NULL) {
                        cpl_table_set_array(slit_curv[det_nr-1], col_name, j, 
                                curv_array) ;
                        cpl_array_delete(curv_array) ;
                    }
                }
            }
            cpl_free(col_name) ;
            cpl_free(curvatures) ; 
            cpl_msg_indent_less() ;
        }
        cpl_msg_indent_less() ;
    }

    /* Save the new SLIT_CURV table */
    out_file=cpl_sprintf("%s_slit_curv.fits", 
            cr2res_get_base_name(cr2res_get_root_name(trace_wave_file)));
    cr2res_io_save_SLIT_CURV(out_file, frameset, parlist, slit_curv, NULL, 
            ext_plist, CR2RES_UTIL_SLIT_CURV_PROCATG, RECIPE_STRING) ;
    cpl_free(out_file);

    /* Save the new TRACE_WAVE table */
    out_file=cpl_sprintf("%s_trace_wave.fits", 
            cr2res_get_base_name(cr2res_get_root_name(trace_wave_file)));
    cr2res_io_save_TRACE_WAVE(out_file, frameset, parlist, trace_wave, NULL, 
            ext_plist, CR2RES_UTIL_SLIT_CURV_TRACE_WAVE_PROCATG,RECIPE_STRING) ;
    cpl_free(out_file);

    /* Free and return */
    for (i=0 ; i<CR2RES_NB_DETECTORS ; i++) {
        if (slit_curv[i] != NULL) cpl_table_delete(slit_curv[i]) ;
        if (trace_wave[i] != NULL) cpl_table_delete(trace_wave[i]) ;
        if (ext_plist[i] != NULL)
            cpl_propertylist_delete(ext_plist[i]) ;
    }

    return (int)cpl_error_get_code();
}

/**@}*/
