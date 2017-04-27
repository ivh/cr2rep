
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
#include "cr2res_extract.h"
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

static char cr2res_util_wave_description[] =
"TODO : Descripe here the recipe in / out / params / basic algo\n"
"trace_wave.fits " CR2RES_COMMAND_LINE "\n"
"extracted.fits " CR2RES_COMMAND_LINE "\n"
"catalog.fits " CR2RES_COMMAND_LINE "\n"
" The recipe produces the following products:\n"
"\n";

/*-----------------------------------------------------------------------------
                                Function code
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @defgroup cr2res_util_wave    Wavelength Calibration
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
                    "cr2res_util_wave",
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
    int                     reduce_det, reduce_order, reduce_trace ;
    cpl_frame           *   fr ;
    const char          *   trace_wave_file ;
    const char          *   extracted_file ;
    const char          *   catalog_file ;

    cpl_table           *   trace_wave_table ;
    cpl_table           *   extracted_table ;
    cpl_table           *   catalog_table ;
    cpl_vector          *   extracted_vec ;
    cpl_bivector        *   template ;
    cpl_polynomial      *   init_guess ;
    const cpl_array     *   init_guess_arr ;
    cpl_polynomial      *   wave_sol ;
    int                     det_nr, nb_traces, trace_id, order, i ;

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
    fr = cpl_frameset_get_position(frameset, 1);
    extracted_file = cpl_frame_get_filename(fr) ;
    fr = cpl_frameset_get_position(frameset, 2);
    catalog_file = cpl_frame_get_filename(fr) ;
    if (trace_wave_file==NULL || extracted_file==NULL || catalog_file==NULL) {
        cpl_msg_error(__func__, "The utility needs 3 files");
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Loop over the detectors */
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {

        /* Initialise */

        /* Compute only one detector */
        if (reduce_det != 0 && det_nr != reduce_det) continue ;

        cpl_msg_info(__func__, "Process detector number %d", det_nr) ;
        cpl_msg_indent_more() ;

        /* Load the TRACE_WAVE table of this detector */
        cpl_msg_info(__func__, "Load the TRACE_WAVE table") ;
        if ((trace_wave_table = cr2res_io_load_TRACE_WAVE(trace_wave_file,
                        det_nr)) == NULL) {
            cpl_msg_error(__func__,"Failed to load table - skip detector");
            cpl_error_reset() ;
            cpl_msg_indent_less() ;
            continue ;
        }
        nb_traces = cpl_table_get_nrow(trace_wave_table) ;

        /* Load the EXTRACT1D table of this detector */
        cpl_msg_info(__func__, "Load the EXTRACT1D table") ;
        if ((extracted_table = cr2res_io_load_EXTRACT_1D(extracted_file,
                        det_nr)) == NULL) {
            cpl_msg_error(__func__,"Failed to load table - skip detector");
            cpl_error_reset() ;
            cpl_msg_indent_less() ;
            continue ;
        }
 
        /* Allocate Data containers */

        /* Loop over the traces spectra */
        for (i=0 ; i<nb_traces ; i++) {
            /* Initialise */

            /* Get Order and trace id */
            order = cpl_table_get(trace_wave_table, "Order", i, NULL) ;
            trace_id = cpl_table_get(trace_wave_table, "TraceNb", i, NULL) ;

            /* Check if this order needs to be skipped */
            if (reduce_order > -1 && order != reduce_order) {
                continue ;
            }

            /* Check if this trace needs to be skipped */
            if (reduce_trace > -1 && trace_id != reduce_trace) {
                continue ;
            }

            cpl_msg_info(__func__, "Process Order %d/Trace %d",order,trace_id) ;
            cpl_msg_indent_more() ;

            /* Get the extracted spectrum */
            if ((extracted_vec = cr2res_extract_EXTRACT1D_get_spectrum(
                    extracted_table, order, trace_id)) == NULL) {
                cpl_msg_error(__func__, "Cannot get the extracted spectrum") ;
                cpl_msg_indent_less() ;
                continue ;
            }

            /* Get the initial guess */
            init_guess_arr = cpl_table_get_array(trace_wave_table, 
                    "Wavelength", i) ;
            init_guess = cr2res_convert_array_to_poly(init_guess_arr) ;

            /* Call the Wavelength Calibration */
            if ((wave_sol = cr2res_wave(extracted_vec, init_guess,
                            catalog_table, template)) == NULL) {
                cpl_msg_error(__func__, "Cannot calibrate in Wavelength") ;
                cpl_polynomial_delete(init_guess) ;
                cpl_vector_delete(extracted_vec) ;
                cpl_error_reset() ;
                cpl_msg_indent_less() ;
                continue ;
            }
            cpl_vector_delete(extracted_vec) ;
            cpl_polynomial_delete(init_guess) ;
            cpl_msg_indent_less() ;
        }
        cpl_table_delete(trace_wave_table) ;
        cpl_table_delete(extracted_table) ;

        /* Deallocate */
        cpl_msg_indent_less() ;
    }

    /* Save the Products */

    /* Free and return */
    return (int)cpl_error_get_code();
}

/**@}*/
