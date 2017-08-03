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

static int cr2res_util_trace_create(cpl_plugin *);
static int cr2res_util_trace_exec(cpl_plugin *);
static int cr2res_util_trace_destroy(cpl_plugin *);
static int cr2res_util_trace(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_util_trace_description[] =
"TODO : Descripe here the recipe in / out / params / basic algo\n"
"science.fits " CR2RES_COMMAND_LINE "\n"
" The recipe produces the following products:\n"
"\n";

/*-----------------------------------------------------------------------------
                                Function code
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @defgroup cr2res_util_trace   Trace Utility
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
                    "cr2res_util_trace",
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
            "cr2res.cr2res_util_trace", 10000);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "min_cluster");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_trace.smooth",
            CPL_TYPE_DOUBLE, "Length of the smoothing kernel",
            "cr2res.cr2res_util_trace", 5.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "smooth");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_trace.opening",
            CPL_TYPE_BOOL, "Use a morphological opening to rejoin clusters",
            "cr2res.cr2res_util_trace", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "opening");
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
    int                     min_cluster, degree, opening, reduce_det ;
    double                  smooth, y_pos ;
    int                 *   orders ;
    const char          *   flat_file ;
    char                *   out_file;
    cpl_image           *   flat_ima ;
    cpl_image           *   debug_ima ;
    int                     det_nr, order ;
    cpl_array           *   wl_array ;
    cpl_table           *   traces[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   ext_plist[CR2RES_NB_DETECTORS] ;
    int                     i, j, nb_orders, trace_nb, trace_id ;

    /* RETRIEVE INPUT PARAMETERS */
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_trace.min_cluster");
    min_cluster = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_trace.degree");
    degree = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_trace.smooth");
    smooth = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_trace.opening");
    opening = cpl_parameter_get_bool(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_trace.detector");
    reduce_det = cpl_parameter_get_int(param);

    /* Check Parameters */
    /* TODO */

    /* Identify the RAW and CALIB frames in the input frameset */
    if (cr2res_dfs_set_groups(frameset) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot identify RAW and CALIB frames") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Get Inputs */
    flat_file = cr2res_extract_filename(frameset, CR2RES_COMMAND_LINE) ;
    if (flat_file == NULL) {
        cpl_msg_error(__func__, "The utility needs a science file");
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Loop over the detectors */
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {

        /* Initialise */
        traces[det_nr-1] = NULL ;
        ext_plist[det_nr-1] = NULL ;

        /* Store the extenѕion header for product saving */
        ext_plist[det_nr-1] = cpl_propertylist_load(flat_file,
                cr2res_io_get_ext_idx(flat_file, det_nr)) ;

        /* Compute only one detector */
        if (reduce_det != 0 && det_nr != reduce_det) continue ;

        cpl_msg_info(__func__, "Process detector number %d", det_nr) ;
        cpl_msg_indent_more() ;

        /* Load the image in which the orders are to extract*/
        cpl_msg_info(__func__, "Load the Image") ;
        if ((flat_ima = cpl_image_load(flat_file, CPL_TYPE_FLOAT, 0, 
                        cr2res_io_get_ext_idx(flat_file, det_nr))) == NULL) {
            cpl_msg_warning(__func__,
                    "Cannot load the image - skip detector");
            cpl_error_reset() ;
            cpl_msg_indent_less() ;
            continue ;
        }

        /* Get the traces */
        cpl_msg_info(__func__, "Compute the traces") ;
        cpl_msg_indent_more() ;
        if ((traces[det_nr-1] = cr2res_trace(flat_ima,
                smooth, opening, degree, min_cluster)) == NULL) {
            cpl_msg_warning(__func__,
                    "Cannot compute the trace - skip detector");
            cpl_error_reset() ;
            cpl_image_delete(flat_ima) ;
            cpl_msg_indent_less() ;
            cpl_msg_indent_less() ;
            continue ;
        }
        cpl_msg_indent_less() ;

        /* Add The Order column using the header */
        cpl_table_new_column(traces[det_nr-1], CR2RES_COL_ORDER, CPL_TYPE_INT) ;

        /* Loop on the traces */
        for (i=0 ; i<cpl_table_get_nrow(traces[det_nr-1]) ; i++) {
            /* Get the current trace Y position */
            y_pos = cr2res_trace_get_trace_ypos(traces[det_nr-1], i) ;

            /* Compute the trace order from the header */
            order = cr2res_pfits_get_order(ext_plist[det_nr-1], y_pos) ;

            /* Store the Order in the table */
            cpl_table_set(traces[det_nr-1], CR2RES_COL_ORDER, i, order);
        }

        /* Add The TraceNb column */
        cpl_table_new_column(traces[det_nr-1], CR2RES_COL_TRACENB,CPL_TYPE_INT);

        orders = cr2res_trace_get_order_numbers(traces[det_nr-1], &nb_orders) ;
        for (i=0 ; i<nb_orders ; i++) {
            /* Initialise */
            trace_nb = 1 ;
            /* Loop on the traces */
            for (j=0 ; j<cpl_table_get_nrow(traces[det_nr-1]) ; j++) {
                if (cpl_table_get(traces[det_nr-1], CR2RES_COL_ORDER, j, 
                            NULL) == orders[i]) {
                    cpl_table_set(traces[det_nr-1], CR2RES_COL_TRACENB, j, 
                            trace_nb);
                    trace_nb ++ ;
                }
            }
        }
        cpl_free(orders) ;

        /* Add The Wavelength column using the header */
        cpl_table_new_column_array(traces[det_nr-1], CR2RES_COL_WAVELENGTH,
                CPL_TYPE_DOUBLE, 2) ;

        /* Loop on the traces */
        for (i=0 ; i<cpl_table_get_nrow(traces[det_nr-1]) ; i++) {
            /* Get the Order number */
            order = cpl_table_get(traces[det_nr-1], CR2RES_COL_ORDER, i, NULL) ;
            trace_id = cpl_table_get(traces[det_nr-1], CR2RES_COL_TRACENB, i, 
                    NULL) ;

            /* Get the Wavelength estimates from the header */
            if ((wl_array = cr2res_wave_get_estimate(flat_file, det_nr,
                            order)) == NULL) {
                cpl_msg_warning(__func__, 
                        "No Wavelength estimate for Detector %d / order %d",
                        det_nr, order) ;
                cpl_error_reset() ;
                cpl_table_set_array(traces[det_nr-1], CR2RES_COL_WAVELENGTH, i,
                        NULL);
            } else {
                /* Store the Wavelength in the table */
                cpl_table_set_array(traces[det_nr-1], CR2RES_COL_WAVELENGTH, i,
                        wl_array);
                cpl_array_delete(wl_array) ;
            }
        }

        /* Debug Image */
        if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
            debug_ima = cr2res_trace_gen_image(traces[det_nr-1],
                    cpl_image_get_size_x(flat_ima),
                    cpl_image_get_size_y(flat_ima)) ;
            cpl_image_save(debug_ima, "debug_trace_image.fits",
                    CPL_BPP_IEEE_FLOAT, NULL, CPL_IO_CREATE) ;
            cpl_image_delete(debug_ima) ;
        }
        cpl_image_delete(flat_ima) ;
        cpl_msg_indent_less() ;
    }

    /* Save the Products */
    out_file = cpl_sprintf("%s_trace.fits", 
            cr2res_get_base_name(cr2res_get_root_name(flat_file)));
    cpl_msg_debug(__func__, "Writing to %s",out_file);
    cr2res_io_save_TRACE_WAVE(out_file, frameset,
            parlist, traces, NULL, ext_plist, RECIPE_STRING) ;
    cpl_free(out_file);

    /* Free and return */
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
        if (ext_plist[det_nr-1] != NULL)
            cpl_propertylist_delete(ext_plist[det_nr-1]) ;
        if (traces[det_nr-1] != NULL)
            cpl_table_delete(traces[det_nr-1]) ;
    }

    return (int)cpl_error_get_code();
}

/**@}*/
