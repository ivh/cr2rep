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
#include "cr2res_slitdec.h"

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

#define RECIPE_STRING "cr2res_util_extract"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/
static cpl_table * cr2res_extract_tab_create(
        cpl_vector      **  spectrum,
        int             *   orders,
        int                 nborders) ;
static cpl_table * cr2res_slit_func_tab_create(
        cpl_vector      **  slit_func,
        int             *   orders,
        int                 nborders) ;
static int cr2res_util_extract_create(cpl_plugin *);
static int cr2res_util_extract_exec(cpl_plugin *);
static int cr2res_util_extract_destroy(cpl_plugin *);
static int cr2res_util_extract(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_util_extract_description[] =
"TODO : Descripe here the recipe in / out / params / basic algo\n"
"science.fits " CR2RES_SCI_1D_RAW "\n"
"trace.fits " CR2RES_TRACE_OPEN_PROCATG "\n"
" The recipe produces the following products:\n"
"\n";

/*-----------------------------------------------------------------------------
                                Function code
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @defgroup cr2res_util_extract 	Optimal Extraction Utility
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
                    "cr2res_util_extract",
                    "Extraction utility",
                    cr2res_util_extract_description,
                    "Thomas Marquart, Yves Jung",
                    PACKAGE_BUGREPORT,
                    cr2res_get_license(),
                    cr2res_util_extract_create,
                    cr2res_util_extract_exec,
                    cr2res_util_extract_destroy)) {
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
static int cr2res_util_extract_create(cpl_plugin * plugin)
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
    p = cpl_parameter_new_value("cr2res.cr2res_util_extract.oversample",
            CPL_TYPE_INT, "factor by which to oversample the extraction",
            "cr2res.cr2res_util_extract", 10);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "oversample");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_extract.swath_width",
            CPL_TYPE_INT, "The swath width", "cr2res.cr2res_util_extract", 256);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "swath_width");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_extract.smooth_slit",
            CPL_TYPE_DOUBLE, "Smoothing along the slit",
            "cr2res.cr2res_util_extract", 0.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "smooth_slit");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_extract.sum_only",
            CPL_TYPE_BOOL, "Flag to only sum along detector",
            "cr2res.cr2res_util_extract", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "sum_only");
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
static int cr2res_util_extract_exec(cpl_plugin * plugin)
{
    cpl_recipe  *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1;

    return cr2res_util_extract(recipe->frames, recipe->parameters);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_util_extract_destroy(cpl_plugin * plugin)
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
static int cr2res_util_extract(
        cpl_frameset            *   frameset,
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter *   param;
    int                     oversample, swath_width, sum_only ;
    int                     extr_height;
    double                  smooth_slit ;
    const char          *   science_file ;
    const char          *   trace_file ;
    cpl_table           *   trace_table ;
    int                 *   orders ;
    int                     nb_orders[CR2RES_NB_DETECTORS] ;
    cpl_image           *   science_ima ;
    cpl_polynomial      **  traces ;
    cpl_vector          *   y_center ;
    hdrl_image          *   model_master[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   model_tmp ;
    cpl_vector          **  spectrum[CR2RES_NB_DETECTORS] ;
    cpl_vector          **  spectrum_error[CR2RES_NB_DETECTORS] ;
    cpl_vector          **  slit_func[CR2RES_NB_DETECTORS] ;
    cpl_table           *   slit_func_tab[CR2RES_NB_DETECTORS] ;
    cpl_table           *   extract_tab[CR2RES_NB_DETECTORS] ;
    int                     det_nr, i ;

    /* RETRIEVE INPUT PARAMETERS */
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_extract.oversample");
    oversample = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_extract.swath_width");
    swath_width = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_extract.smooth_slit");
    smooth_slit = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_extract.sum_only");
    sum_only = cpl_parameter_get_bool(param);

    /* Check Parameters */
    /* TODO */

    /* Identify the RAW and CALIB frames in the input frameset */
    if (cr2res_dfs_set_groups(frameset) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot identify RAW and CALIB frames") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Get Inputs */
    science_file = cr2res_extract_filename(frameset, CR2RES_SCI_1D_RAW) ;
    trace_file = cr2res_extract_filename(frameset, CR2RES_TRACE_OPEN_PROCATG);
    if (science_file == NULL || trace_file == NULL) {
        cpl_msg_error(__func__, "The utility needs a science file and a trace");
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Loop over the detectors */
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
        cpl_msg_info(__func__, "Process detector number %d", det_nr) ;
        cpl_msg_indent_more() ;

        /* Load the trace table of this detector */
        trace_table = cr2res_io_load_TRACE_OPEN(trace_file, det_nr) ;
        if (trace_table == NULL) {
            cpl_msg_error(__func__, "Failed to get trace table for detector #%d", det_nr);
            cpl_error_set(__func__, CPL_ERROR_CONTINUE) ;
            continue ;
        }

        /* Get the list of orders in the trace table */
        orders = cr2res_trace_get_order_numbers(trace_table,
                &(nb_orders[det_nr-1])) ;

        /* Load the image in which the orders are to extract*/
        science_ima = cpl_image_load(science_file, CPL_TYPE_FLOAT, 0, det_nr) ;

        /* Allocate Data containers */
        spectrum[det_nr-1] = cpl_malloc(nb_orders[det_nr-1] *
                sizeof(cpl_vector *)) ;
        slit_func[det_nr-1] = cpl_malloc(nb_orders[det_nr-1] *
                sizeof(cpl_vector *)) ;
        model_master[det_nr-1] = hdrl_image_create(science_ima, NULL) ;
        hdrl_image_mul_scalar(model_master[det_nr-1], (hdrl_value){0.0, 0.0}) ;

        /* Loop over the orders and extract them */
        for (i=0 ; i<nb_orders[det_nr-1] ; i++) {
            cpl_msg_info(__func__, "Process order number %d", orders[i]) ;
            cpl_msg_indent_more() ;

            /* Get the 2 Traces for the current order */
            traces = cr2es_trace_open_get_polynomials(trace_table, orders[i]) ;

            /* Get the values between the 2 traces and the height */
            y_center = cr2res_trace_compute_middle(traces[0], traces[1],
                    cpl_image_get_size_x(science_ima)) ;
            extr_height = cr2res_trace_compute_height(traces[0], traces[1],
                    cpl_image_get_size_x(science_ima)) ;
            cpl_polynomial_delete(traces[0]) ;
            cpl_polynomial_delete(traces[1]) ;
            cpl_free(traces) ;

            if (cr2res_slitdec_vert(science_ima, y_center, extr_height,
                    swath_width, oversample, smooth_slit,
                    &(slit_func[det_nr-1][i]),
                    &(spectrum[det_nr-1][i]),
                    &model_tmp) != 0) {
                cpl_msg_error(__func__,
                        "Cannot extract order %d on detector %d",
                        orders[i], det_nr) ;
                slit_func[det_nr-1][i] = NULL ;
                spectrum[det_nr-1][i] = NULL ;
                model_tmp = NULL ;
            }
            cpl_vector_delete(y_center) ;

            /* Update the model global image */
            if (model_tmp != NULL) {
                hdrl_image_add_image(model_master[det_nr-1], model_tmp) ;
                hdrl_image_delete(model_tmp) ;
            }
            cpl_msg_indent_less() ;
        }
        cpl_image_delete(science_ima) ;
        cpl_table_delete(trace_table) ;

        /* Create the slit_func_tab for the current detector */
        slit_func_tab[det_nr-1] = cr2res_slit_func_tab_create(
                slit_func[det_nr-1], orders, nb_orders[det_nr-1]) ;

        /* Create the extracted_tab for the current detector */
        extract_tab[det_nr-1] = cr2res_extract_tab_create(
                spectrum[det_nr-1], orders, nb_orders[det_nr-1]) ;

		/* Deallocate Vectors */
        for (i=0 ; i<nb_orders[det_nr-1] ; i++) {
            if (slit_func[det_nr-1][i] != NULL)
                cpl_vector_delete(slit_func[det_nr-1][i]) ;
            if (spectrum[det_nr-1][i] != NULL)
                cpl_vector_delete(spectrum[det_nr-1][i]) ;
        }
        cpl_free(spectrum[det_nr-1]) ;
        cpl_free(slit_func[det_nr-1]) ;
        cpl_free(orders) ;
        cpl_msg_indent_less() ;
    }

    /* Save the Products */
    cr2res_io_save_SLIT_MODEL("cr2res_util_extract_model.fits", frameset,
            parlist, model_master, NULL, RECIPE_STRING) ;
    cr2res_io_save_SLIT_FUNC("cr2res_util_extract_slit_func.fits", frameset,
            parlist, slit_func_tab, NULL, RECIPE_STRING) ;
    cr2res_io_save_EXTRACT_1D("cr2res_util_extract_extract_1D.fits", frameset,
            parlist, extract_tab, NULL, RECIPE_STRING) ;

    /* Free and return */
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
        cpl_table_delete(slit_func_tab[det_nr-1]) ;
        cpl_table_delete(extract_tab[det_nr-1]) ;
        hdrl_image_delete(model_master[det_nr-1]) ;
    }
    return (int)cpl_error_get_code();
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Create the extract 1D table to be saved
  @param    spectrum   	The extracted spectra of the different orders
  @param    orders      The orders numbers
  @param    nborders    The number of orders
  @return   the extract_1D table or NULL
 */
/*----------------------------------------------------------------------------*/
static cpl_table * cr2res_extract_tab_create(
        cpl_vector      **  spectrum,
        int             *   orders,
        int                 nborders)
{
    cpl_table       *   out ;
    char            *   col_name ;
    const double    *   pspec ;
    int                 nrows, all_null, i ;

    /* Check entries */
    if (spectrum == NULL) return NULL ;
    if (nborders < 1) return NULL ;

    /* Check the all vectorѕ are not null */
    all_null = 1 ;
    for (i=0 ; i<nborders ; i++)
        if (spectrum[i] != NULL) {
            nrows = cpl_vector_get_size(spectrum[i]) ;
            all_null = 0 ;
        }
    if (all_null == 1) return NULL ;

    /* Check the sizes */
    for (i=0 ; i<nborders ; i++)
        if (spectrum[i] != NULL && cpl_vector_get_size(spectrum[i]) != nrows)
            return NULL ;

    /* Create the table */
    out = cpl_table_new(nrows);
    for (i=0 ; i<nborders ; i++) {
        col_name = cpl_sprintf("%02d_SPEC", orders[i]) ;
        cpl_table_new_column(out, col_name, CPL_TYPE_DOUBLE);
        cpl_free(col_name) ;
    }

    /* Fill the table */
    for (i=0 ; i<nborders ; i++) {
        pspec = cpl_vector_get_data_const(spectrum[i]) ;
        col_name = cpl_sprintf("%02d_SPEC", orders[i]) ;
        cpl_table_copy_data_double(out, col_name, pspec) ;
        cpl_free(col_name) ;
    }
    return out ;
}
/*----------------------------------------------------------------------------*/
/**
  @brief    Create the slit functions table to be saved
  @param    slit_func   The slit functions of the different orders
  @param    orders      The orders numbers
  @param    nborders    The number of orders
  @return   the slit_func table or NULL
 */
/*----------------------------------------------------------------------------*/
static cpl_table * cr2res_slit_func_tab_create(
        cpl_vector      **  slit_func,
        int             *   orders,
        int                 nborders)
{
    cpl_table       *   out ;
    char            *   col_name ;
    const double    *   pslit ;
    int                 nrows, all_null, i ;

    /* Check entries */
    if (slit_func == NULL) return NULL ;
    if (nborders < 1) return NULL ;

    /* Check the all vectorѕ are not null */
    all_null = 1 ;
    for (i=0 ; i<nborders ; i++)
        if (slit_func[i] != NULL) {
            nrows = cpl_vector_get_size(slit_func[i]) ;
            all_null = 0 ;
        }
    if (all_null == 1) return NULL ;

    /* Check the sizes */
    for (i=0 ; i<nborders ; i++)
        if (slit_func[i] != NULL && cpl_vector_get_size(slit_func[i]) != nrows)
            return NULL ;

    /* Create the table */
    out = cpl_table_new(nrows);
    for (i=0 ; i<nborders ; i++) {
        col_name = cpl_sprintf("%02d_SLIT_FUNC", orders[i]) ;
        cpl_table_new_column(out, col_name, CPL_TYPE_DOUBLE);
        cpl_free(col_name) ;
    }

    /* Fill the table */
    for (i=0 ; i<nborders ; i++) {
        pslit = cpl_vector_get_data_const(slit_func[i]) ;
        col_name = cpl_sprintf("%02d_SLIT_FUNC", orders[i]) ;
        cpl_table_copy_data_double(out, col_name, pslit) ;
        cpl_free(col_name) ;
    }
    return out ;
}
