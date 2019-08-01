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

#include "cr2res_pfits.h"
#include "cr2res_dfs.h"
#include "cr2res_io.h"
#include "cr2res_splice.h"

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

#define RECIPE_STRING "cr2res_util_splice"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_util_splice_create(cpl_plugin *);
static int cr2res_util_splice_exec(cpl_plugin *);
static int cr2res_util_splice_destroy(cpl_plugin *);
static int cr2res_util_splice(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_util_splice_description[] =
"TODO : Descripe here the recipe in / out / params / basic algo\n"
"blaze.fits "CR2RES_CAL_FLAT_EXTRACT_1D_PROCATG  "\n"
"trace.fits "CR2RES_TW_PROTYPE  "\n"
"extracted.fits "CR2RES_UTIL_EXTRACT_1D_PROCATG  "\n"
" The recipe produces the following products:\n"
"\n";

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
                    "Splicing utility",
                    cr2res_util_splice_description,
                    "Ansgar Wehrhahn, Thomas Marquart, Yves Jung",
                    PACKAGE_BUGREPORT,
                    cr2res_get_license(),
                    cr2res_util_splice_create,
                    cr2res_util_splice_exec,
                    cr2res_util_splice_destroy)) {
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
static int cr2res_util_splice_create(cpl_plugin * plugin)
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
    p = cpl_parameter_new_value("cr2res.cr2res_util_splice.detector",
            CPL_TYPE_INT, "Only reduce the specified detector",
            "cr2res.cr2res_util_splice", 0);
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
static int cr2res_util_splice_exec(cpl_plugin * plugin)
{
    cpl_recipe  *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1;

    return cr2res_util_splice(recipe->frames, recipe->parameters);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_util_splice_destroy(cpl_plugin * plugin)
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
static int cr2res_util_splice(
        cpl_frameset            *   frameset,
        const cpl_parameterlist *   parlist)
{
    cpl_frameset        *   trace_fset ;
    cpl_frameset        *   extracted_fset ;
    cpl_frameset        *   blaze_fset ;
    const char          *   trace_file ;
    const char          *   extracted_file ;
    const char          *   blaze_file ;
    cpl_table           *   trace_table ;        
    cpl_table           *   extracted_table ;        
    cpl_table           *   blaze_table ;        
    cpl_table           **  trace_tables ;
    cpl_table           **  blaze_tables ;
    cpl_table           **  extracted_tables ;
    cpl_bivector        *   spliced ;
    cpl_bivector        *   spliced_err ;
    cpl_table           *   spliced_table ;
    cpl_error_code          ret ;
    cpl_propertylist    *   ext_plist ;
    char                *   out_file;
    cpl_size                nframes, i ;
    int                     ext_nr_trace, ext_nr_blaze, ext_nr_extracted, 
                            det_nr, nb_det ;

    /* RETRIEVE INPUT PARAMETERS */

    /* Check Parameters */

    /* Identify the RAW and CALIB frames in the input frameset */
    if (cr2res_dfs_set_groups(frameset) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot identify RAW and CALIB frames") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Get Inputs */
    trace_fset = cr2res_extract_frameset(frameset,
            CR2RES_TW_PROTYPE) ;
    extracted_fset = cr2res_extract_frameset(frameset,
            CR2RES_UTIL_EXTRACT_1D_PROCATG) ;
    blaze_fset = cr2res_extract_frameset(frameset, 
            CR2RES_CAL_FLAT_EXTRACT_1D_PROCATG) ;

    /* Tests the Inputs */
    if (trace_fset==NULL || extracted_fset==NULL || blaze_fset==NULL) {
        cpl_msg_error(__func__, "Missing Inputs") ;
        if (trace_fset != NULL) cpl_frameset_delete(trace_fset) ;
        if (extracted_fset != NULL) cpl_frameset_delete(extracted_fset) ;
        if (blaze_fset != NULL) cpl_frameset_delete(blaze_fset) ;
        return -1 ;
    }
    nframes = cpl_frameset_get_size(trace_fset) ;
    if (cpl_frameset_get_size(extracted_fset) != nframes ||
            cpl_frameset_get_size(blaze_fset) != nframes) {
        cpl_msg_error(__func__, "Inconsistent Inputs") ;
        cpl_frameset_delete(trace_fset) ;
        cpl_frameset_delete(extracted_fset) ;
        cpl_frameset_delete(blaze_fset) ;
        return -1 ;
    }

    if (nframes < 1) {
        cpl_msg_error(__func__, "No valid frames could be found") ;
        cpl_frameset_delete(trace_fset) ;
        cpl_frameset_delete(extracted_fset) ;
        cpl_frameset_delete(blaze_fset) ;
        return -1 ;
    }

    cpl_msg_info(__func__, 
            "Detected %"CPL_SIZE_FORMAT" BLAZE/EXTRACTED/TRACE triplet(s)",
            nframes) ;

    /* Count the number of detector to process */

    /* Loop over the Frames */
    nb_det = 0 ;
    for (i=0 ; i<nframes; i++) {

        /* Get the File names */
        trace_file = cpl_frame_get_filename(
                cpl_frameset_get_position(trace_fset, i)) ;
        extracted_file = cpl_frame_get_filename(
                cpl_frameset_get_position(extracted_fset, i)) ;
        blaze_file = cpl_frame_get_filename(
                cpl_frameset_get_position(blaze_fset, i)) ;

        /* Loop over the detectors */
        for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {

            /* Get Extension Numbers */
            ext_nr_trace = cr2res_io_get_ext_idx(trace_file, det_nr, 1) ;
            ext_nr_extracted = cr2res_io_get_ext_idx(extracted_file, det_nr, 1);
            ext_nr_blaze = cr2res_io_get_ext_idx(blaze_file, det_nr, 1) ;


            /* Skіp if any detector is missing */
            if (ext_nr_trace<0 || ext_nr_extracted<0 || ext_nr_blaze<0) 
                continue ;

            /* Try loading */
            trace_table = cr2res_io_load_TRACE_WAVE(trace_file, det_nr);
            blaze_table = cr2res_io_load_EXTRACT_1D(blaze_file, det_nr);
            extracted_table = cr2res_io_load_EXTRACT_1D(extracted_file, det_nr);

            /* Increase counter if ok */
            if (trace_table != NULL && blaze_table != NULL && 
                    extracted_table != NULL) {
                nb_det++ ;
            } else {
                cpl_error_reset() ;
            }

            if (trace_table!= NULL) cpl_table_delete(trace_table);
            if (blaze_table!= NULL) cpl_table_delete(blaze_table);
            if (extracted_table!= NULL) cpl_table_delete(extracted_table);
        }
    }

    if (nb_det == 0) {
        cpl_msg_error(__func__, "No valid data could be found") ;
        cpl_frameset_delete(trace_fset) ;
        cpl_frameset_delete(extracted_fset) ;
        cpl_frameset_delete(blaze_fset) ;
        return -1 ;
    }

    cpl_msg_info(__func__, "Loading %d Detector(s) for splicing", nb_det) ;

    /* Allocate Data containers */
    blaze_tables = cpl_malloc(nb_det * sizeof(cpl_table*)) ;
    extracted_tables = cpl_malloc(nb_det * sizeof(cpl_table*)) ;
    trace_tables = cpl_malloc(nb_det * sizeof(cpl_table*)) ;

    /* Load the inputs */
    /* Loop over the Frames */
    nb_det = 0 ;
    for (i=0 ; i<nframes; i++) {
        /* Get the File names */
        trace_file = cpl_frame_get_filename(
                cpl_frameset_get_position(trace_fset, i)) ;
        extracted_file = cpl_frame_get_filename(
                cpl_frameset_get_position(extracted_fset, i)) ;
        blaze_file = cpl_frame_get_filename(
                cpl_frameset_get_position(blaze_fset, i)) ;

        /* Loop over the detectors */
        for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {

            /* Get Extension Numbers */
            ext_nr_trace = cr2res_io_get_ext_idx(trace_file, det_nr, 1) ;
            ext_nr_extracted = cr2res_io_get_ext_idx(extracted_file, det_nr, 1);
            ext_nr_blaze = cr2res_io_get_ext_idx(blaze_file, det_nr, 1) ;

            /* Skіp if any detector is missing */
            if (ext_nr_trace<0 || ext_nr_extracted<0 || ext_nr_blaze<0) 
                continue ;

            /* Load */
            trace_table = cr2res_io_load_TRACE_WAVE(trace_file, det_nr);
            blaze_table = cr2res_io_load_EXTRACT_1D(blaze_file, det_nr);
            extracted_table = cr2res_io_load_EXTRACT_1D(extracted_file, det_nr);

            /* Increase counter if ok */
            if (trace_table!=NULL && blaze_table!=NULL && 
                    extracted_table!=NULL) {
                trace_tables[nb_det] = trace_table;
                extracted_tables[nb_det] = extracted_table;
                blaze_tables[nb_det] = blaze_table;
                nb_det++ ;
            } else {
                cpl_error_reset() ;
                if (trace_table!= NULL) cpl_table_delete(trace_table);
                if (blaze_table!= NULL) cpl_table_delete(blaze_table);
                if (extracted_table!= NULL) cpl_table_delete(extracted_table);
            }
        }
    }
    cpl_frameset_delete(trace_fset) ;
    cpl_frameset_delete(blaze_fset) ;

    /* Call the splicing */

    ret = cr2res_splice(extracted_tables, blaze_tables, trace_tables,
            nb_det, &spliced, &spliced_err) ;


    /* Deallocate */
    for (i=0 ; i<nb_det ; i++) {
        cpl_table_delete(blaze_tables[i]) ;
        cpl_table_delete(extracted_tables[i]) ;
        cpl_table_delete(trace_tables[i]) ;
    }
    cpl_free(blaze_tables) ;
    cpl_free(extracted_tables) ;
    cpl_free(trace_tables) ;

    if (ret) {
        cpl_msg_error(__func__, "Splicing Error"); 
        spliced_table = NULL ;
    } else {
        /* Store the spliced in a SPLICED table */
        spliced_table = cr2res_splice_SPLICED_1D_create(spliced, spliced_err) ;

        /* Free the bivectors */
        cpl_bivector_delete(spliced) ;
        cpl_bivector_delete(spliced_err) ;

        /* Save the table */
        if (spliced_table == NULL) {
            cpl_msg_error(__func__, "Failed to create the SPLICED table");
        } else {
            extracted_file = cpl_frame_get_filename(
                    cpl_frameset_get_position(extracted_fset, 0)) ;
            ext_plist = cpl_propertylist_load(extracted_file, 1) ; 
            out_file = cpl_sprintf("%s_spliced.fits",
                    cr2res_get_base_name(cr2res_get_root_name(extracted_file)));
            cr2res_io_save_SPLICED_1D(out_file, frameset,  frameset, parlist, 
                    spliced_table, NULL, ext_plist, 
                    CR2RES_UTIL_SPLICE_SPLICED_1D_PROCATG, RECIPE_STRING) ;
            cpl_free(out_file);
            cpl_propertylist_delete(ext_plist);
            cpl_table_delete(spliced_table) ;
        }
    }

    /* Free and return */
    cpl_frameset_delete(extracted_fset) ;
    return (int)cpl_error_get_code();
}

