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

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

#define RECIPE_STRING "cr2res_util_genlines"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_util_genlines_save(cpl_table *, const cpl_parameterlist *,
                cpl_frameset *);
static int cr2res_util_genlines_create(cpl_plugin *);
static int cr2res_util_genlines_exec(cpl_plugin *);
static int cr2res_util_genlines_destroy(cpl_plugin *);
static int cr2res_util_genlines(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_util_genlines_description[] =
"This recipe is used to generate spectrum calibration tables.\n"
"The sof file contains the names of the input ASCII file\n"
"tagged with "CR2RES_COMMAND_LINE".\n"
"The ASCII file must contain two columns:\n"
"1st: Wavelengths in increasing order (the unit is corrected by\n"
"     the factor option to obtain nanometers).\n"
"2nd: The atmospheric emission.\n"
"The ASCII files are in the catalogs/ directory of the CR2RES distribution.\n"
"This recipe produces 1 file:\n"
"First product:     the table with the lines.\n"
"                   (PRO TYPE = "CR2RES_PROTYPE_CATALOG")\n" 
"                   (PRO CATG = "CR2RES_EMMISION_LINES_PROCATG")\n" ;

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
                    "cr2res_util_genlines",
                    "Generate spectrum calibration FITS tables",
                    cr2res_util_genlines_description,
                    "Thomas Marquart, Yves Jung",
                    PACKAGE_BUGREPORT,
                    cr2res_get_license(),
                    cr2res_util_genlines_create,
                    cr2res_util_genlines_exec,
                    cr2res_util_genlines_destroy)) {    
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
static int cr2res_util_genlines_create(cpl_plugin * plugin)
{
    cpl_recipe          *   recipe ;
    cpl_parameter       *   p ;
    cpl_parameterlist   *   collapse_par ;
    hdrl_parameter      *   sigclip_def ;
    hdrl_parameter      *   minmax_def ;

    /* Check that the plugin is part of a valid recipe */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else
        return -1;

    /* Create the parameters list in the cpl_recipe object */
    recipe->parameters = cpl_parameterlist_new();

    /* Fill the parameters list */
    p = cpl_parameter_new_value("cr2res_util_genlines.wl_factor", 
            CPL_TYPE_DOUBLE, "The factor used to multiply the wl",
            "cr2res_util_genlines", 1.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "wl_factor");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res_util_genlines.display", 
            CPL_TYPE_BOOL, "Flag to plot", "cr2res_util_genlines", FALSE);
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
static int cr2res_util_genlines_exec(cpl_plugin * plugin)
{
    cpl_recipe  *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1;

    return cr2res_util_genlines(recipe->frames, recipe->parameters);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_util_genlines_destroy(cpl_plugin * plugin)
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
static int cr2res_util_genlines(
        cpl_frameset            *   frameset,
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter *   par ;
    double                  wl_fac ;
    int                     display ;
    cpl_frame           *   cur_frame ;
    cpl_bivector        *   bivec ;
    double              *   pbivec_x ;
    double              *   pbivec_y ;
    cpl_bivector        *   bivec_fill ;
    double              *   pbivec_fill_x ;
    double              *   pbivec_fill_y ;
    int                     nvals ;
    double                  wavel ;
    cpl_table           *   tab ;
    int                     i ;

    /* Retrieve input parameters */
    par=cpl_parameterlist_find_const(parlist, "cr2res_util_genlines.display");
    display = cpl_parameter_get_bool(par);
    par=cpl_parameterlist_find_const(parlist, "cr2res_util_genlines.wl_factor");
    wl_fac = cpl_parameter_get_double(par);

    /* Identify the RAW and CALIB frames in the input frameset */
    if (cr2res_dfs_set_groups(frameset)) {
        cpl_msg_error(__func__, "Cannot identify RAW and CALIB frames") ;
        return -1 ;
    }

    /* Load the file */
    cur_frame = cpl_frameset_get_position(frameset, 0) ;
    if ((bivec=cpl_bivector_read(cpl_frame_get_filename(cur_frame)))==NULL) {
        cpl_msg_error(__func__, "Cannot load the file in the bivector") ;
        return -1 ;
    }
    nvals = cpl_bivector_get_size(bivec) ;

    /* Use wl_factor */
    cpl_vector_multiply_scalar(cpl_bivector_get_x(bivec), wl_fac) ;

    /* Display if requested */
    if (display) {
        cpl_plot_bivector(
                "set grid;set xlabel 'Wavelength (nm)';set ylabel 'Emission';",
                "t 'Catalog lines' w lines", "", bivec);
    }

    /* Allocate the data container */
    tab = cpl_table_new(nvals) ;
    cpl_table_wrap_double(tab, cpl_bivector_get_x_data(bivec),
            CR2RES_COL_WAVELENGTH) ;
    cpl_table_wrap_double(tab, cpl_bivector_get_y_data(bivec),
            CR2RES_COL_EMISSION) ;

    /* Save the table */
    cpl_msg_info(__func__, "Saving the table with %d rows", nvals) ;
    if (cr2res_util_genlines_save(tab, parlist, frameset) == -1) {
        cpl_msg_error(__func__, "Cannot write the table") ;
        cpl_bivector_delete(bivec) ;
        cpl_table_unwrap(tab, CR2RES_COL_WAVELENGTH) ;
        cpl_table_unwrap(tab, CR2RES_COL_EMISSION) ;
        cpl_table_delete(tab) ;
        return -1 ;
    }
    cpl_bivector_delete(bivec) ;
    cpl_table_unwrap(tab, CR2RES_COL_WAVELENGTH) ;
    cpl_table_unwrap(tab, CR2RES_COL_EMISSION) ;
    cpl_table_delete(tab) ;
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save the product of the recipe
  @param    out_table   the table 
  @param    parlist     the input list of parameters
  @param    set         the input frame set
  @return   0 if everything is ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
static int cr2res_util_genlines_save(
        cpl_table               *   out_table,
        const cpl_parameterlist *   parlist,
        cpl_frameset            *   set)
{
    cpl_propertylist    *   plist ;

    plist = cpl_propertylist_new();
    cpl_propertylist_append_string(plist, "INSTRUME", "CR2RES") ;
    cpl_propertylist_append_string(plist, CPL_DFS_PRO_CATG, 
            CR2RES_EMMISION_LINES_PROCATG) ;
    cpl_propertylist_append_string(plist, CPL_DFS_PRO_TYPE,
            CR2RES_PROTYPE_CATALOG) ;

    if (cpl_dfs_save_table(set, NULL, parlist, set, NULL, out_table,
                NULL, "cr2res_util_genlines", plist, NULL,
                PACKAGE "/" PACKAGE_VERSION,
                "cr2res_util_genlines.fits") != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot save the table") ;
        return -1 ;
    }
    cpl_propertylist_delete(plist) ;

    /* Return */
    return 0 ;
}

