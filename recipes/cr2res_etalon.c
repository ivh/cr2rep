/*
 * This file is part of the CR2RE Pipeline
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

#include "cr2re_utils.h"
#include "cr2re_etalon.h"
#include "cr2re_pfits.h"
#include "cr2re_dfs.h"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_etalon_create(cpl_plugin *);
static int cr2res_etalon_exec(cpl_plugin *);
static int cr2res_etalon_destroy(cpl_plugin *);
static int cr2res_etalon(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_etalon_description[] =
"This example text is used to describe the recipe.\n"
"The description should include the required FITS-files and\n"
"their associated tags, e.g.\n"
"raw-file.fits " CR2RE_ETALON_RAW "\n"
"\n"
"Additionally, it should describe functionality of the expected output."
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
                    CR2RE_BINARY_VERSION,
                    CPL_PLUGIN_TYPE_RECIPE,
                    "cr2res_etalon",
                    "Etalon Test Programm",
                    cr2res_etalon_description,
                    "Thomas Marquart, Yves Jung",
                    PACKAGE_BUGREPORT,
                    cr2re_get_license(),
                    cr2res_etalon_create,
                    cr2res_etalon_exec,
                    cr2res_etalon_destroy)) {    
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
static int cr2res_etalon_create(cpl_plugin * plugin)
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


 
    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Execute the plugin instance given by the interface
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_etalon_exec(cpl_plugin * plugin)
{
    cpl_recipe  *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1;

    return cr2res_etalon(recipe->frames, recipe->parameters);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_etalon_destroy(cpl_plugin * plugin)
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
static int cr2res_etalon(
        cpl_frameset            *   frameset,
        const cpl_parameterlist *   parlist)
{
    cpl_propertylist    *   plist;
    cpl_frame           *   rawframe ;
    cpl_propertylist    *   applist;
    cpl_image           *   bin_image;
    cpl_image           *   etalon_im ;

    /* RETRIEVE INPUT PARAMETERS */
  
    
    /* Identify the RAW and CALIB frames in the input frameset */
    cr2re_dfs_set_groups(frameset) ;
 
    /* Get Data */
    rawframe = cpl_frameset_get_position(frameset, 0);
    plist = cpl_propertylist_load(cpl_frame_get_filename(rawframe), 0);
    etalon_im = cpl_image_load(cpl_frame_get_filename(rawframe), 
            CPL_TYPE_DOUBLE,0,0);
    if (etalon_im == NULL) {
        cpl_propertylist_delete(plist) ;
        return -1 ;
    }

    /* NOW PERFORMING THE DATA REDUCTION */
    bin_image = cr2res_etalon_computation(etalon_im);
    if (bin_image == NULL) {
        cpl_propertylist_delete(plist) ;
        cpl_image_delete(etalon_im) ;
        return -1 ;
    }

    /* Add the product category  */
    applist = cpl_propertylist_duplicate(plist);
    cpl_propertylist_append_string(applist, CPL_DFS_PRO_CATG, "ETALON_PROCATG");

    /* Save Product */
    cpl_dfs_save_image(frameset, plist, parlist, frameset, NULL, bin_image,
                CPL_BPP_IEEE_FLOAT, "cr2res_etalon", applist, NULL, 
                PACKAGE "/" PACKAGE_VERSION, "cr2res_etalon.fits") ;
    
    /* Free and return */
    cpl_propertylist_delete(plist) ;
    cpl_image_delete(etalon_im) ;
    cpl_image_delete(bin_image) ;
    cpl_propertylist_delete(applist) ;
    return (int)cpl_error_get_code();
}

