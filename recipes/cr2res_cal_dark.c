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
#include "cr2re_pfits.h"
#include "cr2re_dfs.h"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_cal_dark_create(cpl_plugin *);
static int cr2res_cal_dark_exec(cpl_plugin *);
static int cr2res_cal_dark_destroy(cpl_plugin *);
static int cr2res_cal_dark(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_cal_dark_description[] =
"TODO : Descripe here the recipe in / out / params / basic algo\n"
"raw-file.fits " CR2RES_DARK_RAW "\n"
"detlin.fits " CR2RES_DETLIN_BPM_PROCATG "\n"
" The recipe produces the following products:\n"
"master_dark.fits " CR2RES_MASTER_DARK_PROCATG "\n"
"master_bpm.fits " CR2RES_MASTER_BPM_PROCATG "\n"
"dark_bpm.fits " CR2RES_DARK_BPM_PROCATG "\n"
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
                    "cr2res_cal_dark",
                    "Dark recipe",
                    cr2res_cal_dark_description,
                    "Thomas Marquart, Yves Jung",
                    PACKAGE_BUGREPORT,
                    cr2re_get_license(),
                    cr2res_cal_dark_create,
                    cr2res_cal_dark_exec,
                    cr2res_cal_dark_destroy)) {    
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
static int cr2res_cal_dark_create(cpl_plugin * plugin)
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
static int cr2res_cal_dark_exec(cpl_plugin * plugin)
{
    cpl_recipe  *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1;

    return cr2res_cal_dark(recipe->frames, recipe->parameters);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_cal_dark_destroy(cpl_plugin * plugin)
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
static int cr2res_cal_dark(
        cpl_frameset            *   frameset,
        const cpl_parameterlist *   parlist)
{
    cpl_propertylist    *   plist;
    cpl_frame           *   rawframe ;
    cpl_propertylist    *   applist;
    cpl_image           *   in ;
    cpl_image           *   master_dark ;

    /* RETRIEVE INPUT PARAMETERS */
  
    
    /* Identify the RAW and CALIB frames in the input frameset */
    cr2re_dfs_set_groups(frameset) ;
 
    /* Get Data */
    rawframe = cpl_frameset_get_position(frameset, 0);
    plist = cpl_propertylist_load(cpl_frame_get_filename(rawframe), 0);
    in = cpl_image_load(cpl_frame_get_filename(rawframe), CPL_TYPE_DOUBLE,0,0);
    if (in == NULL) {
        cpl_propertylist_delete(plist) ;
        return -1 ;
    }

    /* NOW PERFORMING THE DATA REDUCTION */
    master_dark = cpl_image_duplicate(in);
    if (master_dark == NULL) {
        cpl_propertylist_delete(plist) ;
        cpl_image_delete(in) ;
        return -1 ;
    }
    cpl_image_delete(in) ;

    /* Add the product category  */
    applist = cpl_propertylist_duplicate(plist);
    cpl_propertylist_append_string(applist, CPL_DFS_PRO_CATG, 
            CR2RES_MASTER_DARK_PROCATG);

    /* Save Product */
    cpl_dfs_save_image(frameset, plist, parlist, frameset, NULL,
            master_dark, CPL_BPP_IEEE_FLOAT, "cr2res_cal_dark", applist,
            NULL, PACKAGE "/" PACKAGE_VERSION, "cr2res_cal_dark.fits") ;
    
    /* Free and return */
    cpl_propertylist_delete(plist) ;
    cpl_propertylist_delete(applist) ;
    cpl_image_delete(master_dark) ;
    return (int)cpl_error_get_code();
}

