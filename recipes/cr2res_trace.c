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

static int cr2res_trace_create(cpl_plugin *);
static int cr2res_trace_exec(cpl_plugin *);
static int cr2res_trace_destroy(cpl_plugin *);
static int cr2res_trace(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_trace_description[] =
"This example text is used to describe the recipe.\n"
"The description should include the required FITS-files and\n"
"their associated tags, e.g.\n"
"raw-file.fits " CR2RE_TRACE_RAW "\n"
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
                    "cr2res_trace",
                    "Tracing programm",
                    cr2res_trace_description,
                    "Thomas Marquart",
                    PACKAGE_BUGREPORT,
                    cr2re_get_license(),
                    cr2res_trace_create,
                    cr2res_trace_exec,
                    cr2res_trace_destroy)) {    
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
static int cr2res_trace_create(cpl_plugin * plugin)
{
    cpl_recipe    * recipe;                                               
    cpl_parameter * p;
                                                                       
    /* Do not create the recipe if an error code is already set */     
    if (cpl_error_get_code() != CPL_ERROR_NONE) {                      
        cpl_msg_error(cpl_func, "%s():%d: An error is already set: %s",
                      cpl_func, __LINE__, cpl_error_get_where());      
        return (int)cpl_error_get_code();                              
    }                                                                  
                                                                       
    if (plugin == NULL) {                                              
        cpl_msg_error(cpl_func, "Null plugin");                        
        cpl_ensure_code(0, (int)CPL_ERROR_NULL_INPUT);                 
    }                                                                  
                                                                       
    /* Verify plugin type */                                           
    if (cpl_plugin_get_type(plugin) != CPL_PLUGIN_TYPE_RECIPE) {       
        cpl_msg_error(cpl_func, "Plugin is not a recipe");             
        cpl_ensure_code(0, (int)CPL_ERROR_TYPE_MISMATCH);              
    }                                                                  
                                                                       
    /* Get the recipe */                                               
    recipe = (cpl_recipe *)plugin;                                     
                                                                       
    /* Create the parameters list in the cpl_recipe object */          
    recipe->parameters = cpl_parameterlist_new();                      
    if (recipe->parameters == NULL) {                                  
        cpl_msg_error(cpl_func, "Parameter list allocation failed");   
        cpl_ensure_code(0, (int)CPL_ERROR_ILLEGAL_OUTPUT);             
    }                                                                  

    /* Fill the parameters list */
    /* --stropt */
    p = cpl_parameter_new_value("cr2res.cr2res_trace.str_option", 
            CPL_TYPE_STRING, "the string option", "cr2res.cr2res_trace",NULL);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "stropt");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /* --boolopt */
    p = cpl_parameter_new_value("cr2res.cr2res_trace.bool_option", 
            CPL_TYPE_BOOL, "a flag", "cr2res.cr2res_trace", TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "boolopt");
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
static int cr2res_trace_exec(cpl_plugin * plugin)
{

    cpl_recipe * recipe;                                                   
    int recipe_status;                                                     
    cpl_errorstate initial_errorstate = cpl_errorstate_get();              
                                                                           
    /* Return immediately if an error code is already set */               
    if (cpl_error_get_code() != CPL_ERROR_NONE) {                          
        cpl_msg_error(cpl_func, "%s():%d: An error is already set: %s",    
                      cpl_func, __LINE__, cpl_error_get_where());          
        return (int)cpl_error_get_code();                                  
    }                                                                      
                                                                           
    if (plugin == NULL) {                                                  
        cpl_msg_error(cpl_func, "Null plugin");                            
        cpl_ensure_code(0, (int)CPL_ERROR_NULL_INPUT);                     
    }                                                                      
                                                                           
    /* Verify plugin type */                                               
    if (cpl_plugin_get_type(plugin) != CPL_PLUGIN_TYPE_RECIPE) {           
        cpl_msg_error(cpl_func, "Plugin is not a recipe");                 
        cpl_ensure_code(0, (int)CPL_ERROR_TYPE_MISMATCH);                  
    }                                                                      
                                                                           
    /* Get the recipe */                                                   
    recipe = (cpl_recipe *)plugin;                                         
                                                                           
    /* Verify parameter and frame lists */                                 
    if (recipe->parameters == NULL) {                                      
        cpl_msg_error(cpl_func, "Recipe invoked with NULL parameter list");
        cpl_ensure_code(0, (int)CPL_ERROR_NULL_INPUT);                     
    }                                                                      
    if (recipe->frames == NULL) {                                          
        cpl_msg_error(cpl_func, "Recipe invoked with NULL frame set");     
        cpl_ensure_code(0, (int)CPL_ERROR_NULL_INPUT);                     
    }                                                                      
                                                                           
    /* Invoke the recipe */                                                
    recipe_status = cr2res_trace(recipe->frames, recipe->parameters);
                                                                           
    /* Ensure DFS-compliance of the products */                            
    if (cpl_dfs_update_product_header(recipe->frames)) {                   
        if (!recipe_status) recipe_status = (int)cpl_error_get_code();                         
    }                                                                      
                                                                           
    if (!cpl_errorstate_is_equal(initial_errorstate)) {                    
        /* Dump the error history since recipe execution start.            
           At this point the recipe cannot recover from the error */       
        cpl_errorstate_dump(initial_errorstate, CPL_FALSE, NULL);          
    }                                                                      
                                                                           
    return recipe_status;                                                  
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_trace_destroy(cpl_plugin * plugin)
{
    cpl_recipe * recipe;                                          
                                                                  
    if (plugin == NULL) {                                         
        cpl_msg_error(cpl_func, "Null plugin");                   
        cpl_ensure_code(0, (int)CPL_ERROR_NULL_INPUT);            
    }                                                             
                                                                  
    /* Verify plugin type */                                      
    if (cpl_plugin_get_type(plugin) != CPL_PLUGIN_TYPE_RECIPE) {  
        cpl_msg_error(cpl_func, "Plugin is not a recipe");        
        cpl_ensure_code(0, (int)CPL_ERROR_TYPE_MISMATCH);         
    }                                                             
                                                                  
    /* Get the recipe */                                          
    recipe = (cpl_recipe *)plugin;                                
                                                                  
    cpl_parameterlist_delete(recipe->parameters);             
                                                                  
    return 0;                                                    
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Interpret the command line options and execute the data processing
  @param    frameset   the frames list
  @param    parlist    the parameters list
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_trace(
        cpl_frameset            * frameset,
        const cpl_parameterlist * parlist)
{
    const cpl_parameter *   param;
    const char          *   str_option;
    int                     bool_option;
    const cpl_frame     *   rawframe;
    double                  qc_param;
    cpl_propertylist    *   plist;
    cpl_propertylist    *   applist;
    cpl_image           *   image;
    cpl_imagelist       *   biaslist;
    /* Use the errorstate to detect an error in a function that does not
       return an error code. */
    cpl_errorstate          prestate = cpl_errorstate_get();

    /* HOW TO RETRIEVE INPUT PARAMETERS */
    /* --stropt */
    param = cpl_parameterlist_find_const(parlist,
                                         "cr2res.cr2res_trace.str_option");
    str_option = cpl_parameter_get_string(param);

    /* --boolopt */
    param = cpl_parameterlist_find_const(parlist,
                                         "cr2res.cr2res_trace.bool_option");
    bool_option = cpl_parameter_get_bool(param);
  
    if (!cpl_errorstate_is_equal(prestate)) {
        return (int)cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                          "Could not retrieve the input "
                                          "parameters");
    }
    
    /* Identify the RAW and CALIB frames in the input frameset */
    cpl_ensure_code(cr2re_dfs_set_groups(frameset) == CPL_ERROR_NONE,
                    cpl_error_get_code());
 
    /* HOW TO ACCESS INPUT DATA */
    /*  - A required file */
    biaslist = cpl_imagelist_load_frameset(frameset, CPL_TYPE_DOUBLE,0,0);
    if (biaslist== NULL) {
        /* cpl_frameset_find_const() does not set an error code, when a frame
           is not found, so we will set one here. */
        return (int)cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                                          "SOF does not have any file tagged "
                                          "with %s", CR2RE_TRACE_RAW);
    }
    
    /* HOW TO GET THE VALUE OF A FITS KEYWORD */
    /*  - Load only DETector related keys */
    rawframe=cpl_frameset_get_first(frameset);
    plist = cpl_propertylist_load(cpl_frame_get_filename(rawframe),
                                          0);
    if (plist == NULL) {
        /* In this case an error message is added to the error propagation */
        return (int)cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                          "Could not read the FITS header");
    }

    qc_param = cr2re_pfits_get_dit(plist);

    /* Check for a change in the CPL error state */
    /* - if it did change then propagate the error and return */
    cpl_ensure_code(cpl_errorstate_is_equal(prestate), cpl_error_get_code());
    
    /* NOW PERFORMING THE DATA REDUCTION */
    image = cpl_imagelist_collapse_create(biaslist);
    if (image == NULL) {
        return (int)cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                     "Average failed");
    }

    applist = cpl_propertylist_duplicate(plist);

    /* Add the product category  */
    cpl_propertylist_append_string(applist, CPL_DFS_PRO_CATG,
                                   CR2RE_TRACE_PROCATG);

    /* Add a QC parameter  */
    cpl_propertylist_append_double(applist, "ESO QC QCPARAM", qc_param);
    
    /* HOW TO SAVE A DFS-COMPLIANT PRODUCT TO DISK  */
    if (cpl_dfs_save_image(frameset, plist, parlist, frameset, NULL, image,
                           CPL_BPP_IEEE_FLOAT, "cr2res_trace", applist,
                           NULL, PACKAGE "/" PACKAGE_VERSION,
                           "cr2res_trace.fits")) {
        /* Propagate the error */
        (void)cpl_error_set_where(cpl_func);
    }
    

    cpl_propertylist_delete(plist);
    cpl_imagelist_delete(biaslist);
    cpl_propertylist_delete(applist);
    cpl_image_delete(image);

    return (int)cpl_error_get_code();
}
