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
    p = cpl_parameter_new_value("cr2res.cr2res_trace.poly_order",
            CPL_TYPE_INT,
            "polynomial order for the fit to the orders",
            "cr2res.cr2res_trace",4);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "polyorder");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_trace.min_cluster",
            CPL_TYPE_INT,
            "size (number of pixels) of the smallest allowed cluster",
            "cr2res.cr2res_trace",40);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "mincluster");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_trace.smooth",
            CPL_TYPE_DOUBLE,
            "Length of the smoothing kernel, relative to inter-order separation",
            "cr2res.cr2res_trace", 1.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "smooth");
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
    int                     polyorder;
    int                     mincluster;
    double                   smoothfactor;
    const cpl_frame     *   rawframe;
    double                  qc_param;
    cpl_propertylist    *   plist;
    cpl_propertylist    *   applist;
    cpl_image           *   image;
    cpl_image           *   smimage;
    cpl_imagelist       *   imlist;
    cpl_mask            *   mask;
    cpl_size                npix;
    cpl_matrix          *   kernel;
    int                 *   xs;
    int                 *   ys;
    int                 *   clusters;
    int i,j,nx,ny,nclusters;
    int count=0;
    int ordersep=180; //this needs to come from a static calibration, each band

    /* Use the errorstate to detect an error in a function that does not
       return an error code. */
    cpl_errorstate          prestate = cpl_errorstate_get();

    /* HOW TO RETRIEVE INPUT PARAMETERS */
    param = cpl_parameterlist_find_const(parlist,
                                         "cr2res.cr2res_trace.min_cluster");
    mincluster = cpl_parameter_get_int(param);

    param = cpl_parameterlist_find_const(parlist,
                                         "cr2res.cr2res_trace.poly_order");
    polyorder = cpl_parameter_get_int(param);

    param = cpl_parameterlist_find_const(parlist,
                                         "cr2res.cr2res_trace.smooth");
    smoothfactor = cpl_parameter_get_double(param);



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
    imlist = cpl_imagelist_load_frameset(frameset, CPL_TYPE_DOUBLE,0,0);
    if (imlist== NULL) {
        /* cpl_frameset_find_const() does not set an error code, when a frame
           is not found, so we will set one here. */
        return (int)cpl_error_set_message(cpl_func, CPL_ERROR_DATA_NOT_FOUND,
                                          "SOF does not have any file tagged "
                                          "with %s", CR2RE_TRACE_RAW);
    }

    /* HOW TO GET THE VALUE OF A FITS KEYWORD */
    /*  - Load only DETector related keys */
    rawframe=cpl_frameset_get_position(frameset,0);
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
    image = cpl_imagelist_get(imlist,0);
    if (image == NULL) {
        return (int)cpl_error_set_message(cpl_func, cpl_error_get_code(),
                                     "Could not get image out of imagelist");
    }

    applist = cpl_propertylist_duplicate(plist);

    /* Add the product category  */
    cpl_propertylist_append_string(applist, CPL_DFS_PRO_CATG,
                                   CR2RE_TRACE_PROCATG);

    /* Add a QC parameter  */
    cpl_propertylist_append_double(applist, "ESO QC QCPARAM", qc_param);


    /* find the pixels with signal*/
    smimage = cpl_image_duplicate(image);
    ordersep = (int) (ordersep*smoothfactor);
    if (ordersep % 2 == 0) ordersep +=1;
    cpl_msg_debug(cpl_func,cpl_sprintf("ordersep: %d",ordersep));
    kernel = cpl_matrix_new( ordersep ,1);
    cpl_matrix_add_scalar(kernel,0.01);
    if (cpl_image_filter(smimage, image, kernel,
        CPL_FILTER_LINEAR, CPL_BORDER_FILTER) != CPL_ERROR_NONE) {
         return cpl_error_set_message(cpl_func, CPL_ERROR_ILLEGAL_INPUT,
         "The filter went bad");
    }


    mask = cpl_mask_threshold_image_create(image,100,1000);
    cpl_mask_save(mask,"mask.fits",plist,CPL_IO_CREATE);

    npix = cpl_mask_count(mask);
    nx = cpl_mask_get_size_x(mask);
    ny = cpl_mask_get_size_y(mask);
    cpl_msg_debug(cpl_func,cpl_sprintf("mask: %d %d, %d",nx,ny, npix));

    xs=(int *)cpl_malloc(npix*sizeof(int));
    ys=(int *)cpl_malloc(npix*sizeof(int));
    clusters=(int *)cpl_malloc(npix*sizeof(int));


    /* make the arrays of x and y indices that only contail the ones with signal*/
    for(i=1;i<=nx;i++){
        for(j=1;j<=ny;j++){
            if (cpl_mask_get(mask,i,j) == CPL_BINARY_1) {
                xs[count]=i;
                ys[count]=j;
                count++;
            }
        }
    }

    nclusters = cluster(xs,ys,npix,nx,ny,mincluster,clusters);

    /* put the results bac into 2d image form */
    for(i=0;i<npix;i++){
            cpl_image_set(image,xs[i],ys[i],clusters[i]);
    }

    /* HOW TO SAVE A DFS-COMPLIANT PRODUCT TO DISK  */
    if (cpl_dfs_save_image(frameset, plist, parlist, frameset, NULL, smimage,
                           CPL_BPP_IEEE_FLOAT, "cr2res_trace", applist,
                           NULL, PACKAGE "/" PACKAGE_VERSION,
                           "cr2res_trace.fits")) {
        /* Propagate the error */
        (void)cpl_error_set_where(cpl_func);
    }


    cpl_propertylist_delete(plist);
    cpl_propertylist_delete(applist);
    cpl_imagelist_delete(imlist);
    cpl_image_delete(smimage);
    cpl_mask_delete(mask);
    cpl_matrix_delete(kernel);
    cpl_free(xs);
    cpl_free(ys);
    cpl_free(clusters);

    return (int)cpl_error_get_code();
}
