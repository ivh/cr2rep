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

#define RECIPE_STRING "cr2res_cal_dark"

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
"cr2res_cal_dark_master.fits " CR2RES_MASTER_DARK_PROCATG "\n"
"cr2res_cal_dark_bpm.fits " CR2RES_MASTER_BPM_PROCATG "\n"
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
                    "cr2res_cal_dark",
                    "Dark recipe",
                    cr2res_cal_dark_description,
                    "Thomas Marquart, Yves Jung",
                    PACKAGE_BUGREPORT,
                    cr2res_get_license(),
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
	/* --gain */
    p = cpl_parameter_new_value("cr2res_cal_dark.gain", CPL_TYPE_DOUBLE,
       "Gain in [e- / ADU]", "cr2res_cal_dark", 2.5);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "gain");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /* Collapsing related parameters */
    sigclip_def = hdrl_collapse_sigclip_parameter_create(3., 3., 5);
    minmax_def = hdrl_collapse_minmax_parameter_create(1., 1.);
    collapse_par = hdrl_collapse_parameter_create_parlist("cr2res_cal_dark", 
            "", "MEDIAN", sigclip_def, minmax_def) ;
    hdrl_parameter_delete(sigclip_def);
    hdrl_parameter_delete(minmax_def);
    for (p = cpl_parameterlist_get_first(collapse_par) ;
            p != NULL; p = cpl_parameterlist_get_next(collapse_par))
        cpl_parameterlist_append(recipe->parameters,cpl_parameter_duplicate(p));
    cpl_parameterlist_delete(collapse_par);

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
    const cpl_parameter *   par ;
    double                  gain ;
    hdrl_parameter      *   collapse_params ;
    cpl_frameset        *   rawframes ;
    hdrl_imagelist      *   dark_cube ;
    const char          *   fname ;
    cpl_image           *   ima_data ;
    cpl_image           *   ima_err ;
    hdrl_image 			*   master_darks[CR2RES_NB_DETECTORS] ;
    hdrl_image 			* 	hdrl_ima ;
    hdrl_image 			* 	master;
    cpl_image 			*	contrib_map;
    int                     nb_frames, i, ext ;

    /* RETRIEVE INPUT PARAMETERS */
    /* --gain */
    par = cpl_parameterlist_find_const(parlist, "cr2res_cal_dark.gain");
    gain = cpl_parameter_get_double(par);
    /* Collapse parameters */
    collapse_params = hdrl_collapse_parameter_parse_parlist(parlist,
            "cr2res_cal_dark") ;
   
    /* Identify the RAW and CALIB frames in the input frameset */
    if (cr2res_dfs_set_groups(frameset)) {
        hdrl_parameter_destroy(collapse_params) ;
        cpl_msg_error(__func__, "Cannot identify RAW and CALIB frames") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }
	
    /* Extract RAW frames */
    rawframes = cr2res_extract_frameset(frameset, CR2RES_DARK_RAW) ;
    if (rawframes==NULL || cpl_frameset_get_size(rawframes) <= 0) {
        hdrl_parameter_destroy(collapse_params) ;
        cpl_msg_error(__func__, "Cannot find any RAW file") ;
        cpl_error_set(__func__, CPL_ERROR_DATA_NOT_FOUND) ;
        return -1 ;
    }
    nb_frames = cpl_frameset_get_size(rawframes) ;

	/* Loop on the extensions */
    for (ext=1 ; ext<=CR2RES_NB_DETECTORS ; ext++) {
        cpl_msg_info(__func__, "Process Detector nb %i", ext) ;
        cpl_msg_indent_more() ;

        /* Loop on the frames */
        dark_cube = hdrl_imagelist_new();
        for (i=0; i<nb_frames ; i++) {
            /* Identify current file */
            fname=cpl_frame_get_filename(
                    cpl_frameset_get_position(rawframes, i)) ; 
            cpl_msg_info(__func__, "Load Image from File %s / Detector %i", 
                    fname, ext) ;
            
            /* Load the image */
            if ((ima_data=cpl_image_load(fname,CPL_TYPE_DOUBLE,0,ext))==NULL) {
                hdrl_parameter_destroy(collapse_params) ;
                cpl_frameset_delete(rawframes) ;
                hdrl_imagelist_delete(dark_cube) ;
                cpl_msg_error(__func__, 
                        "Cannot load image from File %s / Detector %d", 
                        fname, ext) ;
                cpl_error_set(__func__, CPL_ERROR_DATA_NOT_FOUND) ;
                cpl_msg_indent_less() ;
                return -1 ;
            }

            /* Create the noise image */
            cpl_msg_info(__func__, "Create the associated Noise image");
            if (cr2res_detector_shotnoise_model(ima_data, gain, 10.,
                        &ima_err) != CPL_ERROR_NONE) {
                hdrl_parameter_destroy(collapse_params) ;
                cpl_frameset_delete(rawframes) ;
                hdrl_imagelist_delete(dark_cube) ;
                cpl_image_delete(ima_data); 
                cpl_msg_error(__func__, "Cannot create the Noise image") ;
                cpl_error_set(__func__, CPL_ERROR_DATA_NOT_FOUND) ;
                cpl_msg_indent_less() ;
                return -1 ;
            }

            /* Store Data and Error together in an hdrl image */
            hdrl_ima = hdrl_image_create(ima_data, ima_err);
            cpl_image_delete(ima_data);
            cpl_image_delete(ima_err);
            
            /* Store the hdrl image in the dark_cube */
            hdrl_imagelist_set(dark_cube, hdrl_ima, i);
        }

        /* Get the proper collapsing function and perform frames combination */
        if (hdrl_imagelist_collapse(dark_cube, collapse_params,
                &(master_darks[ext-1]), &contrib_map) != CPL_ERROR_NONE) {
            cpl_msg_warning(__func__, "Cannot collapse Detector %d", ext) ;
            master_darks[ext-1] = NULL ;
            contrib_map = NULL ;
        }
        cpl_image_delete(contrib_map);
        hdrl_imagelist_delete(dark_cube);
    
        cpl_msg_indent_less() ;
    }
    hdrl_parameter_delete(collapse_params);

	/* Save the results */
	if (cr2res_io_save_MASTER_DARK(frameset, "cr2res_cal_dark_master.fits", 
                rawframes, parlist, master_darks, NULL, RECIPE_STRING) != 0) {
        cpl_frameset_delete(rawframes) ;
        for (ext=1 ; ext<=CR2RES_NB_DETECTORS ; ext++) {
            if (master_darks[ext-1] != NULL) 
                hdrl_image_delete(master_darks[ext-1]);
        }
        cpl_msg_error(__func__, "Cannot save the MASTER DARK") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        cpl_msg_indent_less() ;
        return -1 ;

    }
    cpl_frameset_delete(rawframes) ;
    for (ext=1 ; ext<=CR2RES_NB_DETECTORS ; ext++) {
        if (master_darks[ext-1] != NULL) hdrl_image_delete(master_darks[ext-1]);
    }

    return (int)cpl_error_get_code();
}
