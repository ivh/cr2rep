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

#include "cr2res_utils.h"
#include "cr2res_bpm.h"
#include "cr2res_dfs.h"
#include "cr2res_io.h"

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

#define RECIPE_STRING "cr2res_util_bpm_split"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_util_bpm_split_create(cpl_plugin *);
static int cr2res_util_bpm_split_exec(cpl_plugin *);
static int cr2res_util_bpm_split_destroy(cpl_plugin *);
static int cr2res_util_bpm_split(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_util_bpm_split_description[] = "\
BPM splitting                                                           \n\
  Each input BPM is splitted into several BPMs                          \n\
                                                                        \n\
  Inputs                                                                \n\
   	raw.fits " CR2RES_CAL_DARK_BPM_PROCATG " [1 to n]                   \n\
          or " CR2RES_CAL_FLAT_BPM_PROCATG "                            \n\
          or " CR2RES_CAL_DETLIN_BPM_PROCATG "                          \n\
          or " CR2RES_UTIL_BPM_SPLIT_PROCATG "                          \n\
          or " CR2RES_UTIL_NORM_BPM_PROCATG "                           \n\
                                                                        \n\
  Outputs                                                               \n\
    <input_name>_splitted_<bpm_code>.fits " 
    CR2RES_UTIL_BPM_SPLIT_PROCATG "\n\
                                                                        \n\
  Algorithm                                                             \n\
    loop on input raw frames f:                                         \n\
      loop on detectors d:                                              \n\
        loop on bpm types t:                                            \n\
          call cr2res_bpm_from_mask()                                   \n\
           -> bpm_single_type(t, d, f)                                  \n\
      loop on bpm types t:                                              \n\
        Save bpm_single_type(f, t) (UTIL_BPM_SPLIT)                     \n\
                                                                        \n\
  Library functions uѕed:                                               \n\
    cr2res_io_load_BPM()                                                \n\
    cr2res_bpm_from_mask()                                              \n\
    cr2res_io_save_BPM()                                                \n\
";

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
                    "BPM splitting utility",
                    cr2res_util_bpm_split_description,
                    "Thomas Marquart, Yves Jung",
                    PACKAGE_BUGREPORT,
                    cr2res_get_license(),
                    cr2res_util_bpm_split_create,
                    cr2res_util_bpm_split_exec,
                    cr2res_util_bpm_split_destroy)) {    
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
static int cr2res_util_bpm_split_create(cpl_plugin * plugin)
{
    cpl_recipe          *   recipe ;
    cpl_parameter       *   p ;

    /* Check that the plugin is part of a valid recipe */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else
        return -1;

    /* Create the parameters list in the cpl_recipe object */
    recipe->parameters = cpl_parameterlist_new();

    /* Fill the parameters list */
    p = cpl_parameter_new_value("cr2res.cr2res_util_bpm_split.detector",
            CPL_TYPE_INT, "Only reduce the specified detector",
            "cr2res.cr2res_util_bpm_split", 0);
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
static int cr2res_util_bpm_split_exec(cpl_plugin * plugin)
{
    cpl_recipe  *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1;

    return cr2res_util_bpm_split(recipe->frames, recipe->parameters);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_util_bpm_split_destroy(cpl_plugin * plugin)
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
static int cr2res_util_bpm_split(
        cpl_frameset            *   frameset,
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter *   param ;
    int                     reduce_det ;
    cpl_frameset        *   rawframes ;
    const cpl_frame     *   cur_frame ;
    const char          *   cur_fname ;
    cpl_frameset        *   cur_fset ;
    cpl_image           *   ima ;
    cpl_image   *   splitted_bpms[CR2RES_NB_BPM_TYPES][CR2RES_NB_DETECTORS] ;
    cpl_mask            *   my_mask ;
    cpl_propertylist    *   ext_plist[CR2RES_NB_DETECTORS] ;
    char                *   out_file;
    int                     i, j, det_nr, wished_ext_nb; 

    /* Initialise */

    /* RETRIEVE INPUT PARAMETERS */
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_bpm_split.detector");
    reduce_det = cpl_parameter_get_int(param);

    /* Identify the RAW and CALIB frames in the input frameset */
    if (cr2res_dfs_set_groups(frameset)) {
        cpl_msg_error(__func__, "Cannot identify RAW and CALIB frames") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Get Calibration frames */

    /* Get the rawframes */
    rawframes = cr2res_io_find_BPM_all(frameset) ;
    if (rawframes==NULL || cpl_frameset_get_size(rawframes) <= 0) {
        cpl_msg_error(__func__, "Cannot find any RAW file") ;
        cpl_error_set(__func__, CPL_ERROR_DATA_NOT_FOUND) ;
        return -1 ;
    }

    /* Loop on the RAW frames */
    for (i=0 ; i<cpl_frameset_get_size(rawframes) ; i++) {
        /* Get the Current Frame */
        cur_frame = cpl_frameset_get_position(rawframes, i) ;
        cur_fname = cpl_frame_get_filename(cur_frame) ;
        cpl_msg_info(__func__, "Reduce Frame %s", cur_fname) ;
        cpl_msg_indent_more() ;

        /* Loop on the detectors */
        for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
            /* Initialise */
            for (j=0 ; j<CR2RES_NB_BPM_TYPES ; j++)
                splitted_bpms[j][det_nr-1] = NULL ;
            ext_plist[det_nr-1] = NULL ;

            /* Compute only one detector */
            if (reduce_det != 0 && det_nr != reduce_det) continue ;
        
            cpl_msg_info(__func__, "Process Detector %d", det_nr) ;
            cpl_msg_indent_more() ;

            /* Load the image to calibrate */
            ima = cr2res_io_load_BPM(cur_fname, det_nr, 1);

            /* Loop on the BPM types */
            for (j=0 ; j<CR2RES_NB_BPM_TYPES ; j++) {
                my_mask = cr2res_bpm_extract_mask(ima, bpm_types[j]) ;
                splitted_bpms[j][det_nr-1] = 
                    cr2res_bpm_from_mask(my_mask, bpm_types[j]) ;
                cpl_mask_delete(my_mask) ;
            }
            cpl_image_delete(ima); 
            cpl_msg_indent_less() ;

            /* Create the header */
            wished_ext_nb = cr2res_io_get_ext_idx(cur_fname, det_nr, 1) ;
            ext_plist[det_nr-1]=cpl_propertylist_load(cur_fname,wished_ext_nb);
        }

        /* Ѕave Products */

        /* SPLITTED_BPM */
        for (j=0 ; j<CR2RES_NB_BPM_TYPES ; j++) {
            out_file=cpl_sprintf("%s_splitted_%d.fits", 
                    cr2res_get_root_name(cur_fname), bpm_types[j]) ;
            cur_fset = cpl_frameset_new() ;
            cpl_frameset_insert(cur_fset, cpl_frame_duplicate(cur_frame)) ;
            cr2res_io_save_BPM(out_file, frameset, cur_fset, parlist,
                    splitted_bpms[j], NULL, ext_plist, 
                    CR2RES_UTIL_BPM_SPLIT_PROCATG, RECIPE_STRING) ;
            cpl_frameset_delete(cur_fset) ;
            cpl_free(out_file);
        }
        /* Free */
        for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
            for (j=0 ; j<CR2RES_NB_BPM_TYPES ; j++) {
                if (splitted_bpms[j][det_nr-1] != NULL)
                    cpl_image_delete(splitted_bpms[j][det_nr-1]) ;
            }
            if (ext_plist[det_nr-1] != NULL) 
                cpl_propertylist_delete(ext_plist[det_nr-1]) ;
        }
        cpl_msg_indent_less() ;
    }
    cpl_frameset_delete(rawframes) ;
    return (int)cpl_error_get_code();
}

