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
#include "cr2res_pfits.h"
#include "cr2res_dfs.h"
#include "cr2res_flat.h"
#include "cr2res_trace.h"
#include "cr2res_extract.h"
#include "cr2res_io.h"

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

#define RECIPE_STRING "cr2res_util_normflat"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_util_normflat_reduce(
        const cpl_frameset  *   rawframes,
        const cpl_frame     *   detlin_frame,
        const cpl_frame     *   master_dark_frame,
        const cpl_frame     *   master_bpm_frame,
        int                     reduce_det,
        hdrl_image          **  master_flat,
        hdrl_image          **  bpm,
        cpl_propertylist    **  plist) ;
static int cr2res_util_normflat_create(cpl_plugin *);
static int cr2res_util_normflat_exec(cpl_plugin *);
static int cr2res_util_normflat_destroy(cpl_plugin *);
static int cr2res_util_normflat(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_util_normflat_description[] =
"TODO : Descripe here the recipe in / out / params / basic algo\n"
"raw-file.fits " CR2RES_FLAT_RAW "\n"
"detlin.fits " CR2RES_DETLIN_COEFFS_PROCATG "\n"
"master_dark.fits " CR2RES_MASTER_DARK_PROCATG "\n"
"master_bpm.fits " CR2RES_MASTER_BPM_PROCATG "\n"
" The recipe produces the following products:\n"
"cr2res_util_normflat_master.fits " CR2RES_MASTER_FLAT_PROCATG  "\n"
"cr2res_util_normflat_bpm.fits " CR2RES_MASTER_BPM_PROCATG "\n"
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
                    "cr2res_util_normflat",
                    "Flat utility",
                    cr2res_util_normflat_description,
                    "Thomas Marquart, Yves Jung",
                    PACKAGE_BUGREPORT,
                    cr2res_get_license(),
                    cr2res_util_normflat_create,
                    cr2res_util_normflat_exec,
                    cr2res_util_normflat_destroy)) {    
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
static int cr2res_util_normflat_create(cpl_plugin * plugin)
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
    p = cpl_parameter_new_value("cr2res.cr2res_util_normflat.detector",
            CPL_TYPE_INT, "Only reduce the specified detector",
            "cr2res.cr2res_util_normflat", 0);
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
static int cr2res_util_normflat_exec(cpl_plugin * plugin)
{
    cpl_recipe  *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1;

    return cr2res_util_normflat(recipe->frames, recipe->parameters);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_util_normflat_destroy(cpl_plugin * plugin)
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
static int cr2res_util_normflat(
        cpl_frameset            *   frameset,
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter *   param ;
    int                     reduce_det ;
    cpl_frameset        *   rawframes ;
    const cpl_frame     *   detlin_frame ;
    const cpl_frame     *   master_dark_frame ;
    const cpl_frame     *   master_bpm_frame ;
    hdrl_image          *   master_flat[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   bpm[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   ext_plist[CR2RES_NB_DETECTORS] ;
    char                *   out_file;
    int                     i, det_nr; 

    /* Initiaise */
    cr2res_decker decker_values[CR2RES_NB_DECKER_POSITIONS] = 
        {CR2RES_DECKER_NONE, CR2RES_DECKER_1_3, CR2RES_DECKER_2_4} ; 
    char * decker_desc[CR2RES_NB_DECKER_POSITIONS] =
        {"Open", "Decker1", "Decker2"} ;

    /* RETRIEVE INPUT PARAMETERS */
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_normflat.detector");
    reduce_det = cpl_parameter_get_int(param);

    /* Identify the RAW and CALIB frames in the input frameset */
    if (cr2res_dfs_set_groups(frameset)) {
        cpl_msg_error(__func__, "Cannot identify RAW and CALIB frames") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }
	
    /* Get Calibration frames */
    detlin_frame = cpl_frameset_find_const(frameset,
            CR2RES_DETLIN_COEFFS_PROCATG);
    master_dark_frame = cpl_frameset_find_const(frameset,
            CR2RES_MASTER_DARK_PROCATG) ; 
    master_bpm_frame = cpl_frameset_find_const(frameset,
            CR2RES_MASTER_BPM_PROCATG) ;

    /* Loop on the decker positions */
    for (i=0 ; i<CR2RES_NB_DECKER_POSITIONS ; i++) {
        /* Get the Frames for the current decker position */
        rawframes = cr2res_extract_decker_frameset(frameset,
                CR2RES_FLAT_RAW, decker_values[i]) ;
        if (rawframes == NULL) continue ;
        cpl_msg_info(__func__, "Reduce %s Frames", decker_desc[i]) ;
        cpl_msg_indent_more() ;

        /* Loop on the detectors */
        for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
            /* Initialise */
            master_flat[det_nr-1] = NULL ;
            bpm[det_nr-1] = NULL ;
            ext_plist[det_nr-1] = NULL ;

            /* Compute only one detector */
            if (reduce_det != 0 && det_nr != reduce_det) continue ;
        
            cpl_msg_info(__func__, "Process Detector %d", det_nr) ;
            cpl_msg_indent_more() ;

            /* Call the reduction function */
            if (cr2res_util_normflat_reduce(rawframes, detlin_frame, 
                        master_dark_frame, master_bpm_frame, det_nr,
                        &(master_flat[det_nr-1]),
                        &(bpm[det_nr-1]),
                        &(ext_plist[det_nr-1])) == -1) {
                cpl_msg_warning(__func__, 
                        "Failed to reduce detector %d of %s Frames", 
                        det_nr, decker_desc[i]);
            }
            cpl_msg_indent_less() ;
        }
        cpl_frameset_delete(rawframes) ;

        /* Ѕave Products */

        /* MASTER_FLAT */
		out_file = cpl_sprintf("%s_%s_master_flat.fits", RECIPE_STRING,
                decker_desc[i]) ;
        cr2res_io_save_MASTER_FLAT(out_file, frameset, parlist,
                master_flat, NULL, ext_plist, RECIPE_STRING) ;
		cpl_free(out_file);

        /* BPM */
		out_file = cpl_sprintf("%s_%s_master_bpm.fits", RECIPE_STRING,
                decker_desc[i]) ;
        cr2res_io_save_MASTER_BPM(out_file, frameset, parlist,
                bpm, NULL, ext_plist, RECIPE_STRING) ;
		cpl_free(out_file);

        /* Free */
        for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
            if (master_flat[det_nr-1] != NULL)
                hdrl_image_delete(master_flat[det_nr-1]) ;
            if (bpm[det_nr-1] != NULL) 
                hdrl_image_delete(bpm[det_nr-1]) ;
            if (ext_plist[det_nr-1] != NULL) 
                cpl_propertylist_delete(ext_plist[det_nr-1]) ;
        }
        cpl_msg_indent_less() ;
    }
    return (int)cpl_error_get_code();
}

/*----------------------------------------------------------------------------*/
/**
  @brief  
  @param 
  @return  
 */
/*----------------------------------------------------------------------------*/
static int cr2res_util_normflat_reduce(
        const cpl_frameset  *   rawframes,
        const cpl_frame     *   detlin_frame,
        const cpl_frame     *   master_dark_frame,
        const cpl_frame     *   master_bpm_frame,
        int                     reduce_det,
        hdrl_image          **  master_flat,
        hdrl_image          **  bpm,
        cpl_propertylist    **  ext_plist)
{
    const char          *   first_file ;
    cpl_imagelist       *   detlin_coeffs ;
    cpl_image           *   master_dark ;
    cpl_imagelist       *   imlist ;
    hdrl_image          *   master_flat_loc ;
    cpl_image           *   bpm_in ;
    cpl_image           *   bpm_loc ;
    cpl_table           *   traces ;
    cpl_vector          **  spectrum ;
    cpl_vector          **  slit_func ;
    hdrl_image          *   model_master;
    cpl_table           *   slit_func_tab ;
    cpl_table           *   extract_tab ;
    hdrl_image          *   model_tmp ;
    cpl_propertylist    *   plist ;
    int                     i, ext_nr, nb_traces, order, trace_id ;
    
    /* Check Inputs */
    if (rawframes == NULL) return -1 ;

    /* Get the Extension number */
    first_file = cpl_frame_get_filename(
            cpl_frameset_get_position_const(rawframes, 0)) ;
    ext_nr = cr2res_io_get_ext_idx(first_file, reduce_det, 1) ;

    /* Load the extension header for saving */
    plist = cpl_propertylist_load(first_file, ext_nr) ;
    if (plist == NULL) return -1 ;

    /* Load the image list */
    imlist = cpl_imagelist_load_frameset(rawframes, CPL_TYPE_FLOAT, 1, ext_nr) ;
    if (imlist == NULL) {
        cpl_msg_error(__func__, "Failed to Load the images") ;
        cpl_propertylist_delete(plist);
        return -1 ;
    }

    /* Load MASTER DARK */
    if (master_dark_frame != NULL) {
        if ((master_dark = cr2res_io_load_MASTER_DARK(
                        cpl_frame_get_filename(master_dark_frame), 
                        reduce_det, 1)) == NULL) {
            cpl_msg_warning(__func__, "Failed to Load the Master Dark") ;
        }
    } else {
        master_dark = NULL ;
    }

    /* Load DETLIN */
    if (detlin_frame != NULL) {
        if ((detlin_coeffs = cr2res_io_load_DETLIN_COEFFS(
                        cpl_frame_get_filename(detlin_frame), 
                        reduce_det)) == NULL) {
            cpl_msg_warning(__func__, "Failed to Load the Detlin Coeffs") ;
        }
    } else {
        detlin_coeffs = NULL ;
    }

    /* Compute the Master flat */
    if ((master_flat_loc = cr2res_flat(imlist, master_dark, 
                    detlin_coeffs)) == NULL) {
        cpl_msg_error(__func__, "Failed to Compute the master flat") ;
        cpl_propertylist_delete(plist);
        cpl_imagelist_delete(imlist) ;
        if (detlin_coeffs != NULL) cpl_imagelist_delete(detlin_coeffs) ;
        if (master_dark != NULL) cpl_image_delete(master_dark) ;
        return -1 ;
    }
    if (detlin_coeffs != NULL) cpl_imagelist_delete(detlin_coeffs) ;
    if (master_dark != NULL) cpl_image_delete(master_dark) ;
    cpl_imagelist_delete(imlist) ;

    /* Compute BPM */
    if ((bpm_loc = cr2res_bpm_from_master_flat(master_flat_loc,
                    0.5, 2.0, 0.5)) == NULL) {
        cpl_msg_warning(__func__, "Failed to Compute the BPM") ;
    } else {
        /* Merge BPM with the input one */
        if (master_bpm_frame != NULL) {
            if ((bpm_in = cr2res_io_load_MASTER_BPM(
                            cpl_frame_get_filename(master_bpm_frame), 
                            reduce_det)) == NULL) {
                cpl_msg_warning(__func__, "Failed to Load the Master BPM") ;
            }
        } else {
            bpm_in = NULL ;
        }
        /* Combine the 2 BPMs */
        /* TODO */
        if (bpm_in != NULL) 
            cpl_image_add(bpm_loc, bpm_in) ;
    }

    /* Return the results */
    *master_flat = master_flat_loc ;
    if (bpm_loc != NULL) {
        *bpm = hdrl_image_create(bpm_loc, NULL) ;
        cpl_image_delete(bpm_loc) ;
    } else {
        *bpm = NULL ;
    }
    *ext_plist = plist ;
    return 0 ;
}

