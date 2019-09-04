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

#include "cr2res_utils.h"
#include "cr2res_calib.h"
#include "cr2res_pfits.h"
#include "cr2res_dfs.h"
#include "cr2res_flat.h"
#include "cr2res_trace.h"
#include "cr2res_extract.h"
#include "cr2res_io.h"

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

#define RECIPE_STRING "cr2res_util_calib"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static cpl_frameset * cr2res_util_calib_find_RAW(const cpl_frameset * in) ;
static int cr2res_util_calib_create(cpl_plugin *);
static int cr2res_util_calib_exec(cpl_plugin *);
static int cr2res_util_calib_destroy(cpl_plugin *);
static int cr2res_util_calib(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_util_calib_description[] = "\
Frames Calibration                                                      \n\
  Each input file is corrected with BPM / Dark / Flat / Det.Lin. / Cosmics\n\
                                                                        \n\
  Inputs                                                                \n\
    raw.fits " CR2RES_FLAT_RAW" [1 to n]                                \n\
          or " CR2RES_WAVE_RAW"                                         \n\
          or " CR2RES_OBS_NODDING_RAW"                                  \n\
          or " CR2RES_OBS_2D_RAW"                                       \n\
          or " CR2RES_OBS_POL_RAW"                                      \n\
    detlin.fits " CR2RES_CAL_DETLIN_COEFFS_PROCATG " [0 to 1]           \n\
    bpm.fits " CR2RES_CAL_DARK_BPM_PROCATG " [0 to 1]                   \n\
          or " CR2RES_CAL_FLAT_BPM_PROCATG "                            \n\
          or " CR2RES_CAL_DETLIN_BPM_PROCATG "                          \n\
          or " CR2RES_UTIL_BPM_SPLIT_PROCATG "                          \n\
    master_dark.fits " CR2RES_CAL_DARK_MASTER_PROCATG " [0 to 1]        \n\
    master_flat.fits " CR2RES_CAL_FLAT_MASTER_PROCATG " [0 to 1]        \n\
                                                                        \n\
  Outputs                                                               \n\
    <input_name>_calibrated.fits " CR2RES_UTIL_CALIB_PROCATG "          \n\
    or                                                                  \n\
    cr2res_util_calib_calibrated_collapsed.fits " 
    CR2RES_UTIL_CALIB_PROCATG "\n\
                                                                        \n\
  Algorithm                                                             \n\
    loop on detectors d:                                                \n\
      Call cr2res_calib_imagelist()() to calibrate --calib_cosmics_corr,\n\
               detlin, bpm, dark, flat                                  \n\
        -> calibrated(d)                                                \n\
      Collapse the calibrated image list to collapsed(d)                \n\
      if (collapse) save collapsed                                      \n\
      else save every individual calibrated frame                       \n\
  Library functions uѕed:                                               \n\
    cr2res_io_load_image()                                              \n\
    cr2res_calib_imagelist()                                            \n\
    cr2res_io_save_CALIBRATED()                                         \n\
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
                    "Calibration utility",
                    cr2res_util_calib_description,
                    "Thomas Marquart, Yves Jung",
                    PACKAGE_BUGREPORT,
                    cr2res_get_license(),
                    cr2res_util_calib_create,
                    cr2res_util_calib_exec,
                    cr2res_util_calib_destroy)) {    
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
static int cr2res_util_calib_create(cpl_plugin * plugin)
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
    p = cpl_parameter_new_value("cr2res.cr2res_util_calib.detector",
            CPL_TYPE_INT, "Only reduce the specified detector",
            "cr2res.cr2res_util_calib", 0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "detector");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_calib.clean_bad",
            CPL_TYPE_BOOL, "Apply the cleaning to the bad pixels",
            "cr2res.cr2res_util_calib", TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "clean_bad");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_calib.calib_cosmics_corr",
            CPL_TYPE_BOOL, "Correct the Cosmics",
            "cr2res.cr2res_util_calib", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "calib_cosmics_corr");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_calib.collapse",
            CPL_TYPE_STRING, "Collapse the input images (NONE, MEAN or MEDIAN)",
            "cr2res.cr2res_util_calib", "NONE");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "collapse");
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
static int cr2res_util_calib_exec(cpl_plugin * plugin)
{
    cpl_recipe  *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1;

    return cr2res_util_calib(recipe->frames, recipe->parameters);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_util_calib_destroy(cpl_plugin * plugin)
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
static int cr2res_util_calib(
        cpl_frameset            *   frameset,
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter *   param ;
    int                     clean_bad, calib_cosmics_corr, reduce_det ;
    cr2res_collapse         collapse ;
    const char          *   sval ;
    cpl_frameset        *   rawframes ;
    const cpl_frame     *   detlin_frame ;
    const cpl_frame     *   master_dark_frame ;
    const cpl_frame     *   master_flat_frame ;
    const cpl_frame     *   bpm_frame ;
    cpl_vector          *   dits ;
    hdrl_imagelist      *   in ;
    cpl_image           *   contrib ;
    const cpl_frame     *   cur_frame ;
    const char          *   cur_fname ;
    cpl_frameset        *   cur_fset ;
    hdrl_image          *   collapsed_ima[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   calibrated_one[CR2RES_NB_DETECTORS] ;
    hdrl_imagelist      *   calibrated[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   ext_plist[CR2RES_NB_DETECTORS] ;
    char                *   out_file;
    int                     i, det_nr ; 

    /* Needed for sscanf() */
    setlocale(LC_NUMERIC, "C");

    /* Initialise */
    contrib = NULL ;
    collapse = CR2RES_COLLAPSE_UNSPECIFIED ;

    /* RETRIEVE INPUT PARAMETERS */
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_calib.detector");
    reduce_det = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_calib.clean_bad");
    clean_bad = cpl_parameter_get_bool(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_calib.calib_cosmics_corr");
    calib_cosmics_corr = cpl_parameter_get_bool(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_calib.collapse");
    sval = cpl_parameter_get_string(param);
    if (!strcmp(sval, "NONE"))          collapse = CR2RES_COLLAPSE_NONE ;
    else if (!strcmp(sval, "MEAN"))     collapse = CR2RES_COLLAPSE_MEAN ;
    else if (!strcmp(sval, "MEDIAN"))   collapse = CR2RES_COLLAPSE_MEDIAN ;
    if (collapse == CR2RES_COLLAPSE_UNSPECIFIED) {
        cpl_msg_error(__func__, "Cannot understand the collapse method") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Identify the RAW and CALIB frames in the input frameset */
    if (cr2res_dfs_set_groups(frameset)) {
        cpl_msg_error(__func__, "Cannot identify RAW and CALIB frames") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Get Calibration frames */
    detlin_frame = cpl_frameset_find_const(frameset,
            CR2RES_CAL_DETLIN_COEFFS_PROCATG);
    master_dark_frame = cpl_frameset_find_const(frameset,
            CR2RES_CAL_DARK_MASTER_PROCATG) ; 
    master_flat_frame = cpl_frameset_find_const(frameset,
            CR2RES_CAL_FLAT_MASTER_PROCATG) ; 
    bpm_frame = cr2res_io_find_BPM(frameset) ;

    /* Get the rawframes */
    rawframes = cr2res_util_calib_find_RAW(frameset) ;
    if (rawframes==NULL || cpl_frameset_get_size(rawframes) <= 0) {
        cpl_msg_error(__func__, "Cannot find any RAW file") ;
        cpl_error_set(__func__, CPL_ERROR_DATA_NOT_FOUND) ;
        if (rawframes!= NULL) cpl_frameset_delete(rawframes) ;
        return -1 ;
    }

    /* Loop on the detectors */
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
        /* Initialise */
        collapsed_ima[det_nr-1] = NULL ;
        calibrated[det_nr-1] = NULL ;

        /* Compute only one detector */
        if (reduce_det != 0 && det_nr != reduce_det) continue ;
    
        cpl_msg_info(__func__, "Process Detector %d", det_nr) ;
        cpl_msg_indent_more() ;

        /* Load the DITs if necessary */
        if (master_dark_frame != NULL)  dits = cr2res_io_read_dits(rawframes) ;
        else                            dits = NULL ;
        if (cpl_msg_get_level() == CPL_MSG_DEBUG && dits != NULL)
            cpl_vector_dump(dits, stdout) ;

        /* Load image list */
        cpl_msg_info(__func__, "Load the input frames") ;
        if ((in = cr2res_io_load_image_list_from_set(rawframes,
                        det_nr)) == NULL) {
            cpl_msg_warning(__func__, "Cannot load images") ;
            if (dits != NULL) cpl_vector_delete(dits) ;
            cpl_msg_indent_less() ;
            continue ; 
        }
        if (hdrl_imagelist_get_size(in) != cpl_frameset_get_size(rawframes)) {
            cpl_msg_error(__func__, "Inconsistent number of loaded images") ;
            if (dits != NULL) cpl_vector_delete(dits) ;
            hdrl_imagelist_delete(in) ;
            cpl_msg_indent_less() ;
            continue ; 
        }

        /* Calibrate the images */
        cpl_msg_info(__func__, "Calibrate the input images") ;
        if ((calibrated[det_nr-1] = cr2res_calib_imagelist(in, det_nr,
                        clean_bad, 0, master_flat_frame, master_dark_frame, 
                        bpm_frame, detlin_frame, dits)) == NULL) {
            cpl_msg_warning(__func__, "Failed to apply the calibrations") ;
            if (dits != NULL) cpl_vector_delete(dits) ;
            hdrl_imagelist_delete(in) ;
            cpl_msg_indent_less() ;
            continue ;
        }
        hdrl_imagelist_delete(in) ;
        if (dits != NULL) cpl_vector_delete(dits) ;

        /* Collapse */
        if (collapse == CR2RES_COLLAPSE_MEAN) {
            cpl_msg_info(__func__, "Collapse (Mean) the calibrated images") ;
            cpl_msg_indent_more() ;
            hdrl_imagelist_collapse_mean(calibrated[det_nr-1],
                    &(collapsed_ima[det_nr-1]), &contrib) ;
        } else if (collapse == CR2RES_COLLAPSE_MEDIAN) {
            cpl_msg_info(__func__, "Collapse (Median) the calibrated images") ;
            cpl_msg_indent_more() ;
            hdrl_imagelist_collapse_median(calibrated[det_nr-1],
                    &(collapsed_ima[det_nr-1]), &contrib) ;
        }
        if (contrib != NULL) {
            cpl_image_delete(contrib) ;
            contrib = NULL ;
        }
        if (cpl_error_get_code() != CPL_ERROR_NONE) {
            collapsed_ima[det_nr-1] = NULL ;
            cpl_msg_warning(__func__, "Failed to Collapse") ;
            cpl_msg_indent_less() ;
            cpl_msg_indent_less() ;
            continue ;
        }
        cpl_msg_indent_less() ;
        cpl_msg_indent_less() ;
    }

    /* Ѕave Products */
    if (collapse == CR2RES_COLLAPSE_NONE) {
        /* Save individual calibrated images */
        /* Loop on the RAW frames */
        for (i=0 ; i<cpl_frameset_get_size(rawframes) ; i++) {
            /* Get the Current Frame */
            cur_frame = cpl_frameset_get_position(rawframes, i) ;
            cur_fname = cpl_frame_get_filename(cur_frame) ;
        
            /* Save CALIBRATED */
            out_file=cpl_sprintf("%s_calibrated.fits", 
                    cr2res_get_root_name(cr2res_get_base_name(cur_fname))) ;
            cur_fset = cpl_frameset_new() ;
            cpl_frameset_insert(cur_fset, cpl_frame_duplicate(cur_frame)) ;

            for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
                ext_plist[det_nr-1] = cpl_propertylist_load(cur_fname, det_nr) ;
                calibrated_one[det_nr-1] = hdrl_image_duplicate(
                        hdrl_imagelist_get_const(calibrated[det_nr-1], i)) ;
            }
            cr2res_io_save_CALIBRATED(out_file, frameset, cur_fset, parlist,
                    calibrated_one, NULL, ext_plist, CR2RES_UTIL_CALIB_PROCATG, 
                    RECIPE_STRING) ;
            for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
                if (ext_plist[det_nr-1] != NULL) 
                    cpl_propertylist_delete(ext_plist[det_nr-1]) ;
                if (calibrated_one[det_nr-1] != NULL)
                    hdrl_image_delete(calibrated_one[det_nr-1]) ;
            }
            cpl_frameset_delete(cur_fset) ;
            cpl_free(out_file);
        }
    } else {
        /* Save COLLAPSED calibrated image */
        cur_frame = cpl_frameset_get_position(rawframes, 0) ;
        cur_fname = cpl_frame_get_filename(cur_frame) ;
        for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
            ext_plist[det_nr-1] = cpl_propertylist_load(cur_fname, det_nr) ;
        }
        out_file=cpl_sprintf("%s_calibrated_collapsed.fits", RECIPE_STRING) ;
        cr2res_io_save_CALIBRATED(out_file, frameset, rawframes, parlist,
                collapsed_ima, NULL, ext_plist, CR2RES_UTIL_CALIB_PROCATG, 
                RECIPE_STRING) ;
        for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
            if (ext_plist[det_nr-1] != NULL) 
                cpl_propertylist_delete(ext_plist[det_nr-1]) ;
        }
        cpl_free(out_file);
    }

    /* Free */
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
        if (collapsed_ima[det_nr-1] != NULL)
            hdrl_image_delete(collapsed_ima[det_nr-1]) ;
        if (calibrated[det_nr-1] != NULL)
            hdrl_imagelist_delete(calibrated[det_nr-1]) ;
    }
    cpl_frameset_delete(rawframes) ;
    return (int)cpl_error_get_code();
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the RAW frames from a frameset
  @param    set     Input frame set
  @return   the RAW frameset or NULL in error case or if it is missing
    Allowed RAW types : CR2RES_FLAT_RAW
                        CR2RES_WAVE_RAW
                        CR2RES_OBS_NODDING_RAW
                        CR2RES_OBS_2D_RAW
                        CR2RES_OBS_POL_RAW
 */
/*----------------------------------------------------------------------------*/
static cpl_frameset * cr2res_util_calib_find_RAW(const cpl_frameset * in)
{
    cpl_frameset    *   out ;

    /* Check entries */
    if (in == NULL) return NULL ;

    out = cr2res_extract_frameset(in, CR2RES_FLAT_RAW) ;
    if (out == NULL)    
        out = cr2res_extract_frameset(in, CR2RES_WAVE_RAW) ;
    if (out == NULL)    
        out = cr2res_extract_frameset(in, CR2RES_OBS_NODDING_RAW) ;
    if (out == NULL)
        out = cr2res_extract_frameset(in, CR2RES_OBS_2D_RAW) ;
    if (out == NULL)
        out = cr2res_extract_frameset(in, CR2RES_OBS_POL_RAW) ;
    return out ;
}




