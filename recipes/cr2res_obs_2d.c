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
#include "cr2res_calib.h"
#include "cr2res_pfits.h"
#include "cr2res_dfs.h"
#include "cr2res_bpm.h"
#include "cr2res_trace.h"
#include "cr2res_extract.h"
#include "cr2res_io.h"

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

#define RECIPE_STRING "cr2res_obs_2d"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_obs_2d_check_inputs_validity(
        const cpl_frame *   rawframe) ;
static int cr2res_obs_2d_reduce(
        const cpl_frame     *   rawframe,
        const cpl_frame     *   trace_wave_frame,
        const cpl_frame     *   detlin_frame,
        const cpl_frame     *   master_dark_frame,
        const cpl_frame     *   master_flat_frame,
        const cpl_frame     *   bpm_frame,
        int                     calib_cosmics_corr,
        int                     reduce_det,
        int                     reduce_order,
        int                     reduce_trace,
        cpl_table           **  extract,
        cpl_propertylist    **  ext_plist) ;

static int cr2res_obs_2d_create(cpl_plugin *);
static int cr2res_obs_2d_exec(cpl_plugin *);
static int cr2res_obs_2d_destroy(cpl_plugin *);
static int cr2res_obs_2d(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_obs_2d_description[] = "\
2D Observation                                                          \n\
  This recipe is meant for extended objects. In each trace, the pixels  \n\
  are calibrated, and stored in the output table.                       \n\
                                                                        \n\
  Inputs                                                                \n\
    raw.fits " CR2RES_OBS_2D_RAW"                                       \n\
    trace.fits " CR2RES_CAL_FLAT_TW_PROCATG " [1]                       \n\
            or " CR2RES_CAL_FLAT_TW_MERGED_PROCATG "                    \n\
            or " CR2RES_UTIL_TRACE_TW_PROCATG "                         \n\
            or " CR2RES_UTIL_WAVE_TW_PROCATG "                          \n\
            or " CR2RES_CAL_WAVE_TW_PROCATG "                           \n\
            or " CR2RES_UTIL_SLIT_CURV_TW_PROCATG "                     \n\
    detlin.fits " CR2RES_CAL_DETLIN_COEFFS_PROCATG " [0 to 1]           \n\
    bpm.fits " CR2RES_CAL_DARK_BPM_PROCATG " [0 to 1]                   \n\
          or " CR2RES_CAL_FLAT_BPM_PROCATG "                            \n\
          or " CR2RES_CAL_DETLIN_BPM_PROCATG "                          \n\
          or " CR2RES_UTIL_BPM_SPLIT_PROCATG "                          \n\
    master_dark.fits " CR2RES_CAL_DARK_MASTER_PROCATG " [0 to 1]        \n\
    master_flat.fits " CR2RES_CAL_FLAT_MASTER_PROCATG " [0 to 1]        \n\
                                                                        \n\
  Outputs                                                               \n\
    cr2res_obs_2d_extract.fits " CR2RES_OBS_2D_EXTRACT_PROCATG "        \n\
                                                                        \n\
  Algorithm                                                             \n\
    loop on raw frames f:                                               \n\
      loop on detectors d:                                              \n\
        cr2res_obs_2d_reduce()                                          \n\
          -> extract(d)                                                 \n\
      Save extract                                                      \n\
                                                                        \n\
    cr2res_obs_2d_reduce()                                              \n\
      Load the input image                                              \n\
      Apply the calibrations to the image                               \n\
      Load the trace wave table                                         \n\
      Call cr2res_extract2d_traces()                                    \n\
        -> extract                                                      \n\
                                                                        \n\
  Library Functions uѕed                                                \n\
    cr2res_io_find_TRACE_WAVE()                                         \n\
    cr2res_io_find_BPM()                                                \n\
    cr2res_obs_2d_reduce()                                              \n\
    cr2res_obs_2d_check_inputs_validity()                               \n\
    cr2res_pfits_get_dit()                                              \n\
    cr2res_io_load_image()                                              \n\
    cr2res_calib_image()                                                \n\
    cr2res_io_load_TRACE_WAVE()                                         \n\
    cr2res_extract2d_traces()                                           \n\
    cr2res_io_save_EXTRACT_2D()                                         \n\
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
                    "2D Observation recipe",
                    cr2res_obs_2d_description,
                    "Thomas Marquart, Yves Jung",
                    PACKAGE_BUGREPORT,
                    cr2res_get_license(),
                    cr2res_obs_2d_create,
                    cr2res_obs_2d_exec,
                    cr2res_obs_2d_destroy)) {    
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
static int cr2res_obs_2d_create(cpl_plugin * plugin)
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
    p = cpl_parameter_new_value("cr2res.cr2res_obs_2d.detector",
            CPL_TYPE_INT, "Only reduce the specified detector",
            "cr2res.cr2res_obs_2d", 0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "detector");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_2d.order",
            CPL_TYPE_INT, "Only reduce the specified order",
            "cr2res.cr2res_obs_2d", -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "order");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_2d.trace_nb",
            CPL_TYPE_INT, "Only reduce the specified trace number",
            "cr2res.cr2res_obs_2d", -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "trace_nb");
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
static int cr2res_obs_2d_exec(cpl_plugin * plugin)
{
    cpl_recipe  *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1;

    return cr2res_obs_2d(recipe->frames, recipe->parameters);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_obs_2d_destroy(cpl_plugin * plugin)
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
static int cr2res_obs_2d(
        cpl_frameset            *   frameset,
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter *   param ;
    int                     reduce_det, reduce_order, 
                            reduce_trace ;
    cpl_frameset        *   rawframes ;
    cpl_frame           *   rawframe ;
    const cpl_frame     *   detlin_frame ;
    const cpl_frame     *   master_dark_frame ;
    const cpl_frame     *   master_flat_frame ;
    const cpl_frame     *   bpm_frame ;
    const cpl_frame     *   trace_wave_frame ;
    cpl_propertylist    *   ext_plist[CR2RES_NB_DETECTORS] ;
    cpl_table           *   extract[CR2RES_NB_DETECTORS] ;
    char                *   out_file;
    int                     i, det_nr; 

    /* Initialise */

    /* RETRIEVE INPUT PARAMETERS */
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_2d.detector");
    reduce_det = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_2d.order");
    reduce_order = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_2d.trace_nb");
    reduce_trace = cpl_parameter_get_int(param);

    /* Identify the RAW and CALIB frames in the input frameset */
    if (cr2res_dfs_set_groups(frameset)) {
        cpl_msg_error(__func__, "Cannot identify RAW and CALIB frames") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Get Calibration frames */
    trace_wave_frame = cr2res_io_find_TRACE_WAVE(frameset) ;
    if (trace_wave_frame == NULL) {
        cpl_msg_error(__func__, "Could not find TRACE_WAVE frame") ;
        return -1 ;
    }
    detlin_frame = cpl_frameset_find_const(frameset,
            CR2RES_CAL_DETLIN_COEFFS_PROCATG);
    master_dark_frame = cpl_frameset_find_const(frameset,
            CR2RES_CAL_DARK_MASTER_PROCATG) ; 
    master_flat_frame = cpl_frameset_find_const(frameset,
            CR2RES_CAL_FLAT_MASTER_PROCATG) ; 
    bpm_frame = cr2res_io_find_BPM(frameset) ;

    /* Get the Frames for the current decker position */
    rawframes = cr2res_extract_frameset(frameset, CR2RES_OBS_2D_RAW) ;
    if (rawframes == NULL) {
        cpl_msg_error(__func__, "Could not find RAW frames") ;
        return -1 ;
    }
    cpl_msg_indent_more() ;

    /* Loop on the RAW files */
    for (i=0 ; i<cpl_frameset_get_size(rawframes) ; i++) {
        
        /* Current frame */
        rawframe = cpl_frameset_get_position(rawframes, i);

        /* Loop on the detectors */
        for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
            /* Initialise */
            extract[det_nr-1] = NULL ;
            ext_plist[det_nr-1] = NULL ;

            /* Compute only one detector */
            if (reduce_det != 0 && det_nr != reduce_det) continue ;
        
            cpl_msg_info(__func__, "Process Detector %d", det_nr) ;
            cpl_msg_indent_more() ;
            
            /* Call the reduction function */
            if (cr2res_obs_2d_reduce(rawframe, trace_wave_frame, detlin_frame, 
                        master_dark_frame, master_flat_frame, bpm_frame,
                        0, det_nr, reduce_order, reduce_trace,
                        &(extract[det_nr-1]),
                        &(ext_plist[det_nr-1])) == -1) {
                cpl_msg_warning(__func__, "Failed to reduce detector %d", 
                        det_nr);
            }
            cpl_msg_indent_less() ;
        }

        /* Ѕave Products */

        /* Extracted */
        out_file = cpl_sprintf("%s_frame_%d_extracted.fits", 
                RECIPE_STRING, i+1) ;
        cr2res_io_save_EXTRACT_2D(out_file, frameset, rawframes, parlist, 
                extract, NULL, ext_plist, CR2RES_OBS_2D_EXTRACT_PROCATG, 
                RECIPE_STRING) ;
        cpl_free(out_file);

        /* Free */
        for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
            if (extract[det_nr-1] != NULL) 
                cpl_table_delete(extract[det_nr-1]) ;
            if (ext_plist[det_nr-1] != NULL) 
                cpl_propertylist_delete(ext_plist[det_nr-1]) ;
        }
        cpl_msg_indent_less() ;
    }
    cpl_frameset_delete(rawframes) ;
    return (int)cpl_error_get_code();
}
 
/*----------------------------------------------------------------------------*/
/**
  @brief  Execute the 2d observation on one detector
  @param rawframe               Raw science frame
  @param trace_wave_frame       Trace Wave file
  @param detlin_frame           Associated detlin coefficients
  @param master_dark_frame      Associated master dark
  @param master_flat_frame      Associated master flat
  @param bpm_frame              Associated BPM
  @param calib_cosmics_corr     Flag to correct for cosmics
  @param reduce_det             The detector to compute
  @param reduce_order           The order to reduce (-1 for all)
  @param reduce_trace           The trace to reduce (-1 for all)
  @param extract                [out] extracted spectrum
  @param ext_plist              [out] the header for saving the products
  @return  0 if ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
static int cr2res_obs_2d_reduce(
        const cpl_frame     *   rawframe,
        const cpl_frame     *   trace_wave_frame,
        const cpl_frame     *   detlin_frame,
        const cpl_frame     *   master_dark_frame,
        const cpl_frame     *   master_flat_frame,
        const cpl_frame     *   bpm_frame,
        int                     calib_cosmics_corr,
        int                     reduce_det,
        int                     reduce_order,
        int                     reduce_trace,
        cpl_table           **  extract,
        cpl_propertylist    **  ext_plist) 
{
    hdrl_image          *   in ;
    hdrl_image          *   in_calib ;
    cpl_table           *   trace_wave ;
    cpl_propertylist    *   plist ;
    double                  dit ;
    cpl_table           *   extracted ;
    cpl_size                i ;
    int                     det_nr ;

    /* Check Inputs */
    if (extract == NULL || ext_plist == NULL || rawframe == NULL || 
            trace_wave_frame == NULL) return -1 ;

    /* Check raw frames consistency */
    if (cr2res_obs_2d_check_inputs_validity(rawframe) != 1) {
        cpl_msg_error(__func__, "Invalid Inputs") ;
        return -1 ;
    }

    /* Get the DIT */
    plist = cpl_propertylist_load(cpl_frame_get_filename(rawframe), 0) ;
    dit = cr2res_pfits_get_dit(plist) ;
    cpl_propertylist_delete(plist); 
    if (cpl_error_get_code()) {
        cpl_msg_error(__func__, "Cannot read the DIT") ;
        return -1 ;
    }
    cpl_msg_debug(__func__, "DIT value : %g", dit) ;

    /* Load the input image */
    if ((in = cr2res_io_load_image(cpl_frame_get_filename(rawframe),
                    reduce_det)) == NULL) {
        cpl_msg_error(__func__, "Cannot load image") ;
        return -1 ;
    }

    /* Calibrate the image */
    if ((in_calib = cr2res_calib_image(in, reduce_det, 0, master_flat_frame,
            master_dark_frame, bpm_frame, detlin_frame, dit)) == NULL) {
        cpl_msg_error(__func__, "Failed to apply the calibrations") ;
        hdrl_image_delete(in) ;
        return -1 ;
    }
    hdrl_image_delete(in) ;

    /* Load the trace wave */
    cpl_msg_info(__func__, "Load the TRACE WAVE") ;
    if ((trace_wave = cr2res_io_load_TRACE_WAVE(
                    cpl_frame_get_filename(trace_wave_frame), 
                    reduce_det)) == NULL) {
        cpl_msg_error(__func__, "Failed to Load the traces file") ;
        hdrl_image_delete(in_calib) ;
        return -1 ;
    }

    /* Execute the extraction */
    cpl_msg_info(__func__, "Spectra Extraction 2D") ;
    if (cr2res_extract2d_traces(in_calib, trace_wave, reduce_order,
                reduce_trace, &extracted) == -1) {
        cpl_msg_error(__func__, "Failed to extract");
        hdrl_image_delete(in_calib) ;
        cpl_table_delete(trace_wave) ;
        return -1 ;
    }
    cpl_table_delete(trace_wave) ;
    hdrl_image_delete(in_calib) ;

    /* QC parameters */
    plist = cpl_propertylist_new() ;

    /* Compute the QC parameters */
    /* TODO */

    /* Return */
    *extract = extracted ;
    *ext_plist = plist ;

    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief  Run basic checks for the rawframe consistency
  @param    rawframe   The input rawframe
  @return   1 if valid, 0 if not, -1 in error case
 */
/*----------------------------------------------------------------------------*/
static int cr2res_obs_2d_check_inputs_validity(
        const cpl_frame *   rawframe)
{
    /* TODO */

    /* Check Inputs */
    if (rawframe == NULL) return -1 ;
    return 1 ;
}


