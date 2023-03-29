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
#include "cr2res_idp.h"
#include "cr2res_calib.h"
#include "cr2res_pfits.h"
#include "cr2res_dfs.h"
#include "cr2res_bpm.h"
#include "cr2res_trace.h"
#include "cr2res_extract.h"
#include "cr2res_io.h"
#include "cr2res_qc.h"

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

#define RECIPE_STRING "cr2res_obs_staring"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static cpl_frameset * cr2res_obs_staring_find_RAW(
        const cpl_frameset  *   in) ;
static int cr2res_obs_staring_check_inputs_validity(
        const cpl_frameset  *   rawframes) ;
static int cr2res_obs_staring_reduce(
        const cpl_frameset  *   rawframes,
        const cpl_frame     *   trace_wave_frame,
        const cpl_frame     *   detlin_frame,
        const cpl_frame     *   master_dark_frame,
        const cpl_frame     *   master_flat_frame,
        const cpl_frame     *   bpm_frame,
        const cpl_frame     *   blaze_frame,
        const cpl_array     *   slit_frac,
        int                     subtract_nolight_rows,
        int                     subtract_interorder_column,
        int                     calib_cosmics_corr,
        int                     extract_oversample,
        int                     extract_swath_width,
        int                     extract_height,
        double                  extract_smooth_slit,
        double                  extract_smooth_spec,
        int                     reduce_det,
        cpl_table           **  extract,
        cpl_table           **  slitfunc,
        hdrl_image          **  model,
        cpl_propertylist    **  ext_plist) ;
static int cr2res_obs_staring_create(cpl_plugin *);
static int cr2res_obs_staring_exec(cpl_plugin *);
static int cr2res_obs_staring_destroy(cpl_plugin *);
static int cr2res_obs_staring(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_obs_staring_description[] = "\
Staring Observation                                                     \n\
  This recipe handles staring observations.                             \n\
                                                                        \n\
  Inputs                                                                \n\
    raw.fits " CR2RES_OBS_STARING_OTHER_RAW" [1 to n]                   \n\
          or " CR2RES_OBS_STARING_JITTER_RAW" [1 to n]                  \n\
          or " CR2RES_OBS_STARING_WAVE_SKY_RAW" [1 to n]                \n\
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
          or " CR2RES_UTIL_BPM_MERGE_PROCATG "                          \n\
          or " CR2RES_UTIL_BPM_SPLIT_PROCATG "                          \n\
    master_dark.fits " CR2RES_CAL_DARK_MASTER_PROCATG " [0 to 1]        \n\
    master_flat.fits " CR2RES_CAL_FLAT_MASTER_PROCATG " [0 to 1]        \n\
    blaze.fits " CR2RES_CAL_FLAT_EXTRACT_1D_PROCATG " [0 to 1]          \n\
                                                                        \n\
  Outputs                                                               \n\
	cr2res_obs_staring_slitfunc.fits "
	CR2RES_OBS_STARING_SLITFUNC_PROCATG "\n\
	cr2res_obs_staring_model.fits "
	CR2RES_OBS_STARING_SLITMODEL_PROCATG "\n\
	cr2res_obs_staring_extracted.fits "
	CR2RES_OBS_STARING_EXTRACT_PROCATG "\n\
                                                                        \n\
  Algorithm                                                             \n\
    loop on detectors d:                                                \n\
      call cr2res_obs_staring_reduce()                                  \n\
        -> extract(d)                                                   \n\
        -> slitfunc(d)                                                  \n\
        -> model(d)                                                     \n\
    Save extract                                                        \n\
    Save slitfunc                                                       \n\
    Save model                                                          \n\
                                                                        \n\
    cr2res_obs_staring_reduce()                                         \n\
      Load the input raw frames in an image list                        \n\
      Apply the calibrations to the image list                          \n\
      Collapse the image list                                           \n\
      Load the input trace wave                                         \n\
      Recompute a new trace wave with the specified slit fraction       \n\
             (--slit_frac) if needed                                    \n\
      Extract the spectra from the collapsed image                      \n\
        -> extracted                                                    \n\
        -> slit_func                                                    \n\
        -> model_master                                                 \n\
      Compute QC parameters                                             \n\
                                                                        \n\
  Library functions used                                                \n\
    cr2res_io_find_TRACE_WAVE()                                         \n\
    cr2res_io_find_BPM()                                                \n\
    cr2res_obs_staring_reduce()                                         \n\
    cr2res_io_read_dits()                                               \n\
    cr2res_io_load_image_list_from_set()                                \n\
    cr2res_calib_imagelist()                                            \n\
    cr2res_io_load_TRACE_WAVE()                                         \n\
    cr2res_extract_traces()                                             \n\
    cr2res_io_save_COMBINED()                                           \n\
    cr2res_io_save_EXTRACT_1D()                                         \n\
    cr2res_io_save_SLIT_FUNC()                                          \n\
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
                    "Staring Observation recipe",
                    cr2res_obs_staring_description,
                    CR2RES_PIPELINE_AUTHORS,
                    PACKAGE_BUGREPORT,
                    cr2res_get_license(),
                    cr2res_obs_staring_create,
                    cr2res_obs_staring_exec,
                    cr2res_obs_staring_destroy)) {    
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
static int cr2res_obs_staring_create(cpl_plugin * plugin)
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

    p = cpl_parameter_new_value(
            "cr2res.cr2res_obs_staring.subtract_nolight_rows",
            CPL_TYPE_BOOL,
            "Subtract median row from baffled region at detector bottom",
            "cr2res.cr2res_obs_staring", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "subtract_nolight_rows");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value(
            "cr2res.cr2res_obs_staring.subtract_interorder_column",
            CPL_TYPE_BOOL,
            "Subtract column-by-column fit to the pixel values between"
            " spectral orders",
            "cr2res.cr2res_obs_staring", TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
                                                "subtract_interorder_column");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_staring.slit_frac",
            CPL_TYPE_STRING, "Wished slit fraction",
            "cr2res.cr2res_obs_staring", "-1.0, -1.0");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "slit_frac");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_staring.extract_oversample",
            CPL_TYPE_INT, "factor by which to oversample the extraction",
            "cr2res.cr2res_obs_staring", 5);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_oversample");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_staring.extract_swath_width",
            CPL_TYPE_INT, "The swath width", "cr2res.cr2res_obs_staring", 800);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_swath_width");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_staring.extract_height",
            CPL_TYPE_INT, "Extraction height",
            "cr2res.cr2res_obs_staring", -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_height");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_staring.extract_smooth_slit",
            CPL_TYPE_DOUBLE,
            "Smoothing along the slit (1 for high S/N, 5 for low)",
            "cr2res.cr2res_obs_staring", 2.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_smooth_slit");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_staring.extract_smooth_spec",
            CPL_TYPE_DOUBLE, "Smoothing along the spectrum",
            "cr2res.cr2res_obs_staring", 0.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_smooth_spec");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_staring.detector",
            CPL_TYPE_INT, "Only reduce the specified detector",
            "cr2res.cr2res_obs_staring", 0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "detector");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_staring.create_idp",
            CPL_TYPE_BOOL, "Flag to produce  IDP files",
            "cr2res.cr2res_obs_staring", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "idp");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_staring.display_order",
            CPL_TYPE_INT, "Apply the display for the specified order",
            "cr2res.cr2res_obs_staring", 0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "display_order");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_staring.display_trace",
            CPL_TYPE_INT, "Apply the display for the specified trace",
            "cr2res.cr2res_obs_staring", 0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "display_trace");
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
static int cr2res_obs_staring_exec(cpl_plugin * plugin)
{
    cpl_recipe  *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1;

    return cr2res_obs_staring(recipe->frames, recipe->parameters);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_obs_staring_destroy(cpl_plugin * plugin)
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
static int cr2res_obs_staring(
        cpl_frameset            *   frameset,
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter *   param ;
    int                     subtract_nolight_rows, subtract_interorder_column,
                            extract_oversample, create_idp,
                            extract_swath_width, extract_height, reduce_det, 
                            disp_order, disp_trace ;
    double                  extract_smooth_slit, extract_smooth_spec, 
                            slit_low, slit_up ;
    cpl_array           *   slit_frac ;
    cpl_frameset        *   rawframes ;
    const cpl_frame     *   trace_wave_frame ;
    const cpl_frame     *   detlin_frame ;
    const cpl_frame     *   master_dark_frame ;
    const cpl_frame     *   master_flat_frame ;
    const cpl_frame     *   bpm_frame ;
    const cpl_frame     *   blaze_frame ;
    const char          *   sval ;
    cpl_table           *   extract[CR2RES_NB_DETECTORS] ;
    cpl_table           *   slitfunc[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   model[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   qc_main ;
    cpl_propertylist    *   ext_plist[CR2RES_NB_DETECTORS] ;
    char                *   out_file;
    int                     det_nr; 


    /* RETRIEVE INPUT PARAMETERS */
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_staring.subtract_nolight_rows");
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_staring.subtract_interorder_column");
    subtract_interorder_column = cpl_parameter_get_bool(param);
    subtract_nolight_rows = cpl_parameter_get_bool(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_staring.extract_oversample");
    extract_oversample = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_staring.extract_swath_width");
    extract_swath_width = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_staring.extract_height");
    extract_height = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_staring.extract_smooth_slit");
    extract_smooth_slit = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_staring.extract_smooth_spec");
    extract_smooth_spec = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_staring.detector");
    reduce_det = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_staring.create_idp");
    create_idp = cpl_parameter_get_bool(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_staring.display_order");
    disp_order = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_staring.display_trace");
    disp_trace = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_staring.slit_frac");
    sval = cpl_parameter_get_string(param) ;
    if (sscanf(sval, "%lg,%lg", &slit_low, &slit_up) != 2) {
        cpl_msg_error(__func__, "Invalid Slit Fraction specified");
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Check Parameters */
    if (slit_low >= 0.0 && slit_up >= 0.0 && slit_low <= 1.0 && slit_up <= 1.0
            && slit_up > slit_low) {
        slit_frac = cpl_array_new(3, CPL_TYPE_DOUBLE) ;
        cpl_array_set(slit_frac, 0, slit_low) ;
        cpl_array_set(slit_frac, 1, (slit_low+slit_up)/2.0) ;
        cpl_array_set(slit_frac, 2, slit_up) ;
    } else {
        slit_frac = NULL ;
    }

    /* Identify the RAW and CALIB frames in the input frameset */
    if (cr2res_dfs_set_groups(frameset)) {
        if (slit_frac != NULL) cpl_array_delete(slit_frac) ;
        cpl_msg_error(__func__, "Cannot identify RAW and CALIB frames") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Get Calibration frames */
    trace_wave_frame = cr2res_io_find_TRACE_WAVE(frameset) ;
    if (trace_wave_frame == NULL) {
		if (slit_frac != NULL) cpl_array_delete(slit_frac) ;
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
    blaze_frame = cpl_frameset_find_const(frameset,
            CR2RES_CAL_FLAT_EXTRACT_1D_PROCATG) ;

    /* Get the RAW Frames */
    rawframes = cr2res_obs_staring_find_RAW(frameset) ;
    if (rawframes == NULL) {
		if (slit_frac != NULL) cpl_array_delete(slit_frac) ;
        cpl_msg_error(__func__, "Could not find RAW frames") ;
        return -1 ;
    }
      
    /* Loop on the detectors */
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
        /* Initialise */
        extract[det_nr-1] = NULL ;
        slitfunc[det_nr-1] = NULL ;
        model[det_nr-1] = NULL ;
        ext_plist[det_nr-1] = NULL ;

        /* Compute only one detector */
        if (reduce_det != 0 && det_nr != reduce_det) continue ;
    
        cpl_msg_info(__func__, "Process Detector %d", det_nr) ;
        cpl_msg_indent_more() ;

        /* Call the reduction function */
        if (cr2res_obs_staring_reduce(rawframes, 
                    trace_wave_frame, detlin_frame, master_dark_frame, 
                    master_flat_frame, bpm_frame, blaze_frame, slit_frac, 
                    subtract_nolight_rows, subtract_interorder_column,
                    0, extract_oversample, 
                    extract_swath_width, extract_height, extract_smooth_slit, 
                    extract_smooth_spec, det_nr,
                    &(extract[det_nr-1]),
                    &(slitfunc[det_nr-1]),
                    &(model[det_nr-1]),
                    &(ext_plist[det_nr-1])) == -1) {
            cpl_msg_warning(__func__, "Failed to reduce detector %d", det_nr);
            cpl_error_reset() ;
            cpl_msg_indent_less() ;
        }
        cpl_msg_indent_less() ;
    }
    if (slit_frac != NULL) cpl_array_delete(slit_frac) ;

    /* Save Products */

    /* Add ESO.DRS.TMID in the Main Header */
    qc_main = cpl_propertylist_new();
    cpl_propertylist_append_double(qc_main,
            CR2RES_HEADER_DRS_TMID,
            cr2res_utils_get_center_mjd(rawframes)) ;

    out_file = cpl_sprintf("%s_slitfunc.fits", RECIPE_STRING) ;
    cr2res_io_save_SLIT_FUNC(out_file, frameset, frameset, parlist,
            slitfunc, qc_main, ext_plist, CR2RES_OBS_STARING_SLITFUNC_PROCATG,
            RECIPE_STRING) ;
    cpl_free(out_file);

    out_file = cpl_sprintf("%s_model.fits", RECIPE_STRING) ;
    cr2res_io_save_SLIT_MODEL(out_file, frameset, frameset, parlist,
            model, qc_main, ext_plist, CR2RES_OBS_STARING_SLITMODEL_PROCATG,
            RECIPE_STRING) ;
    cpl_free(out_file);

    out_file = cpl_sprintf("%s_extracted.fits", RECIPE_STRING) ;
    cr2res_io_save_EXTRACT_1D(out_file, frameset, frameset, parlist, extract,
            qc_main, ext_plist, CR2RES_OBS_STARING_EXTRACT_PROCATG,
            RECIPE_STRING);
	if (create_idp) {
        cr2res_idp_save(out_file, frameset, rawframes, parlist, 
                extract, ext_plist, RECIPE_STRING) ;
		}
    cpl_free(out_file);
    cpl_frameset_delete(rawframes) ;

    /* Free */
    cpl_propertylist_delete(qc_main) ;
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
        if (extract[det_nr-1] != NULL) 
            cpl_table_delete(extract[det_nr-1]) ;
        if (slitfunc[det_nr-1] != NULL) 
            cpl_table_delete(slitfunc[det_nr-1]) ;
        if (model[det_nr-1] != NULL)
            hdrl_image_delete(model[det_nr-1]) ;
        if (ext_plist[det_nr-1] != NULL) 
            cpl_propertylist_delete(ext_plist[det_nr-1]) ;
    }

    return (int)cpl_error_get_code();
}
 
/*----------------------------------------------------------------------------*/
/**
  @brief    Execute the science recipe on a specific detector
  @param rawframes              Raw science frames
  @param trace_wave_frame       Trace Wave file
  @param detlin_frame           Associated detlin coefficients
  @param master_dark_frame      Associated master dark
  @param master_flat_frame      Associated master flat
  @param bpm_frame              Associated BPM
  @param blaze_frame            Associated Blaze
  @param slit_frac              Specified slit fraction or NULL
  @param subtract_nolight_rows
  @param calib_cosmics_corr     Flag to correct for cosmics
  @param extract_oversample     Extraction related
  @param extract_swath_width    Extraction related
  @param extract_height         Extraction related
  @param extract_smooth_slit    Extraction related
  @param extract_smooth_spec    Extraction related
  @param reduce_det             The detector to compute
  @param extract                [out] extracted spectrum 
  @param slitfunc               [out] slit function
  @param model                  [out] slit model
  @param ext_plist              [out] the header for saving the products
  @return  0 if ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
static int cr2res_obs_staring_reduce(
        const cpl_frameset  *   rawframes,
        const cpl_frame     *   trace_wave_frame,
        const cpl_frame     *   detlin_frame,
        const cpl_frame     *   master_dark_frame,
        const cpl_frame     *   master_flat_frame,
        const cpl_frame     *   bpm_frame,
        const cpl_frame     *   blaze_frame,
        const cpl_array     *   slit_frac,
        int                     subtract_nolight_rows,
        int                     subtract_interorder_column,
        int                     calib_cosmics_corr,
        int                     extract_oversample,
        int                     extract_swath_width,
        int                     extract_height,
        double                  extract_smooth_slit,
        double                  extract_smooth_spec,
        int                     reduce_det,
        cpl_table           **  extract,
        cpl_table           **  slitfunc,
        hdrl_image          **  model,
        cpl_propertylist    **  ext_plist)
{
    hdrl_imagelist      *   in ;
    hdrl_imagelist      *   in_calib ;
    cpl_vector          *   dits ;
    cpl_vector          *   ndits ;
    cpl_table           *   blaze_table ;
    cpl_table           *   trace_wave ;
    cpl_table           *   trace_wave_new ;
    hdrl_image          *   collapsed ;
    cpl_image           *   contrib ;
    cpl_propertylist    *   plist ;
    cpl_size                i ;
    hdrl_image          *   model_master ;
    cpl_table           *   slit_func ;
    cpl_table           *   extracted ;
    int                 *   order_idx_values ;
    double              *   qc_snrs ;
    cpl_array           *   fwhm_array ;
    char                *   key_name ;
    const char          *   first_fname ;
    double                  qc_signal, qc_fwhm ;
    int                     order_zp, nb_order_idx_values,
                            order_real, order_idx, order_idxp ;

    /* Check Inputs */
    if (extract == NULL || ext_plist == NULL || rawframes == NULL
            || trace_wave_frame == NULL) return -1 ;

    /* Check raw frames consistency */
    if (cr2res_obs_staring_check_inputs_validity(rawframes) != 1) {
        cpl_msg_error(__func__, "Invalid Inputs") ;
        return -1 ;
    }

    /* Initialise */
    first_fname = cpl_frame_get_filename(
            cpl_frameset_get_position_const(rawframes, 0)) ;

    /* Get the order zeropoint */
    if ((plist = cpl_propertylist_load(cpl_frame_get_filename(trace_wave_frame),
                    0)) == NULL) {
        cpl_msg_error(__func__, "Cannot read the ORDER_ZP from the input TW") ;
        return -1 ;
    }
    order_zp = cr2res_pfits_get_order_zp(plist) ;
    cpl_propertylist_delete(plist) ;
    if (cpl_error_get_code()) {
        cpl_msg_error(__func__, "Missing ORDER_ZP in the header - Skip") ;
        cpl_error_reset() ;
        /* Negative Zerop to log the fact that it is missing */
        order_zp = -100 ;
    }

    /* Load the DITs if necessary */
    if (master_dark_frame != NULL)  dits = cr2res_io_read_dits(rawframes) ;
    else                            dits = NULL ;
    if (cpl_msg_get_level() == CPL_MSG_DEBUG && dits != NULL) 
        cpl_vector_dump(dits, stdout) ;

    /* Load NDITs */
    ndits = cr2res_io_read_ndits(rawframes) ;

    /* Load image list */
    if ((in = cr2res_io_load_image_list_from_set(rawframes, 
                    reduce_det)) == NULL) {
        cpl_msg_error(__func__, "Cannot load images") ;
        if (dits != NULL) cpl_vector_delete(dits) ;
        return -1 ;
    }
    if (hdrl_imagelist_get_size(in) != cpl_frameset_get_size(rawframes)) {
        cpl_msg_error(__func__, "Inconsistent number of loaded images") ;
        if (dits != NULL) cpl_vector_delete(dits) ;
        hdrl_imagelist_delete(in) ;
        return -1 ;
    }

    /* Calibrate the images */
    if ((in_calib = cr2res_calib_imagelist(in, reduce_det, 0,
            subtract_nolight_rows, subtract_interorder_column, 0, master_flat_frame, 
            master_dark_frame, bpm_frame, detlin_frame, dits, ndits))==NULL) {
        cpl_msg_error(__func__, "Failed to apply the calibrations") ;
        if (dits != NULL) cpl_vector_delete(dits) ;
        if (ndits != NULL) cpl_vector_delete(ndits) ;
        hdrl_imagelist_delete(in) ;
        return -1 ;
    }
    hdrl_imagelist_delete(in) ;
    if (dits != NULL) cpl_vector_delete(dits) ;
    if (ndits != NULL) cpl_vector_delete(ndits) ;

    /* Collapse the image list */
    cpl_msg_info(__func__, "Collapse") ;
    cpl_msg_indent_more() ;
    if (hdrl_imagelist_collapse_mean(in_calib, &collapsed, &contrib) !=
            CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Failed to Collapse") ;
        hdrl_imagelist_delete(in_calib) ;
        cpl_msg_indent_less() ;
        return -1 ;
    }
    cpl_image_delete(contrib) ;
    hdrl_imagelist_delete(in_calib) ;
    cpl_msg_indent_less() ;

    /* Load the trace wave */
    cpl_msg_info(__func__, "Load the TRACE WAVE") ;
    if ((trace_wave = cr2res_io_load_TRACE_WAVE(cpl_frame_get_filename(
                        trace_wave_frame), reduce_det)) == NULL) {
        cpl_msg_error(__func__, "Failed to Load the traces file") ;
        hdrl_image_delete(collapsed) ;
        return -1 ;
    }

	/* Extract at the specified slit fraction */
	if (slit_frac != NULL) {
		if ((trace_wave_new = cr2res_trace_new_slit_fraction(
						trace_wave, slit_frac)) == NULL) {
			cpl_msg_warning(__func__,
	"Failed to compute the traces for user specified slit fraction") ;
			cpl_error_reset() ;
		} else {
			cpl_table_delete(trace_wave) ;
			trace_wave = trace_wave_new ;
			trace_wave_new = NULL ;
		}
	}

    /* Load Blaze */
    blaze_table = NULL ;
    if (blaze_frame != NULL) {
        cpl_msg_info(__func__, "Load the BLAZE") ;
        if ((blaze_table = cr2res_io_load_EXTRACT_1D(cpl_frame_get_filename(
                            blaze_frame), reduce_det)) == NULL) {
            cpl_msg_error(__func__, "Failed to Load the Blaze file") ;
            hdrl_image_delete(collapsed) ;
            cpl_table_delete(trace_wave) ;
            return -1 ;
        }
    }

    /* TODO, make parameters */
    int extract_niter = 10;
    double extract_kappa = 10;

    /* Execute the extraction */
    cpl_msg_info(__func__, "Spectra Extraction") ;
    if (cr2res_extract_traces(collapsed, trace_wave, NULL, blaze_table, -1, -1,
                CR2RES_EXTR_OPT_CURV, extract_height, extract_swath_width, 
                extract_oversample, extract_smooth_slit, extract_smooth_spec,
                extract_niter, extract_kappa, 
                0, 0, 0, &extracted, &slit_func, &model_master) == -1) {
        cpl_msg_error(__func__, "Failed to extract");
        hdrl_image_delete(collapsed) ;
        cpl_table_delete(trace_wave) ;
        return -1 ;
    }
	hdrl_image_delete(collapsed) ;
    if (blaze_table != NULL) cpl_table_delete(blaze_table) ;

    /* Store the extension header for product saving */
    plist = cpl_propertylist_load(first_fname,
            cr2res_io_get_ext_idx(first_fname, reduce_det, 1)) ;

    /* QC - Signal and FWHM */
    qc_signal = cr2res_qc_obs_nodding_signal(extracted) ;
    cpl_propertylist_append_double(plist, CR2RES_HEADER_QC_SIGNAL, qc_signal) ;

    /* QC - SNR */
    qc_snrs = cr2res_qc_snr(trace_wave, extracted, &order_idx_values,
            &nb_order_idx_values) ;
    for (i=0 ; i<nb_order_idx_values ; i++) {
        order_idx = order_idx_values[i] ;
        order_idxp = cr2res_io_convert_order_idx_to_idxp(order_idx) ;
        key_name = cpl_sprintf(CR2RES_HEADER_QC_SNR, order_idxp) ;
        cpl_propertylist_append_double(plist, key_name, qc_snrs[i]) ;
        cpl_free(key_name) ;
    }
    cpl_free(order_idx_values) ;
    cpl_free(qc_snrs) ;

    /* Get the order numbers from the TW rows */
    order_idx_values = cr2res_trace_get_order_idx_values(trace_wave,
            &nb_order_idx_values);

    /* QC - SLIT FWHM */
    fwhm_array = cpl_array_new(nb_order_idx_values, CPL_TYPE_DOUBLE);
    for (i=0 ; i<nb_order_idx_values ; i++) {
        order_idx = order_idx_values[i] ;
        order_idxp = cr2res_io_convert_order_idx_to_idxp(order_idx) ;
        qc_fwhm = cr2res_qc_obs_slit_psf(slit_func, order_idxp,
            extract_oversample);

        key_name = cpl_sprintf(CR2RES_HEADER_QC_SLITFWHM_ORDER, order_idxp) ;
        cpl_propertylist_append_double(plist, key_name, qc_fwhm) ;
        cpl_free(key_name) ;
        cpl_array_set(fwhm_array, i, qc_fwhm) ;
    }
    cpl_free(order_idx_values) ;
    qc_fwhm = cpl_array_get_median(fwhm_array);
    cpl_propertylist_append_double(plist, CR2RES_HEADER_QC_SLITFWHM_MED, 
            qc_fwhm) ;
    if (qc_fwhm < 3.5) {
        cpl_msg_warning(__func__, "Median FWMH of the PSF along the slit "
            "is %gpix, i.e. below the slit width. This means the slit "
            "is likely not evenly filled with light "
            "in the spectral direction. This can result in a "
            "wavelength offset between different postitions along the slit,"
            " and with respect to calibrations."
            , qc_fwhm);
    }

    cpl_array_delete(fwhm_array) ;

    /* QC - Real Orders */
    if (order_zp > 0) {
        /* Get the order numbers from the TW rows */
        order_idx_values = cr2res_trace_get_order_idx_values(trace_wave,
                &nb_order_idx_values);

        /* Compute the Real Order numbers and store them in QCs */
        for (i=0 ; i<nb_order_idx_values ; i++) {
            order_idx = order_idx_values[i] ;
            order_idxp = cr2res_io_convert_order_idx_to_idxp(order_idx) ;
            order_real = cr2res_order_idx_to_real(order_idx, order_zp) ;
            key_name = cpl_sprintf(CR2RES_HEADER_QC_REAL_ORDER, order_idxp) ;
            cpl_propertylist_append_int(plist, key_name, order_real) ;
            cpl_free(key_name) ;
        }
        cpl_free(order_idx_values) ;
    }
	cpl_table_delete(trace_wave) ;

    /* Return */
    *extract = extracted ;
    *slitfunc = slit_func ;
    *model = model_master ;
    *ext_plist = plist ;

    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief  Run basic checks for the rawframes consistency
  @param    rawframes   The input rawframes
  @return   1 if valid, 0 if not, -1 in error case
 */
/*----------------------------------------------------------------------------*/
static int cr2res_obs_staring_check_inputs_validity(
        const cpl_frameset  *   rawframes)
{
    /* Check Inputs */
    if (rawframes == NULL) return -1 ;
    return 1 ;
}
 
/*----------------------------------------------------------------------------*/
/**
  @brief    Get the RAW frames from a frameset
  @param    set     Input frame set
  @return   the RAW frameset or NULL in error case or if it is missing
    Allowed RAW types : CR2RES_OBS_STARING_OTHER_RAW
                        CR2RES_OBS_STARING_JITTER_RAW
                        CR2RES_OBS_STARING_WAVE_SKY_RAW
 */
/*----------------------------------------------------------------------------*/
static cpl_frameset * cr2res_obs_staring_find_RAW(
        const cpl_frameset  *   in)
{
    cpl_frameset    *   out ;

    /* Check entries */
    if (in == NULL) return NULL ;

    out = cr2res_extract_frameset(in, CR2RES_OBS_STARING_OTHER_RAW) ;
    if (out == NULL) {
        out = cr2res_extract_frameset(in, CR2RES_OBS_STARING_JITTER_RAW) ;
    }
    if (out == NULL) {
        out = cr2res_extract_frameset(in, CR2RES_OBS_STARING_WAVE_SKY_RAW) ;
    }
    return out ;
}

