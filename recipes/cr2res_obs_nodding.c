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
#include "cr2res_nodding.h"
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

#define RECIPE_STRING "cr2res_obs_nodding"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_obs_nodding_check_inputs_validity(
        const cpl_frameset  *   rawframes) ;
static int cr2res_obs_nodding_reduce(
        const cpl_frameset  *   rawframes,
        const cpl_frameset  *   raw_flat_frames,
        const cpl_frame     *   trace_wave_frame,
        const cpl_frame     *   detlin_frame,
        const cpl_frame     *   master_dark_frame,
        const cpl_frame     *   master_flat_frame,
        const cpl_frame     *   bpm_frame,
        int                     calib_cosmics_corr,
        int                     extract_oversample,
        int                     extract_swath_width,
        int                     extract_height,
        double                  extract_smooth,
        int                     reduce_det,
        hdrl_image          **  combineda,
        cpl_table           **  extracta,
        cpl_table           **  slitfunca,
        hdrl_image          **  modela,
        hdrl_image          **  combinedb,
        cpl_table           **  extractb,
        cpl_table           **  slitfuncb,
        hdrl_image          **  modelb,
        cpl_propertylist    **  ext_plist) ;

static int cr2res_obs_nodding_create(cpl_plugin *);
static int cr2res_obs_nodding_exec(cpl_plugin *);
static int cr2res_obs_nodding_destroy(cpl_plugin *);
static int cr2res_obs_nodding(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_obs_nodding_description[] = "\
Nodding Observation                                                     \n\
  This recipe handles nodding observations. It expects an even number   \n\
  of rawframes in input, and as many A positions as B positions         \n\
  Inputs                                                                \n\
    raw.fits " CR2RES_OBS_NODDING_RAW" [2 to 2n]                        \n\
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
    cr2res_obs_nodding_extractA.fits " 
    CR2RES_OBS_NODDING_EXTRACTA_PROCATG "\n\
    cr2res_obs_nodding_extractB.fits " 
    CR2RES_OBS_NODDING_EXTRACTB_PROCATG "\n\
    cr2res_obs_nodding_combinedA.fits " 
    CR2RES_OBS_NODDING_COMBINEDA_PROCATG "\n\
    cr2res_obs_nodding_combinedB.fits " 
    CR2RES_OBS_NODDING_COMBINEDB_PROCATG "\n\
    cr2res_obs_nodding_modelA.fits " 
    CR2RES_OBS_NODDING_SLITMODELA_PROCATG "\n\
    cr2res_obs_nodding_modelB.fits " 
    CR2RES_OBS_NODDING_SLITMODELB_PROCATG "\n\
    cr2res_obs_nodding_slitfuncA.fits " 
    CR2RES_OBS_NODDING_SLITFUNCA_PROCATG "\n\
    cr2res_obs_nodding_slitfuncB.fits " 
    CR2RES_OBS_NODDING_SLITFUNCB_PROCATG "\n\
                                                                        \n\
  Algorithm                                                             \n\
    loop on detectors d:                                                \n\
      call cr2res_obs_nodding_reduce()                                  \n\
        -> combined[a|b](d)                                             \n\
        -> extract[a|b](d)                                              \n\
        -> slitfunc[a|b](d)                                             \n\
        -> model[a|b](d)                                                \n\
    Save combineda and combinedb                                        \n\
    Save extracta and extractb                                          \n\
    Save slitfunca and slitfuncb                                        \n\
    Save modela and modelb                                              \n\
                                                                        \n\
    cr2res_obs_nodding_reduce()                                         \n\
      Load the input raw frames in an image list                        \n\
      Apply the calibrations to the image list                          \n\
      Split the list in listA and listB image lists                     \n\
      Compute diffA=listA-listB and diffB=listB-listA                   \n\
      Collapse diffA and diffB                                          \n\
        -> combined[a|b]                                                \n\
      Load the input trace wave                                         \n\
      Compute the slit fractions for A and B by using the nodthrow      \n\
              and the assumption that A and B are at equal distances    \n\
              from the slit center.                                     \n\
      Compute 2 new trace_wave files  with these 2 computed slit fractions\n\
      Extract the spectra for A and for B                               \n\
        -> extracted[a|b]                                               \n\
        -> slit_func[a|b]                                               \n\
        -> model_master[a|b]                                            \n\
      Compute QC parameters                                             \n\
                                                                        \n\
  Library functions uѕed                                                \n\
    cr2res_io_find_TRACE_WAVE()                                         \n\
    cr2res_io_find_BPM()                                                \n\
    cr2res_obs_nodding_reduce()                                         \n\
    cr2res_nodding_read_positions()                                     \n\
    cr2res_io_read_dits()                                               \n\
    cr2res_io_load_image_list_from_set()                                \n\
    cr2res_calib_imagelist()                                            \n\
    cr2res_combine_nodding_split()                                      \n\
    cr2res_io_load_TRACE_WAVE()                                         \n\
    cr2res_pfits_get_nodthrow()                                         \n\
    cr2res_trace_new_slit_fraction()                                    \n\
    cr2res_extract_traces()                                             \n\
    cr2res_qc_obs_nodding_signal()                                      \n\
    cr2res_qc_obs_nodding_transmission()                                \n\
    cr2res_io_save_COMBINED()                                           \n\
    cr2res_io_save_EXTRACT_1D()                                         \n\
    cr2res_io_save_SLIT_FUNC()                                          \n\
    cr2res_io_save_SLIT_MODEL()                                         \n\
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
                    "Nodding Observation recipe",
                    cr2res_obs_nodding_description,
                    "Thomas Marquart, Yves Jung",
                    PACKAGE_BUGREPORT,
                    cr2res_get_license(),
                    cr2res_obs_nodding_create,
                    cr2res_obs_nodding_exec,
                    cr2res_obs_nodding_destroy)) {    
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
static int cr2res_obs_nodding_create(cpl_plugin * plugin)
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
    p = cpl_parameter_new_value("cr2res.cr2res_obs_nodding.extract_oversample",
            CPL_TYPE_INT, "factor by which to oversample the extraction",
            "cr2res.cr2res_obs_nodding", 2);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_oversample");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_nodding.extract_swath_width",
            CPL_TYPE_INT, "The swath width", "cr2res.cr2res_obs_nodding", 24);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_swath_width");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_nodding.extract_height",
            CPL_TYPE_INT, "Extraction height",
            "cr2res.cr2res_obs_nodding", -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_height");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_nodding.extract_smooth",
            CPL_TYPE_DOUBLE,
            "Smoothing along the slit (1 for high S/N, 5 for low)",
            "cr2res.cr2res_obs_nodding", 1.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_smooth");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_nodding.detector",
            CPL_TYPE_INT, "Only reduce the specified detector",
            "cr2res.cr2res_obs_nodding", 0);
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
static int cr2res_obs_nodding_exec(cpl_plugin * plugin)
{
    cpl_recipe  *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1;

    return cr2res_obs_nodding(recipe->frames, recipe->parameters);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_obs_nodding_destroy(cpl_plugin * plugin)
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
static int cr2res_obs_nodding(
        cpl_frameset            *   frameset,
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter *   param ;
    int                     extract_oversample, extract_swath_width,
                            extract_height, reduce_det ;
    double                  extract_smooth ;
    cpl_frameset        *   rawframes ;
    cpl_frameset        *   raw_flat_frames ;
    const cpl_frame     *   trace_wave_frame ;
    const cpl_frame     *   detlin_frame ;
    const cpl_frame     *   master_dark_frame ;
    const cpl_frame     *   master_flat_frame ;
    const cpl_frame     *   bpm_frame ;
    hdrl_image          *   combineda[CR2RES_NB_DETECTORS] ;
    cpl_table           *   extracta[CR2RES_NB_DETECTORS] ;
    cpl_table           *   slitfunca[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   modela[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   combinedb[CR2RES_NB_DETECTORS] ;
    cpl_table           *   extractb[CR2RES_NB_DETECTORS] ;
    cpl_table           *   slitfuncb[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   modelb[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   ext_plist[CR2RES_NB_DETECTORS] ;
    char                *   out_file;
    int                     i, det_nr; 

    /* Initialise */

    /* RETRIEVE INPUT PARAMETERS */
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_nodding.extract_oversample");
    extract_oversample = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_nodding.extract_swath_width");
    extract_swath_width = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_nodding.extract_height");
    extract_height = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_nodding.extract_smooth");
    extract_smooth = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_nodding.detector");
    reduce_det = cpl_parameter_get_int(param);

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

    /* Get the RAW Frames */
    rawframes = cr2res_extract_frameset(frameset, CR2RES_OBS_NODDING_RAW) ;
    if (rawframes == NULL) {
        cpl_msg_error(__func__, "Could not find RAW frames") ;
        return -1 ;
    }
   
    /* Get the RAW flat frames */
    raw_flat_frames = cr2res_extract_frameset(frameset, CR2RES_FLAT_RAW) ;

    /* Loop on the detectors */
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
        /* Initialise */
        combineda[det_nr-1] = NULL ;
        extracta[det_nr-1] = NULL ;
        slitfunca[det_nr-1] = NULL ;
        modela[det_nr-1] = NULL ;
        combinedb[det_nr-1] = NULL ;
        extractb[det_nr-1] = NULL ;
        slitfuncb[det_nr-1] = NULL ;
        modelb[det_nr-1] = NULL ;
        ext_plist[det_nr-1] = NULL ;

        /* Compute only one detector */
        if (reduce_det != 0 && det_nr != reduce_det) continue ;
    
        cpl_msg_info(__func__, "Process Detector %d", det_nr) ;
        cpl_msg_indent_more() ;

        /* Call the reduction function */
        if (cr2res_obs_nodding_reduce(rawframes, raw_flat_frames, 
                    trace_wave_frame, detlin_frame, master_dark_frame, 
                    master_flat_frame, bpm_frame, 0, extract_oversample, 
                    extract_swath_width, extract_height, extract_smooth, det_nr,
                    &(combineda[det_nr-1]),
                    &(extracta[det_nr-1]),
                    &(slitfunca[det_nr-1]),
                    &(modela[det_nr-1]),
                    &(combinedb[det_nr-1]),
                    &(extractb[det_nr-1]),
                    &(slitfuncb[det_nr-1]),
                    &(modelb[det_nr-1]),
                    &(ext_plist[det_nr-1])) == -1) {
            cpl_msg_warning(__func__, "Failed to reduce detector %d", det_nr);
        }
        cpl_msg_indent_less() ;
    }

    /* Ѕave Products */
    out_file = cpl_sprintf("%s_combinedA.fits", RECIPE_STRING) ;
    cr2res_io_save_COMBINED(out_file, frameset, rawframes, parlist,
            combineda, NULL, ext_plist, CR2RES_OBS_NODDING_COMBINEDA_PROCATG, 
            RECIPE_STRING) ;
    cpl_free(out_file);

    out_file = cpl_sprintf("%s_extractedA.fits", RECIPE_STRING) ;
    cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, extracta,
            NULL, ext_plist, CR2RES_OBS_NODDING_EXTRACTA_PROCATG,RECIPE_STRING);
    cpl_free(out_file);

    out_file = cpl_sprintf("%s_slitfuncA.fits", RECIPE_STRING) ;
    cr2res_io_save_SLIT_FUNC(out_file, frameset, rawframes, parlist,
            slitfunca, NULL, ext_plist, CR2RES_OBS_NODDING_SLITFUNCA_PROCATG,
            RECIPE_STRING) ;
    cpl_free(out_file);

    out_file = cpl_sprintf("%s_modelA.fits", RECIPE_STRING) ;
    cr2res_io_save_SLIT_MODEL(out_file, frameset, rawframes, parlist,
            modela, NULL, ext_plist, CR2RES_OBS_NODDING_SLITMODELA_PROCATG,
            RECIPE_STRING) ;
    cpl_free(out_file);

    out_file = cpl_sprintf("%s_combinedB.fits", RECIPE_STRING) ;
    cr2res_io_save_COMBINED(out_file, frameset, rawframes, parlist,
            combinedb, NULL, ext_plist, CR2RES_OBS_NODDING_COMBINEDB_PROCATG, 
            RECIPE_STRING) ;
    cpl_free(out_file);

    out_file = cpl_sprintf("%s_extractedB.fits", RECIPE_STRING) ;
    cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, extractb,
            NULL, ext_plist, CR2RES_OBS_NODDING_EXTRACTB_PROCATG,RECIPE_STRING);
    cpl_free(out_file);

    out_file = cpl_sprintf("%s_slitfuncB.fits", RECIPE_STRING) ;
    cr2res_io_save_SLIT_FUNC(out_file, frameset, rawframes, parlist,
            slitfuncb, NULL, ext_plist, CR2RES_OBS_NODDING_SLITFUNCB_PROCATG,
            RECIPE_STRING) ;
    cpl_free(out_file);

    out_file = cpl_sprintf("%s_modelB.fits", RECIPE_STRING) ;
    cr2res_io_save_SLIT_MODEL(out_file, frameset, rawframes, parlist,
            modelb, NULL, ext_plist, CR2RES_OBS_NODDING_SLITMODELB_PROCATG,
            RECIPE_STRING) ;
    cpl_free(out_file);

    /* Free */
    cpl_frameset_delete(rawframes) ;
    if (raw_flat_frames != NULL) cpl_frameset_delete(raw_flat_frames) ;
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
        if (combineda[det_nr-1] != NULL)
            hdrl_image_delete(combineda[det_nr-1]) ;
        if (extracta[det_nr-1] != NULL) 
            cpl_table_delete(extracta[det_nr-1]) ;
        if (slitfunca[det_nr-1] != NULL) 
            cpl_table_delete(slitfunca[det_nr-1]) ;
        if (modela[det_nr-1] != NULL)
            hdrl_image_delete(modela[det_nr-1]) ;
        if (combinedb[det_nr-1] != NULL)
            hdrl_image_delete(combinedb[det_nr-1]) ;
        if (extractb[det_nr-1] != NULL) 
            cpl_table_delete(extractb[det_nr-1]) ;
        if (slitfuncb[det_nr-1] != NULL) 
            cpl_table_delete(slitfuncb[det_nr-1]) ;
        if (modelb[det_nr-1] != NULL)
            hdrl_image_delete(modelb[det_nr-1]) ;
        if (ext_plist[det_nr-1] != NULL) 
            cpl_propertylist_delete(ext_plist[det_nr-1]) ;
    }

    return (int)cpl_error_get_code();
}
 
/*----------------------------------------------------------------------------*/
/**
  @brief    Execute the science recipe on a specific detector
  @param rawframes              Raw science frames
  @param raw_flat_frames        Raw flat frames 
  @param trace_wave_frame       Trace Wave file
  @param detlin_frame           Associated detlin coefficients
  @param master_dark_frame      Associated master dark
  @param master_flat_frame      Associated master flat
  @param bpm_frame              Associated BPM
  @param calib_cosmics_corr     Flag to correct for cosmics
  @param extract_oversample     Extraction related
  @param extract_swath_width    Extraction related
  @param extract_height         Extraction related
  @param extract_smooth         Extraction related
  @param reduce_det             The detector to compute
  @param combineda              [out] Combined image (A)
  @param extracta               [out] extracted spectrum (A)
  @param slitfunca              [out] slit function (A)
  @param modela                 [out] slit model (A)
  @param combinedb              [out] Combined image (B)
  @param extractb               [out] extracted spectrum (B)
  @param slitfuncb              [out] slit function (B)
  @param modelb                 [out] slit model (B)
  @param ext_plist              [out] the header for saving the products
  @return  0 if ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
static int cr2res_obs_nodding_reduce(
        const cpl_frameset  *   rawframes,
        const cpl_frameset  *   raw_flat_frames,
        const cpl_frame     *   trace_wave_frame,
        const cpl_frame     *   detlin_frame,
        const cpl_frame     *   master_dark_frame,
        const cpl_frame     *   master_flat_frame,
        const cpl_frame     *   bpm_frame,
        int                     calib_cosmics_corr,
        int                     extract_oversample,
        int                     extract_swath_width,
        int                     extract_height,
        double                  extract_smooth,
        int                     reduce_det,
        hdrl_image          **  combineda,
        cpl_table           **  extracta,
        cpl_table           **  slitfunca,
        hdrl_image          **  modela,
        hdrl_image          **  combinedb,
        cpl_table           **  extractb,
        cpl_table           **  slitfuncb,
        hdrl_image          **  modelb,
        cpl_propertylist    **  ext_plist)
{
    hdrl_imagelist      *   in ;
    hdrl_imagelist      *   in_calib ;
    hdrl_imagelist      *   in_a ;
    hdrl_imagelist      *   in_b ;
    hdrl_imagelist      *   diff_a ;
    hdrl_imagelist      *   diff_b ;
    hdrl_image          *   collapsed_a ;
    hdrl_image          *   collapsed_b ;
    cpl_image           *   contrib_a ;
    cpl_image           *   contrib_b ;
    cr2res_nodding_pos  *   nod_positions ;
    cpl_vector          *   dits ;
    cpl_table           *   trace_wave ;
    cpl_table           *   trace_wave_corrected ;
    cpl_table           *   trace_wave_a ;
    cpl_table           *   trace_wave_b ;
    cpl_array           *   slit_frac_a ;
    cpl_array           *   slit_frac_b ;
    cpl_table           *   extracted_a ;
    cpl_table           *   extracted_b ;
    cpl_table           *   slit_func_a ;
    cpl_table           *   slit_func_b ;
    hdrl_image          *   model_master_a ;
    hdrl_image          *   model_master_b ;
    cpl_propertylist    *   plist ;
    cpl_size                nframes, i ;
    double                  slit_length, extr_width_frac, slit_frac_a_bot, 
                            slit_frac_a_mid, slit_frac_a_top, slit_frac_b_bot, 
                            slit_frac_b_mid, slit_frac_b_top, nod_throw ;
    double                  qc_signal_a, qc_signal_b, qc_transm, qc_fwhm_a, qc_fwhm_b ;
    int                     det_nr ;

    /* Check Inputs */
    if (combineda == NULL || combinedb == NULL || extracta == NULL ||
            extractb == NULL || ext_plist == NULL || rawframes == NULL
            || trace_wave_frame == NULL) return -1 ;

    /* Check raw frames consistency */
    if (cr2res_obs_nodding_check_inputs_validity(rawframes) != 1) {
        cpl_msg_error(__func__, "Invalid Inputs") ;
        return -1 ;
    }

    /* Initialise */
    nframes = cpl_frameset_get_size(rawframes) ;

    /* Get the Nodding positions */
    nod_positions = cr2res_nodding_read_positions(rawframes) ;

    if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
        for (i=0 ; i<nframes ; i++) {
            cpl_msg_debug(__func__, "Frame %s - Nodding %c", 
                    cpl_frame_get_filename(
                        cpl_frameset_get_position_const(rawframes, i)), 
                cr2res_nodding_position_char(nod_positions[i])) ;
        }
    }

    /* Load the DITs if necessary */
    if (master_dark_frame != NULL)  dits = cr2res_io_read_dits(rawframes) ;
    else                            dits = NULL ;
    if (cpl_msg_get_level() == CPL_MSG_DEBUG && dits != NULL) 
        cpl_vector_dump(dits, stdout) ;

    /* Load image list */
    if ((in = cr2res_io_load_image_list_from_set(rawframes, 
                    reduce_det)) == NULL) {
        cpl_msg_error(__func__, "Cannot load images") ;
        cpl_free(nod_positions) ;    
        if (dits != NULL) cpl_vector_delete(dits) ;
        return -1 ;
    }
    if (hdrl_imagelist_get_size(in) != cpl_frameset_get_size(rawframes)) {
        cpl_msg_error(__func__, "Inconsistent number of loaded images") ;
        cpl_free(nod_positions) ;    
        if (dits != NULL) cpl_vector_delete(dits) ;
        hdrl_imagelist_delete(in) ;
        return -1 ;
    }

    /* Calibrate the images */
    if ((in_calib = cr2res_calib_imagelist(in, reduce_det, 0, master_flat_frame,
            master_dark_frame, bpm_frame, detlin_frame, dits)) == NULL) {
        cpl_msg_error(__func__, "Failed to apply the calibrations") ;
        cpl_free(nod_positions) ;    
        if (dits != NULL) cpl_vector_delete(dits) ;
        hdrl_imagelist_delete(in) ;
        return -1 ;
    }
    hdrl_imagelist_delete(in) ;
    if (dits != NULL) cpl_vector_delete(dits) ;

    /* Split the image lists */
    if (cr2res_combine_nodding_split(in_calib, nod_positions, &in_a, 
            &in_b)) {
        cpl_msg_error(__func__, "Failed to split the nodding positions") ;
        cpl_free(nod_positions) ;    
        hdrl_imagelist_delete(in_calib) ;
        return -1 ;
    }
    cpl_free(nod_positions) ;    
    hdrl_imagelist_delete(in_calib) ;

    /* Check the sizes of A/B image lists */
    if (hdrl_imagelist_get_size(in_a) != hdrl_imagelist_get_size(in_b)
            || hdrl_imagelist_get_size(in_a) == 0) {
        cpl_msg_error(__func__, "Іnconsistent A / B number of images") ;
        hdrl_imagelist_delete(in_a) ;
        hdrl_imagelist_delete(in_b) ;
        return -1 ;
    }

    /* Compute diff_a = in_a - in_b and diff_b = in_b - in_a */
    diff_a = hdrl_imagelist_duplicate(in_a) ;
    hdrl_imagelist_sub_imagelist(diff_a, in_b) ;
    diff_b = hdrl_imagelist_duplicate(in_b) ;
    hdrl_imagelist_sub_imagelist(diff_b, in_a) ;

    hdrl_imagelist_delete(in_a) ;
    hdrl_imagelist_delete(in_b) ;
    
    /* Collapse A-B and B-A */
    cpl_msg_info(__func__, "Collapse A-B and B-A") ;
    cpl_msg_indent_more() ;
    if (hdrl_imagelist_collapse_mean(diff_a, &collapsed_a, &contrib_a) !=
            CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Failed to Collapse A-B") ;
        hdrl_imagelist_delete(diff_a) ;
        hdrl_imagelist_delete(diff_b) ;
        cpl_msg_indent_less() ;
        return -1 ;
    }
    cpl_image_delete(contrib_a) ;
    hdrl_imagelist_delete(diff_a) ;
    if (hdrl_imagelist_collapse_mean(diff_b, &collapsed_b, &contrib_b) !=
            CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Failed to Collapse B-A") ;
        hdrl_imagelist_delete(diff_b) ;
        hdrl_image_delete(collapsed_a) ;
        cpl_msg_indent_less() ;
        return -1 ;
    }
    cpl_image_delete(contrib_b) ;
    hdrl_imagelist_delete(diff_b) ;
    cpl_msg_indent_less() ;

    /* Load the trace wave */
    cpl_msg_info(__func__, "Load the TRACE WAVE") ;
    if ((trace_wave = cr2res_io_load_TRACE_WAVE(cpl_frame_get_filename(
                        trace_wave_frame), reduce_det)) == NULL) {
        cpl_msg_error(__func__, "Failed to Load the traces file") ;
        hdrl_image_delete(collapsed_a) ;
        hdrl_image_delete(collapsed_b) ;
        return -1 ;
    }

    /* Correct trace_wave with some provided raw flats */
    if (raw_flat_frames != NULL) {
        cpl_msg_info(__func__, "Try to correct the reproducibility error") ;
        cpl_msg_indent_more() ;
        trace_wave_corrected = cr2res_trace_adjust(trace_wave, raw_flat_frames, 
                reduce_det) ;
        if (trace_wave_corrected != NULL) {
            cpl_table_delete(trace_wave) ;
            trace_wave = trace_wave_corrected ;
            trace_wave_corrected = NULL ;
        }
        cpl_msg_indent_less() ;
    }

    /* Compute the slit fractions for A and B positions extraction */   
    /*
    The assumption is made here that :  
        - The slit center is exactly in the middle of A and B poѕitions
        - The nodthrow iѕ the distance in arcseconds between A and B
        - The slit size is 10 arcseconds
        - The A position is above the B position
    */
    slit_length = 10 ;
    if ((plist = cpl_propertylist_load(cpl_frame_get_filename(
                        cpl_frameset_get_position_const(rawframes, 0)),
                    0)) == NULL) {
        cpl_msg_error(__func__, "Cannot read the NODTHROW in the input files") ;
        hdrl_image_delete(collapsed_a) ;
        hdrl_image_delete(collapsed_b) ;
        cpl_table_delete(trace_wave) ;
        return -1 ;
    }
    nod_throw = cr2res_pfits_get_nodthrow(plist) ;
    cpl_propertylist_delete(plist) ;
    if (nod_throw < slit_length / 2.0) {
        extr_width_frac = nod_throw/slit_length ;  
        slit_frac_a_bot = 0.5 ;
        slit_frac_a_mid = slit_frac_a_bot + extr_width_frac/2.0 ;
        slit_frac_a_top = slit_frac_a_bot + extr_width_frac ;
        slit_frac_b_top = 0.5 ;
        slit_frac_b_mid = slit_frac_b_top - extr_width_frac/2.0 ;
        slit_frac_b_bot = slit_frac_b_top - extr_width_frac ;
    } else if (nod_throw <= slit_length) {
        extr_width_frac = (slit_length - nod_throw)/slit_length ;  
        slit_frac_a_top = 1.0 ;
        slit_frac_a_mid = slit_frac_a_top - extr_width_frac/2.0 ;
        slit_frac_a_bot = slit_frac_a_top - extr_width_frac ;
        slit_frac_b_bot = 0.0 ;
        slit_frac_b_mid = slit_frac_b_bot + extr_width_frac/2.0 ;
        slit_frac_b_top = slit_frac_b_bot + extr_width_frac ;
    } else {
        cpl_msg_error(__func__, "NODTHROW > slit length (%g>%g)- abort", 
                nod_throw, slit_length) ;
        hdrl_image_delete(collapsed_a) ;
        hdrl_image_delete(collapsed_b) ;
        cpl_table_delete(trace_wave) ;
        return -1 ;
    }
    cpl_msg_info(__func__, "Nodding A extraction: Slit fraction %g - %g",
            slit_frac_a_bot, slit_frac_a_top) ;
    cpl_msg_info(__func__, "Nodding B extraction: Slit fraction %g - %g",
            slit_frac_b_bot, slit_frac_b_top) ;

    slit_frac_a = cpl_array_new(3, CPL_TYPE_DOUBLE) ;
    cpl_array_set(slit_frac_a, 0, slit_frac_a_bot) ;
    cpl_array_set(slit_frac_a, 1, slit_frac_a_mid) ;
    cpl_array_set(slit_frac_a, 2, slit_frac_a_top) ;

    slit_frac_b = cpl_array_new(3, CPL_TYPE_DOUBLE) ;
    cpl_array_set(slit_frac_b, 0, slit_frac_b_bot) ;
    cpl_array_set(slit_frac_b, 1, slit_frac_b_mid) ;
    cpl_array_set(slit_frac_b, 2, slit_frac_b_top) ;

    /* Recompute the traces for the new slit fractions */
    if ((trace_wave_a = cr2res_trace_new_slit_fraction(trace_wave,
            slit_frac_a)) == NULL) {
        cpl_msg_error(__func__, 
                "Failed to compute the traces for extraction of A") ;
        hdrl_image_delete(collapsed_a) ;
        hdrl_image_delete(collapsed_b) ;
        cpl_table_delete(trace_wave) ;
        cpl_array_delete(slit_frac_a) ;
        cpl_array_delete(slit_frac_b) ;
        return -1 ;
    }

    cpl_array_delete(slit_frac_a) ;
    if ((trace_wave_b = cr2res_trace_new_slit_fraction(trace_wave,
            slit_frac_b)) == NULL) {
        cpl_msg_error(__func__, 
                "Failed to compute the traces for extraction of B") ;
        hdrl_image_delete(collapsed_a) ;
        hdrl_image_delete(collapsed_b) ;
        cpl_table_delete(trace_wave) ;
        cpl_table_delete(trace_wave_a) ;
        cpl_array_delete(slit_frac_b) ;
        return -1 ;
    }
    cpl_array_delete(slit_frac_b) ;
    cpl_table_delete(trace_wave) ;

    /* Execute the extraction */
    cpl_msg_info(__func__, "Spectra Extraction") ;
    if (cr2res_extract_traces(collapsed_a, trace_wave_a, -1, -1,
                CR2RES_EXTR_OPT_CURV, extract_height, extract_swath_width, 
                extract_oversample, extract_smooth,
                &extracted_a, &slit_func_a, &model_master_a) == -1) {
        cpl_msg_error(__func__, "Failed to extract A");
        hdrl_image_delete(collapsed_a) ;
        hdrl_image_delete(collapsed_b) ;
        cpl_table_delete(trace_wave_a) ;
        cpl_table_delete(trace_wave_b) ;
        return -1 ;
    }
/* TODO : Save trace_wave_a and b as products */
    cpl_table_delete(trace_wave_a) ;
    if (cr2res_extract_traces(collapsed_b, trace_wave_b, -1, -1,
                CR2RES_EXTR_OPT_CURV, extract_height, extract_swath_width, 
                extract_oversample, extract_smooth,
                &extracted_b, &slit_func_b, &model_master_b) == -1) {
        cpl_msg_error(__func__, "Failed to extract B");
        cpl_table_delete(extracted_a) ;
        cpl_table_delete(slit_func_a) ;
        hdrl_image_delete(model_master_a) ;
        hdrl_image_delete(collapsed_a) ;
        hdrl_image_delete(collapsed_b) ;
        cpl_table_delete(trace_wave_b) ;
        return -1 ;
    }
    cpl_table_delete(trace_wave_b) ;

    /* QC parameters */
    plist = cpl_propertylist_new() ;

    /* Compute the QC parameters */
    /* TODO : pass the proper inputs */
    qc_signal_a = cr2res_qc_obs_nodding_signal(extracted_a) ;
    qc_signal_b = cr2res_qc_obs_nodding_signal(extracted_b) ;
    qc_transm = cr2res_qc_obs_nodding_transmission(NULL) ;
    qc_fwhm_a = cr2res_qc_obs_nodding_slit_psf(slit_func_a);
    qc_fwhm_b = cr2res_qc_obs_nodding_slit_psf(slit_func_b);

    /* Store the QC parameters in the plist */
    cpl_propertylist_append_double(plist, CR2RES_HEADER_QC_SIGNAL, 
            (qc_signal_a+qc_signal_b)/2.0) ;
    cpl_propertylist_append_double(plist, CR2RES_HEADER_QC_TRANSM,
            qc_transm) ;
    cpl_propertylist_append_double(plist, CR2RES_HEADER_QC_SLITFWHM,
            (qc_fwhm_a+qc_fwhm_b)/2.0) ;

    /* Return */
    *combineda = collapsed_a ;
    *extracta = extracted_a ;
    *slitfunca = slit_func_a ;
    *modela = model_master_a ;

    *combinedb = collapsed_b ;
    *extractb = extracted_b ;
    *slitfuncb = slit_func_b ;
    *modelb = model_master_b ;

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
static int cr2res_obs_nodding_check_inputs_validity(
        const cpl_frameset  *   rawframes)
{
    cpl_propertylist        *   plist ;
    cr2res_nodding_pos      *   nod_positions ;
    cpl_size                    nframes, i ;
    double                      nodthrow, nodthrow_cur ;
    int                         nb_a, nb_b ;

    /* Check Inputs */
    if (rawframes == NULL) return -1 ;

    /* Need even number of frames */
    nframes = cpl_frameset_get_size(rawframes) ;
    if (nframes % 2) {
        cpl_msg_error(__func__, "Require an even number of raw frames") ;
        return 0 ;
    }

    /* Need same number of A and B positions */
    nb_a = nb_b = 0 ;
    nod_positions = cr2res_nodding_read_positions(rawframes) ;
    if (nod_positions == NULL) return -1 ;
    for (i=0 ; i<nframes ; i++) {
        if (nod_positions[i] == CR2RES_NODDING_A) nb_a++ ;
        if (nod_positions[i] == CR2RES_NODDING_B) nb_b++ ;
    }
    cpl_free(nod_positions) ;    

    if (nb_a == 0 || nb_b == 0 || nb_a != nb_b) {
        cpl_msg_error(__func__, "Require as many A and B positions") ;
        return 0 ;
    }

    /* Need the same nod throw in all frames */
    if ((plist = cpl_propertylist_load(cpl_frame_get_filename(
                        cpl_frameset_get_position_const(rawframes, 0)),
                    0)) == NULL) {
        return -1;
    } 
    nodthrow = cr2res_pfits_get_nodthrow(plist);
    cpl_propertylist_delete(plist) ;
    for (i=1 ; i<nframes ; i++) {
        if ((plist = cpl_propertylist_load(cpl_frame_get_filename(
                            cpl_frameset_get_position_const(rawframes, i)),
                        0)) == NULL) {
            return -1;
        } 
        nodthrow_cur = cr2res_pfits_get_nodthrow(plist);
        cpl_propertylist_delete(plist) ;

        if (fabs(nodthrow_cur-nodthrow) > 1e-3) {
            cpl_msg_error(__func__, 
                    "Require constant NOD THROW in the raw frames") ;
            return 0 ;
        }
    }
    return 1 ;
}

