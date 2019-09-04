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

#include <string.h>
#include <cpl.h>

#include "cr2res_utils.h"
#include "cr2res_pol.h"
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

#define RECIPE_STRING "cr2res_obs_pol"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int * cr2res_obs_pol_get_order_numbers(
        const cpl_table **  extracted,
        int                 next,
        int             *   norders) ;
static int cr2res_obs_pol_check_inputs_validity(
        const cpl_frameset  *   rawframes) ;
static int cr2res_obs_pol_reduce(
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
        cpl_table           **  pol_spec_a,
        cpl_table           **  pol_spec_b,
        cpl_propertylist    **  ext_plista,
        cpl_propertylist    **  ext_plistb) ;
static int cr2res_obs_pol_reduce_one(
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
        cpl_table           **  pol_spec,
        cpl_propertylist    **  ext_plist) ;
static int cr2res_obs_pol_create(cpl_plugin *);
static int cr2res_obs_pol_exec(cpl_plugin *);
static int cr2res_obs_pol_destroy(cpl_plugin *);
static int cr2res_obs_pol(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_obs_pol_description[] = "\
Polarimetry Observation                                                 \n\
  The input raw frames are separated in A and B nodding positions. A    \n\
  and B frames are reduced separately. The frames are grouped by blocks \n\
  of 4 frames. Each block generateѕ 8 extractions. For each order       \n\
  found, 8 spectra are passed to the demodulation functions.            \n\
  The results are stored in polarimetry tables (1 per group of 4        \n\
  frames). The polarimetry tables are then merged together if there are \n\
  several group of 4 frames available.                                  \n\
                                                                        \n\
  Inputs                                                                \n\
    raw.fits " CR2RES_OBS_POL_RAW " [4 to 4n]                           \n\
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
    cr2res_obs_pol_specA.fits " CR2RES_OBS_POL_SPECA_PROCATG "          \n\
    cr2res_obs_pol_specB.fits " CR2RES_OBS_POL_SPECB_PROCATG "          \n\
                                                                        \n\
  Algorithm                                                             \n\
    loop on detectors d:                                                \n\
      cr2res_obs_pol_reduce()                                           \n\
        -> pol_speca(d)                                                 \n\
        -> pol_specb(d)                                                 \n\
    Save pol_speca                                                      \n\
    Save pol_specb                                                      \n\
                                                                        \n\
    cr2res_obs_pol_reduce():                                            \n\
      Read the nodding positions of the raw frames                      \n\
      Split the rawframes in A and B positions                          \n\
      Call cr2res_obs_pol_reduce_one() successively on A and B          \n\
        -> pol_speca                                                    \n\
        -> pol_specb                                                    \n\
                                                                        \n\
    cr2res_obs_pol_reduce_one():                                        \n\
      Read the DITs                                                     \n\
      Read the decker positions                                         \n\
      Load the image list                                               \n\
      Apply the calibrations to the image list                          \n\
      Load the trace_wave                                               \n\
      For each group of 4 images g:                                     \n\
        Compute 8 extracted tables (2 per image):                       \n\
            1u, 1d, 2u, 2d, 3u, 3d, 4u, 4d                              \n\
            where u/d are for up and down                               \n\
                  1->4 iѕ derived with cr2res_pol_sort_frames()         \n\
                  The decker info is used to derive the 8 slit fractions\n\
        Count norders the number of different orders in those 8 tables  \n\
            [Note : 1 extracted table has 1 spectrum per order]         \n\
        loop on orders o:                                               \n\
          Get the 8 spectra/wl/error for this order from                \n\
                1u, 1d, 2u, 2d, 3u, 3d, 4u, 4d                          \n\
          Call the cr2res_pol_demod_stokes()                            \n\
            -> demod_stokes(o, g)                                       \n\
          Call cr2res_pol_demod_null()                                  \n\
            -> demod_null(o, g)                                         \n\
          Call cr2res_pol_demod_intens()                                \n\
            -> demod_intenѕ(o, g)                                       \n\
        Create pol_spec(g) from demod_stokes(g),                        \n\
                                demod_null(g),                          \n\
                                demod_intenѕ(g)                         \n\
      Merge pol_spec(g) into pol_spec                                   \n\
                                                                        \n\
  Library functions uѕed                                                \n\
    cr2res_io_find_TRACE_WAVE()                                         \n\
    cr2res_io_find_BPM()                                                \n\
    cr2res_obs_pol_reduce()                                             \n\
    cr2res_nodding_read_positions()                                     \n\
    cr2res_combine_nodding_split_frames()                               \n\
    cr2res_obs_pol_reduce_one()                                         \n\
    cr2res_io_read_dits()                                               \n\
    cr2res_io_read_decker_positions()                                   \n\
    cr2res_io_load_image_list_from_set()                                \n\
    cr2res_calib_imagelist()                                            \n\
    cr2res_io_load_TRACE_WAVE()                                         \n\
    cr2res_pol_sort_frames()                                            \n\
    cr2res_trace_slit_fraction_create()                                 \n\
    cr2res_trace_new_slit_fraction()                                    \n\
    cr2res_extract_traces()                                             \n\
    cr2res_obs_pol_get_order_numbers()                                  \n\
    cr2res_pol_demod_stokes()                                           \n\
    cr2res_pol_demod_null()                                             \n\
    cr2res_pol_demod_intens()                                           \n\
    cr2res_pol_POL_SPEC_create()                                        \n\
    cr2res_pol_spec_pol_merge()                                         \n\
    cr2res_io_save_POL_SPEC()                                           \n\
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
                    "Polarimetry Observation recipe",
                    cr2res_obs_pol_description,
                    "Thomas Marquart, Yves Jung",
                    PACKAGE_BUGREPORT,
                    cr2res_get_license(),
                    cr2res_obs_pol_create,
                    cr2res_obs_pol_exec,
                    cr2res_obs_pol_destroy)) {    
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
static int cr2res_obs_pol_create(cpl_plugin * plugin)
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
    p = cpl_parameter_new_value("cr2res.cr2res_obs_pol.extract_oversample",
            CPL_TYPE_INT, "factor by which to oversample the extraction",
            "cr2res.cr2res_obs_pol", 2);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_oversample");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_pol.extract_swath_width",
            CPL_TYPE_INT, "The swath width", "cr2res.cr2res_obs_pol", 90);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_swath_width");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_pol.extract_height",
            CPL_TYPE_INT, "Extraction height",
            "cr2res.cr2res_obs_pol", -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_height");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_pol.extract_smooth",
            CPL_TYPE_DOUBLE,
            "Smoothing along the slit (1 for high S/N, 5 for low)",
            "cr2res.cr2res_obs_pol", 1.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_smooth");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_pol.detector",
            CPL_TYPE_INT, "Only reduce the specified detector",
            "cr2res.cr2res_obs_pol", 0);
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
static int cr2res_obs_pol_exec(cpl_plugin * plugin)
{
    cpl_recipe  *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1;

    return cr2res_obs_pol(recipe->frames, recipe->parameters);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_obs_pol_destroy(cpl_plugin * plugin)
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
static int cr2res_obs_pol(
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
    cpl_table           *   pol_speca[CR2RES_NB_DETECTORS] ;
    cpl_table           *   pol_specb[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   ext_plista[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   ext_plistb[CR2RES_NB_DETECTORS] ;
    char                *   out_file;
    int                     i, det_nr; 

    /* Initialise */

    /* RETRIEVE INPUT PARAMETERS */
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_pol.extract_oversample");
    extract_oversample = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_pol.extract_swath_width");
    extract_swath_width = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_pol.extract_height");
    extract_height = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_pol.extract_smooth");
    extract_smooth = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_pol.detector");
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
    rawframes = cr2res_extract_frameset(frameset, CR2RES_OBS_POL_RAW) ;
    if (rawframes == NULL) {
        cpl_msg_error(__func__, "Could not find RAW frames") ;
        return -1 ;
    }
   
    /* Get the RAW flat frames */
    raw_flat_frames = cr2res_extract_frameset(frameset, CR2RES_FLAT_RAW) ;

    /* Loop on the detectors */
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
        /* Initialise */
        pol_speca[det_nr-1] = NULL ;
        pol_specb[det_nr-1] = NULL ;
        ext_plista[det_nr-1] = NULL ;
        ext_plistb[det_nr-1] = NULL ;

        /* Compute only one detector */
        if (reduce_det != 0 && det_nr != reduce_det) continue ;
    
        cpl_msg_info(__func__, "Process Detector %d", det_nr) ;
        cpl_msg_indent_more() ;

        /* Call the reduction function */
        if (cr2res_obs_pol_reduce(rawframes, raw_flat_frames, trace_wave_frame, 
                    detlin_frame, master_dark_frame, master_flat_frame, 
                    bpm_frame, 0, extract_oversample, extract_swath_width, 
                    extract_height, extract_smooth, det_nr,
                    &(pol_speca[det_nr-1]),
                    &(pol_specb[det_nr-1]),
                    &(ext_plista[det_nr-1]),
                    &(ext_plistb[det_nr-1])) == -1) {
            cpl_msg_error(__func__, "Failed to reduce detector %d", det_nr);
        }
        cpl_msg_indent_less() ;
    }

    /* Ѕave Products */
    out_file = cpl_sprintf("%s_pol_specA.fits", RECIPE_STRING) ;
    cr2res_io_save_POL_SPEC(out_file, frameset, rawframes, parlist,
            pol_speca, NULL, ext_plista, CR2RES_OBS_POL_SPECA_PROCATG, 
            RECIPE_STRING) ;
    cpl_free(out_file);

    out_file = cpl_sprintf("%s_pol_specB.fits", RECIPE_STRING) ;
    cr2res_io_save_POL_SPEC(out_file, frameset, rawframes, parlist,
            pol_specb, NULL, ext_plistb, CR2RES_OBS_POL_SPECB_PROCATG, 
            RECIPE_STRING) ;
    cpl_free(out_file);

    /* Free */
    cpl_frameset_delete(rawframes) ;
    if (raw_flat_frames != NULL) cpl_frameset_delete(raw_flat_frames) ;
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
        if (pol_speca[det_nr-1] != NULL) 
            cpl_table_delete(pol_speca[det_nr-1]) ;
        if (pol_specb[det_nr-1] != NULL) 
            cpl_table_delete(pol_specb[det_nr-1]) ;
        if (ext_plista[det_nr-1] != NULL) 
            cpl_propertylist_delete(ext_plista[det_nr-1]) ;
        if (ext_plistb[det_nr-1] != NULL) 
            cpl_propertylist_delete(ext_plistb[det_nr-1]) ;
    }

    return (int)cpl_error_get_code();
}
 
/*----------------------------------------------------------------------------*/
/**
  @brief    Execute the polarimetry recipe on a specific detector
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
  @param pol_speca              [out] polarimetry spectrum (A)
  @param pol_specb              [out] polarimetry spectrum (B)
  @param ext_plista             [out] the header for saving the products (A)
  @param ext_plistb             [out] the header for saving the products (B)
  @return  0 if ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
static int cr2res_obs_pol_reduce(
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
        cpl_table           **  pol_speca,
        cpl_table           **  pol_specb,
        cpl_propertylist    **  ext_plista,
        cpl_propertylist    **  ext_plistb)
{
    cpl_frameset        *   rawframes_a ;
    cpl_frameset        *   rawframes_b ;
    cr2res_nodding_pos  *   nod_positions ;
    cpl_table           *   pol_speca_loc ;
    cpl_table           *   pol_specb_loc ;
    cpl_propertylist    *   ext_plista_loc ;
    cpl_propertylist    *   ext_plistb_loc ;
    int                     i ;

    /* Check Inputs */
    if (pol_speca == NULL || pol_specb == NULL || ext_plista == NULL || 
            ext_plistb == NULL || rawframes == NULL || 
            trace_wave_frame == NULL) return -1 ;

    /* Get the Nodding positions */
    nod_positions = cr2res_nodding_read_positions(rawframes) ;
    if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
        for (i=0 ; i<cpl_frameset_get_size(rawframes) ; i++) {
            cpl_msg_debug(__func__, "Frame %s - Nodding %c",
                    cpl_frame_get_filename(
                        cpl_frameset_get_position_const(rawframes,i)),
                    cr2res_nodding_position_char(nod_positions[i])) ;
        }
    }

    /* Split the frames */
    if (cr2res_combine_nodding_split_frames(rawframes, nod_positions, 
                &rawframes_a, &rawframes_b)) {
        cpl_msg_error(__func__, "Failed to split the nodding positions") ;
        cpl_free(nod_positions) ;    
        return -1 ;
    }
    cpl_free(nod_positions) ;    

    /* Reduce A position */
    cpl_msg_info(__func__, "Compute Polarimetry for nodding A position") ;
    cpl_msg_indent_more() ;
    if (cr2res_obs_pol_reduce_one(rawframes_a, raw_flat_frames, 
                trace_wave_frame, detlin_frame, master_dark_frame, 
                master_flat_frame, bpm_frame, 0, extract_oversample, 
                extract_swath_width, extract_height, extract_smooth, reduce_det,
                &pol_speca_loc, &ext_plista_loc) == -1) {
        cpl_msg_error(__func__, "Failed to Reduce A nodding frames") ;
    }
    cpl_msg_indent_less() ;
    if (rawframes_a != NULL) cpl_frameset_delete(rawframes_a);

    /* Reduce B position */
    cpl_msg_info(__func__, "Compute Polarimetry for nodding B position") ;
    cpl_msg_indent_more() ;
    if (cr2res_obs_pol_reduce_one(rawframes_b, raw_flat_frames, 
                trace_wave_frame, detlin_frame, master_dark_frame, 
                master_flat_frame, bpm_frame, 0, extract_oversample, 
                extract_swath_width, extract_height, extract_smooth, reduce_det,
                &pol_specb_loc, &ext_plistb_loc) == -1) {
        cpl_msg_error(__func__, "Failed to Reduce B nodding frames") ;
    }
    cpl_msg_indent_less() ;
    if (rawframes_b != NULL) cpl_frameset_delete(rawframes_b);

    *pol_speca = pol_speca_loc ;
    *pol_specb = pol_specb_loc ;
    *ext_plista = ext_plista_loc ;
    *ext_plistb = ext_plistb_loc ;
    return 0 ;
}
 
/*----------------------------------------------------------------------------*/
/**
  @brief Execute the polarimetry computation for 1 nodding position
  @param rawframes              Raw science frames for 1 position
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
  @param pol_spec               [out] polarimetry spectrum
  @param ext_plist              [out] the header for saving the products
  @return  0 if ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
static int cr2res_obs_pol_reduce_one(
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
        cpl_table           **  pol_spec,
        cpl_propertylist    **  ext_plist)
{
    cpl_vector          *   dits ;
    cr2res_decker       *   decker_positions ;
    hdrl_imagelist      *   in ;
    hdrl_imagelist      *   in_calib ;
    int                 *   pol_sorting ;
    cpl_table           *   trace_wave ;
    cpl_table           *   trace_wave_corrected ;
    cpl_table           *   trace_wave_loc ;
    cpl_array           *   slit_frac ;
    const char          *   fname ;
    char                *   decker_name ;
    cpl_table           *   slit_func ;
    hdrl_image          *   model_master ;
    cpl_table           **  pol_spec_one_group ;
    cpl_table           **  extract_1d ;
    char                *   colname ;
    const double        *   pcol_data ;
    double              *   pvec_data ;
    cpl_vector          **  intens ;
    cpl_vector          **  wl ;
    cpl_vector          **  errors ;
    cpl_vector          **  demod_wl ;
    cpl_bivector        **  demod_stokes ;
    cpl_bivector        **  demod_null ;
    cpl_bivector        **  demod_intens ;
    int                 *   orders ;
    cpl_table           *   pol_spec_merged ;
    cpl_propertylist    *   ext_plist_loc ;
    char                *   out_file;
    cpl_size                nframes, nspec_group, spec_size ;
    int                     ngroups, i, j, k, l, o, norders, frame_idx ;

    /* Initialise */
    *pol_spec = NULL ;
    *ext_plist = NULL ;

    /* Check Inputs */
    if (pol_spec == NULL || ext_plist == NULL || rawframes == NULL || 
            trace_wave_frame == NULL) return -1 ;

    /* Check number of frames */
    nframes = cpl_frameset_get_size(rawframes) ;
    if (nframes == 0 || nframes % CR2RES_POLARIMETRY_GROUP_SIZE) {
        cpl_msg_error(__func__, 
            "Input number of frames is %"CPL_SIZE_FORMAT" and should be multiple of %d",
            nframes, CR2RES_POLARIMETRY_GROUP_SIZE) ;
        return -1 ;
    }

    /* Load the DITs if necessary */
    if (master_dark_frame != NULL)  dits = cr2res_io_read_dits(rawframes) ;
    else                            dits = NULL ;
    if (cpl_msg_get_level() == CPL_MSG_DEBUG && dits != NULL) 
        cpl_vector_dump(dits, stdout) ;

    /* Get the decker positions */
    if ((decker_positions = cr2res_io_read_decker_positions(rawframes))==NULL) {
        cpl_msg_error(__func__, "Cannot read the Decker positions") ;
        if (dits != NULL) cpl_vector_delete(dits) ;
        return -1 ;
    }

    /* Load image list */
    cpl_msg_info(__func__, "Load the images") ;
    if ((in = cr2res_io_load_image_list_from_set(rawframes, 
                    reduce_det)) == NULL) {
        cpl_msg_error(__func__, "Cannot load images") ;
        if (dits != NULL) cpl_vector_delete(dits) ;
        cpl_free(decker_positions) ;
        return -1 ;
    }
    if (hdrl_imagelist_get_size(in) != cpl_frameset_get_size(rawframes)) {
        cpl_msg_error(__func__, "Inconsistent number of loaded images") ;
        if (dits != NULL) cpl_vector_delete(dits) ;
        cpl_free(decker_positions) ;
        hdrl_imagelist_delete(in) ;
        return -1 ;
    }

    /* Calibrate the images */
    cpl_msg_info(__func__, "Apply the calibrations") ;
    if ((in_calib = cr2res_calib_imagelist(in, reduce_det, 0, 0, 
                    master_flat_frame, master_dark_frame, bpm_frame, 
                    detlin_frame, dits)) == NULL) {
        cpl_msg_error(__func__, "Failed to apply the calibrations") ;
        if (dits != NULL) cpl_vector_delete(dits) ;
        cpl_free(decker_positions) ;
        hdrl_imagelist_delete(in) ;
        return -1 ;
    }
    hdrl_imagelist_delete(in) ;
    if (dits != NULL) cpl_vector_delete(dits) ;

    /* Load the trace wave */
    cpl_msg_info(__func__, "Load the TRACE WAVE") ;
    if ((trace_wave = cr2res_io_load_TRACE_WAVE(cpl_frame_get_filename(
                        trace_wave_frame), reduce_det)) == NULL) {
        cpl_msg_error(__func__, "Failed to Load the traces file") ;
        cpl_free(decker_positions) ;
        hdrl_imagelist_delete(in_calib) ;
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

    /* Compute the number of groups */
    ngroups = nframes/CR2RES_POLARIMETRY_GROUP_SIZE ;
    nspec_group = 2*CR2RES_POLARIMETRY_GROUP_SIZE ;

    /* Allocate pol_spec_group containers */
    pol_spec_one_group = cpl_malloc(ngroups * sizeof(cpl_table*)) ;
    for (i = 0; i < ngroups; i++) pol_spec_one_group[i] = NULL;

    /* Loop on the groups */
    for (i=0 ; i<ngroups ; i++) {
        cpl_msg_info(__func__, "Process %d-group number %d/%d", 
                CR2RES_POLARIMETRY_GROUP_SIZE, i+1, ngroups) ;
        cpl_msg_indent_more() ;    

        /* Compute the proper order of the frames group */
        if ((pol_sorting = cr2res_pol_sort_frames(
                cpl_frameset_get_position_const(rawframes, 
                    i * CR2RES_POLARIMETRY_GROUP_SIZE),
                cpl_frameset_get_position_const(rawframes, 
                    i * CR2RES_POLARIMETRY_GROUP_SIZE + 1),
                cpl_frameset_get_position_const(rawframes, 
                    i * CR2RES_POLARIMETRY_GROUP_SIZE + 2),
                cpl_frameset_get_position_const(rawframes, 
                    i * CR2RES_POLARIMETRY_GROUP_SIZE + 3))) == NULL) {
            cpl_msg_warning(__func__, "Failed to sort the files") ;
            continue ;
        } 
        if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
            for (j=0 ; j<CR2RES_POLARIMETRY_GROUP_SIZE ; j++) {
                if (pol_sorting[j] != j) 
                    cpl_msg_warning(__func__, 
                            "Frame #%d moved to position #%d", 
                            j+1, pol_sorting[j]+1) ;
            }
        }

        /* Container for extracted tables: 1u, 1d, 2u, 2d, 3u, 3d, 4u, 4d */
        extract_1d = cpl_malloc(nspec_group * sizeof(cpl_table*));
        for (j = 0; j < nspec_group; j++) extract_1d[j] = NULL;

        /* Loop on the frames */
        for (j=0 ; j<CR2RES_POLARIMETRY_GROUP_SIZE ; j++) {
            frame_idx = i * CR2RES_POLARIMETRY_GROUP_SIZE + pol_sorting[j] ;
            fname = cpl_frame_get_filename(cpl_frameset_get_position_const(
                        rawframes,frame_idx)) ;

            /* Extract up */
            decker_name = cr2res_decker_print_position(
                    decker_positions[frame_idx]) ;
            cpl_msg_info(__func__, 
                    "Extract Up Spectrum from %s (Det %d / Decker %s)", 
                    fname, reduce_det, decker_name) ;
            cpl_free(decker_name) ;
            cpl_msg_indent_more() ;
           
            /* Get slit fraction for the upper trace */
            slit_frac = cr2res_trace_slit_fraction_create(
                    decker_positions[frame_idx], 1) ;

            /* Compute the new trace_wave for the extraction */
            trace_wave_loc = cr2res_trace_new_slit_fraction(trace_wave,
                            slit_frac) ;
            cpl_array_delete(slit_frac) ;

            /* Execute the extraction */
            cpl_msg_info(__func__, "Spectra Extraction") ;
            if (cr2res_extract_traces(
                        hdrl_imagelist_get_const(in_calib, frame_idx),
                        trace_wave_loc, -1, -1, CR2RES_EXTR_OPT_CURV, 
                        extract_height, extract_swath_width, extract_oversample,
                        extract_smooth, &(extract_1d[2*j]), &slit_func, 
                        &model_master) == -1) {
                cpl_msg_error(__func__, "Failed Extraction") ;
                extract_1d[2*j] = NULL ;
            } else {
                cpl_table_delete(slit_func) ;
                hdrl_image_delete(model_master) ;
                
                /* Save the table and the trace for debug */
                if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
                    out_file = cpl_sprintf(
                            "debug_ext_spec_det%d_group%d_%du.fits", 
                            reduce_det, i+1, j+1) ;
                    cpl_table_save(extract_1d[2*j], NULL, NULL,
                            out_file, CPL_IO_CREATE) ;
                    cpl_free(out_file) ; 
                    out_file = cpl_sprintf(
                            "debug_ext_trace_det%d_group%d_%du.fits", 
                            reduce_det, i+1, j+1) ;
                    cpl_table_save(trace_wave_loc, NULL, NULL,
                            out_file, CPL_IO_CREATE) ;
                    cpl_free(out_file) ; 
                }
            }
            cpl_table_delete(trace_wave_loc) ;
            cpl_msg_indent_less() ;

            /* Extract Down */
            decker_name = cr2res_decker_print_position(
                    decker_positions[frame_idx]) ;
            cpl_msg_info(__func__, 
                    "Extract Down Spectrum from %s (Det %d / Decker %s)", 
                    fname, reduce_det, decker_name) ;
            cpl_free(decker_name) ;
            cpl_msg_indent_more() ;
           
            /* Get slit fraction for the lower trace */
            slit_frac = cr2res_trace_slit_fraction_create(
                    decker_positions[frame_idx], 2) ;

            /* Compute the new trace_wave for the extraction */
            trace_wave_loc = cr2res_trace_new_slit_fraction(trace_wave,
                            slit_frac) ;
            cpl_array_delete(slit_frac) ;

            /* Execute the extraction */
            cpl_msg_info(__func__, "Spectra Extraction") ;
            if (cr2res_extract_traces(
                        hdrl_imagelist_get_const(in_calib, frame_idx),
                        trace_wave_loc, -1, -1, CR2RES_EXTR_OPT_CURV, 
                        extract_height, extract_swath_width, extract_oversample,
                        extract_smooth, &(extract_1d[2*j+1]), &slit_func, 
                        &model_master) == -1) {
                cpl_msg_error(__func__, "Failed Extraction") ;
                extract_1d[2*j+1] = NULL ;
            } else {
                cpl_table_delete(slit_func) ;
                hdrl_image_delete(model_master) ;

                /* Save the table and the trace for debug */
                if (cpl_msg_get_level() == CPL_MSG_DEBUG) {
                    out_file = cpl_sprintf(
                            "debug_ext_spec_det%d_group%d_%dd.fits", 
                            reduce_det, i+1, j+1) ;
                    cpl_table_save(extract_1d[2*j+1], NULL, NULL,
                            out_file, CPL_IO_CREATE) ;
                    cpl_free(out_file) ; 
                    out_file = cpl_sprintf(
                            "debug_ext_trace_det%d_group%d_%dd.fits", 
                            reduce_det, i+1, j+1) ;
                    cpl_table_save(trace_wave_loc, NULL, NULL,
                            out_file, CPL_IO_CREATE) ;
                    cpl_free(out_file) ; 
                }
            }
            cpl_table_delete(trace_wave_loc) ;
            cpl_msg_indent_less() ;
        }
        cpl_free(pol_sorting) ;

        /* How many orders */
        if ((orders = cr2res_obs_pol_get_order_numbers(
                        (const cpl_table **)extract_1d, 
                        nspec_group, 
                        &norders)) == NULL) {
            cpl_msg_warning(__func__, "No Order found in the extracted tables");
            for (j=0 ; j<nspec_group ; j++) 
                if (extract_1d[j] != NULL) cpl_table_delete(extract_1d[j]) ;
            cpl_free(extract_1d) ;
            continue ;
        }
        cpl_msg_debug(__func__, "%d different orders found", norders) ;

        /* Allocate data containerѕ */
        demod_wl = cpl_malloc(norders * sizeof(cpl_vector*)) ;
        demod_stokes = cpl_malloc(norders * sizeof(cpl_bivector*)) ;
        demod_null = cpl_malloc(norders * sizeof(cpl_bivector*)) ;
        demod_intens = cpl_malloc(norders * sizeof(cpl_bivector*)) ;
        for (o = 0; o < norders; o++){
            demod_wl[o] = NULL;
            demod_stokes[o] = NULL;
            demod_null[o] = NULL;
            demod_intens[o] = NULL;
        }

        /* Loop on the orders */
        for (o=0 ; o<norders ; o++) {
            cpl_msg_info(__func__, "Compute Polarimetry for order %d",
                    orders[o]) ;
            /* Get the inputs for the demod functions calls */
            intens = cpl_malloc(nspec_group * sizeof(cpl_vector*));
            wl = cpl_malloc(nspec_group * sizeof(cpl_vector*));
            errors = cpl_malloc(nspec_group * sizeof(cpl_vector*));
            for (k = 0; k < nspec_group; k++){
                intens[k] = NULL;
                wl[k] = NULL;
                errors[k] = NULL;
            }

            for (k=0 ; k<nspec_group ; k++) {
                spec_size = cpl_table_get_nrow(extract_1d[k]) ;
                /* Get the SPEC for this order/trace 1 */
                colname = cr2res_dfs_SPEC_colname(orders[o], 1) ;
                pcol_data = cpl_table_get_data_double_const(extract_1d[k], 
                        colname) ;
                cpl_free(colname) ;
                if (pcol_data == NULL) {
                    cpl_error_reset() ;
                } else {
                    intens[k] = cpl_vector_new(spec_size) ;
                    pvec_data = cpl_vector_get_data(intens[k]) ;
                    for (l=0 ; l<spec_size ; l++) 
                        pvec_data[l] = pcol_data[l] ;
                }

                /* Get the WAVELENGTH for this order/trace 1 */
                colname = cr2res_dfs_WAVELENGTH_colname(orders[o], 1) ;
                pcol_data = cpl_table_get_data_double_const(extract_1d[k], 
                        colname) ;
                cpl_free(colname) ;
                if (pcol_data == NULL) {
                    cpl_error_reset() ;
                } else {
                    wl[k] = cpl_vector_new(spec_size) ;
                    pvec_data = cpl_vector_get_data(wl[k]) ;
                    for (l=0 ; l<spec_size ; l++) 
                        pvec_data[l] = pcol_data[l] ;
                }

                /* Get the ERROR for this order/trace 1 */
                colname = cr2res_dfs_SPEC_ERR_colname(orders[o], 1) ;
                pcol_data = cpl_table_get_data_double_const(extract_1d[k], 
                        colname) ;
                cpl_free(colname) ;
                if (pcol_data == NULL) {
                    cpl_error_reset() ;
                } else {
                    errors[k] = cpl_vector_new(spec_size) ;
                    pvec_data = cpl_vector_get_data(errors[k]) ;
                    for (l=0 ; l<spec_size ; l++) 
                        pvec_data[l] = pcol_data[l] ;
                }
            }

            /* Keep the 1st wl of the 8 as reference for the output file*/
            if (wl[0] != NULL)  demod_wl[o] = cpl_vector_duplicate(wl[0]) ;
            else                demod_wl[o] = NULL ;

            /* Call Library Demodulation functions */
            demod_stokes[o] =   cr2res_pol_demod_stokes(intens, wl, errors,
                    nspec_group) ;
            demod_null[o] =     cr2res_pol_demod_null(intens, wl, errors,
                    nspec_group) ;
            demod_intens[o] =   cr2res_pol_demod_intens(intens, wl, errors,
                    nspec_group) ;

            /* Free */
            for (k=0 ; k<nspec_group ; k++) {
                cpl_vector_delete(intens[k]) ;
                cpl_vector_delete(wl[k]) ;
                cpl_vector_delete(errors[k]) ;
            }
            cpl_free(intens) ;
            cpl_free(wl) ;
            cpl_free(errors) ;
        }

        /* Free the extracted tables */
        for (j=0 ; j<nspec_group ; j++) 
            if (extract_1d[j] != NULL) cpl_table_delete(extract_1d[j]) ;
        cpl_free(extract_1d) ;

        /* Create the pol_spec table */
        cpl_msg_info(__func__, "Create the POL_SPEC table for this group");
        pol_spec_one_group[i] = cr2res_pol_POL_SPEC_create(orders, demod_wl, 
                demod_stokes, demod_null, demod_intens, norders) ;

        /* Deallocate */
        for (o=0 ; o<norders ; o++) {
            cpl_vector_delete(demod_wl[o]) ;
            cpl_bivector_delete(demod_stokes[o]) ;
            cpl_bivector_delete(demod_null[o]) ;
            cpl_bivector_delete(demod_intens[o]) ;
        }
        cpl_free(demod_wl) ;
        cpl_free(demod_stokes) ;
        cpl_free(demod_null) ;
        cpl_free(demod_intens) ;
        cpl_free(orders) ;
        cpl_msg_indent_less() ;
    }
    cpl_free(decker_positions) ;
    cpl_table_delete(trace_wave) ;
    hdrl_imagelist_delete(in_calib) ;

    /* Merge the groups together */
    if (ngroups > 1) {
        cpl_msg_info(__func__, "Merge the %d groups POL_SPEC tables into one", 
                ngroups);
        pol_spec_merged = cr2res_pol_spec_pol_merge(
                (const cpl_table **)pol_spec_one_group, ngroups) ;
    } else {
        if (pol_spec_one_group[0] == NULL) {
            pol_spec_merged = NULL ;
        } else {
            pol_spec_merged = cpl_table_duplicate(pol_spec_one_group[0]) ;
        }
    }

    /* Deallocate */
    for (i=0 ; i<ngroups ; i++) cpl_table_delete(pol_spec_one_group[i]) ;
    cpl_free(pol_spec_one_group) ;

    /* Check */
    if (pol_spec_merged == NULL) {
        cpl_msg_error(__func__, "Cannot create the POL_SPEC table");
        return -1 ;
    }

    /* QCs */
    cpl_msg_info(__func__, "Store the QC parameters") ;
    ext_plist_loc = cpl_propertylist_new() ;
    /* TODO */

    /* Return */
    *pol_spec = pol_spec_merged ;
    *ext_plist = ext_plist_loc ;
    return 0 ;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Count and return the order numbers from extracted tables
  @param    extracted       Array of extracted tables
  @param    next            Number of extracted tables
  @param    norders         [out] Number of orders in the output table
  @return   newly allocated int array

  The int array will need to be freed by the caller. Its size iѕ
  norders. It contains the list of orders found in all the input ext tables
 */
/*----------------------------------------------------------------------------*/
static int * cr2res_obs_pol_get_order_numbers(
        const cpl_table **  extracted,
        int                 next,
        int             *   norders)
{
    cpl_array   *   col_names ;
    const char  *   col_name ;
    char        *   col_type ;
    int         *   tmp_orders_list ;
    int             count_orders ;
    int         *   out_orders ;
    cpl_size        j, ncols ;
    int             i, k, order, trace_nb, max_possible_spectra, new_order ;

    /* Check entries */
    if (extracted == NULL || norders == NULL) return NULL ;
    for (i=0 ; i<next ; i++) 
        if (extracted[i] == NULL) return NULL ;

    /* Initialize */
    max_possible_spectra = 1000 ;

    /* Allocate orders list */
    tmp_orders_list = cpl_malloc(max_possible_spectra * sizeof(int)) ;

    /* Count the different orders */
    count_orders = 0 ;
    /* Loop over all columns */
    for (i=0 ; i<next ; i++) {
        col_names = cpl_table_get_column_names(extracted[i]);
        ncols = cpl_table_get_ncol(extracted[i]) ;
        for (j=0 ; j<ncols ; j++) {
            col_name = cpl_array_get_string(col_names, j);
            col_type = cr2res_dfs_SPEC_colname_parse(col_name, &order,
                    &trace_nb) ;
            if (col_type != NULL && !strcmp(col_type, "SPEC")) {
                /* Is the current order a new one ? */
                new_order = 1 ;
                for (k=0 ; k<count_orders ; k++)
                    if (tmp_orders_list[k] == order)
                        new_order = 0 ;

                /* Current order not yet stored */
                if (new_order) {
                    tmp_orders_list[count_orders] = order ;
                    count_orders ++ ;
                }
            } 
            if (col_type != NULL) cpl_free(col_type) ;
        }
        cpl_array_delete(col_names) ;
    }

    /* Allocate and fill output array */
    out_orders = cpl_malloc(count_orders * sizeof(int)) ;
    for (i=0 ; i<count_orders ; i++) out_orders[i] = tmp_orders_list[i] ;

    /* Free and return */
    cpl_free(tmp_orders_list) ;
    *norders = count_orders ;
    return out_orders ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief  Run basic checks for the rawframes consistency
  @param    rawframes   The input rawframes
  @return   1 if valid, 0 if not, -1 in error case
 */
/*----------------------------------------------------------------------------*/
static int cr2res_obs_pol_check_inputs_validity(
        const cpl_frameset  *   rawframes)
{
    cpl_size                    nframes ;

    /* Check Inputs */
    if (rawframes == NULL) return -1 ;

    /* number of frames must be modulo 4 */
    nframes = cpl_frameset_get_size(rawframes) ;
    if (nframes % 4) {
        cpl_msg_error(__func__, "Require a multiple of 4 raw frames") ;
        return 0 ;
    }

    /* TODO : use cr2res_combine_nodding_split_frames() */

    /* TODO : Check that all frames are decker 1_3 or 2_4 */

    return 1 ;
}

