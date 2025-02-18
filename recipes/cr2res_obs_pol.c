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
#include "cr2res_idp.h"
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
        const cpl_frameset  *   rawframes,
        cpl_size            *   ngroups) ;
static int cr2res_obs_pol_reduce(
        const cpl_frameset  *   rawframes,
        const cpl_frame     *   trace_wave_frame,
        const cpl_frame     *   detlin_frame,
        const cpl_frame     *   master_dark_frame,
        const cpl_frame     *   master_flat_frame,
        const cpl_frame     *   bpm_frame,
        const cpl_frame     *   blaze_frame,
        int                     subtract_nolight_rows,
        int                     subtract_interorder_column,
        int                     cosmics,
        int                     extract_oversample,
        int                     extract_swath_width,
        int                     extract_height,
        double                  extract_smooth_slit,
        double                  extract_smooth_spec,
        int                     save_group,
        int                     reduce_det,
        hdrl_image          **  in_calib_1_a,
        hdrl_image          **  in_calib_2_a,
        hdrl_image          **  in_calib_3_a,
        hdrl_image          **  in_calib_4_a,
        cpl_table           **  trace_wave_1u_a,
        cpl_table           **  trace_wave_1d_a,
        cpl_table           **  trace_wave_2u_a,
        cpl_table           **  trace_wave_2d_a,
        cpl_table           **  trace_wave_3u_a,
        cpl_table           **  trace_wave_3d_a,
        cpl_table           **  trace_wave_4u_a,
        cpl_table           **  trace_wave_4d_a,
        cpl_table           **  extract1D_1u_a,
        cpl_table           **  extract1D_1d_a,
        cpl_table           **  extract1D_2u_a,
        cpl_table           **  extract1D_2d_a,
        cpl_table           **  extract1D_3u_a,
        cpl_table           **  extract1D_3d_a,
        cpl_table           **  extract1D_4u_a,
        cpl_table           **  extract1D_4d_a,
        cpl_table           **  pol_speca,
        hdrl_image          **  in_calib_1_b,
        hdrl_image          **  in_calib_2_b,
        hdrl_image          **  in_calib_3_b,
        hdrl_image          **  in_calib_4_b,
        cpl_table           **  trace_wave_1u_b,
        cpl_table           **  trace_wave_1d_b,
        cpl_table           **  trace_wave_2u_b,
        cpl_table           **  trace_wave_2d_b,
        cpl_table           **  trace_wave_3u_b,
        cpl_table           **  trace_wave_3d_b,
        cpl_table           **  trace_wave_4u_b,
        cpl_table           **  trace_wave_4d_b,
        cpl_table           **  extract1D_1u_b,
        cpl_table           **  extract1D_1d_b,
        cpl_table           **  extract1D_2u_b,
        cpl_table           **  extract1D_2d_b,
        cpl_table           **  extract1D_3u_b,
        cpl_table           **  extract1D_3d_b,
        cpl_table           **  extract1D_4u_b,
        cpl_table           **  extract1D_4d_b,
        cpl_table           **  pol_specb,
        cpl_propertylist    **  ext_plista,
        cpl_propertylist    **  ext_plistb) ;
static int cr2res_obs_pol_reduce_one(
        const cpl_frameset  *   rawframes,
        const cpl_frameset  *   raw_background_frames,
        const cpl_frame     *   trace_wave_frame,
        const cpl_frame     *   detlin_frame,
        const cpl_frame     *   master_dark_frame,
        const cpl_frame     *   master_flat_frame,
        const cpl_frame     *   bpm_frame,
        const cpl_frame     *   blaze_frame,
        int                     subtract_nolight_rows,
        int                     subtract_interorder_column,
        int                     cosmics,
        int                     extract_oversample,
        int                     extract_swath_width,
        int                     extract_height,
        double                  extract_smooth_slit,
        double                  extract_smooth_spec,
        int                     save_group,
        int                     reduce_det,
        hdrl_image          *** in_calib,
        cpl_table           *** tw,
        cpl_table           *** extract1D,
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
  of 4 frames. Each block generates 8 extractions. For each order       \n\
  found, 8 spectra are passed to the demodulation functions.            \n\
  The results are stored in polarimetry tables (1 per group of 4        \n\
  frames). The polarimetry tables are then merged together if there are \n\
  several group of 4 frames available.                                  \n\
                                                                        \n\
  Inputs                                                                \n\
    raw.fits " CR2RES_OBS_POLARIMETRY_OTHER_RAW " [4 to 4n]             \n\
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
                  1->4 is derived with cr2res_pol_sort_frames()         \n\
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
            -> demod_intens(o, g)                                       \n\
        Create pol_spec(g) from demod_stokes(g),                        \n\
                                demod_null(g),                          \n\
                                demod_intens(g)                         \n\
      Merge pol_spec(g) into pol_spec                                   \n\
                                                                        \n\
  Library functions used                                                \n\
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
                    CR2RES_PIPELINE_AUTHORS,
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
    p = cpl_parameter_new_value("cr2res.cr2res_obs_pol.subtract_nolight_rows",
            CPL_TYPE_BOOL, 
            "Subtract median row from baffled region at detector bottom",
            "cr2res.cr2res_obs_pol", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "subtract_nolight_rows");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value(
            "cr2res.cr2res_obs_pol.subtract_interorder_column", CPL_TYPE_BOOL,
            "Subtract column-by-column fit to the pixel values between orders",
            "cr2res.cr2res_obs_pol", TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
            "subtract_interorder_column");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_pol.cosmics",
            CPL_TYPE_BOOL, "Find and mark cosmic rays hits as bad",
            "cr2res.cr2res_obs_pol", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "cosmics");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_pol.extract_oversample",
            CPL_TYPE_INT, "factor by which to oversample the extraction",
            "cr2res.cr2res_obs_pol", 5);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_oversample");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_pol.extract_swath_width",
            CPL_TYPE_INT, "The swath width", "cr2res.cr2res_obs_pol", 2048);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_swath_width");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_pol.extract_height",
            CPL_TYPE_INT, "Extraction height", "cr2res.cr2res_obs_pol", -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_height");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_pol.extract_smooth_spec",
            CPL_TYPE_DOUBLE, "Smoothing along the spectrum",
            "cr2res.cr2res_obs_pol", 0.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_smooth_spec");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_pol.extract_smooth_slit",
            CPL_TYPE_DOUBLE, "Smoothing along the slit",
            "cr2res.cr2res_obs_pol", 2.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_smooth_slit");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_pol.detector", CPL_TYPE_INT,
            "Only reduce the specified detector", "cr2res.cr2res_obs_pol", 0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "detector");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_nodding.create_idp",
                                CPL_TYPE_BOOL, "Flag to produce  IDP files",
                                "cr2res.cr2res_obs_nodding", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "idp");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_pol.save_group",
            CPL_TYPE_INT, 
            "Save extra files for the specified group number (0: no save)",
            "cr2res.cr2res_obs_pol", 0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "save_group");
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
    int                     extract_oversample, extract_swath_width, cosmics,
                            extract_height, reduce_det, create_idp, 
                            subtract_nolight_rows,
                            subtract_interorder_column, save_group ;
    double                  extract_smooth_slit, extract_smooth_spec ;
    double                  barycorr;
    cpl_frameset        *   rawframes ;
    cpl_frameset        *   raw_flat_frames ;
    const cpl_frame     *   trace_wave_frame ;
    const cpl_frame     *   detlin_frame ;
    const cpl_frame     *   master_dark_frame ;
    const cpl_frame     *   master_flat_frame ;
    const cpl_frame     *   bpm_frame ;
    const cpl_frame     *   blaze_frame ;
    cpl_propertylist    *   qc_main ;
    hdrl_image          *   in_calib_1_a[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   in_calib_2_a[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   in_calib_3_a[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   in_calib_4_a[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   in_calib_1_b[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   in_calib_2_b[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   in_calib_3_b[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   in_calib_4_b[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	tw_1u_a[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	tw_1d_a[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	tw_2u_a[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	tw_2d_a[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	tw_3u_a[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	tw_3d_a[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	tw_4u_a[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	tw_4d_a[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	tw_1u_b[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	tw_1d_b[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	tw_2u_b[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	tw_2d_b[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	tw_3u_b[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	tw_3d_b[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	tw_4u_b[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	tw_4d_b[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	extract1D_1u_a[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	extract1D_1d_a[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	extract1D_2u_a[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	extract1D_2d_a[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	extract1D_3u_a[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	extract1D_3d_a[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	extract1D_4u_a[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	extract1D_4d_a[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	extract1D_1u_b[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	extract1D_1d_b[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	extract1D_2u_b[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	extract1D_2d_b[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	extract1D_3u_b[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	extract1D_3d_b[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	extract1D_4u_b[CR2RES_NB_DETECTORS] ;
    cpl_table           *  	extract1D_4d_b[CR2RES_NB_DETECTORS] ;
    cpl_table           *   pol_speca[CR2RES_NB_DETECTORS] ;
    cpl_table           *   pol_specb[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   ext_plista[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   ext_plistb[CR2RES_NB_DETECTORS] ;
    char                *   out_file;
    cpl_table           *   eop_table ;
    int                     det_nr, save_products ; 

    /* Initialise */

    /* RETRIEVE INPUT PARAMETERS */
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_pol.subtract_nolight_rows");
    subtract_nolight_rows = cpl_parameter_get_bool(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_pol.subtract_interorder_column");
    subtract_interorder_column = cpl_parameter_get_bool(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_pol.cosmics");
    cosmics = cpl_parameter_get_bool(param);
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
            "cr2res.cr2res_obs_pol.extract_smooth_slit");
    extract_smooth_slit = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_pol.extract_smooth_spec");
    extract_smooth_spec = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_pol.detector");
    reduce_det = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_nodding.create_idp");
    create_idp = cpl_parameter_get_bool(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_pol.save_group");
    save_group = cpl_parameter_get_int(param);

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
    blaze_frame = cpl_frameset_find_const(frameset,
            CR2RES_CAL_FLAT_EXTRACT_1D_PROCATG) ;

    /* Get the RAW Frames */
    rawframes = cr2res_extract_frameset(frameset, 
            CR2RES_OBS_POLARIMETRY_OTHER_RAW) ;
    if (rawframes == NULL) {
        cpl_msg_error(__func__, "Could not find RAW frames") ;
        return -1 ;
    }
   
    /* Get the RAW flat frames */
    raw_flat_frames = cr2res_extract_frameset(frameset, CR2RES_FLAT_RAW) ;

    /* Loop on the detectors */
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
        /* Initialise */
        in_calib_1_a[det_nr-1] = NULL ;
        in_calib_2_a[det_nr-1] = NULL ;
        in_calib_3_a[det_nr-1] = NULL ;
        in_calib_4_a[det_nr-1] = NULL ;
        in_calib_1_b[det_nr-1] = NULL ;
        in_calib_2_b[det_nr-1] = NULL ;
        in_calib_3_b[det_nr-1] = NULL ;
        in_calib_4_b[det_nr-1] = NULL ;
        tw_1u_a[det_nr-1] = NULL ;
        tw_1d_a[det_nr-1] = NULL ;
        tw_2u_a[det_nr-1] = NULL ;
        tw_2d_a[det_nr-1] = NULL ;
        tw_3u_a[det_nr-1] = NULL ;
        tw_3d_a[det_nr-1] = NULL ;
        tw_4u_a[det_nr-1] = NULL ;
        tw_4d_a[det_nr-1] = NULL ;
        tw_1u_b[det_nr-1] = NULL ;
        tw_1d_b[det_nr-1] = NULL ;
        tw_2u_b[det_nr-1] = NULL ;
        tw_2d_b[det_nr-1] = NULL ;
        tw_3u_b[det_nr-1] = NULL ;
        tw_3d_b[det_nr-1] = NULL ;
        tw_4u_b[det_nr-1] = NULL ;
        tw_4d_b[det_nr-1] = NULL ;
        extract1D_1u_a[det_nr-1] = NULL ;
        extract1D_1d_a[det_nr-1] = NULL ;
        extract1D_2u_a[det_nr-1] = NULL ;
        extract1D_2d_a[det_nr-1] = NULL ;
        extract1D_3u_a[det_nr-1] = NULL ;
        extract1D_3d_a[det_nr-1] = NULL ;
        extract1D_4u_a[det_nr-1] = NULL ;
        extract1D_4d_a[det_nr-1] = NULL ;
        extract1D_1u_b[det_nr-1] = NULL ;
        extract1D_1d_b[det_nr-1] = NULL ;
        extract1D_2u_b[det_nr-1] = NULL ;
        extract1D_2d_b[det_nr-1] = NULL ;
        extract1D_3u_b[det_nr-1] = NULL ;
        extract1D_3d_b[det_nr-1] = NULL ;
        extract1D_4u_b[det_nr-1] = NULL ;
        extract1D_4d_b[det_nr-1] = NULL ;
        pol_speca[det_nr-1] = NULL ;
        pol_specb[det_nr-1] = NULL ;
        ext_plista[det_nr-1] = NULL ;
        ext_plistb[det_nr-1] = NULL ;

        /* Compute only one detector */
        if (reduce_det != 0 && det_nr != reduce_det) continue ;
    
        cpl_msg_info(__func__, "Process Detector %d", det_nr) ;
        cpl_msg_indent_more() ;

        /* Call the reduction function */
        if (cr2res_obs_pol_reduce(rawframes, 
                    trace_wave_frame, detlin_frame, master_dark_frame, 
                    master_flat_frame, bpm_frame, blaze_frame, 
                    subtract_nolight_rows, subtract_interorder_column, 
                    cosmics, extract_oversample, extract_swath_width, 
                    extract_height, extract_smooth_slit, 
                    extract_smooth_spec, save_group, det_nr,
                    &(in_calib_1_a[det_nr-1]),
                    &(in_calib_2_a[det_nr-1]),
                    &(in_calib_3_a[det_nr-1]),
                    &(in_calib_4_a[det_nr-1]),
                    &(tw_1u_a[det_nr-1]),
                    &(tw_1d_a[det_nr-1]),
                    &(tw_2u_a[det_nr-1]),
                    &(tw_2d_a[det_nr-1]),
                    &(tw_3u_a[det_nr-1]),
                    &(tw_3d_a[det_nr-1]),
                    &(tw_4u_a[det_nr-1]),
                    &(tw_4d_a[det_nr-1]),
                    &(extract1D_1u_a[det_nr-1]),
                    &(extract1D_1d_a[det_nr-1]),
                    &(extract1D_2u_a[det_nr-1]),
                    &(extract1D_2d_a[det_nr-1]),
                    &(extract1D_3u_a[det_nr-1]),
                    &(extract1D_3d_a[det_nr-1]),
                    &(extract1D_4u_a[det_nr-1]),
                    &(extract1D_4d_a[det_nr-1]),
                    &(pol_speca[det_nr-1]),
                    &(in_calib_1_b[det_nr-1]),
                    &(in_calib_2_b[det_nr-1]),
                    &(in_calib_3_b[det_nr-1]),
                    &(in_calib_4_b[det_nr-1]),
                    &(tw_1u_b[det_nr-1]),
                    &(tw_1d_b[det_nr-1]),
                    &(tw_2u_b[det_nr-1]),
                    &(tw_2d_b[det_nr-1]),
                    &(tw_3u_b[det_nr-1]),
                    &(tw_3d_b[det_nr-1]),
                    &(tw_4u_b[det_nr-1]),
                    &(tw_4d_b[det_nr-1]),
                    &(extract1D_1u_b[det_nr-1]),
                    &(extract1D_1d_b[det_nr-1]),
                    &(extract1D_2u_b[det_nr-1]),
                    &(extract1D_2d_b[det_nr-1]),
                    &(extract1D_3u_b[det_nr-1]),
                    &(extract1D_3d_b[det_nr-1]),
                    &(extract1D_4u_b[det_nr-1]),
                    &(extract1D_4d_b[det_nr-1]),
                    &(pol_specb[det_nr-1]),
                    &(ext_plista[det_nr-1]),
                    &(ext_plistb[det_nr-1])) == -1) {
            cpl_msg_warning(__func__, "Failed to reduce detector %d", 
                    det_nr);
            cpl_error_reset() ;
        }
        cpl_msg_indent_less() ;
    }

	/* Add ESO.DRS.TMID in the Main Header */
	qc_main = cpl_propertylist_new();
	cpl_propertylist_append_double(qc_main,
			CR2RES_HEADER_DRS_TMID,
			cr2res_utils_get_center_mjd(rawframes)) ;

    /* Add barycentric correction */
    eop_table = cr2res_io_get_eop_table() ;
    if (eop_table != NULL) {
        double ra, dec, mjd_obs, geolon, geolat, geoelev;
        cpl_propertylist    *   plist ;
        plist=cpl_propertylist_load(cpl_frame_get_filename(
                    cpl_frameset_get_position_const(rawframes, 0)), 0) ;

        ra = cpl_propertylist_get_double(plist, "RA") ;
        dec = cpl_propertylist_get_double(plist, "DEC") ;
        mjd_obs = cpl_propertylist_get_double(plist, "MJD-OBS") ;
        geolon = cpl_propertylist_get_double(plist, "ESO TEL GEOLON") ;
        geolat = cpl_propertylist_get_double(plist, "ESO TEL GEOLAT") ;
        geoelev = cpl_propertylist_get_double(plist, "ESO TEL GEOELEV") ;

        cpl_propertylist_delete(plist) ;

        barycorr = 0.0 ;
        if (!cpl_error_get_code()) {
            double mjd_cen;
            mjd_cen = cr2res_utils_get_center_mjd(rawframes) ;
            hdrl_barycorr_compute(ra, dec, eop_table, mjd_obs,
                    (mjd_cen-mjd_obs)*24*3600, geolon, geolat, geoelev,
                    0.0, 0.0, 0.0, 0.0, &barycorr);

            cpl_msg_info(__func__, "Barycentric correction: %g m/s",
                    barycorr);
        } else {
            cpl_msg_info(__func__, "Cannot derive Barycentric correction");
            cpl_error_reset() ;
        }
        cpl_table_delete(eop_table) ;
        cpl_propertylist_append_double(qc_main, CR2RES_HEADER_DRS_BARYCORR,
                barycorr);
    }

    /* Add QC NUMSAT */
	cpl_propertylist_append_int(qc_main,
			CR2RES_HEADER_QC_NUMSAT,
			cr2res_qc_numsat(rawframes)) ;

    /* Save Extra Products if at least one is returned */
    save_products = 0 ;
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
        if (in_calib_1_a[det_nr-1] != NULL) save_products = 1 ;
    }

    /* Save only the used RAW - fill rawframes with CALIBS */
	if (trace_wave_frame != NULL)
		cpl_frameset_insert(rawframes,
				cpl_frame_duplicate(trace_wave_frame)) ;
	if (detlin_frame != NULL)
		cpl_frameset_insert(rawframes,
				cpl_frame_duplicate(detlin_frame)) ;
	if (master_dark_frame != NULL)
		cpl_frameset_insert(rawframes,
				cpl_frame_duplicate(master_dark_frame)) ;
	if (master_flat_frame!= NULL)
		cpl_frameset_insert(rawframes,
				cpl_frame_duplicate(master_flat_frame)) ;
	if (bpm_frame!= NULL)
		cpl_frameset_insert(rawframes,
				cpl_frame_duplicate(bpm_frame)) ;
	if (blaze_frame!= NULL)
		cpl_frameset_insert(rawframes,
				cpl_frame_duplicate(blaze_frame)) ;

    if (save_products == 1) {
        /* Calibrated images a */
        out_file = cpl_sprintf("%s_in_calibA_1.fits", RECIPE_STRING) ;
        cr2res_io_save_CALIBRATED(out_file, frameset, rawframes, parlist,
                in_calib_1_a, NULL, ext_plista, 
                CR2RES_OBS_POL_CALIB_A_PROCATG, RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_in_calibA_2.fits", RECIPE_STRING) ;
        cr2res_io_save_CALIBRATED(out_file, frameset, rawframes, parlist,
                in_calib_2_a, NULL, ext_plista, 
                CR2RES_OBS_POL_CALIB_A_PROCATG, RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_in_calibA_3.fits", RECIPE_STRING) ;
        cr2res_io_save_CALIBRATED(out_file, frameset, rawframes, parlist,
                in_calib_3_a, NULL, ext_plista, 
                CR2RES_OBS_POL_CALIB_A_PROCATG, RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_in_calibA_4.fits", RECIPE_STRING) ;
        cr2res_io_save_CALIBRATED(out_file, frameset, rawframes, parlist,
                in_calib_4_a, NULL, ext_plista, 
                CR2RES_OBS_POL_CALIB_A_PROCATG, RECIPE_STRING);
        cpl_free(out_file);

        /* TW Nodding a */
        out_file = cpl_sprintf("%s_twA_1u.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                tw_1u_a, NULL, ext_plista, CR2RES_OBS_POL_TWA_PROCATG, 
                RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_twA_1d.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                tw_1d_a, NULL, ext_plista, CR2RES_OBS_POL_TWA_PROCATG, 
                RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_twA_2u.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                tw_2u_a, NULL, ext_plista, CR2RES_OBS_POL_TWA_PROCATG, 
                RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_twA_2d.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                tw_2d_a, NULL, ext_plista, CR2RES_OBS_POL_TWA_PROCATG, 
                RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_twA_3u.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                tw_3u_a, NULL, ext_plista, CR2RES_OBS_POL_TWA_PROCATG, 
                RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_twA_3d.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                tw_3d_a, NULL, ext_plista, CR2RES_OBS_POL_TWA_PROCATG, 
                RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_twA_4u.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                tw_4u_a, NULL, ext_plista, CR2RES_OBS_POL_TWA_PROCATG, 
                RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_twA_4d.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                tw_4d_a, NULL, ext_plista, CR2RES_OBS_POL_TWA_PROCATG, 
                RECIPE_STRING);
        cpl_free(out_file);

        /* Extracted Nodding A */
        out_file = cpl_sprintf("%s_extractedA_1u.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                extract1D_1u_a, NULL, ext_plista, 
                CR2RES_OBS_POL_EXTRACTA_PROCATG, RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_extractedA_1d.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                extract1D_1d_a, NULL, ext_plista, 
                CR2RES_OBS_POL_EXTRACTA_PROCATG, RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_extractedA_2u.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                extract1D_2u_a, NULL, ext_plista, 
                CR2RES_OBS_POL_EXTRACTA_PROCATG, RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_extractedA_2d.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                extract1D_2d_a, NULL, ext_plista, 
                CR2RES_OBS_POL_EXTRACTA_PROCATG, RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_extractedA_3u.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                extract1D_3u_a, NULL, ext_plista, 
                CR2RES_OBS_POL_EXTRACTA_PROCATG, RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_extractedA_3d.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                extract1D_3d_a, NULL, ext_plista, 
                CR2RES_OBS_POL_EXTRACTA_PROCATG, RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_extractedA_4u.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                extract1D_4u_a, NULL, ext_plista, 
                CR2RES_OBS_POL_EXTRACTA_PROCATG, RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_extractedA_4d.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                extract1D_4d_a, NULL, ext_plista, 
                CR2RES_OBS_POL_EXTRACTA_PROCATG, RECIPE_STRING);
        cpl_free(out_file);

        /* Calibrated images b */
        out_file = cpl_sprintf("%s_in_calibB_1.fits", RECIPE_STRING) ;
        cr2res_io_save_CALIBRATED(out_file, frameset, rawframes, parlist,
                in_calib_1_b, NULL, ext_plista, 
                CR2RES_OBS_POL_CALIB_B_PROCATG, RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_in_calibB_2.fits", RECIPE_STRING) ;
        cr2res_io_save_CALIBRATED(out_file, frameset, rawframes, parlist,
                in_calib_2_b, NULL, ext_plista, 
                CR2RES_OBS_POL_CALIB_B_PROCATG, RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_in_calibB_3.fits", RECIPE_STRING) ;
        cr2res_io_save_CALIBRATED(out_file, frameset, rawframes, parlist,
                in_calib_3_b, NULL, ext_plista, 
                CR2RES_OBS_POL_CALIB_B_PROCATG, RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_in_calibB_4.fits", RECIPE_STRING) ;
        cr2res_io_save_CALIBRATED(out_file, frameset, rawframes, parlist,
                in_calib_4_b, NULL, ext_plista, 
                CR2RES_OBS_POL_CALIB_B_PROCATG, RECIPE_STRING);
        cpl_free(out_file);

        /* TW Nodding b */
        out_file = cpl_sprintf("%s_twB_1u.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                tw_1u_b, NULL, ext_plistb, CR2RES_OBS_POL_TWB_PROCATG, 
                RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_twB_1d.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                tw_1d_b, NULL, ext_plistb, CR2RES_OBS_POL_TWB_PROCATG, 
                RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_twB_2u.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                tw_2u_b, NULL, ext_plistb, CR2RES_OBS_POL_TWB_PROCATG, 
                RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_twB_2d.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                tw_2d_b, NULL, ext_plistb, CR2RES_OBS_POL_TWB_PROCATG, 
                RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_twB_3u.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                tw_3u_b, NULL, ext_plistb, CR2RES_OBS_POL_TWB_PROCATG, 
                RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_twB_3d.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                tw_3d_b, NULL, ext_plistb, CR2RES_OBS_POL_TWB_PROCATG, 
                RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_twB_4u.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                tw_4u_b, NULL, ext_plistb, CR2RES_OBS_POL_TWB_PROCATG, 
                RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_twB_4d.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                tw_4d_b, NULL, ext_plistb, CR2RES_OBS_POL_TWB_PROCATG, 
                RECIPE_STRING);
        cpl_free(out_file);

        /* Extracted Nodding B */
        out_file = cpl_sprintf("%s_extractedB_1u.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                extract1D_1u_b, NULL, ext_plistb, 
                CR2RES_OBS_POL_EXTRACTB_PROCATG, RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_extractedB_1d.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                extract1D_1d_b, NULL, ext_plistb, 
                CR2RES_OBS_POL_EXTRACTB_PROCATG, RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_extractedB_2u.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                extract1D_2u_b, NULL, ext_plistb, 
                CR2RES_OBS_POL_EXTRACTB_PROCATG, RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_extractedB_2d.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                extract1D_2d_b, NULL, ext_plistb, 
                CR2RES_OBS_POL_EXTRACTB_PROCATG, RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_extractedB_3u.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                extract1D_3u_b, NULL, ext_plistb, 
                CR2RES_OBS_POL_EXTRACTB_PROCATG, RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_extractedB_3d.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                extract1D_3d_b, NULL, ext_plistb, 
                CR2RES_OBS_POL_EXTRACTB_PROCATG, RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_extractedB_4u.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                extract1D_4u_b, NULL, ext_plistb, 
                CR2RES_OBS_POL_EXTRACTB_PROCATG, RECIPE_STRING);
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_extractedB_4d.fits", RECIPE_STRING) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, 
                extract1D_4d_b, NULL, ext_plistb, 
                CR2RES_OBS_POL_EXTRACTB_PROCATG, RECIPE_STRING);
        cpl_free(out_file);
    } 

    /* Polarimetry Spectra */
    out_file = cpl_sprintf("%s_pol_specA.fits", RECIPE_STRING) ;
    cr2res_io_save_POL_SPEC(out_file, frameset, rawframes, parlist,
            pol_speca, qc_main, ext_plista, CR2RES_OBS_POL_SPECA_PROCATG, 
            RECIPE_STRING) ;
    if (create_idp) {
        cr2res_idp_save(out_file, frameset, rawframes, parlist,
                        pol_speca, qc_main, ext_plista,
                        CR2RES_OBS_POL_EXTRACTA_IDP_PROCATG,
                        RECIPE_STRING);
    }
    cpl_free(out_file);

    out_file = cpl_sprintf("%s_pol_specB.fits", RECIPE_STRING) ;
    cr2res_io_save_POL_SPEC(out_file, frameset, rawframes, parlist,
            pol_specb, qc_main, ext_plistb, CR2RES_OBS_POL_SPECB_PROCATG, 
            RECIPE_STRING) ;
    if (create_idp) {
        cr2res_idp_save(out_file, frameset, rawframes, parlist,
                        pol_specb, qc_main, ext_plistb,
                        CR2RES_OBS_POL_EXTRACTB_IDP_PROCATG,
                        RECIPE_STRING);
    }
    //cpl_msg_info(__func__,"ERR: %d", cpl_error_get_code());
    cpl_free(out_file);

    /* Free */
    cpl_frameset_delete(rawframes) ;
    cpl_propertylist_delete(qc_main) ;
    if (raw_flat_frames != NULL) cpl_frameset_delete(raw_flat_frames) ;
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
        if (in_calib_1_a[det_nr-1] != NULL) 
            hdrl_image_delete(in_calib_1_a[det_nr-1]) ;
        if (in_calib_2_a[det_nr-1] != NULL) 
            hdrl_image_delete(in_calib_2_a[det_nr-1]) ;
        if (in_calib_3_a[det_nr-1] != NULL) 
            hdrl_image_delete(in_calib_3_a[det_nr-1]) ;
        if (in_calib_4_a[det_nr-1] != NULL) 
            hdrl_image_delete(in_calib_4_a[det_nr-1]) ;
        if (in_calib_1_b[det_nr-1] != NULL) 
            hdrl_image_delete(in_calib_1_b[det_nr-1]) ;
        if (in_calib_2_b[det_nr-1] != NULL) 
            hdrl_image_delete(in_calib_2_b[det_nr-1]) ;
        if (in_calib_3_b[det_nr-1] != NULL) 
            hdrl_image_delete(in_calib_3_b[det_nr-1]) ;
        if (in_calib_4_b[det_nr-1] != NULL) 
            hdrl_image_delete(in_calib_4_b[det_nr-1]) ;
        if (tw_1u_a[det_nr-1] != NULL) cpl_table_delete(tw_1u_a[det_nr-1]) ;
        if (tw_1d_a[det_nr-1] != NULL) cpl_table_delete(tw_1d_a[det_nr-1]) ;
        if (tw_2u_a[det_nr-1] != NULL) cpl_table_delete(tw_2u_a[det_nr-1]) ;
        if (tw_2d_a[det_nr-1] != NULL) cpl_table_delete(tw_2d_a[det_nr-1]) ;
        if (tw_3u_a[det_nr-1] != NULL) cpl_table_delete(tw_3u_a[det_nr-1]) ;
        if (tw_3d_a[det_nr-1] != NULL) cpl_table_delete(tw_3d_a[det_nr-1]) ;
        if (tw_4u_a[det_nr-1] != NULL) cpl_table_delete(tw_4u_a[det_nr-1]) ;
        if (tw_4d_a[det_nr-1] != NULL) cpl_table_delete(tw_4d_a[det_nr-1]) ;
        if (tw_1u_b[det_nr-1] != NULL) cpl_table_delete(tw_1u_b[det_nr-1]) ;
        if (tw_1d_b[det_nr-1] != NULL) cpl_table_delete(tw_1d_b[det_nr-1]) ;
        if (tw_2u_b[det_nr-1] != NULL) cpl_table_delete(tw_2u_b[det_nr-1]) ;
        if (tw_2d_b[det_nr-1] != NULL) cpl_table_delete(tw_2d_b[det_nr-1]) ;
        if (tw_3u_b[det_nr-1] != NULL) cpl_table_delete(tw_3u_b[det_nr-1]) ;
        if (tw_3d_b[det_nr-1] != NULL) cpl_table_delete(tw_3d_b[det_nr-1]) ;
        if (tw_4u_b[det_nr-1] != NULL) cpl_table_delete(tw_4u_b[det_nr-1]) ;
        if (tw_4d_b[det_nr-1] != NULL) cpl_table_delete(tw_4d_b[det_nr-1]) ;
 
        if (extract1D_1u_a[det_nr-1] != NULL)
            cpl_table_delete(extract1D_1u_a[det_nr-1]) ;
        if (extract1D_1d_a[det_nr-1] != NULL)
            cpl_table_delete(extract1D_1d_a[det_nr-1]) ;
        if (extract1D_2u_a[det_nr-1] != NULL)
            cpl_table_delete(extract1D_2u_a[det_nr-1]) ;
        if (extract1D_2d_a[det_nr-1] != NULL)
            cpl_table_delete(extract1D_2d_a[det_nr-1]) ;
        if (extract1D_3u_a[det_nr-1] != NULL)
            cpl_table_delete(extract1D_3u_a[det_nr-1]) ;
        if (extract1D_3d_a[det_nr-1] != NULL)
            cpl_table_delete(extract1D_3d_a[det_nr-1]) ;
        if (extract1D_4u_a[det_nr-1] != NULL)
            cpl_table_delete(extract1D_4u_a[det_nr-1]) ;
        if (extract1D_4d_a[det_nr-1] != NULL)
            cpl_table_delete(extract1D_4d_a[det_nr-1]) ;
        if (extract1D_1u_b[det_nr-1] != NULL)
            cpl_table_delete(extract1D_1u_b[det_nr-1]) ;
        if (extract1D_1d_b[det_nr-1] != NULL)
            cpl_table_delete(extract1D_1d_b[det_nr-1]) ;
        if (extract1D_2u_b[det_nr-1] != NULL)
            cpl_table_delete(extract1D_2u_b[det_nr-1]) ;
        if (extract1D_2d_b[det_nr-1] != NULL)
            cpl_table_delete(extract1D_2d_b[det_nr-1]) ;
        if (extract1D_3u_b[det_nr-1] != NULL)
            cpl_table_delete(extract1D_3u_b[det_nr-1]) ;
        if (extract1D_3d_b[det_nr-1] != NULL)
            cpl_table_delete(extract1D_3d_b[det_nr-1]) ;
        if (extract1D_4u_b[det_nr-1] != NULL)
            cpl_table_delete(extract1D_4u_b[det_nr-1]) ;
        if (extract1D_4d_b[det_nr-1] != NULL)
            cpl_table_delete(extract1D_4d_b[det_nr-1]) ;
        if (pol_speca[det_nr-1] != NULL) cpl_table_delete(pol_speca[det_nr-1]) ;
        if (pol_specb[det_nr-1] != NULL) cpl_table_delete(pol_specb[det_nr-1]) ;
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
  @param blaze_frame            Associated Blaze
  @param subtract_nolight_rows
  @param cosmics                Flag to correct for cosmics
  @param extract_oversample     Extraction related
  @param extract_swath_width    Extraction related
  @param extract_height         Extraction related
  @param extract_smooth_slit    Extraction related
  @param extract_smooth_spec    Extraction related
  @param save_group             Group for which extra products will be saved
  @param reduce_det             The detector to compute
  @param in_calib_[1-4]_a       [out] the calibrated image (A)
  @param trace_wave_[1u-4d]_a   [out] the TW tables used for extraction (A)
  @param extract1D_[1u-4d]_a    [out] the extracted table (A)
  @param pol_speca              [out] polarimetry spectrum (A)
  @param in_calib_[1-4]_b       [out] the calibrated image (B)
  @param trace_wave_[1u-4d]_b   [out] the TW tables used for extraction (B)
  @param extract1D_[1u-4d]_b    [out] the extracted table (B)
  @param pol_specb              [out] polarimetry spectrum (B)
  @param ext_plista             [out] the header for saving the products (A)
  @param ext_plistb             [out] the header for saving the products (B)
  @return  0 if ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
static int cr2res_obs_pol_reduce(
        const cpl_frameset  *   rawframes,
        const cpl_frame     *   trace_wave_frame,
        const cpl_frame     *   detlin_frame,
        const cpl_frame     *   master_dark_frame,
        const cpl_frame     *   master_flat_frame,
        const cpl_frame     *   bpm_frame,
        const cpl_frame     *   blaze_frame,
        int                     subtract_nolight_rows,
        int                     subtract_interorder_column,
        int                     cosmics,
        int                     extract_oversample,
        int                     extract_swath_width,
        int                     extract_height,
        double                  extract_smooth_slit,
        double                  extract_smooth_spec,
        int                     save_group,
        int                     reduce_det,
        hdrl_image          **  in_calib_1_a,
        hdrl_image          **  in_calib_2_a,
        hdrl_image          **  in_calib_3_a,
        hdrl_image          **  in_calib_4_a,
        cpl_table           **  trace_wave_1u_a,
        cpl_table           **  trace_wave_1d_a,
        cpl_table           **  trace_wave_2u_a,
        cpl_table           **  trace_wave_2d_a,
        cpl_table           **  trace_wave_3u_a,
        cpl_table           **  trace_wave_3d_a,
        cpl_table           **  trace_wave_4u_a,
        cpl_table           **  trace_wave_4d_a,
        cpl_table           **  extract1D_1u_a,
        cpl_table           **  extract1D_1d_a,
        cpl_table           **  extract1D_2u_a,
        cpl_table           **  extract1D_2d_a,
        cpl_table           **  extract1D_3u_a,
        cpl_table           **  extract1D_3d_a,
        cpl_table           **  extract1D_4u_a,
        cpl_table           **  extract1D_4d_a,
        cpl_table           **  pol_speca,
        hdrl_image          **  in_calib_1_b,
        hdrl_image          **  in_calib_2_b,
        hdrl_image          **  in_calib_3_b,
        hdrl_image          **  in_calib_4_b,
        cpl_table           **  trace_wave_1u_b,
        cpl_table           **  trace_wave_1d_b,
        cpl_table           **  trace_wave_2u_b,
        cpl_table           **  trace_wave_2d_b,
        cpl_table           **  trace_wave_3u_b,
        cpl_table           **  trace_wave_3d_b,
        cpl_table           **  trace_wave_4u_b,
        cpl_table           **  trace_wave_4d_b,
        cpl_table           **  extract1D_1u_b,
        cpl_table           **  extract1D_1d_b,
        cpl_table           **  extract1D_2u_b,
        cpl_table           **  extract1D_2d_b,
        cpl_table           **  extract1D_3u_b,
        cpl_table           **  extract1D_3d_b,
        cpl_table           **  extract1D_4u_b,
        cpl_table           **  extract1D_4d_b,
        cpl_table           **  pol_specb,
        cpl_propertylist    **  ext_plista,
        cpl_propertylist    **  ext_plistb)
{
    cpl_frameset        *   rawframes_a ;
    cpl_frameset        *   rawframes_b ;
    cr2res_nodding_pos  *   nod_positions ;
    hdrl_image          **  in_calib_a_loc ;
    hdrl_image          **  in_calib_b_loc ;
    cpl_table           **  trace_wave_a_loc ;
    cpl_table           **  trace_wave_b_loc ;
    cpl_table           **  extract1D_a_loc ;
    cpl_table           **  extract1D_b_loc ;
    cpl_table           *   pol_speca_loc ;
    cpl_table           *   pol_specb_loc ;
    cpl_propertylist    *   ext_plista_loc ;
    cpl_propertylist    *   ext_plistb_loc ;
    cpl_size                ngroups ;

    /* Check Inputs */
    if (pol_speca == NULL || pol_specb == NULL || ext_plista == NULL || 
            ext_plistb == NULL || rawframes == NULL || 
            trace_wave_frame == NULL) return -1 ;

    if (cr2res_obs_pol_check_inputs_validity(rawframes, &ngroups) != 1) {
        cpl_msg_error(__func__, "Validity Checks failed") ;
        return -1 ;
    }

    /* Initialize */
    *in_calib_1_a = *in_calib_2_a = *in_calib_3_a = *in_calib_4_a = NULL ;
    *trace_wave_1u_a = *trace_wave_1d_a = *trace_wave_2u_a = 
        *trace_wave_2d_a = *trace_wave_3u_a = *trace_wave_3d_a =
        *trace_wave_4u_a = *trace_wave_4d_a = NULL ;
    *extract1D_1u_a = *extract1D_1d_a = *extract1D_2u_a =
        *extract1D_2d_a = *extract1D_3u_a = *extract1D_3d_a =
        *extract1D_4u_a = *extract1D_4d_a = NULL ;
    *in_calib_1_b = *in_calib_2_b = *in_calib_3_b = *in_calib_4_b = NULL ;
    *trace_wave_1u_b = *trace_wave_1d_b = *trace_wave_2u_b = 
        *trace_wave_2d_b = *trace_wave_3u_b = *trace_wave_3d_b =
        *trace_wave_4u_b = *trace_wave_4d_b = NULL ;
    *extract1D_1u_b = *extract1D_1d_b = *extract1D_2u_b =
        *extract1D_2d_b = *extract1D_3u_b = *extract1D_3d_b =
        *extract1D_4u_b = *extract1D_4d_b = NULL ;

    /* Get the Nodding positions */
    nod_positions = cr2res_nodding_read_positions(rawframes) ;

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
    if (cr2res_obs_pol_reduce_one(rawframes_a, rawframes_b,
                trace_wave_frame, detlin_frame, master_dark_frame, 
                master_flat_frame, bpm_frame, blaze_frame, 
                subtract_nolight_rows, subtract_interorder_column, cosmics, 
                extract_oversample, extract_swath_width, extract_height, 
                extract_smooth_slit, extract_smooth_spec, save_group,
                reduce_det, 
                &in_calib_a_loc, &trace_wave_a_loc, &extract1D_a_loc, 
                &pol_speca_loc, &ext_plista_loc) == -1) {
        cpl_msg_warning(__func__, "Failed to Reduce A nodding frames") ;
    }
    cpl_msg_indent_less() ;

    /* Reduce B position */
    cpl_msg_info(__func__, "Compute Polarimetry for nodding B position") ;
    cpl_msg_indent_more() ;
    if (cr2res_obs_pol_reduce_one(rawframes_b, rawframes_a,
                trace_wave_frame, detlin_frame, master_dark_frame, 
                master_flat_frame, bpm_frame, blaze_frame,
                subtract_nolight_rows, subtract_interorder_column, cosmics,
                extract_oversample, extract_swath_width, extract_height,
                extract_smooth_slit, extract_smooth_spec, save_group, 
                reduce_det, 
                &in_calib_b_loc, &trace_wave_b_loc, &extract1D_b_loc, 
                &pol_specb_loc, &ext_plistb_loc) == -1) {
        cpl_msg_warning(__func__, "Failed to Reduce B nodding frames") ;
    }
    cpl_msg_indent_less() ;
    if (rawframes_a != NULL) cpl_frameset_delete(rawframes_a);
    if (rawframes_b != NULL) cpl_frameset_delete(rawframes_b);

    /* Return */
    if (in_calib_a_loc != NULL) {
        *in_calib_1_a = in_calib_a_loc[0] ;
        *in_calib_2_a = in_calib_a_loc[1] ;
        *in_calib_3_a = in_calib_a_loc[2] ;
        *in_calib_4_a = in_calib_a_loc[3] ;
        cpl_free(in_calib_a_loc) ;
    }
    if (trace_wave_a_loc != NULL) {
        *trace_wave_1u_a = trace_wave_a_loc[0] ;
        *trace_wave_1d_a = trace_wave_a_loc[1] ;
        *trace_wave_2u_a = trace_wave_a_loc[2] ;
        *trace_wave_2d_a = trace_wave_a_loc[3] ;
        *trace_wave_3u_a = trace_wave_a_loc[4] ;
        *trace_wave_3d_a = trace_wave_a_loc[5] ;
        *trace_wave_4u_a = trace_wave_a_loc[6] ;
        *trace_wave_4d_a = trace_wave_a_loc[7] ;
        cpl_free(trace_wave_a_loc) ;
    } 
    if (extract1D_a_loc != NULL) {
        *extract1D_1u_a = extract1D_a_loc[0] ;
        *extract1D_1d_a = extract1D_a_loc[1] ;
        *extract1D_2u_a = extract1D_a_loc[2] ;
        *extract1D_2d_a = extract1D_a_loc[3] ;
        *extract1D_3u_a = extract1D_a_loc[4] ;
        *extract1D_3d_a = extract1D_a_loc[5] ;
        *extract1D_4u_a = extract1D_a_loc[6] ;
        *extract1D_4d_a = extract1D_a_loc[7] ;
        cpl_free(extract1D_a_loc) ;
    }
    if (in_calib_b_loc != NULL) {
        *in_calib_1_b = in_calib_b_loc[0] ;
        *in_calib_2_b = in_calib_b_loc[1] ;
        *in_calib_3_b = in_calib_b_loc[2] ;
        *in_calib_4_b = in_calib_b_loc[3] ;
        cpl_free(in_calib_b_loc) ;
    }
    if (trace_wave_b_loc != NULL) {
        *trace_wave_1u_b = trace_wave_b_loc[0] ;
        *trace_wave_1d_b = trace_wave_b_loc[1] ;
        *trace_wave_2u_b = trace_wave_b_loc[2] ;
        *trace_wave_2d_b = trace_wave_b_loc[3] ;
        *trace_wave_3u_b = trace_wave_b_loc[4] ;
        *trace_wave_3d_b = trace_wave_b_loc[5] ;
        *trace_wave_4u_b = trace_wave_b_loc[6] ;
        *trace_wave_4d_b = trace_wave_b_loc[7] ;
        cpl_free(trace_wave_b_loc) ;
    }
    if (extract1D_b_loc != NULL) {
        *extract1D_1u_b = extract1D_b_loc[0] ;
        *extract1D_1d_b = extract1D_b_loc[1] ;
        *extract1D_2u_b = extract1D_b_loc[2] ;
        *extract1D_2d_b = extract1D_b_loc[3] ;
        *extract1D_3u_b = extract1D_b_loc[4] ;
        *extract1D_3d_b = extract1D_b_loc[5] ;
        *extract1D_4u_b = extract1D_b_loc[6] ;
        *extract1D_4d_b = extract1D_b_loc[7] ;
        cpl_free(extract1D_b_loc) ;
    }
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
  @param blaze_frame            Associated Blaze frame
  @param subtract_nolight_rows
  @param cosmics                Flag to correct for cosmics
  @param extract_oversample     Extraction related
  @param extract_swath_width    Extraction related
  @param extract_height         Extraction related
  @param extract_smooth_slit    Extraction related
  @param extract_smooth_spec    Extraction related
  @param save_group             Group for which extra products will be saved
  @param reduce_det             The detector to compute
  @param in_calib   [out] the array of the 4 calibrated images of a group 
  @param tw         [out] the array of the 8 TW tables (1u, 1d, ..., 4d)
  @param extract1D  [out] the array of the 8 extracted tables (1u, ... 4d)
  @param pol_spec   [out] polarimetry spectrum
  @param ext_plist  [out] the header for saving the products
  @return  0 if ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
static int cr2res_obs_pol_reduce_one(
        const cpl_frameset  *   rawframes,
        const cpl_frameset  *   raw_background_frames,
        const cpl_frame     *   trace_wave_frame,
        const cpl_frame     *   detlin_frame,
        const cpl_frame     *   master_dark_frame,
        const cpl_frame     *   master_flat_frame,
        const cpl_frame     *   bpm_frame,
        const cpl_frame     *   blaze_frame,
        int                     subtract_nolight_rows,
        int                     subtract_interorder_column,
        int                     cosmics,
        int                     extract_oversample,
        int                     extract_swath_width,
        int                     extract_height,
        double                  extract_smooth_slit,
        double                  extract_smooth_spec,
        int                     save_group,
        int                     reduce_det,
        hdrl_image          *** in_calib,
        cpl_table           *** tw,
        cpl_table           *** extract1D,
        cpl_table           **  pol_spec,
        cpl_propertylist    **  ext_plist)
{
    cpl_vector          *   dits ;
    cpl_vector          *   ndits ;
    cr2res_decker       *   decker_positions ;
    hdrl_imagelist      *   in ;
    hdrl_imagelist      *   in_calib_loc ;
    cpl_image           *   contrib;
    hdrl_image          *   backgr;
    cpl_table           *   blaze_table ;
    cpl_table           *   trace_wave ;
    //cpl_table           *   trace_wave_corrected ;
    cpl_table           **  trace_wave_extract ;
    const char          *   fname ;
    char                *   decker_name ;
    cpl_table           *   slit_func ;
    hdrl_image          *   model_master ;
    hdrl_image          **  input_images ;
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
    cpl_propertylist    *   plist ;
    cpl_propertylist    *   ext_plist_loc ;
    const char          *   first_fname ;
    cpl_size                nframes, nspec_group, spec_size ;
    int                     ngroups, i, j, k, l, o, norders, frame_idx ;
    double                  gain, error_factor, blaze_norm ;
    int                 *   order_idx_values ;
    int                     nb_order_idx_values, order_zp;

    /* TODO, make parameters */
    int extract_niter = 10;
    double extract_kappa = 10;

    /* Check Inputs */
    if (in_calib == NULL || tw == NULL || extract1D == NULL || 
            pol_spec == NULL || ext_plist == NULL || rawframes == NULL || 
            trace_wave_frame == NULL)
        return -1 ;

    /* Initialise */
    *in_calib = NULL ;
    *tw = NULL ;
    *extract1D = NULL ;
    *pol_spec = NULL ;
    *ext_plist = NULL ;
    first_fname = cpl_frame_get_filename(
            cpl_frameset_get_position_const(rawframes, 0)) ;

    /* Get the Gain */
    if (reduce_det == 1) gain = CR2RES_GAIN_CHIP1 ;
    else if (reduce_det == 2) gain = CR2RES_GAIN_CHIP2 ;
    else if (reduce_det == 3) gain = CR2RES_GAIN_CHIP3 ;
    else {
        cpl_msg_error(__func__, "Failed to get the Gain value") ;
        return -1 ;
    }

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
    
    /* Load the NDITs */
    ndits = cr2res_io_read_ndits(rawframes) ;

    /* Set the error factor. Note that we don't multiply by 4 because */
    /* each the beams get extracted from each frame, then errors propagated on.*/
    error_factor = gain * cpl_vector_get(ndits, 0);

    for (i=0; i<cpl_vector_get_size(ndits); i++){
        if (cpl_vector_get(ndits,i) != cpl_vector_get(ndits, 0))
            cpl_msg_warning(__func__, "Raw frames have different NDIT! "
                "Error spectrum will likely be scaled incorrectly.");
    }

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
        if (ndits != NULL) cpl_vector_delete(ndits) ;
        cpl_free(decker_positions) ;
        return -1 ;
    }
    if (hdrl_imagelist_get_size(in) != cpl_frameset_get_size(rawframes)) {
        cpl_msg_error(__func__, "Inconsistent number of loaded images") ;
        if (dits != NULL) cpl_vector_delete(dits) ;
        if (ndits != NULL) cpl_vector_delete(ndits) ;
        cpl_free(decker_positions) ;
        hdrl_imagelist_delete(in) ;
        return -1 ;
    }

    /* Calibrate the images */
    cpl_msg_info(__func__, "Apply the calibrations") ;
    cpl_msg_indent_more() ;
    if ((in_calib_loc = cr2res_calib_imagelist(in, reduce_det, 0,
            subtract_nolight_rows, subtract_interorder_column, cosmics, 
            master_flat_frame, master_dark_frame, bpm_frame, detlin_frame, 
            dits, ndits))==NULL) {
        cpl_msg_error(__func__, "Failed to apply the calibrations") ;
        cpl_msg_indent_less() ;
        if (dits != NULL) cpl_vector_delete(dits) ;
        if (ndits != NULL) cpl_vector_delete(ndits) ;
        cpl_free(decker_positions) ;
        hdrl_imagelist_delete(in) ;
        return -1 ;
    }
    cpl_msg_indent_less() ;
    hdrl_imagelist_delete(in) ;
    if (dits != NULL) {
        cpl_vector_delete(dits) ;
        dits = NULL ;
    }
    if (ndits != NULL) cpl_vector_delete(ndits) ;

    /* Apply Background Correction using the other nodding position */
    if (cpl_frameset_get_size(raw_background_frames) > 1 ) {

        hdrl_imagelist      *   in_backgr ;
        cpl_msg_info(__func__, "Apply the background Correction") ;
        cpl_msg_indent_more() ;

        /* Load image list for BACKGROUND frames */
        cpl_msg_info(__func__, "Load the background images") ;
        if ((in = cr2res_io_load_image_list_from_set(raw_background_frames, 
                        reduce_det)) == NULL) {
            cpl_msg_error(__func__, "Cannot load background images") ;
            cpl_msg_indent_less() ;
            cpl_free(decker_positions) ;
            hdrl_imagelist_delete(in_calib_loc) ;
            return -1 ;
        }
        /* Load the DITs & NDITs */
        if (master_dark_frame != NULL)  dits = cr2res_io_read_dits(rawframes) ;
        else                            dits = NULL ;
        ndits = cr2res_io_read_ndits(rawframes) ;

        /* Calibrate the background same as science images */
        cpl_msg_info(__func__, "Apply the calibrations to background") ;
        cpl_msg_indent_more() ;
        if ((in_backgr = cr2res_calib_imagelist(in, reduce_det, 0, 
                        subtract_nolight_rows, 1, cosmics, master_flat_frame, 
                        master_dark_frame, bpm_frame, detlin_frame, 
                        dits, ndits)) == NULL) {
            cpl_msg_error(__func__,
                            "Failed to apply the calibrations to background") ;
            cpl_msg_indent_less() ;
            cpl_msg_indent_less() ;
            if (dits != NULL) cpl_vector_delete(dits) ;
            if (ndits != NULL) cpl_vector_delete(ndits) ;
            cpl_free(decker_positions) ;
            hdrl_imagelist_delete(in) ;
            hdrl_imagelist_delete(in_calib_loc) ;
            return -1 ;
        }
        cpl_msg_indent_less() ;
        hdrl_imagelist_delete(in) ;

        /* Collapse background images */
        if (hdrl_imagelist_collapse_mean(in_backgr, &backgr, &contrib) != \
                        CPL_ERROR_NONE) {
            cpl_msg_error(__func__, "Failed to collapse background") ;
            cpl_msg_indent_less() ;
            if (dits != NULL) cpl_vector_delete(dits) ;
            if (ndits != NULL) cpl_vector_delete(ndits) ;
            cpl_free(decker_positions) ;
            hdrl_imagelist_delete(in_calib_loc) ;
            hdrl_imagelist_delete(in_backgr) ;
            return -1;
        }
        if (contrib != NULL) cpl_image_delete(contrib) ;
        if (in_backgr != NULL) hdrl_imagelist_delete(in_backgr) ;

        /* Subtract background from calibrated images */
        if (hdrl_imagelist_sub_image(in_calib_loc, backgr) != CPL_ERROR_NONE) {
            cpl_msg_error(__func__, "Failed to subtract background") ;
            cpl_msg_indent_less() ;
            if (dits != NULL) cpl_vector_delete(dits) ;
            if (ndits != NULL) cpl_vector_delete(ndits) ;
            cpl_free(decker_positions) ;
            hdrl_imagelist_delete(in_calib_loc);
            hdrl_image_delete(backgr);
            return -1;
        }
        hdrl_image_delete(backgr);
        cpl_msg_indent_less() ;
    } else {
        cpl_msg_warning(__func__, "No background subtraction");
    }
    if (dits != NULL) cpl_vector_delete(dits) ;
    if (ndits != NULL) cpl_vector_delete(ndits) ;

    /* Load the trace wave */
    cpl_msg_info(__func__, "Load the TRACE WAVE") ;
    if ((trace_wave = cr2res_io_load_TRACE_WAVE(cpl_frame_get_filename(
                        trace_wave_frame), reduce_det)) == NULL) {
        cpl_msg_error(__func__, "Failed to Load the traces file") ;
        cpl_free(decker_positions) ;
        hdrl_imagelist_delete(in_calib_loc) ;
        return -1 ;
    }

    /* Correct trace_wave with some provided raw flats */
    cpl_msg_warning(__func__, "TRACE ADJUST NOT YET IMPLEMENTED") ;
    /*if (raw_flat_frames != NULL) {
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
    }*/

    /* Load Blaze */
    blaze_table = NULL ;
    blaze_norm = 0 ;
    if (blaze_frame != NULL) {
        cpl_msg_info(__func__, "Load the BLAZE") ;
        if ((blaze_table = cr2res_io_load_EXTRACT_1D(cpl_frame_get_filename(
                            blaze_frame), reduce_det)) == NULL) {
            cpl_msg_error(__func__, "Failed to Load the Blaze file") ;
            cpl_free(decker_positions) ;
            hdrl_imagelist_delete(in_calib_loc) ;
            cpl_table_delete(trace_wave) ;
            return -1 ;
        }
        cpl_propertylist    *   blaze_plist ;
        blaze_plist = cpl_propertylist_load_regexp(cpl_frame_get_filename(
                            blaze_frame),0,CR2RES_HEADER_QC_BLAZE_NORM,0);
        if(cpl_propertylist_get_size(blaze_plist)>0){
            blaze_norm = cpl_propertylist_get_double(blaze_plist, CR2RES_HEADER_QC_BLAZE_NORM);
        }
        else {
            blaze_norm = -1;
            cpl_msg_warning(__func__, "QC BLAZE NORM value not found, reverting to per trace normalization") ;
        }
        if(blaze_plist!=NULL){
            cpl_propertylist_delete(blaze_plist);
        }
    }

    /* Compute the number of groups */
    ngroups = nframes/CR2RES_POLARIMETRY_GROUP_SIZE ;
    nspec_group = 2*CR2RES_POLARIMETRY_GROUP_SIZE ;

    /* Allocate pol_spec_group containers */
    pol_spec_one_group = cpl_malloc(ngroups * sizeof(cpl_table*)) ;
    for (i = 0; i < ngroups; i++) pol_spec_one_group[i] = NULL;

    /* Loop on the groups */
    for (i = 0; i < ngroups; i++) {
        int *pol_sorting;
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

        /* Container for the current group calibrated images 1 2 3 4 */ 
        input_images = cpl_malloc(CR2RES_POLARIMETRY_GROUP_SIZE * 
                sizeof(hdrl_image*));
        for (j = 0; j < CR2RES_POLARIMETRY_GROUP_SIZE; j++)
            input_images[j] = NULL ;

        /* Container for extracted tables: 1u, 1d, 2u, 2d, 3u, 3d, 4u, 4d */
        extract_1d = cpl_malloc(nspec_group * sizeof(cpl_table*));
        trace_wave_extract = cpl_malloc(nspec_group * sizeof(cpl_table*));
        for (j = 0; j < nspec_group; j++) {
            extract_1d[j] = NULL;
            trace_wave_extract[j] = NULL;
        }

        /* Loop on the frames */
        for (j=0 ; j<CR2RES_POLARIMETRY_GROUP_SIZE ; j++) {
            frame_idx = i * CR2RES_POLARIMETRY_GROUP_SIZE + pol_sorting[j] ;
            fname = cpl_frame_get_filename(cpl_frameset_get_position_const(
                        rawframes,frame_idx)) ;

            /* Store the input image */
            input_images[j] = hdrl_image_duplicate(
                    hdrl_imagelist_get_const(in_calib_loc, frame_idx)) ;

            /* Extract up */
            decker_name = cr2res_decker_print_position(
                    decker_positions[frame_idx]) ;
            cpl_msg_info(__func__, 
                    "Extract Up Spectrum from %s (Det %d / Decker %s)", 
                    fname, reduce_det, decker_name) ;
            cpl_free(decker_name) ;
            cpl_msg_indent_more();

            trace_wave_extract[2 * j] = cr2res_pol_get_beam_trace(trace_wave, decker_positions[frame_idx], 1);

            /* Execute the extraction */
            cpl_msg_info(__func__, "Spectra Extraction") ;
            if (cr2res_extract_traces(input_images[j], trace_wave_extract[2*j],
                    NULL, blaze_table, blaze_norm, -1, -1, CR2RES_EXTR_OPT_CURV,
                    extract_height, extract_swath_width, extract_oversample, 
                    extract_smooth_slit, extract_smooth_spec,
                    extract_niter, extract_kappa, error_factor, 0, 0, 0, 
                    &(extract_1d[2*j]), &slit_func, &model_master) == -1) {
                cpl_msg_error(__func__, "Failed Extraction") ;
                extract_1d[2*j] = NULL ;
                cpl_table_delete(trace_wave_extract[2*j]) ;
                trace_wave_extract[2*j] = NULL ;
            } else {
                cpl_table_delete(slit_func) ;
                hdrl_image_delete(model_master) ;
            }
            cpl_msg_indent_less() ;

            /* Extract Down */
            decker_name = cr2res_decker_print_position(
                    decker_positions[frame_idx]) ;
            cpl_msg_info(__func__, 
                    "Extract Down Spectrum from %s (Det %d / Decker %s)", 
                    fname, reduce_det, decker_name) ;
            cpl_free(decker_name) ;
            cpl_msg_indent_more() ;
           
            trace_wave_extract[2*j+1] = cr2res_pol_get_beam_trace(trace_wave, decker_positions[frame_idx], 2);

            /* Execute the extraction */
            cpl_msg_info(__func__, "Spectra Extraction") ;
            if (cr2res_extract_traces(input_images[j],
                    trace_wave_extract[2*j+1], NULL, blaze_table, blaze_norm, -1, -1, 
                    CR2RES_EXTR_OPT_CURV, extract_height, 
                    extract_swath_width, extract_oversample,
                    extract_smooth_slit, extract_smooth_spec,
                    extract_niter, extract_kappa, error_factor, 0, 0, 0, 
                    &(extract_1d[2*j+1]), &slit_func, &model_master)== -1) {
                cpl_msg_error(__func__, "Failed Extraction") ;
                extract_1d[2*j+1] = NULL ;
                cpl_table_delete(trace_wave_extract[2*j+1]) ;
                trace_wave_extract[2*j+1] = NULL ;
            } else {
                cpl_table_delete(slit_func) ;
                hdrl_image_delete(model_master) ;
            }
            cpl_msg_indent_less() ;
        }
        cpl_free(pol_sorting) ;

        /* How many orders */
        if ((orders = cr2res_obs_pol_get_order_numbers(
                        (const cpl_table **)extract_1d, 
                        nspec_group, 
                        &norders)) == NULL) {
            cpl_msg_warning(__func__, "No Order found in the extracted tables");
            for (j = 0; j < CR2RES_POLARIMETRY_GROUP_SIZE; j++) {
                hdrl_image_delete(input_images[j]) ;
                input_images[j] = NULL ;
            }
            cpl_free(input_images) ;
            for (j=0 ; j<nspec_group ; j++) {
                if (extract_1d[j] != NULL) 
                    cpl_table_delete(extract_1d[j]) ;
                if (trace_wave_extract[j] != NULL) 
                    cpl_table_delete(trace_wave_extract[j]) ;
            }
            cpl_free(extract_1d) ;
            cpl_free(trace_wave_extract) ;
            continue ;
        }
        cpl_msg_debug(__func__, "%d different orders found", norders) ;

        /* Allocate data containers */
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

        intens = cpl_malloc(nspec_group * sizeof(cpl_vector*));
        wl = cpl_malloc(nspec_group * sizeof(cpl_vector*));
        errors = cpl_malloc(nspec_group * sizeof(cpl_vector*));
        for (k = 0; k < nspec_group; k++){
            intens[k] = NULL;
            wl[k] = NULL;
            errors[k] = NULL;
        }

        /* Loop on the orders */
        for (o=0 ; o<norders ; o++) {
            cpl_msg_info(__func__, "Compute Polarimetry for order %d",
                    orders[o]) ;
            /* Get the inputs for the demod functions calls */
            for (k=0 ; k<nspec_group ; k++) {
                spec_size = cpl_table_get_nrow(extract_1d[k]) ;
                /* Get the SPEC for this order/trace 1 */
                colname = cr2res_dfs_SPEC_colname(orders[o], 1) ;
                if (cpl_table_has_valid(extract_1d[k], colname)){
                    if ((pcol_data = cpl_table_get_data_double_const(
                            extract_1d[k], colname)) == NULL) {
                        cpl_error_reset() ;
                    } else {
                        intens[k] = cpl_vector_new(spec_size) ;
                        pvec_data = cpl_vector_get_data(intens[k]) ;
                        for (l=0 ; l<spec_size ; l++) 
                            pvec_data[l] = pcol_data[l] ;
                    }
                }
                cpl_free(colname) ;

                /* Get the WAVELENGTH for this order/trace 1 */
                colname = cr2res_dfs_WAVELENGTH_colname(orders[o], 1) ;
                if (cpl_table_has_valid(extract_1d[k], colname)){
                    if ((pcol_data = cpl_table_get_data_double_const(
                            extract_1d[k], colname)) == NULL) {
                        cpl_error_reset() ;
                    } else {
                        wl[k] = cpl_vector_new(spec_size) ;
                        pvec_data = cpl_vector_get_data(wl[k]) ;
                        for (l=0 ; l<spec_size ; l++) 
                            pvec_data[l] = pcol_data[l] ;
                    }
                }
                cpl_free(colname) ;

                /* Get the ERROR for this order/trace 1 */
                colname = cr2res_dfs_SPEC_ERR_colname(orders[o], 1) ;
                if (cpl_table_has_valid(extract_1d[k], colname)){
                    if ((pcol_data = cpl_table_get_data_double_const(
                            extract_1d[k], colname)) == NULL) {
                        cpl_error_reset() ;
                    } else {
                        errors[k] = cpl_vector_new(spec_size) ;
                        pvec_data = cpl_vector_get_data(errors[k]) ;
                        for (l=0 ; l<spec_size ; l++) 
                            pvec_data[l] = pcol_data[l] ;
                    }
                }
                cpl_free(colname) ;
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
                if (intens[k] != NULL) {
                    cpl_vector_delete(intens[k]) ;
                    intens[k] = NULL;
                }
                if (wl[k] != NULL) {
                    cpl_vector_delete(wl[k]) ;
                    wl[k] = NULL;
                }
                if (errors[k] != NULL) {
                    cpl_vector_delete(errors[k]) ;
                    errors[k] = NULL;
                }
            }
        }
        cpl_free(intens) ;
        cpl_free(wl) ;
        cpl_free(errors) ;

        /* Save the extra products for this group */
        if (i+1==save_group) {
            *in_calib = input_images ;
            *tw = trace_wave_extract ; 
            *extract1D = extract_1d ; 
        } else {
            /* Free the input images */
            for (j = 0; j < CR2RES_POLARIMETRY_GROUP_SIZE; j++) 
                if (input_images[j] != NULL) 
                    hdrl_image_delete(input_images[j]) ;
            cpl_free(input_images) ;
            /* Free the extracted tables / the TW */
            for (j=0 ; j<nspec_group ; j++) {
                if (trace_wave_extract[j] != NULL) 
                    cpl_table_delete(trace_wave_extract[j]) ;
                if (extract_1d[j] != NULL) cpl_table_delete(extract_1d[j]) ;
            }
            cpl_free(trace_wave_extract) ;
            cpl_free(extract_1d) ;
        }

        /* Create the pol_spec table */
        cpl_msg_info(__func__, "Create the POL_SPEC table for this group");
        pol_spec_one_group[i] = cr2res_pol_POL_SPEC_create(orders, demod_wl, 
                demod_stokes, demod_null, demod_intens, norders) ;

        /* Deallocate */
        for (o=0 ; o<norders ; o++) {
            if (demod_wl[o] != NULL) cpl_vector_delete(demod_wl[o]) ;
            if (demod_stokes[o] != NULL) cpl_bivector_delete(demod_stokes[o]) ;
            if (demod_null[o] != NULL) cpl_bivector_delete(demod_null[o]) ;
            if (demod_intens[o] != NULL) cpl_bivector_delete(demod_intens[o]) ; 
        }
        cpl_free(demod_wl) ;
        cpl_free(demod_stokes) ;
        cpl_free(demod_null) ;
        cpl_free(demod_intens) ;
        cpl_free(orders) ;
        cpl_msg_indent_less() ;
    }
    if (blaze_table != NULL) cpl_table_delete(blaze_table) ;
    cpl_free(decker_positions) ;
    hdrl_imagelist_delete(in_calib_loc) ;

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
        cpl_table_delete(trace_wave) ;
        for (j=0 ; j<nspec_group ; j++) {
        cpl_table **temp_extract1D = extract1D[j];
        if (temp_extract1D != NULL && *temp_extract1D != NULL) 
            cpl_table_delete(*temp_extract1D);
        cpl_table **temp_tw = tw[j];
            if (temp_tw != NULL && *temp_tw != NULL)
                cpl_table_delete(*temp_tw) ;
        }
        cpl_free(*tw) ; *tw = NULL ;
        cpl_free(*extract1D) ; *extract1D = NULL ;
        return -1 ;
    }

    /* Store the extension header for product saving */
    ext_plist_loc = cpl_propertylist_load(first_fname,
            cr2res_io_get_ext_idx(first_fname, reduce_det, 1)) ;

    cpl_msg_info(__func__,"Order zero-point: %d", order_zp );
    /* Real Orders in QCs */
    if (order_zp > 0) {
        /* Get the order numbers from the TW rows */
        order_idx_values = cr2res_trace_get_order_idx_values(trace_wave,
                &nb_order_idx_values);
        /* Compute the Real Order numbers and store them in QCs */
        for (i = 0; i < nb_order_idx_values; i++) {
            char *key_name;
            char *spec_name;
            char *err_name;
            cpl_vector *spec_vec;
            cpl_vector *err_vec;
            double snr;
            int order_real, order_idx, order_idxp;
            order_idx = order_idx_values[i] ;
            order_idxp = cr2res_io_convert_order_idx_to_idxp(order_idx) ;
            order_real = cr2res_order_idx_to_real(order_idx, order_zp) ;
            key_name = cpl_sprintf(CR2RES_HEADER_QC_REAL_ORDER, order_idxp) ;
            cpl_propertylist_append_int(ext_plist_loc, key_name, order_real) ;
            cpl_free(key_name) ;

            /* Compute SNR and store QC*/
            key_name = cpl_sprintf(CR2RES_HEADER_QC_SNR, order_idxp);
            spec_name = cr2res_dfs_POL_INTENS_colname(order_idxp);
            err_name = cr2res_dfs_POL_INTENS_ERROR_colname(order_idxp);

            spec_vec = cpl_vector_wrap(cpl_table_get_nrow(pol_spec_merged),
                cpl_table_get_data_double(pol_spec_merged, spec_name));
            err_vec = cpl_vector_wrap(cpl_table_get_nrow(pol_spec_merged),
                cpl_table_get_data_double(pol_spec_merged, err_name));
            snr = cr2res_qc_compute_snr(spec_vec, err_vec);
            cpl_propertylist_append_double(ext_plist_loc, key_name, snr);

            cpl_vector_unwrap(spec_vec);
            cpl_vector_unwrap(err_vec);
            cpl_free(key_name);
            cpl_free(spec_name) ;
            cpl_free(err_name) ;
        }
        cpl_free(order_idx_values) ;
    }
    cpl_table_delete(trace_wave) ;

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

  The int array will need to be freed by the caller. Its size is
  norders. It contains the list of orders found in all the input ext tables
 */
/*----------------------------------------------------------------------------*/
static int * cr2res_obs_pol_get_order_numbers(
        const cpl_table **  extracted,
        int                 next,
        int             *   norders)
{
    const char  *   col_name ;
    char        *   col_type ;
    int         *   tmp_orders_list ;
    int             count_orders ;
    int         *   out_orders ;
    cpl_size        j;
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
    for (i = 0; i < next; i++) {
        cpl_size ncols;
        cpl_array *col_names;
        col_names = cpl_table_get_column_names(extracted[i]);
        ncols = cpl_table_get_ncol(extracted[i]) ;
        for (j=0 ; j<ncols ; j++) {
            col_name = cpl_array_get_string(col_names, j);
            col_type = cr2res_dfs_SPEC_colname_parse(col_name, &order,
                    &trace_nb) ;
            if (col_type != NULL && !strcmp(col_type, CR2RES_COL_SPEC_SUFFIX)) {
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
  @param    ngroups     The number of groups
  @return   1 if valid, 0 if not, -1 in error case
 */
/*----------------------------------------------------------------------------*/
static int cr2res_obs_pol_check_inputs_validity(
        const cpl_frameset  *   rawframes,
        cpl_size            *   ngroups)
{
    cr2res_nodding_pos  *       nod_positions ;
    cpl_frameset        *       raw_a ;
    cpl_frameset        *       raw_b ;
    cpl_size                    nframes, nframes_a, nframes_b, i ;

    /* Check Inputs */
    if (rawframes == NULL || ngroups == NULL) return -1 ;
    nframes = cpl_frameset_get_size(rawframes) ;
    
    /* Check that the MJD-OBS is increasing (only when several frames) */
    if (nframes > 1) {
        cpl_propertylist *plist;
        cpl_frame *cur_frame;
        const char *cur_fname;
        double mjd_obs_previous;
        /* Get the first MJDOBS */
        cur_frame = cpl_frameset_get_position((cpl_frameset*)rawframes, 0) ;
        cur_fname = cpl_frame_get_filename(cur_frame);
        plist = cpl_propertylist_load(cur_fname, 0) ;
        mjd_obs_previous = cr2res_pfits_get_mjd_obs(plist) ;
        cpl_propertylist_delete(plist) ;
        if (cpl_error_get_code() != CPL_ERROR_NONE) {
            cpl_msg_error(__func__, "Cannot access MJD-OBS") ;
            return -1 ;
        }
        /* Loop on the n-1 last frames */
        for (i=1 ; i<nframes ; i++) {
            double mjd_obs;
            /* Get the current MJDOBS */
            cur_frame = cpl_frameset_get_position((cpl_frameset*)rawframes, i) ;
            cur_fname = cpl_frame_get_filename(cur_frame);
            plist = cpl_propertylist_load(cur_fname, 0) ;
            mjd_obs = cr2res_pfits_get_mjd_obs(plist) ;
            cpl_propertylist_delete(plist) ;
            if (cpl_error_get_code() != CPL_ERROR_NONE) {
                cpl_msg_error(__func__, "Cannot access MJD-OBS") ;
                return -1 ;
            }
            if (mjd_obs <= mjd_obs_previous) {
                cpl_msg_error(__func__, "Rawframes MJD-OBS is not increasing") ;
                return 0 ;
            }
            mjd_obs_previous = mjd_obs ;
        }
    }

    /* Get the nodding positions */
    nod_positions = cr2res_nodding_read_positions(rawframes) ;
    if (nod_positions == NULL) {
        cpl_msg_error(__func__, "Cannot read nodding positions") ;
        return -1 ;
    }
    if (cpl_msg_get_level() == CPL_MSG_DEBUG)
        for (i=0 ; i<nframes ; i++) 
            cpl_msg_debug(__func__, "Frame %s - Nodding %c",
                    cpl_frame_get_filename(
                        cpl_frameset_get_position_const(rawframes,i)),
                    cr2res_nodding_position_char(nod_positions[i])) ;

    /* Split the frames */
    if (cr2res_combine_nodding_split_frames(rawframes, nod_positions, 
                &raw_a, &raw_b)) {
        cpl_msg_error(__func__, "Cannot split the nodding positions") ;
        cpl_free(nod_positions) ;    
        return -1 ;
    }
    cpl_free(nod_positions) ;    

    /* Check the sizes of A and B */
    nframes_a = cpl_frameset_get_size(raw_a) ;
    nframes_b = cpl_frameset_get_size(raw_b) ;
    if (nframes_a != nframes_b) {
        cpl_msg_error(__func__, "Number of A and B positions differ") ;
        cpl_frameset_delete(raw_a);
        cpl_frameset_delete(raw_b);
        return 0 ;
    }
    cpl_frameset_delete(raw_a);
    cpl_frameset_delete(raw_b);

    /* number of frames must be modulo 4 */
    if (nframes_a % CR2RES_POLARIMETRY_GROUP_SIZE) {
        cpl_msg_error(__func__, "Require a multiple of 4 raw frames") ;
        return 0 ;
    }

    /* Return */
    *ngroups = nframes_a/CR2RES_POLARIMETRY_GROUP_SIZE ;
    return 1 ;
}

