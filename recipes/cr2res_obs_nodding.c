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
#include "cr2res_nodding.h"
#include "cr2res_calib.h"
#include "cr2res_pfits.h"
#include "cr2res_dfs.h"
#include "cr2res_bpm.h"
#include "cr2res_trace.h"
#include "cr2res_extract.h"
#include "cr2res_io.h"
#include "cr2res_qc.h"
#include "cr2res_photom.h"

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

static int cr2res_obs_nodding_astrometry_compare(
        const cpl_frame   *   frame1,
        const cpl_frame   *   frame2) ;
static cpl_frameset * cr2res_obs_nodding_find_RAW(
        const cpl_frameset  *   in,
        int                 *   type) ;
static int cr2res_obs_nodding_check_inputs_validity(
        const cpl_frameset  *   rawframes) ;
static int cr2res_obs_nodding_reduce(
        const cpl_frameset  *   rawframes,
        const cpl_frame     *   trace_wave_frame,
        const cpl_frame     *   detlin_frame,
        const cpl_frame     *   master_dark_frame,
        const cpl_frame     *   master_flat_frame,
        const cpl_frame     *   bpm_frame,
        const cpl_frame     *   blaze_frame,
        int                     nodding_invert,
        int                     subtract_nolight_rows,
        int                     subtract_interorder_column,
        int                     cosmics,
        int                     error_method,
        int                     extract_oversample,
        int                     extract_swath_width,
        int                     extract_height,
        double                  extract_smooth_slit,
        double                  extract_smooth_spec,
        int                     extract_niter,
        double                  extract_kappa,
        int                     reduce_det,
        int                     disp_det,
        int                     disp_order_idx,
        int                     disp_trace,
        hdrl_image          **  combineda,
        cpl_table           **  extracta,
        cpl_table           **  slitfunca,
        hdrl_image          **  modela,
        cpl_table           **  twa,
        hdrl_image          **  combinedb,
        cpl_table           **  extractb,
        cpl_table           **  slitfuncb,
        hdrl_image          **  modelb,
        cpl_table           **  twb,
        cpl_table           **  extractc,
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
  If the input are standard stars, it computes a post processing step   \n\
  to determine the throughput.                                          \n\
  If the input are spectro astrometric data, it will apply the nodding  \n\
  on each of the sub-groups                                             \n\
                                                                        \n\
  Inputs                                                                \n\
    raw.fits " CR2RES_CAL_NODDING_OTHER_RAW" [2 to 2n]                  \n\
          or " CR2RES_CAL_NODDING_JITTER_RAW" [2 to 2n]                 \n\
          or " CR2RES_OBS_NODDING_OTHER_RAW" [2 to 2n]                  \n\
          or " CR2RES_OBS_NODDING_JITTER_RAW" [2 to 2n]                 \n\
          or " CR2RES_OBS_ASTROMETRY_OTHER_RAW" [2 to 2n]               \n\
          or " CR2RES_OBS_ASTROMETRY_JITTER_RAW" [2 to 2n]              \n\
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
    photo_flux.fits " CR2RES_PHOTO_FLUX_PROCATG " [0 to 1]              \n\
    blaze.fits " CR2RES_CAL_FLAT_EXTRACT_1D_PROCATG " [0 to 1]          \n\
                                                                        \n\
  Outputs                                                               \n\
    cr2res_obs_nodding_extractedA.fits " 
    CR2RES_OBS_NODDING_EXTRACTA_PROCATG "\n\
    cr2res_obs_nodding_extractedB.fits " 
    CR2RES_OBS_NODDING_EXTRACTB_PROCATG "\n\
    cr2res_obs_nodding_combinedA.fits " 
    CR2RES_OBS_NODDING_COMBINEDA_PROCATG "\n\
    cr2res_obs_nodding_combinedB.fits " 
    CR2RES_OBS_NODDING_COMBINEDB_PROCATG "\n\
    cr2res_obs_nodding_trace_wave_A.fits "
    CR2RES_OBS_NODDING_TWA_PROCATG "\n\
    cr2res_obs_nodding_trace_wave_B.fits "
    CR2RES_OBS_NODDING_TWB_PROCATG "\n\
    cr2res_obs_nodding_modelA.fits " 
    CR2RES_OBS_NODDING_SLITMODELA_PROCATG "\n\
    cr2res_obs_nodding_modelB.fits " 
    CR2RES_OBS_NODDING_SLITMODELB_PROCATG "\n\
    cr2res_obs_nodding_slitfuncA.fits " 
    CR2RES_OBS_NODDING_SLITFUNCA_PROCATG "\n\
    cr2res_obs_nodding_slitfuncB.fits " 
    CR2RES_OBS_NODDING_SLITFUNCB_PROCATG "\n\
    cr2res_obs_nodding_throughput.fits " 
    CR2RES_OBS_NODDING_THROUGHPUT_PROCATG "\n\
                                                                        \n\
  Algorithm                                                             \n\
    loop on detectors d:                                                \n\
      call cr2res_obs_nodding_reduce()                                  \n\
        -> combined[a|b](d)                                             \n\
        -> extract[a|b|c](d)                                            \n\
        -> slitfunc[a|b](d)                                             \n\
        -> model[a|b](d)                                                \n\
        -> tw[a|b](d)                                                   \n\
    Save combineda and combinedb                                        \n\
    Save extracta, extractb, extractc                                   \n\
    Save slitfunca and slitfuncb                                        \n\
    Save modela and modelb                                              \n\
    Save throughput                                                     \n\
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
        -> trace_wave[a|b]                                              \n\
      Compute QC parameters                                             \n\
      If STD star, compute the throughput                               \n\
        -> throughput                                                   \n\
                                                                        \n\
  Library functions used                                                \n\
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
    cr2res_qc_obs_slit_psf()                                    \n\
    cr2res_photom_engine()                                              \n\
    cr2res_io_save_COMBINED()                                           \n\
    cr2res_io_save_EXTRACT_1D()                                         \n\
    cr2res_io_save_SLIT_FUNC()                                          \n\
    cr2res_io_save_SLIT_MODEL()                                         \n\
    cr2res_io_save_THROUGHPUT()                                         \n\
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
                    CR2RES_PIPELINE_AUTHORS,
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
    p = cpl_parameter_new_value(
            "cr2res.cr2res_obs_nodding.subtract_nolight_rows",
            CPL_TYPE_BOOL, 
            "Subtract median row from baffled region at detector bottom",
            "cr2res.cr2res_obs_nodding", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "subtract_nolight_rows");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);
    
    p = cpl_parameter_new_value(
            "cr2res.cr2res_obs_nodding.subtract_interorder_column",
            CPL_TYPE_BOOL,
            "Subtract column-by-column fit to the pixel values between orders",
            "cr2res.cr2res_obs_nodding", TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, 
            "subtract_interorder_column");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_nodding.cosmics",
            CPL_TYPE_BOOL, "Find and mark cosmic rays hits as bad",
            "cr2res.cr2res_obs_nodding", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "cosmics");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_enum("cr2res.cr2res_obs_nodding.error_method", 
            CPL_TYPE_STRING, "The 1d extraction error calculation method",
            "cr2res.cr2res_obs_nodding",
            "Poisson", 2, "Poisson", "Horne");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "error_method");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_nodding.nodding_invert",
            CPL_TYPE_BOOL, "Flag to use when A is above B",
            "cr2res.cr2res_obs_nodding", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "nodding_invert");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_nodding.extract_oversample",
            CPL_TYPE_INT, "factor by which to oversample the extraction",
            "cr2res.cr2res_obs_nodding", 7);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_oversample");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_nodding.extract_swath_width",
            CPL_TYPE_INT, "The swath width", "cr2res.cr2res_obs_nodding", 2048);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_swath_width");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_nodding.extract_height",
            CPL_TYPE_INT, "Extraction height", "cr2res.cr2res_obs_nodding", -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_height");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_nodding.extract_smooth_slit",
            CPL_TYPE_DOUBLE, "Smoothing along the slit",
            "cr2res.cr2res_obs_nodding", 1.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_smooth_slit");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_nodding.extract_smooth_spec",
            CPL_TYPE_DOUBLE, "Smoothing along spectrum",
            "cr2res.cr2res_obs_nodding", 8e-8);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_smooth_spec");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_nodding.detector",
            CPL_TYPE_INT, "Only reduce the specified detector",
            "cr2res.cr2res_obs_nodding", 0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "detector");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_nodding.create_idp",
            CPL_TYPE_BOOL, "Flag to produce  IDP files",
            "cr2res.cr2res_obs_nodding", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "idp");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_nodding.display_detector",
            CPL_TYPE_INT, "Apply the display for the specified detector",
            "cr2res.cr2res_obs_nodding", 0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "display_detector");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_nodding.display_order",
            CPL_TYPE_INT, "Apply the display for the specified order",
            "cr2res.cr2res_obs_nodding", 1000);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "display_order");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_nodding.display_trace",
            CPL_TYPE_INT, "Apply the display for the specified trace",
            "cr2res.cr2res_obs_nodding", 0);
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
                            extract_height, reduce_det, 
                            disp_order_idx, disp_trace, disp_det, 
                            nodding_invert, create_idp, subtract_nolight_rows,
                            subtract_interorder_column, cosmics,
                            error_method;
    double                  extract_smooth_slit, extract_smooth_spec;
    double                  dit, gain;
    double                  ra, dec, mjd_obs, mjd_cen, geolon, geolat, geoelev,
                            barycorr;
    cpl_frameset        *   rawframes ;
    cpl_frameset        *   raw_flat_frames ;
    const cpl_frame     *   trace_wave_frame ;
    const cpl_frame     *   detlin_frame ;
    const cpl_frame     *   master_dark_frame ;
    const cpl_frame     *   photo_flux_frame ;
    const cpl_frame     *   blaze_frame ;
    const cpl_frame     *   master_flat_frame ;
    const cpl_frame     *   bpm_frame ;
    cpl_size            *   labels ;
    cpl_size                nlabels, l ;
    char                *   product_name_addon ;
    hdrl_image          *   combineda[CR2RES_NB_DETECTORS] ;
    cpl_table           *   extracta[CR2RES_NB_DETECTORS] ;
    cpl_table           *   slitfunca[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   modela[CR2RES_NB_DETECTORS] ;
    cpl_table           *   twa[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   combinedb[CR2RES_NB_DETECTORS] ;
    cpl_table           *   extractb[CR2RES_NB_DETECTORS] ;
    cpl_table           *   slitfuncb[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   modelb[CR2RES_NB_DETECTORS] ;
    cpl_table           *   twb[CR2RES_NB_DETECTORS] ;
    cpl_table           *   extractc[CR2RES_NB_DETECTORS] ;
    cpl_table           *   throughput[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   qc_main ;
    cpl_propertylist    *   ext_plist[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   ext_plist_photom[CR2RES_NB_DETECTORS] ;
    char                *   cur_setting ;
    char                *   out_file;
    cpl_table           *   eop_table ;
    int                     det_nr, type; 

    /* Initialise */
    gain = 0.0 ;
    barycorr = 0.0;

    /* RETRIEVE INPUT PARAMETERS */
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_nodding.subtract_nolight_rows");
    subtract_nolight_rows = cpl_parameter_get_bool(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_nodding.subtract_interorder_column");
    subtract_interorder_column = cpl_parameter_get_bool(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_nodding.nodding_invert");
    nodding_invert = cpl_parameter_get_bool(param);
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
            "cr2res.cr2res_obs_nodding.extract_smooth_slit");
    extract_smooth_slit = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_nodding.extract_smooth_spec");
    extract_smooth_spec = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_nodding.detector");
    reduce_det = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_nodding.create_idp");
    create_idp = cpl_parameter_get_bool(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_nodding.display_detector");
    disp_det = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_nodding.display_order");
    disp_order_idx = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_nodding.display_trace");
    disp_trace = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_nodding.cosmics");
    cosmics = cpl_parameter_get_bool(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_nodding.error_method");
    if(strcmp(cpl_parameter_get_string(param),"Horne") == 0)
        error_method = CR2RES_EXTRACT_ERROR_HORNE;
    else
        error_method = CR2RES_EXTRACT_ERROR_POISSON;

    /* TODO, make parameters, maybe */
    int extract_niter = 30;
    double extract_kappa = 10;

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
    if (master_dark_frame != NULL) {
        cpl_msg_warning(__func__,
                "Providing a MASTER DARK is not recommended for this recipe") ;
    }
    blaze_frame = cpl_frameset_find_const(frameset,
            CR2RES_CAL_FLAT_EXTRACT_1D_PROCATG) ; 
    photo_flux_frame = cpl_frameset_find_const(frameset,
            CR2RES_PHOTO_FLUX_PROCATG) ; 
    master_flat_frame = cpl_frameset_find_const(frameset,
            CR2RES_CAL_FLAT_MASTER_PROCATG) ; 
    bpm_frame = cr2res_io_find_BPM(frameset) ;

    /* Get the RAW Frames */
    rawframes = cr2res_obs_nodding_find_RAW(frameset, &type) ;
    if (rawframes == NULL) {
        cpl_msg_error(__func__, "Could not find RAW frames") ;
        return -1 ;
    }
      
    /* Get the RAW flat frames */
    raw_flat_frames = cr2res_extract_frameset(frameset, CR2RES_FLAT_RAW) ;

    /* Label the raw frames with the different angles */
    if ((labels = cpl_frameset_labelise(rawframes, 
                cr2res_obs_nodding_astrometry_compare, &nlabels)) == NULL) {
        cpl_msg_error(__func__, "Cannot label input frames") ;
        cpl_frameset_delete(rawframes) ;
        if (raw_flat_frames != NULL) cpl_frameset_delete(raw_flat_frames) ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Check the number of angles */
    if (type != 3 && nlabels != 1) {
        cpl_msg_error(__func__, "Expect only one DROT POSANG value - abort") ;
        cpl_frameset_delete(rawframes) ;
        if (raw_flat_frames != NULL) cpl_frameset_delete(raw_flat_frames) ;
        cpl_free(labels) ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Loop on the settings */
    for (l = 0; l < (int)nlabels; l++) {
        double drot_posang;
        cpl_frameset *raw_one_angle;
        cpl_propertylist *plist;
        /* Get the frames for the current angle */
        raw_one_angle = cpl_frameset_extract(rawframes, labels, (cpl_size)l) ;

        /* Get the current angle */
        plist = cpl_propertylist_load(cpl_frame_get_filename(
                    cpl_frameset_get_position(raw_one_angle, 0)), 0) ;
        drot_posang = cr2res_pfits_get_drot_posang(plist) ;
        cpl_propertylist_delete(plist) ;

        cpl_msg_info(__func__, "Process Angle %g", drot_posang) ;
        cpl_msg_indent_more() ;

        /* Loop on the detectors */
        for (det_nr = 1; det_nr <= CR2RES_NB_DETECTORS; det_nr++) {
            /* Initialise */
            combineda[det_nr - 1] = NULL;
            extracta[det_nr - 1] = NULL;
            slitfunca[det_nr - 1] = NULL;
            modela[det_nr - 1] = NULL;
            twa[det_nr - 1] = NULL;
            combinedb[det_nr - 1] = NULL;
            extractb[det_nr - 1] = NULL;
            slitfuncb[det_nr - 1] = NULL;
            modelb[det_nr - 1] = NULL;
            twb[det_nr - 1] = NULL;
            extractc[det_nr - 1] = NULL;
            ext_plist[det_nr - 1] = NULL;
            ext_plist_photom[det_nr - 1] = NULL;
            throughput[det_nr - 1] = NULL;

            /* Compute only one detector */
            if (reduce_det != 0 && det_nr != reduce_det)
                continue;

            cpl_msg_info(__func__, "Process Detector %d", det_nr);
            cpl_msg_indent_more();

            /* Call the reduction function */
            if (cr2res_obs_nodding_reduce(
                    raw_one_angle, trace_wave_frame,
                    detlin_frame, master_dark_frame, master_flat_frame,
                    bpm_frame, blaze_frame, nodding_invert,
                    subtract_nolight_rows, subtract_interorder_column, cosmics,
                    error_method,
                    extract_oversample, extract_swath_width, extract_height,
                    extract_smooth_slit, extract_smooth_spec, extract_niter,
                    extract_kappa, det_nr, disp_det, disp_order_idx, disp_trace,
                    &(combineda[det_nr - 1]), &(extracta[det_nr - 1]),
                    &(slitfunca[det_nr - 1]), &(modela[det_nr - 1]),
                    &(twa[det_nr - 1]), &(combinedb[det_nr - 1]),
                    &(extractb[det_nr - 1]), &(slitfuncb[det_nr - 1]),
                    &(modelb[det_nr - 1]), &(twb[det_nr - 1]),
                    &(extractc[det_nr - 1]), &(ext_plist[det_nr - 1])) == -1) {
                cpl_msg_warning(__func__, "Failed to reduce detector %d",
                                det_nr);
                cpl_error_reset();
            }
            else if (type == 2) {
                cpl_msg_info(
                    __func__,
                    "Sensitivity / Conversion / Throughput computation");
                cpl_msg_indent_more();

                /* Define the gain */
                if (det_nr == 1)
                    gain = CR2RES_GAIN_CHIP1;
                if (det_nr == 2)
                    gain = CR2RES_GAIN_CHIP2;
                if (det_nr == 3)
                    gain = CR2RES_GAIN_CHIP3;

                /* Get the RA and DEC observed */
                plist = cpl_propertylist_load(
                    cpl_frame_get_filename(
                        cpl_frameset_get_position_const(raw_one_angle, 0)),
                    0);
                ra = cr2res_pfits_get_ra(plist);
                dec = cr2res_pfits_get_dec(plist);
                dit = cr2res_pfits_get_dit(plist);
                cur_setting = cpl_strdup(cr2res_pfits_get_wlen_id(plist));
                cr2res_format_setting(cur_setting);
                cpl_propertylist_delete(plist);
                if (cpl_error_get_code()) {
                    cpl_msg_indent_less();
                    cpl_error_reset();
                    cpl_msg_warning(__func__, "Missing Header Informations");
                }
                else {
                    /* Compute the photometry */
                    if (cr2res_photom_engine(
                            extracta[det_nr - 1],
                            cpl_frame_get_filename(photo_flux_frame),
                            cur_setting, ra, dec, gain, dit, disp_det == det_nr,
                            disp_order_idx, disp_trace,
                            &(throughput[det_nr - 1]),
                            &(ext_plist_photom[det_nr - 1])) == -1) {
                        cpl_msg_warning(__func__,
                                        "Failed to reduce detector %d", det_nr);
                        cpl_error_reset();
                    }
                }
                cpl_free(cur_setting);
                cpl_msg_indent_less();
            }
            cpl_msg_indent_less();
        }

        /* Save Products */
        if (nlabels == 1)   product_name_addon = cpl_sprintf(".fits") ;
        else                product_name_addon = cpl_sprintf("_%g.fits",
                drot_posang);

        /* Add the photom QC to the std ones */
        for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
            if (ext_plist_photom[det_nr-1] != NULL &&
                    ext_plist[det_nr-1] != NULL) {
                cpl_propertylist_append(ext_plist[det_nr-1],
                        ext_plist_photom[det_nr-1]) ;
            } 
        }

        /* Add ESO.DRS.TMID in the Main Header */
        qc_main = cpl_propertylist_new();
        cpl_propertylist_append_double(qc_main,
                CR2RES_HEADER_DRS_TMID,
                cr2res_utils_get_center_mjd(raw_one_angle)) ;

        /* Add barycentric correction */
        eop_table = cr2res_io_get_eop_table() ;
        if (eop_table != NULL) {
            plist=cpl_propertylist_load(cpl_frame_get_filename(
                        cpl_frameset_get_position_const(raw_one_angle, 0)), 0) ;

            ra = cpl_propertylist_get_double(plist, "RA") ;
            dec = cpl_propertylist_get_double(plist, "DEC") ;
            mjd_obs = cpl_propertylist_get_double(plist, "MJD-OBS") ;
            geolon = cpl_propertylist_get_double(plist, "ESO TEL GEOLON") ;
            geolat = cpl_propertylist_get_double(plist, "ESO TEL GEOLAT") ;
            geoelev = cpl_propertylist_get_double(plist, "ESO TEL GEOELEV") ;

            cpl_propertylist_delete(plist) ;

            barycorr = 0.0 ;
            if (!cpl_error_get_code()) {
                mjd_cen = cr2res_utils_get_center_mjd(raw_one_angle) ;
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
            
        }


        if(error_method == CR2RES_EXTRACT_ERROR_HORNE){
            cpl_propertylist_append_string(qc_main, CR2RES_HEADER_DRS_ERRMETHOD,
                    "Horne") ;
            }
        else {
            cpl_propertylist_append_string(qc_main, CR2RES_HEADER_DRS_ERRMETHOD,
                "Poisson") ;
            }
        /* Add QC NUMSAT */
        cpl_propertylist_append_int(qc_main,
                CR2RES_HEADER_QC_NUMSAT,
                cr2res_qc_numsat(raw_one_angle)) ;


        /* Save only the used RAW - fill raw_one_angle with CALIBS */
        if (trace_wave_frame != NULL) 
            cpl_frameset_insert(raw_one_angle,
                    cpl_frame_duplicate(trace_wave_frame)) ;
        if (detlin_frame != NULL) 
            cpl_frameset_insert(raw_one_angle,
                    cpl_frame_duplicate(detlin_frame)) ;
        if (master_dark_frame != NULL) 
            cpl_frameset_insert(raw_one_angle,
                    cpl_frame_duplicate(master_dark_frame)) ;
        if (master_flat_frame!= NULL) 
            cpl_frameset_insert(raw_one_angle,
                    cpl_frame_duplicate(master_flat_frame)) ;
        if (bpm_frame!= NULL) 
            cpl_frameset_insert(raw_one_angle,
                    cpl_frame_duplicate(bpm_frame)) ;
        if (blaze_frame!= NULL) 
            cpl_frameset_insert(raw_one_angle,
                    cpl_frame_duplicate(blaze_frame)) ;
        if (photo_flux_frame!= NULL) 
            cpl_frameset_insert(raw_one_angle,
                    cpl_frame_duplicate(photo_flux_frame)) ;

        out_file = cpl_sprintf("%s_combinedA%s", RECIPE_STRING, 
                product_name_addon) ;
        cr2res_io_save_COMBINED(out_file, frameset, raw_one_angle, parlist,
                combineda, qc_main, ext_plist, 
                CR2RES_OBS_NODDING_COMBINEDA_PROCATG, RECIPE_STRING) ;
        cpl_free(out_file);

        char *qc_perdet_keywords[] = { CR2RES_HEADER_QC_SIGNAL,
                                       CR2RES_HEADER_QC_STANDARD_FLUX,
                                       CR2RES_HEADER_QC_SLITFWHM_MED };
        
        char *qc_pertrace_keywords[] = { CR2RES_HEADER_QC_SNR_BASE,
                                         CR2RES_HEADER_QC_SLITFWHM_ORDER_BASE };
                                
        int qc_perdet_size = sizeof(qc_perdet_keywords) / sizeof(char *);
        int qc_pertrace_size = sizeof(qc_pertrace_keywords) / sizeof(char *);

        int qc_sizes_ext[] = {CR2RES_NB_DETECTORS};
       cpl_propertylist** qc_plists_ext[] = {ext_plist} ; 
       char * qc_res_avg, * qc_res_rmsd, * qc_ref;

       int qc_i;

       for (qc_i=0; qc_i<qc_perdet_size ; qc_i++) {

                    qc_res_avg = cpl_sprintf("%s %s",
                            qc_perdet_keywords[qc_i], "AVG");
                    qc_res_rmsd = cpl_sprintf("%s %s",
                            qc_perdet_keywords[qc_i], "RMS");

                cr2res_qc_calculate_mean_and_rmsd(qc_plists_ext, 1, qc_sizes_ext, 
                                 qc_perdet_keywords[qc_i], qc_main,
                                 qc_res_avg, qc_res_rmsd);

                cpl_free(qc_res_avg);
                cpl_free(qc_res_rmsd);
       }


       int min_order = INT_MAX;
       int max_order = INT_MIN;
       int current_min, current_max;

       for (qc_i = 0; qc_i < CR2RES_NB_DETECTORS; qc_i++) {
           if (twa[qc_i] == NULL) {
               continue;
           }
           current_min =
               cpl_table_get_column_min(twa[qc_i], "Order");
           if (cpl_error_get_code() != CPL_ERROR_NONE) {
               continue;
           }
           if (current_min < min_order) {
               min_order = current_min;
           }

           current_max =
               cpl_table_get_column_max(twa[qc_i], "Order");
           if (cpl_error_get_code() != CPL_ERROR_NONE) {
               continue;
           }
           if (current_max > max_order) {
               max_order = current_max;
           }
       }

       for (qc_i = 0; qc_i < qc_pertrace_size; qc_i++) {
           for (int order_id = min_order; order_id <= max_order; order_id++) {

               qc_ref = cpl_sprintf("%s%d", qc_pertrace_keywords[qc_i],
                                    order_id);
               qc_res_avg =
                   cpl_sprintf("%s%d %s", qc_pertrace_keywords[qc_i],
                               order_id, "AVG");
               qc_res_rmsd =
                   cpl_sprintf("%s%d %s", qc_pertrace_keywords[qc_i],
                               order_id, "RMS");

               cr2res_qc_calculate_mean_and_rmsd(qc_plists_ext, 1, qc_sizes_ext,
                                                 qc_ref, qc_main, qc_res_avg,
                                                 qc_res_rmsd);

               cpl_free(qc_ref);
               cpl_free(qc_res_avg);
               cpl_free(qc_res_rmsd);
           }
       }

        out_file = cpl_sprintf("%s_extractedA%s", RECIPE_STRING,
                product_name_addon) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, raw_one_angle, parlist, 
                extracta, qc_main, ext_plist, 
                CR2RES_OBS_NODDING_EXTRACTA_PROCATG, RECIPE_STRING);
        if (create_idp) {
            cr2res_idp_save(out_file, frameset, raw_one_angle, parlist,
                            extracta, qc_main, ext_plist,
                            CR2RES_OBS_NODDING_EXTRACTA_IDP_PROCATG,
                            RECIPE_STRING);
        }
        cpl_free(out_file);

        /* for (qc_i=0; qc_i<qc_perdet_size ; qc_i++) {

                    qc_res_avg = cpl_sprintf("%s %s",
                            qc_perdet_keywords[qc_i], "AVG");
                    qc_res_rmsd = cpl_sprintf("%s %s",
                            qc_perdet_keywords[qc_i], "RMS");

            cpl_propertylist_erase(qc_main, qc_res_avg);
            if(cpl_error_get_code() != CPL_ERROR_NONE) {
                cpl_msg_debug(__func__, "Cannot erase %s", qc_perdet_keywords[qc_i]);
                cpl_error_reset();
            }

            cpl_propertylist_erase(qc_main, qc_res_rmsd);
            if(cpl_error_get_code() != CPL_ERROR_NONE) {
                cpl_msg_debug(__func__, "Cannot erase %s", qc_perdet_keywords[qc_i]);
                cpl_error_reset();
            }


                cpl_free(qc_res_avg);
                cpl_free(qc_res_rmsd);
        }
 */
        out_file = cpl_sprintf("%s_slitfuncA%s", RECIPE_STRING,
                product_name_addon) ;
        cr2res_io_save_SLIT_FUNC(out_file, frameset, raw_one_angle, parlist,
                slitfunca, qc_main, ext_plist, 
                CR2RES_OBS_NODDING_SLITFUNCA_PROCATG, RECIPE_STRING) ;
        cpl_free(out_file);

        out_file = cpl_sprintf("%s_modelA%s", RECIPE_STRING,
                product_name_addon) ;
        cr2res_io_save_SLIT_MODEL(out_file, frameset, raw_one_angle, parlist,
                modela, qc_main, ext_plist, 
                CR2RES_OBS_NODDING_SLITMODELA_PROCATG, RECIPE_STRING) ;
        cpl_free(out_file);

        out_file = cpl_sprintf("%s_trace_wave_A%s", RECIPE_STRING,
                product_name_addon) ;
        cr2res_io_save_TRACE_WAVE(out_file, frameset, raw_one_angle, parlist,
                twa, qc_main, ext_plist, CR2RES_OBS_NODDING_TWA_PROCATG,
                RECIPE_STRING) ;
        cpl_free(out_file);

        out_file = cpl_sprintf("%s_combinedB%s", RECIPE_STRING,
                product_name_addon) ;
        cr2res_io_save_COMBINED(out_file, frameset, raw_one_angle, parlist,
                combinedb, qc_main, ext_plist, 
                CR2RES_OBS_NODDING_COMBINEDB_PROCATG, RECIPE_STRING) ;
        cpl_free(out_file);

        out_file = cpl_sprintf("%s_extractedB%s", RECIPE_STRING,
                product_name_addon) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, raw_one_angle, parlist, 
                extractb, qc_main, ext_plist, 
                CR2RES_OBS_NODDING_EXTRACTB_PROCATG, RECIPE_STRING);
        if (create_idp) {
            cr2res_idp_save(out_file, frameset, raw_one_angle, parlist,
                            extractb, qc_main, ext_plist,
                            CR2RES_OBS_NODDING_EXTRACTB_IDP_PROCATG,
                            RECIPE_STRING);
        }
        cpl_free(out_file);

        out_file = cpl_sprintf("%s_slitfuncB%s", RECIPE_STRING,
                product_name_addon) ;
        cr2res_io_save_SLIT_FUNC(out_file, frameset, raw_one_angle, parlist,
                slitfuncb, qc_main, ext_plist, 
                CR2RES_OBS_NODDING_SLITFUNCB_PROCATG, RECIPE_STRING) ;
        cpl_free(out_file);

        out_file = cpl_sprintf("%s_modelB%s", RECIPE_STRING,
                product_name_addon) ;
        cr2res_io_save_SLIT_MODEL(out_file, frameset, raw_one_angle, parlist,
                modelb, qc_main, ext_plist, 
                CR2RES_OBS_NODDING_SLITMODELB_PROCATG, RECIPE_STRING) ;
        cpl_free(out_file);
        
        out_file = cpl_sprintf("%s_trace_wave_B%s", RECIPE_STRING,
                product_name_addon) ;
        cr2res_io_save_TRACE_WAVE(out_file, frameset, raw_one_angle, parlist,
                twb, qc_main, ext_plist, CR2RES_OBS_NODDING_TWB_PROCATG,
                RECIPE_STRING) ;
        cpl_free(out_file);

        out_file = cpl_sprintf("%s_extracted_combined%s", RECIPE_STRING, 
                product_name_addon) ;
        cr2res_io_save_EXTRACT_1D(out_file, frameset, raw_one_angle, parlist, 
                extractc, qc_main, ext_plist, 
                CR2RES_OBS_NODDING_EXTRACTC_PROCATG, RECIPE_STRING);
        if (create_idp) {
            cr2res_idp_save(out_file, frameset, raw_one_angle, parlist,
                            extractc, qc_main, ext_plist,
                            CR2RES_OBS_NODDING_EXTRACTC_IDP_PROCATG,
                            RECIPE_STRING);
        }
        cpl_free(out_file);

        if (type == 2) {

            qc_res_avg = cpl_sprintf("%s %s", CR2RES_HEADER_QC_THROUGHPUT, "AVG");
            qc_res_rmsd = cpl_sprintf("%s %s", CR2RES_HEADER_QC_THROUGHPUT, "RMS");

            cr2res_qc_calculate_mean_and_rmsd(qc_plists_ext, 1, qc_sizes_ext,
                                              CR2RES_HEADER_QC_THROUGHPUT, qc_main,
                                              qc_res_avg, qc_res_rmsd);

            cpl_free(qc_res_avg);
            cpl_free(qc_res_rmsd);

            out_file = cpl_sprintf("%s_throughput%s", RECIPE_STRING,
                    product_name_addon) ;
            cr2res_io_save_THROUGHPUT(out_file, frameset, raw_one_angle, 
                    parlist, throughput, qc_main, ext_plist, 
                    CR2RES_OBS_NODDING_THROUGHPUT_PROCATG, RECIPE_STRING) ;
            cpl_free(out_file);
        }
        cpl_free(product_name_addon) ;

        /* Free */
        cpl_propertylist_delete(qc_main) ;
        for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
            if (combineda[det_nr-1] != NULL)
                hdrl_image_delete(combineda[det_nr-1]) ;
            if (extracta[det_nr-1] != NULL) 
                cpl_table_delete(extracta[det_nr-1]) ;
            if (slitfunca[det_nr-1] != NULL) 
                cpl_table_delete(slitfunca[det_nr-1]) ;
            if (modela[det_nr-1] != NULL)
                hdrl_image_delete(modela[det_nr-1]) ;
            if (twa[det_nr-1] != NULL) 
                cpl_table_delete(twa[det_nr-1]) ;
            if (combinedb[det_nr-1] != NULL)
                hdrl_image_delete(combinedb[det_nr-1]) ;
            if (extractb[det_nr-1] != NULL) 
                cpl_table_delete(extractb[det_nr-1]) ;
            if (slitfuncb[det_nr-1] != NULL) 
                cpl_table_delete(slitfuncb[det_nr-1]) ;
            if (modelb[det_nr-1] != NULL)
                hdrl_image_delete(modelb[det_nr-1]) ;
            if (twb[det_nr-1] != NULL) 
                cpl_table_delete(twb[det_nr-1]) ;
            if (extractc[det_nr-1] != NULL) 
                cpl_table_delete(extractc[det_nr-1]) ;
            if (throughput[det_nr-1] != NULL) 
                cpl_table_delete(throughput[det_nr-1]) ;
            if (ext_plist[det_nr-1] != NULL) 
                cpl_propertylist_delete(ext_plist[det_nr-1]) ;
            if (ext_plist_photom[det_nr-1] != NULL) 
                cpl_propertylist_delete(ext_plist_photom[det_nr-1]) ;
        }
        cpl_frameset_delete(raw_one_angle) ;
        cpl_msg_indent_less() ;
    }
    cpl_free(labels);
    cpl_frameset_delete(rawframes) ;
    if (raw_flat_frames != NULL) cpl_frameset_delete(raw_flat_frames) ;

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
  @param blaze_frame            Associated Blaze
  @param nodding_invert         Flag to use if A is above B
  @param subtract_nolight_rows
  @param cosmics                Flag to correct for cosmics
  @param extract_oversample     Extraction related
  @param extract_swath_width    Extraction related
  @param extract_height         Extraction related
  @param extract_smooth_slit    Extraction: smoothing along slit
  @param extract_smooth_spec    Extraction: smoothing along spectrum
  @param reduce_det             The detector to compute
  @param disp_det               The detector to display
  @param disp_order_idx         The order index to display
  @param disp_trace             The trace number to display
  @param combineda              [out] Combined image (A)
  @param extracta               [out] extracted spectrum (A)
  @param slitfunca              [out] slit function (A)
  @param modela                 [out] slit model (A)
  @param twa                    [out] trace wave (A)
  @param combinedb              [out] Combined image (B)
  @param extractb               [out] extracted spectrum (B)
  @param slitfuncb              [out] slit function (B)
  @param extractc               [out] extracted A and B combined spectrum
  @param modelb                 [out] slit model (B)
  @param twb                    [out] trace wave (B)
  @param ext_plist              [out] the header for saving the products
  @return  0 if ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
static int cr2res_obs_nodding_reduce(
        const cpl_frameset  *   rawframes,
        const cpl_frame     *   trace_wave_frame,
        const cpl_frame     *   detlin_frame,
        const cpl_frame     *   master_dark_frame,
        const cpl_frame     *   master_flat_frame,
        const cpl_frame     *   bpm_frame,
        const cpl_frame     *   blaze_frame,
        int                     nodding_invert,
        int                     subtract_nolight_rows,
        int                     subtract_interorder_column,
        int                     cosmics,
        int                     error_method,
        int                     extract_oversample,
        int                     extract_swath_width,
        int                     extract_height,
        double                  extract_smooth_slit,
        double                  extract_smooth_spec,
        int                     extract_niter,
        double                  extract_kappa,
        int                     reduce_det,
        int                     disp_det,
        int                     disp_order_idx,
        int                     disp_trace,
        hdrl_image          **  combineda,
        cpl_table           **  extracta,
        cpl_table           **  slitfunca,
        hdrl_image          **  modela,
        cpl_table           **  twa,
        hdrl_image          **  combinedb,
        cpl_table           **  extractb,
        cpl_table           **  slitfuncb,
        hdrl_image          **  modelb,
        cpl_table           **  twb,
        cpl_table           **  extractc,
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
    cpl_vector          *   ndits=NULL ;
    cpl_table           *   trace_wave ;
 //   cpl_table           *   trace_wave_corrected ;
    cpl_table           *   trace_wave_a ;
    cpl_table           *   trace_wave_b ;
    cpl_table           *   blaze_table ;
    cpl_array           *   slit_frac_a ;
    cpl_array           *   slit_frac_b ;
    cpl_table           *   extracted_a ;
    cpl_table           *   extracted_b ;
    cpl_table           *   extracted_combined ;
    cpl_table           *   slit_func_a ;
    cpl_table           *   slit_func_b ;
    hdrl_image          *   model_master_a ;
    hdrl_image          *   model_master_b ;
    cpl_propertylist    *   plist ;
    cpl_size                nframes, i ;
    char                *   key_name ;
    const char          *   first_fname ;
    double                  slit_length, extr_width_frac, slit_frac_a_bot, 
                            slit_frac_a_mid, slit_frac_a_top, slit_frac_b_bot, 
                            slit_frac_b_mid, slit_frac_b_top, nod_throw,
                            gain, error_factor ;
    double                  qc_signal_a, qc_signal_b, qc_fwhm_a, 
                            qc_fwhm_b, qc_standard_flux_a,
                            qc_standard_flux_b, qc_fwhm_med, blaze_norm ;
    cpl_array           *   fwhm_a_array ;
    cpl_array           *   fwhm_b_array ;
    char                *   cur_setting ;
    int                 *   order_idx_values ;
    double              *   qc_snrs ;
    double              *   qc_der_snrs ;
    int                     nb_order_idx_values,
                            order_zp, order_idx, order_idxp ;

    /* Check Inputs */
    if (combineda == NULL || combinedb == NULL || 
            twa == NULL || twb == NULL || 
            slitfunca == NULL || slitfuncb == NULL || 
            modela == NULL || modelb == NULL || 
            extracta == NULL || extractb == NULL || extractc == NULL || 
            ext_plist == NULL || rawframes == NULL || trace_wave_frame == NULL)
        return -1 ;

    /* Get the Gain */
    if (reduce_det == 1) gain = CR2RES_GAIN_CHIP1 ;
    else if (reduce_det == 2) gain = CR2RES_GAIN_CHIP2 ;
    else if (reduce_det == 3) gain = CR2RES_GAIN_CHIP3 ;
    else {
        cpl_msg_error(__func__, "Failed to get the Gain value") ;
        return -1 ;
    }

    /* Check raw frames consistency */
    if (cr2res_obs_nodding_check_inputs_validity(rawframes) != 1) {
        cpl_msg_error(__func__, "Invalid Inputs") ;
        return -1 ;
    }

    /* Initialise */
    nframes = cpl_frameset_get_size(rawframes) ;
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

    /* Get the Nodding positions */
    cpl_msg_info(__func__, "Get the Nodding positions") ;
    cpl_msg_indent_more() ;
    nod_positions = cr2res_nodding_read_positions(rawframes) ;
    for (i=0 ; i<nframes ; i++) {
        cpl_msg_info(__func__, "Frame %s - Nodding %c", 
                cpl_frame_get_filename(
                    cpl_frameset_get_position_const(rawframes, i)), 
            cr2res_nodding_position_char(nod_positions[i])) ;
    }
    cpl_msg_indent_less() ;

    /* Load the DITs if necessary */
    dits = cr2res_io_read_dits(rawframes) ;
    //if (master_dark_frame != NULL)  dits = cr2res_io_read_dits(rawframes) ;
    //else                            dits = NULL ;
    if (cpl_msg_get_level() == CPL_MSG_DEBUG && dits != NULL) 
        cpl_vector_dump(dits, stdout) ;

    /*Load the NDITs */
    ndits = cr2res_io_read_ndits(rawframes);
    if (cpl_msg_get_level() == CPL_MSG_DEBUG && ndits != NULL) 
        cpl_vector_dump(ndits, stdout) ;

    /* Load image list */
    cpl_msg_info(__func__, "Load the image list") ;
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
    cpl_msg_info(__func__, "Apply the Calibrations") ;
    cpl_msg_indent_more() ;
    if ((in_calib = cr2res_calib_imagelist(in, reduce_det, 0,
        subtract_nolight_rows, subtract_interorder_column, cosmics, 
        master_flat_frame, master_dark_frame, bpm_frame, detlin_frame, dits, 
        ndits))==NULL) {
        cpl_msg_error(__func__, "Failed to apply the calibrations") ;
        cpl_msg_indent_less() ;
        cpl_free(nod_positions) ;    
        if (dits != NULL) cpl_vector_delete(dits) ;
        if (ndits != NULL) cpl_vector_delete(ndits) ;
        hdrl_imagelist_delete(in) ;
        return -1 ;
    }
    hdrl_imagelist_delete(in) ;
    if (dits != NULL) cpl_vector_delete(dits) ;
    cpl_msg_indent_less() ;

    /* Split the image lists */
    cpl_msg_info(__func__, "Split the images in A and B lists") ;
    cpl_msg_indent_more() ;
    if (cr2res_combine_nodding_split(in_calib, nod_positions, &in_a, 
            &in_b)) {
        cpl_msg_error(__func__, "Failed to split the nodding positions") ;
        cpl_msg_indent_less() ;
        cpl_free(nod_positions) ;    
        hdrl_imagelist_delete(in_calib) ;
        cpl_vector_delete(ndits) ;
        return -1 ;
    }
    
    /* error factor is gain * ndit *nexp */
    if (error_method == CR2RES_EXTRACT_ERROR_HORNE) {
        error_factor = -1.0;
    }
    else {
        error_factor = gain * cpl_vector_get(ndits, 0) * 
                                        hdrl_imagelist_get_size(in_a);
    }
    
    for (i=0; i<cpl_vector_get_size(ndits); i++){
        if (cpl_vector_get(ndits,i) != cpl_vector_get(ndits, 0))
            cpl_msg_warning(__func__, "Raw frames have different NDIT! "
                "Error spectrum will likely be scaled incorrectly.");
    }
    cpl_vector_delete(ndits) ;
    cpl_free(nod_positions) ;    
    hdrl_imagelist_delete(in_calib) ;
    cpl_msg_indent_less() ;

    /* Check the sizes of A/B image lists */
    if ((hdrl_imagelist_get_size(in_a)) != hdrl_imagelist_get_size(in_b)
            || hdrl_imagelist_get_size(in_a) == 0) {
        cpl_msg_error(__func__, "Inconsistent A / B number of images") ;
        hdrl_imagelist_delete(in_a) ;
        hdrl_imagelist_delete(in_b) ;
        return -1 ;
    }

    /* Compute diff_a = in_a - in_b and diff_b = in_b - in_a */
    cpl_msg_info(__func__, "Compute difference image lists A-B and B-A") ;
    diff_a = hdrl_imagelist_duplicate(in_a) ;
    hdrl_imagelist_sub_imagelist(diff_a, in_b) ;
    diff_b = hdrl_imagelist_duplicate(in_b) ;
    hdrl_imagelist_sub_imagelist(diff_b, in_a) ;
    hdrl_imagelist_delete(in_a) ;
    hdrl_imagelist_delete(in_b) ;
    
    /* Collapse A-B and B-A */
    cpl_msg_info(__func__, "Collapse image lists A-B and B-A") ;
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

    /* Compute the slit fractions for A and B positions extraction */   
    cpl_msg_info(__func__, "Compute the slit fractions for A and B positions");
    cpl_msg_indent_more() ;
    /*
    The assumption is made here that :  
        - The slit center is exactly in the middle of A and B positions
        - The nodthrow is the distance in arcseconds between A and B
        - The slit size is 10 arcseconds
        - The B position is above the A position (--nodding-invert=false)
    */
    slit_length = 10 ;
    if ((plist = cpl_propertylist_load(first_fname, 0)) == NULL) {
        cpl_msg_error(__func__, "Cannot read the NODTHROW in the input files") ;
        cpl_msg_indent_less() ;
        hdrl_image_delete(collapsed_a) ;
        hdrl_image_delete(collapsed_b) ;
        cpl_table_delete(trace_wave) ;
        return -1 ;
    }
    nod_throw = cr2res_pfits_get_nodthrow(plist) ;
    cpl_propertylist_delete(plist) ;
    if (nod_throw < slit_length / 2.0) {
        extr_width_frac = nod_throw/slit_length ;  
        if (nodding_invert == 1) {
            slit_frac_a_bot = 0.5 ;
            slit_frac_a_mid = slit_frac_a_bot + extr_width_frac/2.0 ;
            slit_frac_a_top = slit_frac_a_bot + extr_width_frac ;
            slit_frac_b_top = 0.5 ;
            slit_frac_b_mid = slit_frac_b_top - extr_width_frac/2.0 ;
            slit_frac_b_bot = slit_frac_b_top - extr_width_frac ;
        } else {
            slit_frac_a_top = 0.5 ;
            slit_frac_a_mid = slit_frac_a_top - extr_width_frac/2.0 ;
            slit_frac_a_bot = slit_frac_a_top - extr_width_frac ;
            slit_frac_b_bot = 0.5 ;
            slit_frac_b_mid = slit_frac_b_bot + extr_width_frac/2.0 ;
            slit_frac_b_top = slit_frac_b_bot + extr_width_frac ;
        }
    } else if (nod_throw <= slit_length) {
        extr_width_frac = (slit_length - nod_throw)/slit_length ;  
        if (nodding_invert == 1) {
            slit_frac_a_top = 1.0 ;
            slit_frac_a_mid = slit_frac_a_top - extr_width_frac/2.0 ;
            slit_frac_a_bot = slit_frac_a_top - extr_width_frac ;
            slit_frac_b_bot = 0.0 ;
            slit_frac_b_mid = slit_frac_b_bot + extr_width_frac/2.0 ;
            slit_frac_b_top = slit_frac_b_bot + extr_width_frac ;
        } else {
            slit_frac_a_bot = 0.0 ;
            slit_frac_a_mid = slit_frac_a_bot + extr_width_frac/2.0 ;
            slit_frac_a_top = slit_frac_a_bot + extr_width_frac ;
            slit_frac_b_top = 1.0 ;
            slit_frac_b_mid = slit_frac_b_top - extr_width_frac/2.0 ;
            slit_frac_b_bot = slit_frac_b_top - extr_width_frac ;
        }
    } else {
        cpl_msg_error(__func__, "NODTHROW > slit length (%g>%g)- abort", 
                nod_throw, slit_length) ;
        cpl_msg_indent_less() ;
        hdrl_image_delete(collapsed_a) ;
        hdrl_image_delete(collapsed_b) ;
        cpl_table_delete(trace_wave) ;
        return -1 ;
    }
    cpl_msg_info(__func__, "Nod Throw : %g arcsecs", nod_throw) ;
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
    cpl_msg_indent_less() ;

    /* Recompute the traces for the new slit fractions */
    cpl_msg_info(__func__, "Recompute the traces for the A slit fraction") ;
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

    cpl_msg_info(__func__, "Recompute the traces for the B slit fraction") ;
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
   
    /* Load Blaze */
    blaze_table = NULL ;
    blaze_norm = 0;
    if (blaze_frame != NULL) {
        cpl_msg_info(__func__, "Load the BLAZE") ;
        if ((blaze_table = cr2res_io_load_EXTRACT_1D(cpl_frame_get_filename(
                            blaze_frame), reduce_det)) == NULL) {
            cpl_msg_error(__func__, "Failed to Load the Blaze file") ;
            hdrl_image_delete(collapsed_a) ;
            hdrl_image_delete(collapsed_b) ;
            cpl_table_delete(trace_wave) ;
            cpl_table_delete(trace_wave_a) ;
            cpl_table_delete(trace_wave_b) ;
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

    /* Execute the extraction */
    cpl_msg_info(__func__, "A position Spectra Extraction") ;
    cpl_msg_indent_more() ;
    if (cr2res_extract_traces(collapsed_a, trace_wave_a, NULL, blaze_table, blaze_norm, -1,
                -1, CR2RES_EXTR_OPT_CURV, extract_height, extract_swath_width, 
                extract_oversample, extract_smooth_slit, extract_smooth_spec,
                extract_niter, extract_kappa, error_factor,
                disp_det==reduce_det, disp_order_idx, disp_trace,
                &extracted_a, &slit_func_a, &model_master_a) == -1) {
        cpl_msg_error(__func__, "Failed to extract A");
        cpl_msg_indent_less() ;
        hdrl_image_delete(collapsed_a) ;
        hdrl_image_delete(collapsed_b) ;
        cpl_table_delete(trace_wave_a) ;
        cpl_table_delete(trace_wave_b) ;
        cpl_table_delete(trace_wave) ;
        return -1 ;
    }
    cpl_msg_indent_less() ;

    cpl_msg_info(__func__, "B position Spectra Extraction") ;
    cpl_msg_indent_more() ;
    if (cr2res_extract_traces(collapsed_b, trace_wave_b, NULL, blaze_table, blaze_norm, -1,
                -1, CR2RES_EXTR_OPT_CURV, extract_height, extract_swath_width, 
                extract_oversample, extract_smooth_slit, extract_smooth_spec,
                extract_niter, extract_kappa, error_factor,
                disp_det==reduce_det, disp_order_idx, disp_trace,
                &extracted_b, &slit_func_b, &model_master_b) == -1) {
        cpl_msg_error(__func__, "Failed to extract B");
        cpl_msg_indent_less() ;
        cpl_table_delete(extracted_a) ;
        cpl_table_delete(slit_func_a) ;
        hdrl_image_delete(model_master_a) ;
        hdrl_image_delete(collapsed_a) ;
        hdrl_image_delete(collapsed_b) ;
        cpl_table_delete(trace_wave_a) ;
        cpl_table_delete(trace_wave_b) ;
        cpl_table_delete(trace_wave) ;
        return -1 ;
    }
    cpl_msg_indent_less() ;
    
    if (blaze_table != NULL) cpl_table_delete(blaze_table) ;

    /* Combine both a and b extracted spectra together */
    cpl_msg_info(__func__, "A and B spectra combination") ;
    cpl_msg_indent_more() ;
    extracted_combined = cr2res_combine_extracted(extracted_a, extracted_b) ;
    cpl_msg_indent_less() ;

    /* Get the Setting */
    plist = cpl_propertylist_load(first_fname, 0) ;
    cur_setting = cpl_strdup(cr2res_pfits_get_wlen_id(plist)) ;
    cr2res_format_setting(cur_setting) ;
    cpl_propertylist_delete(plist) ;

    /* Store the extension header for product saving */
    plist = cpl_propertylist_load(first_fname,
            cr2res_io_get_ext_idx(first_fname, reduce_det, 1)) ;
                    
    /* QC - Signal and FWHM */
    cpl_msg_info(__func__, "QC parameters computation") ;
    qc_signal_a = cr2res_qc_obs_nodding_signal(extracted_a) ;
    qc_signal_b = cr2res_qc_obs_nodding_signal(extracted_b) ;
    cpl_propertylist_append_double(plist, CR2RES_HEADER_QC_SIGNAL, 
            (qc_signal_a+qc_signal_b)/2.0) ;

    /* QC STANDARD FLUX */
    qc_standard_flux_a =
        cr2res_qc_obs_nodding_standard_flux(extracted_a, cur_setting) ;
    qc_standard_flux_b =
        cr2res_qc_obs_nodding_standard_flux(extracted_b, cur_setting) ;
    cpl_propertylist_append_double(plist, CR2RES_HEADER_QC_STANDARD_FLUX, 
            (qc_standard_flux_a+qc_standard_flux_b)/2.0) ;
    cpl_free(cur_setting) ;

    /* QC - SNR on nodding A position */
    qc_snrs = cr2res_qc_snr(trace_wave_a, extracted_a, &order_idx_values,
            &nb_order_idx_values) ;
    qc_der_snrs = cr2res_qc_der_snr(trace_wave_a, extracted_a, &order_idx_values,
            &nb_order_idx_values) ;
    for (i=0 ; i<nb_order_idx_values ; i++) {
        order_idx = order_idx_values[i] ;
        order_idxp = cr2res_io_convert_order_idx_to_idxp(order_idx) ;
        key_name = cpl_sprintf(CR2RES_HEADER_QC_SNR, order_idxp) ;
        cpl_propertylist_append_double(plist, key_name, qc_snrs[i]) ;
        cpl_free(key_name) ;
        key_name = cpl_sprintf(CR2RES_HEADER_QC_DER_SNR, order_idxp) ;
        cpl_propertylist_append_double(plist, key_name, qc_der_snrs[i]) ;
        cpl_free(key_name) ;
    }
    cpl_free(order_idx_values) ;
    cpl_free(qc_snrs) ;

    /* Get the order numbers from the TW rows */
    order_idx_values = cr2res_trace_get_order_idx_values(trace_wave, 
            &nb_order_idx_values);

    /* QC - SLIT FWHM */
    fwhm_a_array = cpl_array_new(nb_order_idx_values, CPL_TYPE_DOUBLE);
    fwhm_b_array = cpl_array_new(nb_order_idx_values, CPL_TYPE_DOUBLE);
    for (i=0 ; i<nb_order_idx_values ; i++) {
        order_idx = order_idx_values[i] ;
        order_idxp = cr2res_io_convert_order_idx_to_idxp(order_idx) ;
        qc_fwhm_a = cr2res_qc_obs_slit_psf(slit_func_a, order_idxp,
            extract_oversample);
        qc_fwhm_b = cr2res_qc_obs_slit_psf(slit_func_b, order_idxp,
            extract_oversample);

        key_name = cpl_sprintf(CR2RES_HEADER_QC_SLITFWHM_ORDER, order_idxp) ;
        cpl_propertylist_append_double(plist, key_name,
                (qc_fwhm_a+qc_fwhm_b)/2.0) ;
        cpl_free(key_name) ;
        cpl_array_set(fwhm_a_array, i, qc_fwhm_a) ;
        cpl_array_set(fwhm_b_array, i, qc_fwhm_b) ;
    }
    qc_fwhm_a = cpl_array_get_median(fwhm_a_array);
    qc_fwhm_b = cpl_array_get_median(fwhm_b_array);
    qc_fwhm_med = (qc_fwhm_a+qc_fwhm_b)/2.0 ;
    cpl_propertylist_append_double(plist, CR2RES_HEADER_QC_SLITFWHM_MED,
            qc_fwhm_med) ;
    if (qc_fwhm_med < 3.5) {
        cpl_msg_warning(__func__, "Median FWHM of the PSF along the slit "
            "is %gpix, i.e. below the slit width. This means the slit "
            "is likely not evenly filled with light "
            "in the spectral direction. This can result in a "
            "wavelength offset between A and B nodding positions, and with "
            "respect to calibrations."
            , qc_fwhm_med);
    }
    cpl_array_delete(fwhm_a_array) ;
    cpl_array_delete(fwhm_b_array) ;

    /* QC - Real Orders  */
    if (order_zp > 0) {
        /* Compute the Real Order numbers and store them in QCs */
        for (i = 0; i < nb_order_idx_values; i++) {
            int order_real;
            order_idx = order_idx_values[i] ;
            order_idxp = cr2res_io_convert_order_idx_to_idxp(order_idx) ;

            order_real = cr2res_order_idx_to_real(order_idx, order_zp) ;
            key_name = cpl_sprintf(CR2RES_HEADER_QC_REAL_ORDER, order_idxp) ;
            cpl_propertylist_append_int(plist, key_name, order_real) ;
            cpl_free(key_name) ;
        }
    }
    cpl_table_delete(trace_wave) ;
    cpl_free(order_idx_values) ;

    /* Return */
    *combineda = collapsed_a ;
    *extracta = extracted_a ;
    *slitfunca = slit_func_a ;
    *modela = model_master_a ;
    *twa = trace_wave_a ;

    *combinedb = collapsed_b ;
    *extractb = extracted_b ;
    *slitfuncb = slit_func_b ;
    *modelb = model_master_b ;
    *twb = trace_wave_b ;

    *extractc = extracted_combined ;
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
    double                      nodthrow;
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
    for (i = 1; i < nframes; i++) {
        double nodthrow_cur;
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
 
/*----------------------------------------------------------------------------*/
/**
  @brief    Get the RAW frames from a frameset
  @param    set     Input frame set
  @param    type    [out]   0 - none,     
                            1 - OBS_NODDING, 
                            2 - CAL_NODDING, 
                            3 - OBS_ASTROMETRY
  @return   the RAW frameset or NULL in error case or if it is missing
    Allowed RAW types : CR2RES_CAL_NODDING_OTHER_RAW
                        CR2RES_CAL_NODDING_JITTER_RAW
                        CR2RES_OBS_NODDING_OTHER_RAW
                        CR2RES_OBS_NODDING_JITTER_RAW
                        CR2RES_OBS_ASTROMETRY_OTHER_RAW
                        CR2RES_OBS_ASTROMETRY_JITTER_RAW
 */
/*----------------------------------------------------------------------------*/
static cpl_frameset * cr2res_obs_nodding_find_RAW(
        const cpl_frameset  *   in,
        int                 *   type)
{
    cpl_frameset    *   out ;

    /* Check entries */
    if (in == NULL) return NULL ;

    out = cr2res_extract_frameset(in, CR2RES_OBS_NODDING_OTHER_RAW) ;
    if (out == NULL) {
        out = cr2res_extract_frameset(in, CR2RES_OBS_NODDING_JITTER_RAW) ;
    }
    if (out != NULL) {
        if (type != NULL) *type = 1 ;
    } else {
        out = cr2res_extract_frameset(in, CR2RES_CAL_NODDING_OTHER_RAW) ;
        if (out == NULL) {
            out = cr2res_extract_frameset(in, CR2RES_CAL_NODDING_JITTER_RAW) ;
        }
        if (out != NULL) {
            if (type != NULL) *type = 2 ;
        } else {
            out = cr2res_extract_frameset(in, CR2RES_OBS_ASTROMETRY_OTHER_RAW) ;
            if (out == NULL) {
                out = cr2res_extract_frameset(in, 
                        CR2RES_OBS_ASTROMETRY_JITTER_RAW) ;
            }
            if (out != NULL) {
                if (type != NULL) *type = 3 ;
            }
        }
    }
    if (out == NULL) {
        if (type != NULL) *type = 0 ;
    }
    return out ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Comparison function to identify different astrometry groups
  @param    frame1  first frame 
  @param    frame2  second frame 
  @return   0 if frame1!=frame2, 1 if frame1==frame2, -1 in error case
 */
/*----------------------------------------------------------------------------*/
static int cr2res_obs_nodding_astrometry_compare(
        const cpl_frame   *   frame1,
        const cpl_frame   *   frame2)
{
    int                     comparison ;
    cpl_propertylist    *   plist1 ;
    cpl_propertylist    *   plist2 ;
    double                    dval1, dval2 ;

    /* Test entries */
    if (frame1==NULL || frame2==NULL) return -1 ;

    /* Get property lists */
    if ((plist1=cpl_propertylist_load(cpl_frame_get_filename(frame1),0))==NULL){
        cpl_msg_error(__func__, "getting header from reference frame");
        return -1 ;
    }
    if ((plist2=cpl_propertylist_load(cpl_frame_get_filename(frame2),0))==NULL){
        cpl_msg_error(__func__, "getting header from reference frame");
        cpl_propertylist_delete(plist1) ;
        return -1 ;
    }

    /* Test status */
    if (cpl_error_get_code()) {
        cpl_propertylist_delete(plist1) ;
        cpl_propertylist_delete(plist2) ;
        return -1 ;
    }

    comparison = 1 ;

    /* Compare the DROT angle used */
    dval1 = cr2res_pfits_get_drot_posang(plist1) ;
    dval2 = cr2res_pfits_get_drot_posang(plist2) ;
    if (cpl_error_get_code()) {
        cpl_msg_error(__func__, "Cannot get the angle");
        cpl_propertylist_delete(plist1) ;
        cpl_propertylist_delete(plist2) ;
        return -1 ;
    }
    if (fabs(dval1-dval2) > 1e-3) comparison = 0 ;

    cpl_propertylist_delete(plist1) ;
    cpl_propertylist_delete(plist2) ;
    return comparison ;
}


