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
#include "cr2res_calib.h"
#include "cr2res_pfits.h"
#include "cr2res_dfs.h"
#include "cr2res_qc.h"
#include "cr2res_flat.h"
#include "cr2res_bpm.h"
#include "cr2res_trace.h"
#include "cr2res_extract.h"
#include "cr2res_io.h"

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

#define RECIPE_STRING "cr2res_cal_flat"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_cal_flat_compare(
        const cpl_frame   *   frame1,
        const cpl_frame   *   frame2) ;
static int cr2res_cal_flat_reduce(
        const cpl_frameset  *   rawframes,
        const cpl_frame     *   tw_frame,
        const cpl_frame     *   detlin_frame,
        const cpl_frame     *   master_dark_frame,
        const cpl_frame     *   bpm_frame,
        int                     filter_traces,
        int                     subtract_nolight_rows,
        int                     cosmics,
        double                  bpm_low,
        double                  bpm_high,
        double                  bpm_linemax,
        int                     trace_degree,
        int                     trace_min_cluster,
        int                     trace_smooth_x,
        int                     trace_smooth_y,
        double                  trace_threshold,
        int                     trace_opening,
        cr2res_extr_method      extr_method,
        int                     extract_oversample,
        int                     extract_swath_width,
        int                     extract_height,
        double                  extract_smooth_slit,
        double                  extract_smooth_spec,
        int                     reduce_det,
        int                     reduce_order,
        int                     reduce_trace,
        hdrl_image          **  master_flat,
        cpl_table           **  trace_wave,
        cpl_table           **  slit_func,
        cpl_table           **  extract_1d,
        hdrl_image          **  slit_model,
        cpl_image           **  bpm,
        cpl_propertylist    **  ext_plist) ;
static int cr2res_cal_flat_create(cpl_plugin *);
static int cr2res_cal_flat_exec(cpl_plugin *);
static int cr2res_cal_flat_destroy(cpl_plugin *);
static int cr2res_cal_flat(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_cal_flat_description[] = "\
Flat                                                                    \n\
  Compute the master flat and perform the traces detection.             \n\
  If a TW file is provided, its traces and slit curvatures are used     \n\
  for the extraction needed to create the master flat.                  \n\
  If not provided, the measured traces are used, and the slit is        \n\
  considered vertical.                                                  \n\
                                                                        \n\
  Inputs                                                                \n\
    raw.fits " CR2RES_FLAT_RAW " [1 to n]                               \n\
    trace.fits " CR2RES_CAL_FLAT_TW_PROCATG " [0 to 1]                  \n\
            or " CR2RES_CAL_FLAT_TW_MERGED_PROCATG "                    \n\
            or " CR2RES_UTIL_TRACE_TW_PROCATG "                         \n\
            or " CR2RES_UTIL_WAVE_TW_PROCATG "                          \n\
            or " CR2RES_CAL_WAVE_TW_PROCATG "                           \n\
            or " CR2RES_UTIL_SLIT_CURV_TW_PROCATG "                     \n\
    detlin.fits " CR2RES_CAL_DETLIN_COEFFS_PROCATG " [0 to 1]           \n\
    master_dark.fits " CR2RES_CAL_DARK_MASTER_PROCATG " [0 to 1]        \n\
    bpm.fits " CR2RES_CAL_DARK_BPM_PROCATG " [0 to 1]                   \n\
          or " CR2RES_CAL_FLAT_BPM_PROCATG "                            \n\
          or " CR2RES_CAL_DETLIN_BPM_PROCATG "                          \n\
          or " CR2RES_UTIL_BPM_SPLIT_PROCATG "                          \n\
          or " CR2RES_UTIL_NORM_BPM_PROCATG "                           \n\
                                                                        \n\
  Outputs                                                               \n\
    cr2res_cal_flat_[setting]_[Decker]_bpm.fits " 
    CR2RES_CAL_FLAT_BPM_PROCATG "\n\
    cr2res_cal_flat_[setting]_[Decker]_blaze.fits " 
    CR2RES_CAL_FLAT_EXTRACT_1D_PROCATG "\n\
    cr2res_cal_flat_[setting]_[Decker]_slit_model.fits " 
    CR2RES_CAL_FLAT_SLIT_MODEL_PROCATG "\n\
    cr2res_cal_flat_[setting]_[Decker]_slit_func.fits " 
    CR2RES_CAL_FLAT_SLIT_FUNC_PROCATG "\n\
    cr2res_cal_flat_[setting]_[Decker]_master_flat.fits " 
    CR2RES_CAL_FLAT_MASTER_PROCATG "\n\
    cr2res_cal_flat_[setting]_[Decker]_tw.fits " 
    CR2RES_CAL_FLAT_TW_PROCATG "\n\
    cr2res_cal_flat_[setting]_tw_merged.fits " 
    CR2RES_CAL_FLAT_TW_MERGED_PROCATG "\n\
                                                                        \n\
  Algorithm                                                             \n\
    group the input frames by different settings                        \n\
    loop on groups g:                                                   \n\
      loop on decker positions p:                                       \n\
        loop on detectors d:                                            \n\
          cr2res_cal_flat_reduce() computes (master_flat, trace_wave,   \n\
                 slit_func,extract_1d,slit_model,bpm)(g,p,d)            \n\
        Save slit_model(g,p)                                            \n\
        Save extract_1d(g,p)                                            \n\
        Save master_flat(g,p)                                           \n\
        Save trace_wave(g,p)                                            \n\
        Save slit_func(g,p)                                             \n\
        Save bpm(g,p)                                                   \n\
      Merge the trace_wave(g,p,d) into trace_wave(g,d)                  \n\
      Save trace_wave(g)                                                \n\
                                                                        \n\
    cr2res_cal_flat_reduce()                                            \n\
      Load the images list                                              \n\
      Apply the detlin / dark / bpm calibrations                        \n\
      Average the images to avg                                         \n\
      Compute the traces with cr2res_trace(--trace_degree,              \n\
                 --trace_min_cluster, --trace_smooth, --trace_opening)  \n\
                on avg                                                  \n\
      For the following step, the computed TW is used if there is no    \n\
      input TW provided, the provided one is used otherwise.            \n\
      loop on the traces t:                                             \n\
        cr2res_extract_slitdec_curved(--extract_oversample,             \n\
                 --extract_swath_width, --extract_height,               \n\
                 --extract_smooth_slit, --extract_smooth_spec)          \n\
          -> slit_func(t), extract_1d(t), slit_model(t)                 \n\
      Compute the master flat with cr2res_master_flat(avg, slit_model,  \n\
                 --bpm_low, --bpm_high, --bpm_lines_ratio)              \n\
        -> master_flat, bpm                                             \n\
      Merge the bpm with the input bpm                                  \n\
      Compute QCs                                                       \n\
      store the qc parameters in the returned property list             \n\
                                                                        \n\
  Library functions used:                                               \n\
    cr2res_io_extract_decker_frameset()                                 \n\
    cr2res_trace()                                                      \n\
    cr2res_extract_slitdec_curved()                                     \n\
    cr2res_master_flat()                                                \n\
    cr2res_qc_flat_lamp_ints()                                          \n\
    cr2res_qc_flat_mean_level()                                         \n\
    cr2res_qc_flat_mean_med_flux()                                      \n\
    cr2res_qc_flat_med_snr()                                            \n\
    cr2res_qc_overexposed()                                             \n\
    cr2res_qc_flat_trace_center_y()                                     \n\
    cr2res_trace_merge()                                                \n\
    cr2res_io_save_SLIT_MODEL()                                         \n\
    cr2res_io_save_EXTRACT_1D()                                         \n\
    cr2res_io_save_MASTER_FLAT()                                        \n\
    cr2res_io_save_TRACE_WAVE()                                         \n\
    cr2res_io_save_SLIT_FUNC()                                          \n\
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
                    "Flat recipe",
                    cr2res_cal_flat_description,
                    CR2RES_PIPELINE_AUTHORS,
                    PACKAGE_BUGREPORT,
                    cr2res_get_license(),
                    cr2res_cal_flat_create,
                    cr2res_cal_flat_exec,
                    cr2res_cal_flat_destroy)) {
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
static int cr2res_cal_flat_create(cpl_plugin * plugin)
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
    p = cpl_parameter_new_value("cr2res.cr2res_cal_flat.subtract_nolight_rows",
            CPL_TYPE_BOOL, "Subtract the no-light rows",
            "cr2res.cr2res_cal_flat", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "subtract_nolight_rows");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_flat.cosmics",
            CPL_TYPE_BOOL, "Find and mark cosmic rays hits as bad",
            "cr2res.cr2res_cal_flat", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "cosmics");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_flat.bpm_low",
            CPL_TYPE_DOUBLE, "Low threshold for BPM detection",
            "cr2res.cr2res_cal_flat", 0.8);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "bpm_low");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_flat.bpm_high",
            CPL_TYPE_DOUBLE, "High threshold for BPM detection",
            "cr2res.cr2res_cal_flat", 1.2);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "bpm_high");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_flat.bpm_lines_ratio",
            CPL_TYPE_DOUBLE, "Maximum ratio of bad pixels per line",
            "cr2res.cr2res_cal_flat", 0.5);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "bpm_lines_ratio");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_flat.trace_degree",
            CPL_TYPE_INT, "polynomial degree for the fit to the orders",
            "cr2res.cr2res_cal_flat", 2);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "trace_degree");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_flat.trace_min_cluster",
            CPL_TYPE_INT, "size in pixels of the smallest allowed cluster",
            "cr2res.cr2res_cal_flat", 100000);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "trace_min_cluster");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_flat.trace_smooth_x",
            CPL_TYPE_INT, "Length of the smoothing kernel in x",
            "cr2res.cr2res_cal_flat", 111);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "trace_smooth_x");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_flat.trace_smooth_y",
            CPL_TYPE_INT, "Length of the smoothing kernel in y",
            "cr2res.cr2res_cal_flat", 401);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "trace_smooth_y");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_flat.trace_threshold",
            CPL_TYPE_DOUBLE, "Detection Threshold",
            "cr2res.cr2res_cal_flat", 1.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "trace_threshold");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_flat.trace_opening",
            CPL_TYPE_BOOL, "Use a morphological opening to rejoin clusters",
            "cr2res.cr2res_cal_flat", TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "trace_opening");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /*
    p = cpl_parameter_new_value("cr2res.cr2res_cal_flat.trace_filter",
            CPL_TYPE_BOOL, "Only keep the predefined order traces",
            "cr2res.cr2res_cal_flat", TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "trace_filter");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);
    */

    p = cpl_parameter_new_value("cr2res.cr2res_cal_flat.extract_method",
            CPL_TYPE_STRING, "Extraction method (SUM / MEDIAN / TILTSUM / "
            "OPT_CURV )",
            "cr2res.cr2res_util_extract", "OPT_CURV");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_method");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_flat.extract_oversample",
            CPL_TYPE_INT, "factor by which to oversample the extraction",
            "cr2res.cr2res_cal_flat", 5);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_oversample");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_flat.extract_swath_width",
            CPL_TYPE_INT, "The swath width", "cr2res.cr2res_cal_flat", 800);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_swath_width");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_flat.extract_height",
            CPL_TYPE_INT, "Extraction height",
            "cr2res.cr2res_cal_flat", -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_height");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_flat.extract_smooth_slit",
            CPL_TYPE_DOUBLE,
            "Smoothing along the slit",
            "cr2res.cr2res_cal_flat", 3.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_smooth_slit");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_flat.extract_smooth_spec",
            CPL_TYPE_DOUBLE,
            "Smoothing along the spectrum",
            "cr2res.cr2res_cal_flat", 2.0) ;
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_smooth_spec");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_flat.detector",
            CPL_TYPE_INT, "Only reduce the specified detector",
            "cr2res.cr2res_cal_flat", 0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "detector");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_flat.order",
            CPL_TYPE_INT, "Only reduce the specified order",
            "cr2res.cr2res_cal_flat", -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "order");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_flat.trace_nb",
            CPL_TYPE_INT, "Only reduce the specified trace number",
            "cr2res.cr2res_cal_flat", -1);
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
static int cr2res_cal_flat_exec(cpl_plugin * plugin)
{
    cpl_recipe  *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1;

    return cr2res_cal_flat(recipe->frames, recipe->parameters);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_cal_flat_destroy(cpl_plugin * plugin)
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
static int cr2res_cal_flat(
        cpl_frameset            *   frameset,
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter *   param ;
    int                     subtract_nolight_rows, cosmics, 
                            trace_degree, trace_min_cluster,
                            trace_opening, trace_filter,
                            extract_oversample, extract_swath_width,
                            extract_height, reduce_det, reduce_order,
                            reduce_trace, trace_smooth_x, trace_smooth_y ;
    double                  bpm_low, bpm_high, bpm_lines_ratio,
                            trace_threshold, extract_smooth_slit,
                            extract_smooth_spec ;
    cr2res_extr_method      extr_method;
    const char          *   sval ;
    const cpl_frame     *   trace_wave_frame ;
    const cpl_frame     *   detlin_frame ;
    const cpl_frame     *   master_dark_frame ;
    const cpl_frame     *   bpm_frame ;
    cpl_frameset        *   rawframes ;
    const char          *   used_tag ;
    cpl_frameset        *   raw_one_setting_decker ;
    cpl_size            *   labels ;
    cpl_size                nlabels ;
    cpl_propertylist    *   qc_main ;

    hdrl_image          *   master_flat[CR2RES_NB_DETECTORS] ;
    cpl_table * trace_wave[CR2RES_NB_DECKER_POSITIONS][CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *
        ext_plist[CR2RES_NB_DECKER_POSITIONS][CR2RES_NB_DETECTORS] ;
    cpl_table           *   trace_wave_merged[CR2RES_NB_DETECTORS];
    cpl_table           *   merged ;
    cpl_table           *   trace_wave1 ;
    cpl_table           *   trace_wave2 ;
    cpl_table           *   slit_func[CR2RES_NB_DETECTORS] ;
    cpl_table           *   extract_1d[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   slit_model[CR2RES_NB_DETECTORS] ;
    cpl_image           *   bpm[CR2RES_NB_DETECTORS] ;
    char                *   qc_name;
    char                *   out_file;
    int                 *   orders ;
    cpl_vector          *   qc_mean,
                        *   qc_median,
                        *   qc_flux,
                        *   qc_rms,
                        *   qc_s2n,
                        *   qc_nbbad,
                        *   qc_centery,
                        *   qc_bltot,
                        *   qc_blgood,
                        *   qc_orderpos,
                        *   qc_overexposed ;
    int                     l, i, j, det_nr, nb_orders;

    /* Initialise */
    cr2res_decker decker_values[CR2RES_NB_DECKER_POSITIONS] =
        {CR2RES_DECKER_NONE, CR2RES_DECKER_1_3, CR2RES_DECKER_2_4} ;
    char * decker_desc[CR2RES_NB_DECKER_POSITIONS] =
        {"Open", "Decker1", "Decker2"} ;

    /* RETRIEVE INPUT PARAMETERS */
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_flat.subtract_nolight_rows");
    subtract_nolight_rows = cpl_parameter_get_bool(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_flat.cosmics");
    cosmics = cpl_parameter_get_bool(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_flat.bpm_low");
    bpm_low = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_flat.bpm_high");
    bpm_high = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_flat.bpm_lines_ratio");
    bpm_lines_ratio = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_flat.trace_degree");
    trace_degree = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_flat.trace_min_cluster");
    trace_min_cluster = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_flat.trace_smooth_x");
    trace_smooth_x = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_flat.trace_smooth_y");
    trace_smooth_y = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_flat.trace_threshold");
    trace_threshold = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_flat.trace_opening");
    trace_opening = cpl_parameter_get_bool(param);
    /*
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_flat.trace_filter");
    trace_filter = cpl_parameter_get_bool(param);
    */
    trace_filter = 0 ;
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_flat.extract_method");
    sval = cpl_parameter_get_string(param);
    if (!strcmp(sval, ""))              extr_method = CR2RES_EXTR_OPT_CURV;
    else if (!strcmp(sval, "OPT_CURV")) extr_method = CR2RES_EXTR_OPT_CURV;
    else if (!strcmp(sval, "SUM"))      extr_method = CR2RES_EXTR_SUM;
    else if (!strcmp(sval, "MEDIAN"))   extr_method = CR2RES_EXTR_MEDIAN;
    else if (!strcmp(sval, "TILTSUM"))  extr_method = CR2RES_EXTR_TILTSUM;
    else {
        cpl_msg_error(__func__, "Invalid Extraction Method specified");
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1;
    }
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_flat.extract_oversample");
    extract_oversample = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_flat.extract_swath_width");
    extract_swath_width = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_flat.extract_height");
    extract_height = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_flat.extract_smooth_slit");
    extract_smooth_slit = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_flat.extract_smooth_spec");
    extract_smooth_spec = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_flat.detector");
    reduce_det = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_flat.order");
    reduce_order = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_flat.trace_nb");
    reduce_trace = cpl_parameter_get_int(param);

    /* Identify the RAW and CALIB frames in the input frameset */
    if (cr2res_dfs_set_groups(frameset)) {
        cpl_msg_error(__func__, "Cannot identify RAW and CALIB frames") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Get Calibration frames */
    trace_wave_frame = cr2res_io_find_TRACE_WAVE(frameset) ;
    detlin_frame = cpl_frameset_find_const(frameset,
            CR2RES_CAL_DETLIN_COEFFS_PROCATG);
    master_dark_frame = cpl_frameset_find_const(frameset,
            CR2RES_CAL_DARK_MASTER_PROCATG) ;
    bpm_frame = cr2res_io_find_BPM(frameset) ;

    /* Extract RAW frames */
    used_tag = CR2RES_FLAT_RAW ;
    rawframes = cr2res_extract_frameset(frameset, used_tag) ;
    if (rawframes==NULL || cpl_frameset_get_size(rawframes) <= 0) {
        cpl_msg_error(__func__, "Cannot find any RAW file") ;
        cpl_error_set(__func__, CPL_ERROR_DATA_NOT_FOUND) ;
        return -1 ;
    }

    /* Label the raw frames with the different settings */
    if ((labels = cpl_frameset_labelise(rawframes, cr2res_cal_flat_compare,
                &nlabels)) == NULL) {
        cpl_msg_error(__func__, "Cannot label input frames") ;
        cpl_frameset_delete(rawframes) ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Loop on the settings */
    for (l = 0; l < (int)nlabels; l++) {
        cpl_frameset *raw_one_setting;
        cpl_propertylist *plist;
        char *setting_id;
        /* Get the frames for the current setting */
        raw_one_setting = cpl_frameset_extract(rawframes, labels, (cpl_size)l) ;

        /* Get the current setting */
        plist = cpl_propertylist_load(cpl_frame_get_filename(
                    cpl_frameset_get_position(raw_one_setting, 0)), 0) ;
        setting_id = cpl_strdup(cr2res_pfits_get_wlen_id(plist)) ;
        cr2res_format_setting(setting_id) ;
        cpl_propertylist_delete(plist) ;
        
        cpl_msg_info(__func__, "Process SETTING %s", setting_id) ;
        cpl_msg_indent_more() ;

        /* Loop on the decker positions */
        for (i=0 ; i<CR2RES_NB_DECKER_POSITIONS ; i++) {
            /* Initialise */
            for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
                trace_wave[i][det_nr-1] = NULL ;
                ext_plist[i][det_nr-1] = NULL ;
            }

            /* Get the Frames for the current decker position */
            raw_one_setting_decker = cr2res_io_extract_decker_frameset(
                    raw_one_setting, used_tag, decker_values[i]) ;
            if (raw_one_setting_decker == NULL) {
                cpl_msg_info(__func__, "No files for decker: %s",
                        decker_desc[i]) ;
                continue ;
            }
            cpl_msg_info(__func__, "Reduce %s Frames", decker_desc[i]) ;
            cpl_msg_indent_more() ;

            /* Loop on the detectors */
            for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
                /* Initialise */
                master_flat[det_nr-1] = NULL ;
                slit_func[det_nr-1] = NULL ;
                extract_1d[det_nr-1] = NULL ;
                slit_model[det_nr-1] = NULL ;
                bpm[det_nr-1] = NULL ;

                /* Compute only one detector */
                if (reduce_det != 0 && det_nr != reduce_det) continue ;

                cpl_msg_info(__func__, "Process Detector %d", det_nr) ;
                cpl_msg_indent_more() ;

                /* Call the reduction function */
                if (cr2res_cal_flat_reduce(raw_one_setting_decker,
                            trace_wave_frame, detlin_frame, master_dark_frame, 
                            bpm_frame, trace_filter, subtract_nolight_rows, 
                            cosmics, bpm_low, bpm_high, 
                            bpm_lines_ratio, trace_degree, trace_min_cluster, 
                            trace_smooth_x, trace_smooth_y, trace_threshold, 
                            trace_opening, extr_method, extract_oversample, 
                            extract_swath_width, extract_height, 
                            extract_smooth_slit, extract_smooth_spec, det_nr, 
                            reduce_order, reduce_trace,
                            &(master_flat[det_nr-1]),
                            &(trace_wave[i][det_nr-1]),
                            &(slit_func[det_nr-1]),
                            &(extract_1d[det_nr-1]),
                            &(slit_model[det_nr-1]),
                            &(bpm[det_nr-1]),
                            &(ext_plist[i][det_nr-1])) == -1) {
                    cpl_msg_warning(__func__,
                            "Failed to reduce detector %d of %s Frames",
                            det_nr, decker_desc[i]);
                    cpl_error_reset() ;
                }
                cpl_msg_indent_less() ;
            }


            /* Compute QCs for the main header */
            qc_main = NULL ;
            if (reduce_det == 0) {
                qc_mean = cpl_vector_new(CR2RES_NB_DETECTORS) ;
                qc_median = cpl_vector_new(CR2RES_NB_DETECTORS) ;
                qc_flux = cpl_vector_new(CR2RES_NB_DETECTORS) ;
                qc_rms = cpl_vector_new(CR2RES_NB_DETECTORS) ;
                qc_s2n = cpl_vector_new(CR2RES_NB_DETECTORS) ;
                qc_nbbad = cpl_vector_new(CR2RES_NB_DETECTORS) ;
                qc_centery = cpl_vector_new(CR2RES_NB_DETECTORS) ;
                qc_blgood = cpl_vector_new(CR2RES_NB_DETECTORS) ;
                qc_bltot = cpl_vector_new(CR2RES_NB_DETECTORS) ;
                for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
                    cpl_vector_set(qc_mean, det_nr-1, 
                            cpl_propertylist_get_double(ext_plist[i][det_nr-1],
                                CR2RES_HEADER_QC_FLAT_MEAN)) ;
                    cpl_vector_set(qc_median, det_nr-1, 
                            cpl_propertylist_get_double(ext_plist[i][det_nr-1],
                                CR2RES_HEADER_QC_FLAT_MEDIAN)) ;
                    cpl_vector_set(qc_flux, det_nr-1, 
                            cpl_propertylist_get_double(ext_plist[i][det_nr-1],
                                CR2RES_HEADER_QC_FLAT_FLUX)) ;
                    cpl_vector_set(qc_rms, det_nr-1, 
                            cpl_propertylist_get_double(ext_plist[i][det_nr-1],
                                CR2RES_HEADER_QC_FLAT_RMS)) ;
                    cpl_vector_set(qc_s2n, det_nr-1, 
                            cpl_propertylist_get_double(ext_plist[i][det_nr-1],
                                CR2RES_HEADER_QC_FLAT_S2N)) ;
                    cpl_vector_set(qc_nbbad, det_nr-1, 
                        (double)cpl_propertylist_get_int(ext_plist[i][det_nr-1],
                                CR2RES_HEADER_QC_FLAT_NBBAD)) ;
                    cpl_vector_set(qc_centery, det_nr-1, 
                            cpl_propertylist_get_double(ext_plist[i][det_nr-1],
                                CR2RES_HEADER_QC_FLAT_CENTERY)) ;
                    cpl_vector_set(qc_blgood, det_nr-1, 
                        (double)cpl_propertylist_get_int(ext_plist[i][det_nr-1],
                                CR2RES_HEADER_QC_BLAZE_NGOOD)) ;
                    cpl_vector_set(qc_bltot, det_nr-1, 
                            cpl_propertylist_get_double(ext_plist[i][det_nr-1],
                                CR2RES_HEADER_QC_BLAZE_TOT)) ;
                }

                qc_main = cpl_propertylist_new() ;
				cpl_propertylist_append_double(qc_main,
						CR2RES_HEADER_QC_FLAT_MEAN_AVG,
						cpl_vector_get_mean(qc_mean)) ;
				cpl_propertylist_append_double(qc_main,
						CR2RES_HEADER_QC_FLAT_MEAN_RMS,
						cpl_vector_get_stdev(qc_mean)) ;
				cpl_propertylist_append_double(qc_main,
						CR2RES_HEADER_QC_FLAT_MEDIAN_AVG,
						cpl_vector_get_mean(qc_median)) ;
				cpl_propertylist_append_double(qc_main,
						CR2RES_HEADER_QC_FLAT_MEDIAN_RMS,
						cpl_vector_get_stdev(qc_median)) ;
				cpl_propertylist_append_double(qc_main,
						CR2RES_HEADER_QC_FLAT_FLUX_AVG,
						cpl_vector_get_mean(qc_flux)) ;
				cpl_propertylist_append_double(qc_main,
						CR2RES_HEADER_QC_FLAT_FLUX_RMS,
						cpl_vector_get_stdev(qc_flux)) ;
				cpl_propertylist_append_double(qc_main,
						CR2RES_HEADER_QC_FLAT_RMS_AVG,
						cpl_vector_get_mean(qc_rms)) ;
				cpl_propertylist_append_double(qc_main,
						CR2RES_HEADER_QC_FLAT_RMS_RMS,
						cpl_vector_get_stdev(qc_rms)) ;
				cpl_propertylist_append_double(qc_main,
						CR2RES_HEADER_QC_FLAT_S2N_AVG,
						cpl_vector_get_mean(qc_s2n)) ;
				cpl_propertylist_append_double(qc_main,
						CR2RES_HEADER_QC_FLAT_S2N_RMS,
						cpl_vector_get_stdev(qc_s2n)) ;
				cpl_propertylist_append_double(qc_main,
						CR2RES_HEADER_QC_FLAT_NBBAD_AVG,
						cpl_vector_get_mean(qc_nbbad)) ;
				cpl_propertylist_append_double(qc_main,
						CR2RES_HEADER_QC_FLAT_NBBAD_RMS,
						cpl_vector_get_stdev(qc_nbbad)) ;
				cpl_propertylist_append_double(qc_main,
						CR2RES_HEADER_QC_FLAT_CENTERY_AVG,
						cpl_vector_get_mean(qc_centery)) ;
				cpl_propertylist_append_double(qc_main,
						CR2RES_HEADER_QC_FLAT_CENTERY_RMS,
						cpl_vector_get_stdev(qc_centery)) ;
				cpl_propertylist_append_double(qc_main,
						CR2RES_HEADER_QC_BLAZE_NORM,
						cpl_vector_get_sum(qc_bltot)/cpl_vector_get_sum(qc_blgood)) ;
                cpl_vector_delete(qc_mean) ;
                cpl_vector_delete(qc_median) ;
                cpl_vector_delete(qc_flux) ;
                cpl_vector_delete(qc_rms) ;
                cpl_vector_delete(qc_s2n) ;
                cpl_vector_delete(qc_nbbad) ;
                cpl_vector_delete(qc_centery) ;
                cpl_vector_delete(qc_blgood) ;
                cpl_vector_delete(qc_bltot) ;

                if (trace_wave[i][0] != NULL) {
                    /* Special treatment for ORDERPOSn and OVEREXPOSEDn */
                    /* Get all orders - Use the TW from detector 1 */
                    orders = cr2res_trace_get_order_idx_values(
                            trace_wave[i][0], &nb_orders) ;
                    /* Loop on the orders */
                    for (j=0 ; j<nb_orders ; j++) {
                        qc_orderpos = cpl_vector_new(CR2RES_NB_DETECTORS) ;
                        qc_overexposed = cpl_vector_new(CR2RES_NB_DETECTORS) ;
                        cpl_vector_fill(qc_orderpos, 0.0) ;
                        cpl_vector_fill(qc_overexposed, 0.0) ;
                        for (det_nr = 1; det_nr <= CR2RES_NB_DETECTORS; det_nr++)
                        {
                                // OVEREXPOSEDn
                                qc_name = cpl_sprintf("%s%02d",
                                                      CR2RES_HEADER_QC_OVEREXPOSED, orders[j]);
                                if (cpl_propertylist_has(ext_plist[i][det_nr - 1], qc_name))
                                {
                                        cpl_vector_set(qc_overexposed, det_nr - 1,
                                                       cpl_propertylist_get_double(
                                                           ext_plist[i][det_nr - 1], qc_name));
                                }
                                cpl_free(qc_name);
                                // ORDERPOSn
                                qc_name = cpl_sprintf("%s%02d",
                                                      CR2RES_HEADER_QC_FLAT_ORDERPOS, orders[j]);
                                if (cpl_propertylist_has(ext_plist[i][det_nr - 1], qc_name))
                                {
                                        cpl_vector_set(qc_orderpos, det_nr - 1,
                                                       cpl_propertylist_get_double(
                                                           ext_plist[i][det_nr - 1], qc_name));
                                }
                                cpl_free(qc_name);
                        }

                        // OVEREXPOSEDn.AVG/RMS
                        qc_name = cpl_sprintf("%s%02d AVG",
                                CR2RES_HEADER_QC_OVEREXPOSED, orders[j]) ;
                        cpl_propertylist_append_double(qc_main,
                                qc_name, cpl_vector_get_mean(qc_overexposed)) ;
                        cpl_free(qc_name) ;

                        qc_name = cpl_sprintf("%s%02d RMS",
                                CR2RES_HEADER_QC_OVEREXPOSED, orders[j]) ;
                        cpl_propertylist_append_double(qc_main,
                                qc_name, cpl_vector_get_stdev(qc_overexposed)) ;
                        cpl_free(qc_name) ;
                        cpl_vector_delete(qc_overexposed) ;

                        // ORDERPOSn.AVG/RMS
                        qc_name = cpl_sprintf("%s%02d AVG",
                                CR2RES_HEADER_QC_FLAT_ORDERPOS, orders[j]) ;
                        cpl_propertylist_append_double(qc_main,
                                qc_name, cpl_vector_get_mean(qc_orderpos)) ;
                        cpl_free(qc_name) ;

                        qc_name = cpl_sprintf("%s%02d RMS",
                                CR2RES_HEADER_QC_FLAT_ORDERPOS, orders[j]) ;
                        cpl_propertylist_append_double(qc_main,
                                qc_name, cpl_vector_get_stdev(qc_orderpos)) ;
                        cpl_free(qc_name) ;
                        cpl_vector_delete(qc_orderpos) ;
                    }
                    cpl_free(orders) ;
                }
            }

            /* Save Products */

            /* Save only the used RAW ? : raw_one_... instead of 2nd frameset */
            /* Beware that the calibration PRO RECi CAL will be missing */

            /* MASTER_FLAT */
            if (nlabels == 1) {
                out_file = cpl_sprintf("%s_%s_master_flat.fits", 
                        RECIPE_STRING, decker_desc[i]) ;
            } else {
                out_file = cpl_sprintf("%s_%s_%s_master_flat.fits", 
                        RECIPE_STRING, setting_id, decker_desc[i]) ;
            }
            cr2res_io_save_MASTER_FLAT(out_file, frameset,
                    frameset, parlist, master_flat, qc_main, 
                    ext_plist[i], CR2RES_CAL_FLAT_MASTER_PROCATG,
                    RECIPE_STRING);
            cpl_free(out_file);

            /* SLIT_MODEL */
            if (nlabels == 1) {
                out_file = cpl_sprintf("%s_%s_slit_model.fits", 
                        RECIPE_STRING, decker_desc[i]) ;
            } else {
                out_file = cpl_sprintf("%s_%s_%s_slit_model.fits", 
                        RECIPE_STRING, setting_id, decker_desc[i]) ;
            }
            cr2res_io_save_SLIT_MODEL(out_file, frameset,
                    frameset, parlist, slit_model, qc_main, 
                    ext_plist[i], CR2RES_CAL_FLAT_SLIT_MODEL_PROCATG,
                    RECIPE_STRING);
            cpl_free(out_file);

            /* BLAZE */
            if (nlabels == 1) {
                out_file = cpl_sprintf("%s_%s_blaze.fits", 
                        RECIPE_STRING, decker_desc[i]) ;
            } else {
                out_file = cpl_sprintf("%s_%s_%s_blaze.fits", 
                        RECIPE_STRING, setting_id, decker_desc[i]) ;
            }
            cr2res_io_save_EXTRACT_1D(out_file, frameset, 
                    frameset, parlist, extract_1d, qc_main, 
                    ext_plist[i], CR2RES_CAL_FLAT_EXTRACT_1D_PROCATG,
                    RECIPE_STRING);
            cpl_free(out_file);

            /* TRACE_WAVE */
            if (nlabels == 1) {
                out_file = cpl_sprintf("%s_%s_tw.fits", 
                        RECIPE_STRING, decker_desc[i]) ;
            } else {
                out_file = cpl_sprintf("%s_%s_%s_tw.fits", 
                        RECIPE_STRING, setting_id, decker_desc[i]) ;
            }
            cr2res_io_save_TRACE_WAVE(out_file, frameset,
                    frameset, parlist, trace_wave[i], qc_main, 
                    ext_plist[i], CR2RES_CAL_FLAT_TW_PROCATG, RECIPE_STRING);
            cpl_free(out_file);

            /* SLIT_FUNC */
            if (nlabels == 1) {
                out_file = cpl_sprintf("%s_%s_slit_func.fits", 
                        RECIPE_STRING, decker_desc[i]) ;
            } else {
                out_file = cpl_sprintf("%s_%s_%s_slit_func.fits", 
                        RECIPE_STRING, setting_id, decker_desc[i]) ;
            }
            cr2res_io_save_SLIT_FUNC(out_file, frameset,
                    frameset, parlist, slit_func, qc_main, 
                    ext_plist[i], CR2RES_CAL_FLAT_SLIT_FUNC_PROCATG, 
                    RECIPE_STRING);
            cpl_free(out_file);

            /* BPM */
            if (nlabels == 1) {
                out_file = cpl_sprintf("%s_%s_bpm.fits", 
                        RECIPE_STRING, decker_desc[i]) ;
            } else {
                out_file = cpl_sprintf("%s_%s_%s_bpm.fits", 
                        RECIPE_STRING, setting_id, decker_desc[i]) ;
            }
            cr2res_io_save_BPM(out_file, frameset,
                    frameset, parlist, bpm, qc_main, ext_plist[i], 
                    CR2RES_CAL_FLAT_BPM_PROCATG, RECIPE_STRING) ;
            cpl_free(out_file);

            /* Free */
            if (qc_main != NULL) cpl_propertylist_delete(qc_main) ;
            cpl_frameset_delete(raw_one_setting_decker) ;
            for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
                if (master_flat[det_nr-1] != NULL)
                    hdrl_image_delete(master_flat[det_nr-1]) ;
                if (slit_func[det_nr-1] != NULL)
                    cpl_table_delete(slit_func[det_nr-1]) ;
                if (extract_1d[det_nr-1] != NULL)
                    cpl_table_delete(extract_1d[det_nr-1]) ;
                if (slit_model[det_nr-1] != NULL)
                    hdrl_image_delete(slit_model[det_nr-1]) ;
                if (bpm[det_nr-1] != NULL)
                    cpl_image_delete(bpm[det_nr-1]) ;
            }
            cpl_msg_indent_less() ;
        }

        /* Merge the Decker positions TRACE_WAVE files in a single one */
        for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
            /* Initialise */
            if (trace_wave[0][det_nr-1] != NULL)
                trace_wave_merged[det_nr-1] =
                    cpl_table_duplicate(trace_wave[0][det_nr-1]) ;
            else
                trace_wave_merged[det_nr-1] = NULL ;
            /* Loop on the other detectors */
            for (i=1 ; i<CR2RES_NB_DECKER_POSITIONS ; i++) {
                trace_wave1 = trace_wave_merged[det_nr-1] ;
                trace_wave2 = trace_wave[i][det_nr-1] ;
                if (trace_wave1 == NULL && trace_wave2 == NULL) {
                    /* Do nothing - go to next iteration */
                    trace_wave_merged[det_nr-1] = NULL ;
                } else if (trace_wave1 == NULL && trace_wave2 != NULL) {
                    trace_wave_merged[det_nr-1] = 
                        cpl_table_duplicate(trace_wave2) ;
                } else if (trace_wave1 != NULL && trace_wave2 == NULL) {
                    /* Do nothing - go to next iteration */
                } else {
                    merged = cr2res_trace_merge(trace_wave1, trace_wave2) ;
                    if (merged == NULL) {
                        cpl_msg_error(__func__, "Failed merging") ;
                    } else {
                        cpl_table_delete(trace_wave_merged[det_nr-1]) ;
                        trace_wave_merged[det_nr-1] = merged ;
                    }
                }
            }
        }

        /* Free the traces */
        for (i=0 ; i<CR2RES_NB_DECKER_POSITIONS ; i++)
            for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++)
                if (trace_wave[i][det_nr-1] != NULL)
                    cpl_table_delete(trace_wave[i][det_nr-1]) ;


        /* Save TRACE_WAVE_MERGED */
        if (nlabels == 1) {
            out_file = cpl_sprintf("%s_tw_merged.fits", 
                    RECIPE_STRING) ;
        } else {
            out_file = cpl_sprintf("%s_%s_tw_merged.fits", RECIPE_STRING, 
                    setting_id) ;
        }

        /* Save only the used RAW ? : raw_one_... instead of 2nd frameset */
        /* Beware that the calibration PRO RECi CAL will be missing */

        cr2res_io_save_TRACE_WAVE(out_file, frameset, frameset, parlist,
                trace_wave_merged, NULL, ext_plist[0],
                CR2RES_CAL_FLAT_TW_MERGED_PROCATG, RECIPE_STRING) ;
        cpl_free(out_file);
        cpl_free(setting_id) ;

        /* Free extensions */
        for (i=0 ; i<CR2RES_NB_DECKER_POSITIONS ; i++)
            for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++)
                if (ext_plist[i][det_nr-1] != NULL)
                    cpl_propertylist_delete(ext_plist[i][det_nr-1]) ;

        /* Free the merged traces */
        for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++)
            if (trace_wave_merged[det_nr-1] != NULL)
                cpl_table_delete(trace_wave_merged[det_nr-1]) ;

        cpl_frameset_delete(raw_one_setting) ;
        cpl_msg_indent_less() ;
    }
    cpl_free(labels);
    cpl_frameset_delete(rawframes) ;

    return (int)cpl_error_get_code();
}

/*----------------------------------------------------------------------------*/
/**
  @brief Compute the flat for 1 setting, 1 decker position, 1 detector
  @param rawframes          Input raw frames (same setting, same decker)
  @param tw_frame           Optional TW frame or NULL
  @param detlin_frame       Associated detlin coefficients
  @param master_dark_frame  Associated master dark
  @param bpm_frame          Associated BPM
  @param filter_traces      Flag to Filter out traces
  @param subtract_nolight_rows
  @param cosmics            Flag to correct for cosmics
  @param bpm_low            Threshold for BPM detection
  @param bpm_high           Threshold for BPM detection
  @param bpm_linemax        Max fraction of BPM per line
  @param trace_degree       Trace computation related
  @param trace_min_cluster  Trace computation related
  @param trace_smooth_x     Trace computation related
  @param trace_smooth_y     Trace computation related
  @param trace_threshold    Trace computation related
  @param trace_opening      Trace computation related
  @param extr_method        The wished extraction method
  @param extract_oversample Extraction related
  @param extract_swath_width Extraction related
  @param extract_height     Extraction related
  @param extract_smooth_slit     Extraction related
  @param extract_smooth_spec     Extraction related
  @param reduce_det         The detector to compute
  @param reduce_order       The order to compute (-1 for all)
  @param reduce_trace       The trace to compute (-1 for all)
  @param master_flat        [out] Master flat
  @param trace_wave         [out] trace wave
  @param slit_func          [out] slit function
  @param extract_1d         [out] the blaze
  @param slit_model         [out] slit model
  @param bpm                [out] the BPM
  @param ext_plist          [out] the header for saving the products
  @return   0 if ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
static int cr2res_cal_flat_reduce(
        const cpl_frameset  *   rawframes,
        const cpl_frame     *   tw_frame,
        const cpl_frame     *   detlin_frame,
        const cpl_frame     *   master_dark_frame,
        const cpl_frame     *   bpm_frame,
        int                     filter_traces,
        int                     subtract_nolight_rows,
        int                     cosmics,
        double                  bpm_low,
        double                  bpm_high,
        double                  bpm_linemax,
        int                     trace_degree,
        int                     trace_min_cluster,
        int                     trace_smooth_x,
        int                     trace_smooth_y,
        double                  trace_threshold,
        int                     trace_opening,
        cr2res_extr_method      extr_method,
        int                     extract_oversample,
        int                     extract_swath_width,
        int                     extract_height,
        double                  extract_smooth_slit,
        double                  extract_smooth_spec,
        int                     reduce_det,
        int                     reduce_order,
        int                     reduce_trace,
        hdrl_image          **  master_flat,
        cpl_table           **  trace_wave,
        cpl_table           **  slit_func,
        cpl_table           **  extract_1d,
        hdrl_image          **  slit_model,
        cpl_image           **  bpm,
        cpl_propertylist    **  ext_plist)
{
    const char          *   first_file ;
    hdrl_image          *   first_image ;
    hdrl_imagelist      *   imlist ;
    hdrl_imagelist      *   imlist_calibrated ;
    hdrl_image          *   collapsed ;
    cpl_image           *   contrib ;
    cpl_propertylist    *   plist ;
    hdrl_image          *   master_flat_loc ;
    cpl_image           *   my_master_flat ;
    cpl_image           *   bpm_im ;
    cpl_image           *   bpm_flat ;
    cpl_mask            *   bpm_mask ;
    cpl_table           *   computed_traces ;
    cpl_table           *   filtered_traces ;
    cpl_table           *   traces ;
    cpl_bivector        **  spectrum ;
    cpl_vector          **  slit_func_vec ;
    hdrl_image          *   model_master;
    cpl_table           *   slit_func_tab ;
    cpl_table           *   extract_tab ;
    hdrl_image          *   model_tmp ;
    cpl_vector          *   dits ;
    cpl_vector          *   ndits ;
    int                 *   qc_order_nb ;
    double              *   qc_order_pos ;
    int                 *   orders ;
    double                  qc_mean, qc_median, qc_flux, qc_rms, qc_s2n, 
                            qc_trace_centery, dit,  
                            gain, error_factor, bl_tot ;
    int                     i, ext_nr, nb_traces,
                            nb_orders, qc_nbbad, nbvals, ngood, bl_good ;

    /* TODO, make parameters */
    int extract_niter = 10;
    double extract_kappa = 10;

    /* Check Inputs */
    if (rawframes == NULL) return -1 ;
    if (extr_method != CR2RES_EXTR_OPT_CURV && 
            extr_method != CR2RES_EXTR_SUM) {
        cpl_msg_error(__func__, "Failed to read the dits") ;
        return -1 ;
    }

    /* Get the Gain */
    if (reduce_det == 1) gain = CR2RES_GAIN_CHIP1 ;
    else if (reduce_det == 2) gain = CR2RES_GAIN_CHIP2 ;
    else if (reduce_det == 3) gain = CR2RES_GAIN_CHIP3 ;
    else {
        cpl_msg_error(__func__, "Failed to get the Gain value") ;
        return -1 ;
    }

    /* Get the First RAW file  */
    first_file = cpl_frame_get_filename(
            cpl_frameset_get_position_const(rawframes, 0)) ;

    /* Get the DIT for the Dark correction */
    if ((dits = cr2res_io_read_dits(rawframes)) == NULL) {
        cpl_msg_error(__func__, "Failed to read the dits") ;
        return -1 ;
    }
    /*Load the NDITs */
    ndits = cr2res_io_read_ndits(rawframes);

    /* Load the image list */
    imlist = cr2res_io_load_image_list_from_set(rawframes, reduce_det) ;
    if (imlist == NULL) {
        cpl_msg_error(__func__, "Failed to load the images") ;
        cpl_vector_delete(dits) ;
        cpl_vector_delete(ndits) ;
        return -1 ;
    }

    /* Set the error factor. */
    error_factor = gain * cpl_vector_get(ndits, 0) *
                                cpl_frameset_get_size(rawframes) ;

    for (i=0; i<cpl_vector_get_size(ndits); i++){
        if (cpl_vector_get(ndits,i) != cpl_vector_get(ndits, 0))
            cpl_msg_warning(__func__, "Raw frames have different NDIT! "
                "Error spectrum will likely be scaled incorrectly.");
    }

    /* Compute traces */
    cpl_msg_info(__func__, "Compute the traces") ;
    cpl_msg_indent_more() ;
    if ((computed_traces = cr2res_trace(
                    hdrl_image_get_image(hdrl_imagelist_get(imlist, 0)),
                    trace_smooth_x, trace_smooth_y, trace_threshold, 
                    trace_opening, trace_degree, trace_min_cluster)) == NULL) {
        cpl_msg_error(__func__, "Failed compute the traces") ;
        cpl_vector_delete(dits) ;
        cpl_vector_delete(ndits) ;
        hdrl_imagelist_delete(imlist) ;
        cpl_msg_indent_less() ;
        return -1 ;
    }
    cpl_msg_indent_less() ;

    /* Calibrate the Data */
    cpl_msg_info(__func__, "Calibrate the input images") ;
    cpl_msg_indent_more() ;
    if ((imlist_calibrated = cr2res_calib_imagelist(imlist, reduce_det, 
            0, subtract_nolight_rows, 0, cosmics, NULL, 
            master_dark_frame, bpm_frame, detlin_frame, dits, ndits))==NULL) {
        cpl_msg_error(__func__, "Failed to Calibrate the Data") ;
        cpl_vector_delete(dits) ;
        cpl_vector_delete(ndits) ;
        hdrl_imagelist_delete(imlist) ;
        cpl_table_delete(computed_traces) ;
        cpl_msg_indent_less() ;
        return -1 ;
    } else {
        /* Replace the calibrated image in the list */
        hdrl_imagelist_delete(imlist) ;
        imlist = imlist_calibrated ;
    }
    cpl_vector_delete(dits) ;
    cpl_vector_delete(ndits) ;
    cpl_msg_indent_less() ;

    /* Collapse */
    cpl_msg_info(__func__, "Collapse the input images") ;
    cpl_msg_indent_more() ;
    if (hdrl_imagelist_collapse_mean(imlist, &collapsed, &contrib) !=
            CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Failed to Collapse") ;
        hdrl_imagelist_delete(imlist) ;
        cpl_table_delete(computed_traces) ;
        cpl_msg_indent_less() ;
        return -1 ;
    }
    hdrl_imagelist_delete(imlist) ;
    cpl_image_delete(contrib) ;
    cpl_msg_indent_less() ;

    /* Add The remaining Columns to the trace table */
    cr2res_trace_add_extra_columns(computed_traces, first_file, reduce_det) ;

    /* Filter out traces */
    if (filter_traces) {
        char *setting_id;
        int zp_order;
        cpl_msg_info(__func__, "Filter out the traces") ;
        cpl_msg_indent_more() ;

        /* Get the setting and the zp_order */
        plist = cpl_propertylist_load(first_file, 0) ;
        setting_id = cpl_strdup(cr2res_pfits_get_wlen_id(plist)) ;
        zp_order = cr2res_pfits_get_order_zp(plist) ;
        cpl_propertylist_delete(plist) ;
        filtered_traces = cr2res_trace_filter(computed_traces,
                setting_id, zp_order) ;
        cpl_free(setting_id) ;
        cpl_table_delete(computed_traces) ;
        computed_traces = filtered_traces ;
        cpl_msg_indent_less() ;
        if (cpl_table_get_nrow(computed_traces) == 0) {
            cpl_msg_error(__func__, "All traces are filtered out") ;
            hdrl_image_delete(collapsed) ;
            cpl_table_delete(computed_traces) ;
            cpl_msg_indent_less() ;
            return -1 ;
        }
    }

    /* Use the input traces for extraction if they are provided */
    if (tw_frame) {
        traces = cr2res_io_load_TRACE_WAVE(cpl_frame_get_filename(tw_frame), 
                reduce_det) ;
        if (traces == NULL) {
            cpl_msg_error(__func__, "Cannot Load the provided TW file") ;
            hdrl_image_delete(collapsed) ;
            cpl_table_delete(computed_traces) ;
            return -1 ;
        }
    } else {
        traces = cpl_table_duplicate(computed_traces) ;
    }
    cpl_table_delete(computed_traces) ;

    /* Extract */
    nb_traces = cpl_table_get_nrow(traces) ;
    spectrum = cpl_malloc(nb_traces * sizeof(cpl_bivector *)) ;
    slit_func_vec = cpl_malloc(nb_traces * sizeof(cpl_vector *)) ;
    model_master = hdrl_image_new(CR2RES_DETECTOR_SIZE, CR2RES_DETECTOR_SIZE) ;
    hdrl_image_mul_scalar(model_master, (hdrl_value){0.0, 0.0}) ;

    /* Loop over the traces and extract them */
    cpl_msg_info(__func__, "Extract the traces") ;
    cpl_msg_indent_more() ;
    for (i=0 ; i<nb_traces ; i++) {
        int order, trace_id;
        /* Initialise */
        slit_func_vec[i] = NULL ;
        spectrum[i] = NULL ;
        model_tmp = NULL ;

        /* Get Order and trace id */
        order = cpl_table_get(traces, CR2RES_COL_ORDER, i, NULL) ;
        trace_id = cpl_table_get(traces, CR2RES_COL_TRACENB, i, NULL) ;

        /* Check if this order needs to be skipped */
        if (reduce_order > -1 && order != reduce_order) continue ;

        /* Check if this trace needs to be skipped */
        if (reduce_trace > -1 && trace_id != reduce_trace) continue ;

        cpl_msg_info(__func__, "Process Order %d/Trace %d", order, trace_id) ;
        cpl_msg_indent_more() ;

        /* Call the Extraction */
        if (extr_method == CR2RES_EXTR_SUM) {
            if (cr2res_extract_sum_vert(collapsed, traces, order, trace_id, 
                        extract_height, &(slit_func_vec[i]), &(spectrum[i]), 
                        &model_tmp) != 0) {
                cpl_msg_error(__func__, "Cannot (sum-)extract the trace") ;
                slit_func_vec[i] = NULL ;
                spectrum[i] = NULL ;
                model_tmp = NULL ;
                cpl_error_reset() ;
                cpl_msg_indent_less() ;
                continue ;
            }
        } else if (extr_method == CR2RES_EXTR_MEDIAN) {
            if (cr2res_extract_median(collapsed, traces, order, trace_id, 
                        extract_height, &(slit_func_vec[i]), &(spectrum[i]), 
                        &model_tmp) != 0) {
                cpl_msg_error(__func__, "Cannot (median-)extract the trace") ;
                slit_func_vec[i] = NULL ;
                spectrum[i] = NULL ;
                model_tmp = NULL ;
                cpl_error_reset() ;
                cpl_msg_indent_less() ;
                continue ;
            }
        } else if (extr_method == CR2RES_EXTR_TILTSUM) {
            if (cr2res_extract_sum_tilt(collapsed, traces, order, trace_id, 
                        extract_height, &(slit_func_vec[i]), &(spectrum[i]), 
                        &model_tmp) != 0) {
                cpl_msg_error(__func__, "Cannot (tiltsum-)extract the trace") ;
                slit_func_vec[i] = NULL ;
                spectrum[i] = NULL ;
                model_tmp = NULL ;
                cpl_error_reset() ;
                cpl_msg_indent_less() ;
                continue ;
            }
        } else if (extr_method == CR2RES_EXTR_OPT_CURV) {
            if (cr2res_extract_slitdec_curved(collapsed, traces, NULL,order, 
                        trace_id, extract_height, extract_swath_width, 
                        extract_oversample, extract_smooth_slit, 
                        extract_smooth_spec, extract_niter, extract_kappa,
                        error_factor,
                        &(slit_func_vec[i]), &(spectrum[i]), &model_tmp) != 0) {
                cpl_msg_error(__func__, "Cannot (slitdec-) extract the trace") ;
                slit_func_vec[i] = NULL ;
                spectrum[i] = NULL ;
                model_tmp = NULL ;
                cpl_error_reset() ;
                cpl_msg_indent_less() ;
                continue ;
            }
        }

        /* Update the model global image */
        if (model_tmp != NULL) {
            hdrl_image_add_image(model_master, model_tmp) ;
            hdrl_image_delete(model_tmp) ;
        }
        cpl_msg_indent_less() ;
    }
    cpl_msg_indent_less() ;

    /* Create the slit_func_tab for the current detector */
    slit_func_tab = cr2res_extract_SLITFUNC_create(slit_func_vec, traces) ;

    /* Create the extracted_tab for the current detector */
    extract_tab = cr2res_extract_EXTRACT1D_create(spectrum, traces) ;

    /* Deallocate Vectors */
    for (i=0 ; i<nb_traces ; i++) {
        if (slit_func_vec[i] != NULL) cpl_vector_delete(slit_func_vec[i]) ;
        if (spectrum[i] != NULL) cpl_bivector_delete(spectrum[i]) ;
    }
    cpl_free(spectrum) ;
    cpl_free(slit_func_vec) ;

    /* Compute the Master flat */
    cpl_msg_info(__func__, "Compute the master flat") ;
    cpl_msg_indent_more() ;
    if ((master_flat_loc = cr2res_master_flat(collapsed,
                    model_master, bpm_low, bpm_high, bpm_linemax,
                    &bpm_flat)) == NULL) {
        cpl_msg_error(__func__, "Failed compute the Master Flat") ;
        cpl_table_delete(slit_func_tab) ;
        cpl_table_delete(traces) ;
        cpl_table_delete(extract_tab) ;
        hdrl_image_delete(model_master) ;
        hdrl_image_delete(collapsed) ;
        cpl_msg_indent_less() ;
        return -1 ;
    }
    cpl_msg_indent_less() ;
    hdrl_image_delete(collapsed) ;

    /* Create BPM image */
    bpm_im = NULL ;
    if (bpm_frame != NULL) {
        if ((bpm_im = cr2res_io_load_BPM(
                        cpl_frame_get_filename(bpm_frame),
                        reduce_det, 1)) == NULL) {
            cpl_msg_warning(__func__, "Failed to Load the Master BPM") ;
        }
    }
    if (bpm_im == NULL) {
        bpm_im = cpl_image_duplicate(bpm_flat) ;
    } else {
        if (cpl_image_or(bpm_im, NULL, bpm_flat)) {
            cpl_msg_error(__func__, "Failed to add the mask to the BPM") ;
            cpl_table_delete(traces) ;
            cpl_table_delete(slit_func_tab) ;
            cpl_table_delete(extract_tab) ;
            hdrl_image_delete(model_master) ;
            hdrl_image_delete(master_flat_loc) ;
            cpl_image_delete(bpm_im) ;
            cpl_image_delete(bpm_flat) ;
            cpl_msg_indent_less() ;
            return -1 ;
        }
    }
    cpl_image_delete(bpm_flat) ;

    /* Compute the QC parameters */

    /* Load the first RAW image */
    first_image = cr2res_io_load_image(first_file, reduce_det) ;

    /* Set the BPM */
    bpm_mask = cr2res_bpm_extract_mask(bpm_im, CR2RES_BPM_ALL) ;

    cpl_image_reject_from_mask(hdrl_image_get_image(first_image), bpm_mask) ;
    qc_mean = cpl_image_get_mean(hdrl_image_get_image(first_image)) ;
    qc_median = cpl_image_get_median(hdrl_image_get_image(first_image)) ;
    plist = cpl_propertylist_load(first_file, 0) ;
    dit  = cr2res_pfits_get_dit(plist) ;
    cpl_propertylist_delete(plist) ;
    qc_flux = qc_mean / dit ;

    my_master_flat = cpl_image_duplicate(hdrl_image_get_image(master_flat_loc));
    cpl_image_reject_from_mask(my_master_flat, bpm_mask) ;
    ngood = CR2RES_DETECTOR_SIZE* CR2RES_DETECTOR_SIZE - 
        cpl_mask_count(bpm_mask) ;
    cpl_mask_delete(bpm_mask);
    qc_rms = sqrt(cpl_image_get_sqflux(my_master_flat) / ngood) ;
    cpl_image_delete(my_master_flat) ;

    qc_s2n = cr2res_qc_flat_s2n(extract_tab) ;

    cr2res_util_blaze_stat(extract_tab, &bl_good, &bl_tot); 

    qc_trace_centery = cr2res_qc_flat_trace_center_y(traces) ;
    qc_nbbad = cr2res_bpm_count(bpm_im, CR2RES_BPM_FLAT) ;
    cr2res_qc_flat_order_positions(traces, &qc_order_nb, &qc_order_pos,&nbvals);

    /* Load the extension header for saving */
    ext_nr = cr2res_io_get_ext_idx(first_file, reduce_det, 1) ;
    plist = cpl_propertylist_load(first_file, ext_nr) ;
    if (plist == NULL) {
        hdrl_image_delete(first_image) ;
        cpl_table_delete(traces) ;
        cpl_table_delete(slit_func_tab) ;
        cpl_table_delete(extract_tab) ;
        hdrl_image_delete(model_master) ;
        hdrl_image_delete(master_flat_loc) ;
        cpl_image_delete(bpm_im) ;
        cpl_msg_error(__func__, "Failed to load the plist") ;
        return -1 ;
    }

    /* QC.OVEREXPOSED */

    /* Get all orders */
    orders = cr2res_trace_get_order_idx_values(traces, &nb_orders) ;
    /* Loop on the orders */
    for (i=0 ; i<nb_orders ; i++) {
        double qc_overexposed;
        qc_overexposed = cr2res_qc_overexposed(
                hdrl_image_get_image(first_image), traces, orders[i]) ;
        char * qc_name = cpl_sprintf("%s%02d",
                CR2RES_HEADER_QC_OVEREXPOSED, orders[i]) ;
        cpl_propertylist_append_double(plist, qc_name, qc_overexposed);
        cpl_free(qc_name) ;
    }
    cpl_free(orders) ;
    hdrl_image_delete(first_image) ;

    /* Store the QC parameters in the plist */
    if (qc_order_nb != NULL && qc_order_pos != NULL) {
        for (i=0 ; i<nbvals ; i++) {
            char * qc_name = cpl_sprintf("%s%02d",
                    CR2RES_HEADER_QC_FLAT_ORDERPOS, qc_order_nb[i]) ;
            cpl_propertylist_append_double(plist, qc_name, qc_order_pos[i]);
            cpl_free(qc_name) ;
        }
    }
    if (qc_order_nb != NULL) cpl_free(qc_order_nb) ;
    if (qc_order_pos != NULL) cpl_free(qc_order_pos) ;

    cpl_propertylist_append_double(plist, CR2RES_HEADER_QC_FLAT_MEAN, 
            qc_mean) ;
    cpl_propertylist_append_double(plist, CR2RES_HEADER_QC_FLAT_MEDIAN, 
            qc_median);
    cpl_propertylist_append_double(plist, CR2RES_HEADER_QC_FLAT_FLUX, 
            qc_flux) ;
    cpl_propertylist_append_double(plist, CR2RES_HEADER_QC_FLAT_RMS, 
            qc_rms) ;
    cpl_propertylist_append_double(plist, CR2RES_HEADER_QC_FLAT_S2N, 
            qc_s2n) ;
    cpl_propertylist_append_double(plist, CR2RES_HEADER_QC_FLAT_CENTERY,
            qc_trace_centery) ;
    cpl_propertylist_append_int(plist, CR2RES_HEADER_QC_FLAT_NBBAD,
            qc_nbbad) ;
    cpl_propertylist_append_int(plist, CR2RES_HEADER_QC_BLAZE_NGOOD, 
            bl_good) ;
    cpl_propertylist_append_double(plist, CR2RES_HEADER_QC_BLAZE_TOT,
            bl_tot) ;

    /* Return the results */
    *trace_wave = traces ;
    *slit_func = slit_func_tab ;
    *extract_1d = extract_tab ;
    *slit_model = model_master ;
    *master_flat = master_flat_loc ;
    *ext_plist = plist ;
    *bpm = bpm_im ;
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Comparison function to identify different settings
  @param    frame1  first frame 
  @param    frame2  second frame 
  @return   0 if frame1!=frame2, 1 if frame1==frame2, -1 in error case
 */
/*----------------------------------------------------------------------------*/
static int cr2res_cal_flat_compare(
        const cpl_frame   *   frame1,
        const cpl_frame   *   frame2)
{
    int                     comparison ;
    cpl_propertylist    *   plist1 ;
    cpl_propertylist    *   plist2 ;
    const char          *   sval1 ;
    const char          *   sval2 ;

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

    /* Compare the SETTING used */
    sval1 = cr2res_pfits_get_wlen_id(plist1) ;
    sval2 = cr2res_pfits_get_wlen_id(plist2) ;
    if (cpl_error_get_code()) {
        cpl_msg_error(__func__, "Cannot get the reference wavelength");
        cpl_propertylist_delete(plist1) ;
        cpl_propertylist_delete(plist2) ;
        return -1 ;
    }
    if (strcmp(sval1, sval2)) comparison = 0 ;

    cpl_propertylist_delete(plist1) ;
    cpl_propertylist_delete(plist2) ;
    return comparison ;
}


