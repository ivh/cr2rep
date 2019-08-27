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
#include "cr2res_detlin.h"
#include "cr2res_pfits.h"
#include "cr2res_trace.h"
#include "cr2res_qc.h"
#include "cr2res_dfs.h"
#include "cr2res_bpm.h"
#include "cr2res_io.h"

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

#define RECIPE_STRING "cr2res_cal_detlin"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/
static int cr2res_cal_detlin_update(
        cpl_image           *   global_bpm,
        cpl_image           *   new_bpm,
        hdrl_imagelist      *   global_coeffs,
        hdrl_imagelist      *   new_coeffs) ;
static int cr2res_cal_detlin_compare(
        const cpl_frame   *   frame1,
        const cpl_frame   *   frame2) ;
static int cr2res_cal_detlin_reduce(
        const cpl_frameset  *   rawframes,
        double                  bpm_kappa,
        int                     trace_degree,
        int                     trace_min_cluster,
        int                     trace_smooth_x,
        int                     trace_smooth_y,
        double                  trace_threshold,
        int                     trace_opening,
        int                     trace_collapse,
        int                     reduce_det,
        int                     plotx,
        int                     ploty,
        hdrl_imagelist      **  coeffs,
        cpl_image           **  bpm,
        cpl_propertylist    **  ext_plist) ;
static int cr2res_cal_detlin_create(cpl_plugin *);
static int cr2res_cal_detlin_exec(cpl_plugin *);
static int cr2res_cal_detlin_destroy(cpl_plugin *);
static int cr2res_cal_detlin(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_cal_detlin_description[] = "\
Detector Linearity                                                      \n\
  Measure the pixels non-linearity                                      \n\
                                                                        \n\
  Inputs                                                                \n\
    raw.fits " CR2RES_DETLIN_RAW " [3 to n]                             \n\
                                                                        \n\
  Outputs                                                               \n\
    cr2res_cal_detlin_coeffs.fits " CR2RES_CAL_DETLIN_COEFFS_PROCATG "  \n\
    cr2res_cal_detlin_bpm.fits " CR2RES_CAL_DETLIN_BPM_PROCATG "        \n\
                                                                        \n\
  Algorithm                                                             \n\
    group the input raw frames by different settings                    \n\
    loop on groups g:                                                   \n\
      loop on detectors d:                                              \n\
        cr2res_cal_detlin_reduce() computes bpm(g, d) and coeffs(g, d)  \n\
        fill global_bpm(d) with bpm(g, d)                               \n\
        fill global_coeffs(d) with coeffs(g, d)                         \n\
                                                                        \n\
      if (--single_settings)                                            \n\
        save bpm(g) file                                                \n\
        save coeffs(g) file                                             \n\
    save global_bpm file                                                \n\
    save global_coeffs file                                             \n\
                                                                        \n\
    cr2res_cal_detlin_reduce()                                          \n\
      load input imlist and dits                                        \n\
      compute the traces (from 1. image, or collapsed if --trace_collapse)\n\
        use cr2res_trace(--trace_smooth, --trace_degree,                \n\
                         --trace_min_cluster, --trace_opening)          \n\
      loop on the detector pixels pix:                                  \n\
        if the pixel is within a trace:                                 \n\
          cr2res_detlin_compute() computes polynomial(pix) and errors(pix)\n\
      use the coeffs for the bpm computation                            \n\
      store the qc parameters in the returned property list             \n\
                                                                        \n\
  Library Functions used                                                \n\
    cr2res_trace()                                                      \n\
    cr2res_detlin_compute()                                             \n\
    cr2res_qc_detlin_gain()                                             \n\
    cr2res_qc_detlin_median()                                           \n\
    cr2res_qc_detlin_min_max_level()                                    \n\
    cr2res_io_save_BPM()                                                \n\
    cr2res_io_save_DETLIN_COEFFS()                                      \n\
" ;

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
                    "Detector Linearity recipe",
                    cr2res_cal_detlin_description,
                    "Thomas Marquart, Yves Jung",
                    PACKAGE_BUGREPORT,
                    cr2res_get_license(),
                    cr2res_cal_detlin_create,
                    cr2res_cal_detlin_exec,
                    cr2res_cal_detlin_destroy)) {    
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
static int cr2res_cal_detlin_create(cpl_plugin * plugin)
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
    p = cpl_parameter_new_value("cr2res.cr2res_cal_detlin.bpm_kappa",
            CPL_TYPE_DOUBLE, "Kappa threshold for BPM detection",
            "cr2res.cr2res_cal_detlin", 1.8);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "bpm_kappa");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_detlin.trace_degree",
            CPL_TYPE_INT, "polynomial degree for the fit to the orders",
            "cr2res.cr2res_cal_detlin", 5);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "trace_degree");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_detlin.trace_min_cluster",
            CPL_TYPE_INT, "size in pixels of the smallest allowed cluster",
            "cr2res.cr2res_cal_detlin", 10000);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "trace_min_cluster");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_detlin.trace_smooth_x",
            CPL_TYPE_INT, "Length of the smoothing kernel in x",
            "cr2res.cr2res_cal_detlin", 11);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "trace_smooth_x");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_detlin.trace_smooth_y",
            CPL_TYPE_INT, "Length of the smoothing kernel in y",
            "cr2res.cr2res_cal_detlin", 201);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "trace_smooth_y");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_detlin.trace_threshold",
            CPL_TYPE_DOUBLE, "Detection Threshold",
            "cr2res.cr2res_cal_detlin", -200.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "trace_threshold");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_detlin.trace_opening",
            CPL_TYPE_BOOL, "Use a morphological opening to rejoin clusters",
            "cr2res.cr2res_cal_detlin", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "trace_opening");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_detlin.single_settings",
            CPL_TYPE_BOOL, "Create the products for each setting",
            "cr2res.cr2res_cal_detlin", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "single_settings");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_detlin.trace_collapse",
            CPL_TYPE_BOOL, "Collapse the input frames for the trace analysis",
            "cr2res.cr2res_cal_detlin", TRUE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "trace_collapse");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_detlin.detector",
            CPL_TYPE_INT, "Only reduce the specified detector",
            "cr2res.cr2res_cal_detlin", 0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "detector");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_detlin.plot_x",
            CPL_TYPE_INT, "X position for the plot",
            "cr2res.cr2res_cal_detlin", 0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "plot_x");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_detlin.plot_y",
            CPL_TYPE_INT, "Y position for the plot",
            "cr2res.cr2res_cal_detlin", 0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "plot_y");
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
static int cr2res_cal_detlin_exec(cpl_plugin * plugin)
{
    cpl_recipe  *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1;

    return cr2res_cal_detlin(recipe->frames, recipe->parameters);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_cal_detlin_destroy(cpl_plugin * plugin)
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
static int cr2res_cal_detlin(
        cpl_frameset            *   frameset,
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter *   param ;
    int                     trace_degree, trace_min_cluster, trace_collapse,
                            trace_opening, single_settings, reduce_det, 
                            trace_smooth_x, trace_smooth_y, plot_x, plot_y ;
    double                  bpm_kappa, trace_threshold ;
    cpl_frameset        *   rawframes ;
    cpl_size            *   labels ;
    cpl_size                nlabels ;
    cpl_frameset        *   raw_one ;
    char                *   setting_id ;
    hdrl_imagelist      *   coeffs_merged[CR2RES_NB_DETECTORS] ;
    cpl_image           *   bpm_merged[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   ext_plist_merged[CR2RES_NB_DETECTORS] ;
    hdrl_imagelist      *   coeffs[CR2RES_NB_DETECTORS] ;
    cpl_image           *   bpm[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   ext_plist[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   plist ;
    char                *   out_file;
    int                     i, l, det_nr; 

    /* Initialise */

    /* RETRIEVE INPUT PARAMETERS */
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_detlin.bpm_kappa");
    bpm_kappa = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_detlin.trace_degree");
    trace_degree = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_detlin.trace_min_cluster");
    trace_min_cluster = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_detlin.trace_smooth_x");
    trace_smooth_x = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_detlin.trace_smooth_y");
    trace_smooth_y = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_detlin.trace_threshold");
    trace_threshold = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_detlin.trace_opening");
    trace_opening = cpl_parameter_get_bool(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_detlin.single_settings");
    single_settings = cpl_parameter_get_bool(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_detlin.trace_collapse");
    trace_collapse = cpl_parameter_get_bool(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_detlin.detector");
    reduce_det = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_detlin.plot_x");
    plot_x = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_detlin.plot_y");
    plot_y = cpl_parameter_get_int(param);

    /* Identify the RAW and CALIB frames in the input frameset */
    if (cr2res_dfs_set_groups(frameset)) {
        cpl_msg_error(__func__, "Cannot identify RAW and CALIB frames") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Get Calibration frames */

    /* Retrieve raw frames */
    if ((rawframes = cr2res_extract_frameset(frameset, 
                    CR2RES_DETLIN_RAW)) == NULL) {
        cpl_msg_error(__func__, "No raw frame in input") ;
        return -1 ;
    }

    /* Labelise the raw frames with the different settings*/
    if ((labels = cpl_frameset_labelise(rawframes, cr2res_cal_detlin_compare,
                &nlabels)) == NULL) {
        cpl_msg_error(__func__, "Cannot labelise input frames") ;
        cpl_frameset_delete(rawframes) ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Initialise */
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
        coeffs_merged[det_nr-1] = NULL ;
        bpm_merged[det_nr-1] = NULL ;
        ext_plist_merged[det_nr-1] = NULL ;
    }

    /* Loop on the settings */
    for (l=0 ; l<(int)nlabels ; l++) {
        /* Get the frames for the current setting */
        raw_one = cpl_frameset_extract(rawframes, labels, (cpl_size)l) ;

        /* Get the current setting */
        plist = cpl_propertylist_load(cpl_frame_get_filename(
                    cpl_frameset_get_position(raw_one, 0)), 0) ;
        setting_id = cpl_strdup(cr2res_pfits_get_wlen_id(plist)) ;
        cr2res_format_setting(setting_id) ;
        cpl_propertylist_delete(plist) ;

        cpl_msg_info(__func__, "Process SETTING %s", setting_id) ;
        cpl_msg_indent_more() ;

        /* Loop on the detectors */
        for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
            cpl_msg_info(__func__, "Process Detector %d", det_nr) ;

            /* Initialise */
            coeffs[det_nr-1] = NULL ;
            bpm[det_nr-1] = NULL ;
            ext_plist[det_nr-1] = NULL ;

            /* Compute only one detector */
            if (reduce_det != 0 && det_nr != reduce_det) continue ;
    
            /* Call the reduction function */
            cpl_msg_indent_more() ;
            if (cr2res_cal_detlin_reduce(raw_one, bpm_kappa,
                        trace_degree, trace_min_cluster, trace_smooth_x,
                        trace_smooth_y, trace_threshold, trace_opening, 
                        trace_collapse, det_nr, plot_x, plot_y,
                        &(coeffs[det_nr-1]),
                        &(bpm[det_nr-1]),
                        &(ext_plist[det_nr-1])) == -1) {
                cpl_msg_warning(__func__, 
                        "Failed to reduce SETTING %s / det %d", 
                        setting_id, det_nr);
                cpl_error_reset() ;
            } 
            cpl_msg_indent_less() ;

            /* Merge the products */
            if (ext_plist[det_nr-1] != NULL && coeffs[det_nr-1] != NULL
                    && bpm[det_nr-1] != NULL) {
                /* Take the first header as it is */
                if (ext_plist_merged[det_nr-1] == NULL) {
                    ext_plist_merged[det_nr-1] =
                        cpl_propertylist_duplicate(ext_plist[det_nr-1]) ;
                }

                /* Handle the coefficients and bpm */
                if (coeffs_merged[det_nr-1] == NULL) {
                    /* First non-null solution encountered */
                    coeffs_merged[det_nr-1] =
                        hdrl_imagelist_duplicate(coeffs[det_nr-1]) ;
                    bpm_merged[det_nr-1] = cpl_image_duplicate(bpm[det_nr-1]) ;
                } else {
                    /* Merge */
                    cr2res_cal_detlin_update(bpm_merged[det_nr-1],
                            bpm[det_nr-1], coeffs_merged[det_nr-1],
                            coeffs[det_nr-1]) ;
                }
            }
        }

        /* Save the products */
        if (single_settings) {
            /* BPM */
            out_file = cpl_sprintf("%s_%s_bpm.fits", RECIPE_STRING,
                    setting_id) ;
            cr2res_io_save_BPM(out_file, frameset, raw_one, parlist, bpm, NULL, 
                    ext_plist, CR2RES_CAL_DETLIN_BPM_PROCATG, RECIPE_STRING) ;
            cpl_free(out_file);

            /* COEFFS */
            out_file = cpl_sprintf("%s_%s_coeffs.fits", RECIPE_STRING,
                    setting_id) ;
            cr2res_io_save_DETLIN_COEFFS(out_file, frameset, raw_one, parlist, 
                    coeffs, NULL, ext_plist, CR2RES_CAL_DETLIN_COEFFS_PROCATG, 
                    RECIPE_STRING) ;
            cpl_free(out_file);
        }

        /* Free */
        cpl_frameset_delete(raw_one) ;
        for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
            if (coeffs[det_nr-1] != NULL) 
                hdrl_imagelist_delete(coeffs[det_nr-1]) ;
            if (bpm[det_nr-1] != NULL) 
                cpl_image_delete(bpm[det_nr-1]) ;
            if (ext_plist[det_nr-1] != NULL) 
                cpl_propertylist_delete(ext_plist[det_nr-1]) ;
        }
        cpl_free(setting_id) ;
        cpl_msg_indent_less() ;
    }
    cpl_free(labels) ;

    /* Save the merged products */
    /* BPM */
    out_file = cpl_sprintf("%s_bpm.fits", RECIPE_STRING) ;
    cr2res_io_save_BPM(out_file, frameset, rawframes, parlist, bpm_merged, 
            NULL, ext_plist_merged, CR2RES_CAL_DETLIN_BPM_PROCATG, 
            RECIPE_STRING) ;
    cpl_free(out_file);

    /* COEFFS */
    out_file = cpl_sprintf("%s_coeffs.fits", RECIPE_STRING) ;
    cr2res_io_save_DETLIN_COEFFS(out_file, frameset, rawframes, parlist, 
            coeffs_merged, NULL, ext_plist_merged, 
            CR2RES_CAL_DETLIN_COEFFS_PROCATG, RECIPE_STRING) ;
    cpl_free(out_file);

    cpl_frameset_delete(rawframes) ;
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
        if (coeffs_merged[det_nr-1] != NULL) 
            hdrl_imagelist_delete(coeffs_merged[det_nr-1]) ;
        if (bpm_merged[det_nr-1] != NULL) 
            cpl_image_delete(bpm_merged[det_nr-1]) ;
        if (ext_plist_merged[det_nr-1] != NULL) 
            cpl_propertylist_delete(ext_plist_merged[det_nr-1]) ;
    }

    return (int)cpl_error_get_code();
}

/*----------------------------------------------------------------------------*/
/**
  @brief  Compute the non-linearity for a single setting, single detector
  @param rawframes          Raw frames from a single setting
  @param bpm_kappa          Kappa value for BPM detection
  @param trace_degree
  @param trace_min_cluster
  @param trace_smooth_x
  @param trace_smooth_y
  @param trace_threshold
  @param trace_opening      
  @param trace_collapse     Flag to collapse (or not) before tracing
  @param reduce_det         The detector to compute 
  @param plotx              xposition to plot
  @param ploty              yposition to plot
  @param coeffs             [out] the non-linearity coeffs and error
  @param bpm                [out] the BPM
  @param ext_plist          [out] the header for saving the products
  @return  0 if ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
static int cr2res_cal_detlin_reduce(
        const cpl_frameset  *   rawframes,
        double                  bpm_kappa,
        int                     trace_degree,
        int                     trace_min_cluster,
        int                     trace_smooth_x,
        int                     trace_smooth_y,
        double                  trace_threshold,
        int                     trace_opening,
        int                     trace_collapse,
        int                     reduce_det,
        int                     plotx,
        int                     ploty,
        hdrl_imagelist      **  coeffs,
        cpl_image           **  bpm,
        cpl_propertylist    **  ext_plist)
{
    cpl_frameset        *   sorted_frames ;
    const char          *   first_file ;
    hdrl_imagelist      *   imlist ;
    hdrl_image          *   cur_im ;
    double              *   pcur_im ;
    cpl_vector          *   dits ;
    hdrl_image          *   collapsed ;
    cpl_image           *   contrib ;
    hdrl_image          *   master_flat_loc ;
    cpl_table           *   traces ;
    cpl_imagelist       *   coeffs_loc ;
    cpl_image           *   cur_coeffs ;
    double              *   pcur_coeffs ;
    cpl_imagelist       *   errors_loc ;
    cpl_image           *   cur_errors ;
    double              *   pcur_errors ;
    cpl_propertylist    *   plist ;
    cpl_image           *   bpm_loc ;
    int                 *   pbpm_loc ;
    cpl_image           *   trace_image ;
    int                 *   pti ;
    cpl_polynomial      *   fitted_poly ;
    cpl_vector          *   fitted_errors ;
    cpl_vector          *   fitvals ;
    cpl_mask            *   bpm_mask ;
    int                     i, j, k, idx, ext_nr_data, order, trace_id, nx, ny;
    cpl_size                max_degree, l ;
    double                  low_thresh, high_thresh, median, sigma ;
    int                     qc_nb_bad, qc_nbfailed, qc_nbsuccess ;
    double                  qc_median, qc_fitquality, qc_gain,
                            qc_min_level, qc_max_level ;
    
    /* Check Inputs */
    if (rawframes == NULL) return -1 ;

    /* Sort the frames by increasing DIT */
    if ((sorted_frames = cr2res_detlin_sort_frames(rawframes)) == NULL) {
        cpl_msg_error(__func__, "Failed sorting frames by increasing DITs") ;
        return -1 ;
    }

    /* Initialise */
    max_degree = 2 ;

    /* Get the Extension number */
    first_file = cpl_frame_get_filename(
            cpl_frameset_get_position_const(sorted_frames, 0)) ;
    ext_nr_data = cr2res_io_get_ext_idx(first_file, reduce_det, 1) ;

    /* Load the extension header for saving */
    plist = cpl_propertylist_load(first_file, ext_nr_data) ;
    if (plist == NULL) {
        cpl_frameset_delete(sorted_frames) ;
        return -1 ;
    }

    /* Load the image list */
    if ((imlist = cr2res_io_load_image_list_from_set(sorted_frames,
                    reduce_det)) == NULL) {
        cpl_propertylist_delete(plist);
        cpl_frameset_delete(sorted_frames) ;
        cpl_msg_error(__func__, "Failed to Load the images") ;
        return -1 ;
    }

    /* Load the DITs */
    if ((dits = cr2res_io_read_dits(sorted_frames)) == NULL) {
        hdrl_imagelist_delete(imlist) ;
        cpl_propertylist_delete(plist);
        cpl_frameset_delete(sorted_frames) ;
        cpl_msg_error(__func__, "Failed to Load the DIT values") ;
        return -1 ;
    }
    cpl_frameset_delete(sorted_frames) ;

    /* Collapse all input images for the traces detection (only if wished) */
    if (trace_collapse) {
        cpl_msg_info(__func__, "Collapse the input images") ;
        cpl_msg_indent_more() ;
        if (hdrl_imagelist_collapse_mean(imlist, &collapsed, &contrib) !=
                CPL_ERROR_NONE) {
            cpl_msg_error(__func__, "Failed to Collapse") ;
            cpl_propertylist_delete(plist);
            hdrl_imagelist_delete(imlist) ;
            cpl_vector_delete(dits); 
            cpl_msg_indent_less() ;
            return -1 ;
        }
        cpl_msg_indent_less() ;
    } else {
        /* Only use the first image */
        collapsed = hdrl_image_duplicate(hdrl_imagelist_get(imlist, 0)) ;
    }

    /* Compute traces */
    cpl_msg_info(__func__, "Compute the traces") ;
    cpl_msg_indent_more() ;
    nx = hdrl_image_get_size_x(collapsed) ;
    ny = hdrl_image_get_size_y(collapsed) ;
    if ((traces = cr2res_trace(hdrl_image_get_image(collapsed), 
                    trace_smooth_x, trace_smooth_y, trace_threshold, 
                    trace_opening, trace_degree, trace_min_cluster)) == NULL) {
        cpl_msg_error(__func__, "Failed compute the traces") ;
        hdrl_imagelist_delete(imlist) ;
        cpl_vector_delete(dits); 
        cpl_propertylist_delete(plist);
        hdrl_image_delete(collapsed) ;
        cpl_image_delete(contrib);
        cpl_msg_indent_less() ;
        return -1 ;
    }
    hdrl_image_delete(collapsed) ;
    cpl_image_delete(contrib);
    cpl_msg_indent_less() ;

    /* Allocate */
    bpm_loc = cpl_image_new(nx, ny, CPL_TYPE_INT) ;
    pbpm_loc = cpl_image_get_data_int(bpm_loc) ;
    coeffs_loc = cpl_imagelist_new() ;
    errors_loc = cpl_imagelist_new() ;

    /* Initialise the coeffs cube */
    for (l=0 ; l<=max_degree ; l++) {
        cpl_imagelist_set(coeffs_loc, cpl_image_new(nx, ny,CPL_TYPE_DOUBLE),l) ;
        cpl_imagelist_set(errors_loc, cpl_image_new(nx, ny,CPL_TYPE_DOUBLE),l) ;
    }

    /* Initialise the BPM as Out of Order */
    cpl_image_add_scalar(bpm_loc, CR2RES_BPM_OUTOFORDER); 

    /* Create the trace image */
    trace_image = cr2res_trace_gen_image(traces, nx, ny) ;
    pti = cpl_image_get_data_int(trace_image) ;
    cpl_table_delete(traces) ;

    /* Loop over the traces and compute the non-linearity */
    cpl_msg_info(__func__, "Compute Non Linearity") ;
    cpl_msg_indent_more() ;

    /* Loop on the traces pixels */
    qc_nbfailed = 0 ;
    qc_nbsuccess = 0 ;
    for (j=0 ; j<ny ; j++) {
        for (i=0 ; i<nx ; i++) {
            idx = i + j*nx ;
            if (pti[idx] > 0) {
                /* Prepare values to fit */
                fitvals = cpl_vector_new(cpl_vector_get_size(dits)) ;
                for (k=0 ; k<cpl_vector_get_size(dits) ; k++) {
                    cur_im = hdrl_imagelist_get(imlist,  k) ;
                    pcur_im = cpl_image_get_data_double(
                            hdrl_image_get_image(cur_im)) ;
                    cpl_vector_set(fitvals, k, (float)(pcur_im[idx])) ;
                }


                /* We are in a trace, let's compute the linearity */
                if (cr2res_detlin_compute(dits, fitvals, max_degree,
                            &fitted_poly, &fitted_errors) == -1) {
                    qc_nbfailed++ ;
                    /* Store the null Coefficients in the output image list */
                    pbpm_loc[idx] = CR2RES_BPM_DETLIN ;
                    for (l=0 ; l<=max_degree ; l++) {
                        cur_coeffs = cpl_imagelist_get(coeffs_loc, l) ;
                        pcur_coeffs = cpl_image_get_data_double(cur_coeffs) ;
                        pcur_coeffs[idx] = 0.0 ;
                        cur_errors = cpl_imagelist_get(errors_loc, l) ;
                        pcur_errors = cpl_image_get_data_double(cur_errors) ;
                        pcur_errors[idx] = 0.0 ;
                    } 
                } else {
                    qc_nbsuccess++ ;

                    /* Plot the values and the fit */
                    if (plotx==i+1 && ploty==j+1) {
                        cpl_bivector * toplot_measure =
                            cpl_bivector_wrap_vectors(dits, fitvals) ;
                        cpl_vector * poly_eval = cr2res_polynomial_eval_vector(
                                fitted_poly, dits) ;
                        cpl_bivector * toplot_fitted =
                            cpl_bivector_wrap_vectors(dits, poly_eval) ;
                        cpl_plot_bivector(
                        "set grid;set xlabel 'dits (s)';set ylabel 'int';",
                        "t 'Measured Detlin' w lines", "", toplot_measure);
                        cpl_plot_bivector(
                        "set grid;set xlabel 'dits (s)';set ylabel 'int';",
                        "t 'Fitted Detlin' w lines", "", toplot_fitted);
                        cpl_bivector_unwrap_vectors(toplot_fitted) ;
                        cpl_vector_delete(poly_eval) ;
                        cpl_bivector_unwrap_vectors(toplot_measure) ;
                    }

                    /* Store the Coefficients in the output image list */
                    pbpm_loc[idx] = 0 ;
                    for (l=0 ; l<=max_degree ; l++) {
                        cur_coeffs = cpl_imagelist_get(coeffs_loc, l) ;
                        pcur_coeffs = cpl_image_get_data_double(cur_coeffs) ;
                        pcur_coeffs[idx] =
                            cpl_polynomial_get_coeff(fitted_poly, &l) ;
                        cur_errors = cpl_imagelist_get(errors_loc, l) ;
                        pcur_errors = cpl_image_get_data_double(cur_errors) ;

                        /* Store error */
                        if (fitted_errors != NULL) {
                            pcur_errors[idx] = cpl_vector_get(fitted_errors, l);
                        } else {
                            pcur_errors[idx] = 0.0 ;
                        }
                    }
                    if (fitted_errors != NULL) cpl_vector_delete(fitted_errors);
                    cpl_polynomial_delete(fitted_poly) ;
                }
                cpl_vector_delete(fitvals) ;
            }
        }
    }
    cpl_msg_indent_less() ;
    cpl_image_delete(trace_image) ;
    hdrl_imagelist_delete(imlist) ;
    cpl_vector_delete(dits); 

    /* Use the second coefficient stats for the BPM detection */
    cpl_msg_info(__func__, "BPM detection") ;
    cur_coeffs = cpl_imagelist_get(coeffs_loc, 1) ;
    pcur_coeffs = cpl_image_get_data_double(cur_coeffs) ;
    bpm_mask = cpl_mask_new(nx, ny) ;
    cpl_mask_threshold_image(bpm_mask, cur_coeffs,
            (double)CR2RES_BPM_OUTOFORDER-0.5, 
            (double)CR2RES_BPM_OUTOFORDER+0.5,
            CPL_BINARY_1) ;
    cpl_image_reject_from_mask(cur_coeffs, bpm_mask) ;
    cpl_mask_delete(bpm_mask); 
    median = cpl_image_get_median_dev(cur_coeffs, &sigma) ;
    cpl_image_accept_all(cur_coeffs) ;

    low_thresh = median - bpm_kappa * sigma ;
    high_thresh = median + bpm_kappa * sigma ;
    qc_nb_bad = 0;
    for (j=0 ; j<ny ; j++) {
        for (i=0 ; i<nx ; i++) {
            idx = i + j*nx ;
            if (pbpm_loc[idx] != CR2RES_BPM_OUTOFORDER && 
                    (pcur_coeffs[idx] < low_thresh || 
                     pcur_coeffs[idx] > high_thresh)) {
                pbpm_loc[idx] = CR2RES_BPM_DETLIN ;
                qc_nb_bad ++ ;
            }
        }
    }

    /* Compute the QC parameters */
    /* TODO  */
    qc_fitquality = 0.0;
    qc_median = cr2res_qc_detlin_median(coeffs_loc) ;
    qc_gain = cr2res_qc_detlin_gain(coeffs_loc) ;
    /* TODO */
    qc_min_level = qc_max_level = 0 ;
    cr2res_qc_detlin_min_max_level(NULL, &qc_min_level, &qc_max_level) ;

    /* Store the QC parameters in the plist */
    cpl_propertylist_append_int(plist, CR2RES_HEADER_QC_DETLIN_NBAD, 
            qc_nb_bad) ;
    cpl_propertylist_append_int(plist, CR2RES_HEADER_QC_DETLIN_NBFAILED, 
            qc_nbfailed) ;
    cpl_propertylist_append_int(plist, CR2RES_HEADER_QC_DETLIN_NBSUCCESS,
            qc_nbsuccess) ;
    cpl_propertylist_append_double(plist, CR2RES_HEADER_QC_DETLIN_FITQUALITY, 
            qc_fitquality) ;
    cpl_propertylist_append_double(plist, CR2RES_HEADER_QC_DETLIN_MEDIAN, 
            qc_median) ;
    cpl_propertylist_append_double(plist, CR2RES_HEADER_QC_DETLIN_GAIN, 
            qc_gain) ;
    cpl_propertylist_append_double(plist, CR2RES_HEADER_QC_DETLIN_MINLEVEL,
            qc_min_level) ;
    cpl_propertylist_append_double(plist, CR2RES_HEADER_QC_DETLIN_MAXLEVEL,
            qc_max_level) ;

    /* Return the results */
    *coeffs = hdrl_imagelist_create(coeffs_loc, errors_loc) ;
    cpl_imagelist_delete(coeffs_loc) ;
    cpl_imagelist_delete(errors_loc) ;
    *ext_plist = plist ;
    *bpm = bpm_loc ;
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
static int cr2res_cal_detlin_compare(
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

/*----------------------------------------------------------------------------*/
/**
  @brief Only pixels not yet computed (CR2RES_BPM_OUTOFORDER) are updated
  @param 
  @return   0 if ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_cal_detlin_update(
        cpl_image           *   global_bpm,
        cpl_image           *   new_bpm,
        hdrl_imagelist      *   global_coeffs,
        hdrl_imagelist      *   new_coeffs)
{
    cpl_size        i, j, k, nx, ny, ni, idx ;
    int         *   pglobal_bpm ;
    int         *   pnew_bpm ;
    hdrl_image  *   global_coeffs_ima ;
    hdrl_image  *   new_coeffs_ima ;
    
    /* Check Inputs */
    if (global_bpm == NULL || new_bpm == NULL || global_coeffs == NULL
            || new_coeffs == NULL) return 0 ;

    /* Initialise */
    nx = cpl_image_get_size_x(global_bpm) ;
    ny = cpl_image_get_size_y(global_bpm) ;
    ni = hdrl_imagelist_get_size(global_coeffs) ;

    if (cpl_image_get_size_x(new_bpm) != nx ||
            cpl_image_get_size_y(new_bpm) != ny ||
            hdrl_imagelist_get_size(new_coeffs) != ni) {
        return -1 ;
    }

    pglobal_bpm = cpl_image_get_data_int(global_bpm) ;
    pnew_bpm = cpl_image_get_data_int(new_bpm) ;

    /* Loop on the pixels */
    for (j=0 ; j<ny ; j++) {
        for (i=0 ; i<nx ; i++) {
            idx = i+j*nx ;
            if (pglobal_bpm[idx] == CR2RES_BPM_OUTOFORDER &&
                    pnew_bpm[idx] != CR2RES_BPM_OUTOFORDER) {
                /* Put the new value */
                pglobal_bpm[idx] = pnew_bpm[idx] ;
                for (k=0 ; k<ni ; k++) {
                    global_coeffs_ima = hdrl_imagelist_get(global_coeffs, k) ;
                    new_coeffs_ima = hdrl_imagelist_get(new_coeffs, k) ;

                    hdrl_image_set_pixel(global_coeffs_ima, i+1, j+1,
                            hdrl_image_get_pixel(new_coeffs_ima, i+1,
                                j+1, NULL)) ;
                }
            }
        }
    }
    return 0 ;
}
