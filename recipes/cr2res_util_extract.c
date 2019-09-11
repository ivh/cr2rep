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
#include "hdrl.h"

#include "cr2res_utils.h"
#include "cr2res_pfits.h"
#include "cr2res_dfs.h"
#include "cr2res_io.h"
#include "cr2res_trace.h"
#include "cr2res_extract.h"

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

#define RECIPE_STRING "cr2res_util_extract"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static cpl_frameset * cr2res_util_extract_find_RAW(const cpl_frameset * in) ;
static int cr2res_util_extract_create(cpl_plugin *);
static int cr2res_util_extract_exec(cpl_plugin *);
static int cr2res_util_extract_destroy(cpl_plugin *);
static int cr2res_util_extract(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_util_extract_description[] = "\
Spectrum Extraction                                                     \n\
  This utility performs the optimal extraction along precomputed traces \n\
                                                                        \n\
  Inputs                                                                \n\
    raw.fits " CR2RES_FLAT_RAW " [1 to n]                               \n\
          or " CR2RES_CALIBRATED_PROTYPE "                              \n\
          or " CR2RES_WAVE_RAW "                                        \n\
          or " CR2RES_OBS_NODDING_RAW "                                 \n\
          or " CR2RES_OBS_2D_RAW "                                      \n\
          or " CR2RES_OBS_POL_RAW "                                     \n\
    trace.fits " CR2RES_CAL_FLAT_TW_PROCATG " [1]                       \n\
            or " CR2RES_CAL_FLAT_TW_MERGED_PROCATG "                    \n\
            or " CR2RES_UTIL_TRACE_TW_PROCATG "                         \n\
            or " CR2RES_UTIL_WAVE_TW_PROCATG "                          \n\
            or " CR2RES_CAL_WAVE_TW_PROCATG "                           \n\
            or " CR2RES_UTIL_SLIT_CURV_TW_PROCATG "                     \n\
    bpm.fits " CR2RES_CAL_DARK_BPM_PROCATG " [0 to 1]                   \n\
          or " CR2RES_CAL_FLAT_BPM_PROCATG "                            \n\
          or " CR2RES_CAL_DETLIN_BPM_PROCATG "                          \n\
          or " CR2RES_UTIL_BPM_SPLIT_PROCATG "                          \n\
                                                                        \n\
  Outputs                                                               \n\
    <input_name>_extr1D.fits " CR2RES_UTIL_EXTRACT_1D_PROCATG "         \n\
    <input_name>_extrSlitFu.fits " CR2RES_UTIL_SLIT_FUNC_PROCATG "      \n\
    <input_name>_extrModel.fits " CR2RES_UTIL_SLIT_MODEL_PROCATG "      \n\
                                                                        \n\
  Algorithm                                                             \n\
    loop on raw frames f:                                               \n\
      loop on detectors d:                                              \n\
        Load the trace wave                                             \n\
        Recompute a new trace wave with the specified slit fraction     \n\
                 (--slit_frac) if needed                                \n\
        Load the image to extract                                       \n\
        Load the BPM and set them in the image                          \n\
        Run the extraction cr2res_extract_traces(--method,--height,     \n\
                 --swath_width,--oversample,--smooth_slit)              \n\
          -> creates SLIT_MODEL(f,d), SLIT_FUNC(f,d), EXTRACT_1D(f,d)   \n\
      Save SLIT_MODEL(f), SLIT_FUNC(f), EXTRACT_1D(f)                   \n\
                                                                        \n\
  Library functions uѕed                                                \n\
    cr2res_io_load_TRACE_WAVE()                                         \n\
    cr2res_trace_new_slit_fraction()                                    \n\
    cr2res_io_load_image()                                              \n\
    cr2res_io_load_BPM()                                                \n\
    cr2res_extract_traces()                                             \n\
    cr2res_io_save_SLIT_MODEL()                                         \n\
    cr2res_io_save_SLIT_FUNC()                                          \n\
    cr2res_io_save_EXTRACT_1D()                                         \n\
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
                    "Optimal Extraction utility",
                    cr2res_util_extract_description,
                    "Thomas Marquart, Yves Jung",
                    PACKAGE_BUGREPORT,
                    cr2res_get_license(),
                    cr2res_util_extract_create,
                    cr2res_util_extract_exec,
                    cr2res_util_extract_destroy)) {
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
static int cr2res_util_extract_create(cpl_plugin * plugin)
{
    cpl_recipe    * recipe;
    cpl_parameter * p;

    /* Check that the plugin is part of a valid recipe */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else
        return -1;

    /* Create the parameters list in the cpl_recipe object */
    recipe->parameters = cpl_parameterlist_new();

    /* Fill the parameters list */
    p = cpl_parameter_new_value("cr2res.cr2res_util_extract.oversample",
            CPL_TYPE_INT, "factor by which to oversample the extraction",
            "cr2res.cr2res_util_extract", 5);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "oversample");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_extract.swath_width",
            CPL_TYPE_INT, "The swath width", "cr2res.cr2res_util_extract", 90);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "swath_width");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_extract.height",
            CPL_TYPE_INT, "Extraction height",
            "cr2res.cr2res_util_extract", -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "height");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_extract.smooth_slit",
            CPL_TYPE_DOUBLE,
            "Smoothing along the slit (1 for high S/N, 5 for low)",
            "cr2res.cr2res_util_extract", 1.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "smooth_slit");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_extract.method",
            CPL_TYPE_STRING, "Extraction method (SUM / OPT_VERT / OPT_CURV )",
            "cr2res.cr2res_util_extract", "OPT_CURV");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "method");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_extract.slit_frac",
            CPL_TYPE_STRING, "Wished slit fraction",
            "cr2res.cr2res_util_extract", "-1.0, -1.0");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "slit_frac");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_extract.detector",
            CPL_TYPE_INT, "Only reduce the specified detector",
            "cr2res.cr2res_util_extract", 0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "detector");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_extract.order",
            CPL_TYPE_INT, "Only reduce the specified order",
            "cr2res.cr2res_util_extract", -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "order");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_extract.trace_nb",
            CPL_TYPE_INT, "Only reduce the specified trace number",
            "cr2res.cr2res_util_extract", -1);
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
static int cr2res_util_extract_exec(cpl_plugin * plugin)
{
    cpl_recipe  *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1;

    return cr2res_util_extract(recipe->frames, recipe->parameters);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_util_extract_destroy(cpl_plugin * plugin)
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
static int cr2res_util_extract(
        cpl_frameset            *   frameset,
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter *   param;
    int                     oversample, swath_width, extr_height,
                            reduce_det, reduce_order, reduce_trace ;
    double                  smooth_slit, slit_low, slit_up ;
    cpl_array           *   slit_frac ;
    cpl_frameset        *   rawframes ;
    const cpl_frame     *   cur_frame ;
    const char          *   cur_fname ;
    cpl_frameset        *   cur_fset ;
    const cpl_frame     *   trace_frame ;
    const cpl_frame     *   bpm_frame ;
    const char          *   sval ;
    char                *   out_file;
    hdrl_image          *   model_master[CR2RES_NB_DETECTORS] ;
    cpl_table           *   slit_func_tab[CR2RES_NB_DETECTORS] ;
    cpl_table           *   extract_tab[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   ext_plist[CR2RES_NB_DETECTORS] ;
    cpl_table           *   trace_table ;
    cpl_table           *   trace_table_new ;
    hdrl_image          *   science_hdrl;
    cpl_image           *   bpm_img;
    cpl_mask            *   bpm_mask;
    int                     det_nr, ext_nr, order, i ;
    cr2res_extr_method      extr_method;

    /* Needed for sscanf() */
    setlocale(LC_NUMERIC, "C");

    /* RETRIEVE INPUT PARAMETERS */
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_extract.oversample");
    oversample = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_extract.swath_width");
    swath_width = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_extract.smooth_slit");
    smooth_slit = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_extract.detector");
    reduce_det = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_extract.order");
    reduce_order = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_extract.trace_nb");
    reduce_trace = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_extract.height");
    extr_height = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_extract.method");
    sval = cpl_parameter_get_string(param);
    if (!strcmp(sval, ""))              extr_method = CR2RES_EXTR_OPT_CURV;
    else if (!strcmp(sval, "OPT_VERT")) extr_method = CR2RES_EXTR_OPT_VERT;
    else if (!strcmp(sval, "OPT_CURV")) extr_method = CR2RES_EXTR_OPT_CURV;
    else if (!strcmp(sval, "SUM"))      extr_method = CR2RES_EXTR_SUM;
    else {
        cpl_msg_error(__func__, "Invalid Extraction Method specified");
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1;
    }
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_extract.slit_frac");
    sval = cpl_parameter_get_string(param) ;
    if (sscanf(sval, "%lg,%lg", &slit_low, &slit_up) != 2) {
        cpl_msg_error(__func__, "Invalid Slit Fraction specified");
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Check Parameters */
    if (slit_low > 0.0 && slit_up > 0.0 && slit_low < 1.0 && slit_up < 1.0 
            && slit_up > slit_low) {
        slit_frac = cpl_array_new(3, CPL_TYPE_DOUBLE) ;
        cpl_array_set(slit_frac, 0, slit_low) ;
        cpl_array_set(slit_frac, 1, (slit_low+slit_up)/2.0) ;
        cpl_array_set(slit_frac, 2, slit_up) ;
    } else {
        slit_frac = NULL ;
    }

    /* Identify the RAW and CALIB frames in the input frameset */
    if (cr2res_dfs_set_groups(frameset) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot identify RAW and CALIB frames") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        if (slit_frac != NULL) cpl_array_delete(slit_frac) ;
        return -1 ;
    }

    /* Get Calibration frames */
    trace_frame = cr2res_io_find_TRACE_WAVE(frameset) ;
    bpm_frame = cr2res_io_find_BPM(frameset) ;

    /* Get the rawframes */
    rawframes = cr2res_util_extract_find_RAW(frameset) ;
    if (rawframes==NULL || cpl_frameset_get_size(rawframes) <= 0) {
        cpl_msg_error(__func__, "Cannot find any RAW file") ;
        cpl_error_set(__func__, CPL_ERROR_DATA_NOT_FOUND) ;
        if (rawframes!= NULL) cpl_frameset_delete(rawframes) ;
        return -1 ;
    }

    if (trace_frame == NULL) {
        cpl_msg_error(__func__, "The utility needs a trace wave frame");
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        if (slit_frac != NULL) cpl_array_delete(slit_frac) ;
        return -1 ;
    }
   
    /* Loop on the RAW frames */
    for (i=0 ; i<cpl_frameset_get_size(rawframes) ; i++) {
        /* Get the Current Frame */
        cur_frame = cpl_frameset_get_position(rawframes, i) ;
        cur_fname = cpl_frame_get_filename(cur_frame) ;
        cpl_msg_info(__func__, "Reduce Frame %s", cur_fname) ;
        cpl_msg_indent_more() ;

        /* Loop over the detectors */
        for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {

            /* Initialise */
            model_master[det_nr-1] = NULL ;
            slit_func_tab[det_nr-1] = NULL ;
            extract_tab[det_nr-1] = NULL ;
            ext_plist[det_nr-1] = NULL ;

            /* Store the extenѕion header for product saving */
            ext_nr = cr2res_io_get_ext_idx(cur_fname, det_nr, 1) ;
            if (ext_nr < 0) continue ;
            ext_plist[det_nr-1] = cpl_propertylist_load(cur_fname, ext_nr) ;

            /* Compute only one detector */
            if (reduce_det != 0 && det_nr != reduce_det) continue ;

            cpl_msg_info(__func__, "Process detector number %d", det_nr) ;
            cpl_msg_indent_more() ;

            /* Load the trace table of this detector */
            cpl_msg_info(__func__, "Load the trace table") ;
            if ((trace_table = cr2res_io_load_TRACE_WAVE(
                            cpl_frame_get_filename(trace_frame),
                            det_nr)) == NULL) {
                cpl_msg_error(__func__,
                        "Failed to get trace table - skip detector");
                cpl_error_reset() ;
                cpl_msg_indent_less() ;
                continue ;
            }

            /* Extract at the specified slit fraction */
            if (slit_frac != NULL) {
                if ((trace_table_new = cr2res_trace_new_slit_fraction(
                                trace_table, slit_frac)) == NULL) {
                    cpl_msg_warning(__func__,
            "Failed to compute the traces for user specified slit fraction") ;
                    cpl_error_reset() ;
                } else {
                    cpl_table_delete(trace_table) ;
                    trace_table = trace_table_new ;
                    trace_table_new = NULL ;
                }
            }

            /* Load the image in which the traces are to extract */
            cpl_msg_info(__func__, "Load the Image") ;
            if ((science_hdrl = cr2res_io_load_image(cur_fname, det_nr))==NULL){
                cpl_table_delete(trace_table) ;
                cpl_msg_error(__func__, 
                        "Failed to load the image - skip detector");
                cpl_error_reset() ;
                cpl_msg_indent_less() ;
                continue ;
            }

            /* Load the BPM and assign to hdrl-mask*/
            if (bpm_frame != NULL) {
                cpl_msg_info(__func__, "Load and assign the BPM") ;
                if ((bpm_img = cr2res_io_load_BPM(
                                cpl_frame_get_filename(bpm_frame), 
                                det_nr, 1))==NULL) {
                    cpl_table_delete(trace_table) ;
                    hdrl_image_delete(science_hdrl) ;
                    cpl_msg_error(__func__, 
                            "Failed to load BPM - skip detector");
                    cpl_error_reset() ;
                    cpl_msg_indent_less() ;
                    continue ;
                } else {
                    bpm_mask = cpl_mask_threshold_image_create(bpm_img, 0, 
                            INT_MAX);
                    cpl_mask_or(bpm_mask, hdrl_image_get_mask_const(science_hdrl));
                    if (hdrl_image_reject_from_mask(science_hdrl, 
                                bpm_mask) != CPL_ERROR_NONE) {
                        cpl_msg_error(__func__, 
                            "Failed to assign BPM to image - skip detector");
                        cpl_table_delete(trace_table) ;
                        hdrl_image_delete(science_hdrl) ;
                        cpl_mask_delete(bpm_mask);
                        cpl_image_delete(bpm_img);
                        cpl_error_reset() ;
                        cpl_msg_indent_less() ;
                        continue ;
                    }
                }
                cpl_mask_delete(bpm_mask);
                cpl_image_delete(bpm_img);
            }
            /* Compute the extraction */
            cpl_msg_info(__func__, "Spectra Extraction") ;
            if (cr2res_extract_traces(science_hdrl, trace_table, reduce_order,
                    reduce_trace, extr_method, extr_height, swath_width,
                    oversample, smooth_slit, &(extract_tab[det_nr-1]),
                    &(slit_func_tab[det_nr-1]), &(model_master[det_nr-1]))==-1){
                cpl_table_delete(trace_table) ;
                cpl_msg_error(__func__, "Failed to extract - skip detector");
                cpl_error_reset() ;
                cpl_msg_indent_less() ;
                continue ;
            }
            hdrl_image_delete(science_hdrl) ;
            cpl_table_delete(trace_table) ;
            cpl_msg_indent_less() ;
        }
        cpl_array_delete(slit_frac) ;

        /* Generate the currently used frameset */
        /* TODO : add calibrations */
        cur_fset = cpl_frameset_new() ;
        cpl_frameset_insert(cur_fset, cpl_frame_duplicate(cur_frame)) ;

        /* Save the Products */
        out_file = cpl_sprintf("%s_extrModel.fits",
                cr2res_get_base_name(cr2res_get_root_name(cur_fname)));
        cr2res_io_save_SLIT_MODEL(out_file, frameset, cur_fset, parlist, 
                model_master, NULL, ext_plist, CR2RES_UTIL_SLIT_MODEL_PROCATG, 
                RECIPE_STRING) ;
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_extrSlitFu.fits",
                cr2res_get_base_name(cr2res_get_root_name(cur_fname)));
        cr2res_io_save_SLIT_FUNC(out_file, frameset, cur_fset, parlist, 
                slit_func_tab, NULL, ext_plist, CR2RES_UTIL_SLIT_FUNC_PROCATG, 
                RECIPE_STRING) ;
        cpl_free(out_file);
        out_file = cpl_sprintf("%s_extr1D.fits",
                cr2res_get_base_name(cr2res_get_root_name(cur_fname)));
        cr2res_io_save_EXTRACT_1D(out_file, frameset, cur_fset, parlist, 
                extract_tab, NULL, ext_plist, CR2RES_UTIL_EXTRACT_1D_PROCATG, 
                RECIPE_STRING) ;
        cpl_free(out_file);
        cpl_frameset_delete(cur_fset) ;

        /* Free and return */
        for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
            if (ext_plist[det_nr-1] != NULL)
                cpl_propertylist_delete(ext_plist[det_nr-1]) ;
            if (slit_func_tab[det_nr-1] != NULL)
                cpl_table_delete(slit_func_tab[det_nr-1]) ;
            if (extract_tab[det_nr-1] != NULL)
                cpl_table_delete(extract_tab[det_nr-1]) ;
            if (model_master[det_nr-1] != NULL) 
                hdrl_image_delete(model_master[det_nr-1]) ;
        }
        cpl_msg_indent_less() ;
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
                        CR2RES_CALIBRATED_PROTYPE
                        CR2RES_WAVE_RAW
                        CR2RES_OBS_NODDING_RAW
                        CR2RES_OBS_2D_RAW
                        CR2RES_OBS_POL_RAW
 */
/*----------------------------------------------------------------------------*/
static cpl_frameset * cr2res_util_extract_find_RAW(const cpl_frameset * in)
{
    cpl_frameset    *   out ;

    /* Check entries */
    if (in == NULL) return NULL ;

    out = cr2res_extract_frameset(in, CR2RES_FLAT_RAW) ;
    if (out == NULL)
        out = cr2res_extract_frameset(in, CR2RES_CALIBRATED_PROTYPE) ;
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

