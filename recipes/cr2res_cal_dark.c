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

#include <math.h>
#include <string.h>
#include <cpl.h>
#include "hdrl.h"

#include "cr2res_utils.h"
#include "cr2res_pfits.h"
#include "cr2res_dfs.h"
#include "cr2res_io.h"
#include "cr2res_bpm.h"
#include "cr2res_qc.h"
#include "cr2res_calib.h"

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

#define RECIPE_STRING "cr2res_cal_dark"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_cal_dark_compare(
        const cpl_frame   *   frame1,
        const cpl_frame   *   frame2) ;

static int cr2res_cal_dark_create(cpl_plugin *);
static int cr2res_cal_dark_exec(cpl_plugin *);
static int cr2res_cal_dark_destroy(cpl_plugin *);
static int cr2res_cal_dark(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_cal_dark_description[] = "\
Dark                                                                    \n\
  Compute the master dark                                               \n\
                                                                        \n\
  Inputs                                                                \n\
    raw.fits " CR2RES_DARK_RAW " [3 to n]                               \n\
                                                                        \n\
  Outputs                                                               \n\
    cr2res_cal_dark_[setting]_DITxNDIT_master.fits " 
    CR2RES_CAL_DARK_MASTER_PROCATG "\n\
    cr2res_cal_dark_[setting]_DITxNDIT_bpm.fits " 
    CR2RES_CAL_DARK_BPM_PROCATG "\n\
                                                                        \n\
  Algorithm                                                             \n\
    group the input frames by different values of DET SEQ1 DIT          \n\
               or/and DET NDIT or/and WLEN ID                           \n\
    loop on groups g:                                                   \n\
      loop on detectors d:                                              \n\
        Load the images and create the associate error for each of      \n\
               them using cr2res_detector_shotnoise_model(--gain)       \n\
        Collapse the images with hdrl_imagelist_collapse(--collapse.*)  \n\
        Compute BPM from the collapsed master dark using                \n\
               cr2res_bpm_compute(--bpm_kappa, --bpm_lines_ratio)       \n\
        Set the BPM in the master dark                                  \n\
        Compute the QCs with statistics and                             \n\
               cr2res_dark_qc_ron(--ron_hsize, --ron_nsamples)          \n\
      save master dark(g)                                               \n\
      save bpm(g)                                                       \n\
                                                                        \n\
  Library Functions used                                                \n\
    cr2res_detector_shotnoise_model()                                   \n\
    cr2res_bpm_compute()                                                \n\
    cr2res_bpm_from_mask()                                              \n\
    cr2res_dark_qc_ron()                                                \n\
    cr2res_bpm_count()                                                  \n\
    cr2res_io_save_MASTER_DARK()                                        \n\
    cr2res_io_save_BPM()                                                \n\
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
                    "Dark recipe",
                    cr2res_cal_dark_description,
                    CR2RES_PIPELINE_AUTHORS,
                    PACKAGE_BUGREPORT,
                    cr2res_get_license(),
                    cr2res_cal_dark_create,
                    cr2res_cal_dark_exec,
                    cr2res_cal_dark_destroy)) {    
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
static int cr2res_cal_dark_create(cpl_plugin * plugin)
{
    cpl_recipe          *   recipe ;
    cpl_parameter       *   p ;
    cpl_parameterlist   *   collapse_par ;
    hdrl_parameter      *   sigclip_def ;
    hdrl_parameter      *   minmax_def ;
    hdrl_parameter      *   mode_def ;

    /* Check that the plugin is part of a valid recipe */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else
        return -1;

    /* Create the parameters list in the cpl_recipe object */
    recipe->parameters = cpl_parameterlist_new();

    /* Fill the parameters list */

    /* --detector */
    p = cpl_parameter_new_value("cr2res.cr2res_cal_dark.detector",
            CPL_TYPE_INT, "Only reduce the specified detector",
            "cr2res.cr2res_cal_dark", 0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "detector");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /* --bpm_method */
    p = cpl_parameter_new_value("cr2res.cr2res_cal_dark.bpm_method",
            CPL_TYPE_STRING, "Method (DEFAULT, GLOBAL, LOCAL or RUNNING)",
            "cr2res.cr2res_cal_dark", "DEFAULT");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "bpm_method");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /* --bpm_kappa */
    p = cpl_parameter_new_value("cr2res_cal_dark.bpm_kappa", CPL_TYPE_DOUBLE,
           "Kappa Threshold for the BPM", "cr2res_cal_dark", -1.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "bpm_kappa");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /* --bpm_lines_ratio */
    p = cpl_parameter_new_value("cr2res_cal_dark.bpm_lines_ratio", 
            CPL_TYPE_DOUBLE, "Maximum ratio of bad pixels per line", 
            "cr2res_cal_dark", 0.5);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "bpm_lines_ratio");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /* --ron_hsize */
    p = cpl_parameter_new_value("cr2res_cal_dark.ron_hsize", 
            CPL_TYPE_INT, "Half size of the window for RON computation", 
            "cr2res_cal_dark", 6);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ron_hsize");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /* --ron_nsamples */
    p = cpl_parameter_new_value("cr2res_cal_dark.ron_nsamples", 
            CPL_TYPE_INT, "Number of samples for RON computation", 
            "cr2res_cal_dark", 100);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "ron_nsamples");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /* --gain */
    p = cpl_parameter_new_value("cr2res_cal_dark.gain", CPL_TYPE_DOUBLE,
       "Gain in [e- / ADU]", "cr2res_cal_dark", 2.1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "gain");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    /* Collapsing related parameters */
    sigclip_def = hdrl_collapse_sigclip_parameter_create(3., 3., 5);
    minmax_def = hdrl_collapse_minmax_parameter_create(1., 1.);
    mode_def = hdrl_collapse_mode_parameter_create(10., 1., 0., 
            HDRL_MODE_MEDIAN, 0);
    collapse_par = hdrl_collapse_parameter_create_parlist("cr2res_cal_dark", 
            "collapse", "MEAN", sigclip_def, minmax_def, mode_def) ;
    hdrl_parameter_delete(sigclip_def);
    hdrl_parameter_delete(minmax_def);
    hdrl_parameter_delete(mode_def);
    for (p = cpl_parameterlist_get_first(collapse_par) ;
            p != NULL; p = cpl_parameterlist_get_next(collapse_par))
        cpl_parameterlist_append(recipe->parameters,cpl_parameter_duplicate(p));
    cpl_parameterlist_delete(collapse_par);

    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Execute the plugin instance given by the interface
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_cal_dark_exec(cpl_plugin * plugin)
{
    cpl_recipe  *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1;

    return cr2res_cal_dark(recipe->frames, recipe->parameters);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_cal_dark_destroy(cpl_plugin * plugin)
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
static int cr2res_cal_dark(
        cpl_frameset            *   frameset,
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter *   par ;
    int                     reduce_det, ron_hsize, ron_nsamples, ndit ;
    double                  gain, dit, bpm_kappa, bpm_lines_ratio ;
    cr2res_bpm_method       bpm_method, my_bpm_method ;
    const char          *   sval ;
    hdrl_parameter      *   collapse_params ;
    cpl_frameset        *   rawframes ;
    cpl_frameset        *   raw_one ;
    cpl_frameset        *   raw_one_calib ;
    cpl_size            *   labels ;
    cpl_size                nlabels ;
    cpl_propertylist    *   plist ;
    cpl_propertylist    *   qc_main ;
    char                *   setting_id ;
    double                  bpm_kappa_global_default, bpm_kappa_local_default,
                            bpm_kappa_running_default, my_bpm_kappa ;
        
    hdrl_image          *   master_darks[CR2RES_NB_DETECTORS] ;
    cpl_image           *   bpms[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   ext_plist[CR2RES_NB_DETECTORS] ;

    hdrl_imagelist      *   dark_cube ;
    cpl_mask            *   my_bpm ;
    const char          *   fname ;
    hdrl_image          *   ima_data ;
    hdrl_image          *   ima_data_err ;
    cpl_image           *   ima_err ;

    double                  ron1, ron2, med, mean, sigma ;
    cpl_image           *   contrib_map;
    cpl_mask            *   bpm ;
    char                *   filename ;
    cpl_frame           *   frame ;
    int                     i, l, det_nr, nb_bad ;
    int                     single_dit_ndit_setting ;
    int                     original_ndit ;
    double                  original_dit ;
    char                *   original_setting ;

    /* Initialise */
    bpm_kappa_global_default = 0.5 ;
    bpm_kappa_local_default = 6.0 ;
    bpm_kappa_running_default = 5.0 ;

    /* RETRIEVE INPUT PARAMETERS */
    /* --detector */
    par = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_dark.detector");
    reduce_det = cpl_parameter_get_int(par);

    /* --bpm_method */
    par = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_dark.bpm_method");
    sval = cpl_parameter_get_string(par);
    if (!strcmp(sval, "GLOBAL"))        bpm_method = CR2RES_BPM_GLOBAL_STATS ;
    else if (!strcmp(sval, "LOCAL"))    bpm_method = CR2RES_BPM_LOCAL_STATS ;
    else if (!strcmp(sval, "RUNNING"))  bpm_method = CR2RES_BPM_RUNNING_FILTER ;
    else if (!strcmp(sval, "DEFAULT"))  bpm_method = CR2RES_BPM_UNSPECIFIED ;
    else {
        cpl_msg_error(__func__, "Unsupported bpm method") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* --ron_hsize */
    par = cpl_parameterlist_find_const(parlist, "cr2res_cal_dark.ron_hsize");
    ron_hsize = cpl_parameter_get_int(par);

    /* --ron_nsamples */
    par = cpl_parameterlist_find_const(parlist, "cr2res_cal_dark.ron_nsamples");
    ron_nsamples = cpl_parameter_get_int(par);

    /* --bpm_kappa */
    par = cpl_parameterlist_find_const(parlist, "cr2res_cal_dark.bpm_kappa");
    bpm_kappa = cpl_parameter_get_double(par);
	/* Set default depending on the method */
    if (bpm_kappa < 0.0) {
        if (bpm_method == CR2RES_BPM_GLOBAL_STATS) 
            bpm_kappa = bpm_kappa_global_default ;
        if (bpm_method == CR2RES_BPM_LOCAL_STATS) 
            bpm_kappa = bpm_kappa_local_default ;
        if (bpm_method == CR2RES_BPM_RUNNING_FILTER) 
            bpm_kappa = bpm_kappa_running_default ;
    }

    /* --bpm_lines_ratio */
    par = cpl_parameterlist_find_const(parlist,
            "cr2res_cal_dark.bpm_lines_ratio");
    bpm_lines_ratio = cpl_parameter_get_double(par);

    /* --gain */
    par = cpl_parameterlist_find_const(parlist, "cr2res_cal_dark.gain");
    gain = cpl_parameter_get_double(par);

    /* Collapse parameters */
    collapse_params = hdrl_collapse_parameter_parse_parlist(parlist,
            "cr2res_cal_dark.collapse") ;
   
    /* Identify the RAW and CALIB frames in the input frameset */
    if (cr2res_dfs_set_groups(frameset)) {
        hdrl_parameter_destroy(collapse_params) ;
        cpl_msg_error(__func__, "Cannot identify RAW and CALIB frames") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }
    
    /* Retrieve calibration data */

    /* Extract RAW frames */
    rawframes = cr2res_extract_frameset(frameset, CR2RES_DARK_RAW) ;
    if (rawframes==NULL || cpl_frameset_get_size(rawframes) <= 0) {
        hdrl_parameter_destroy(collapse_params) ;
        cpl_msg_error(__func__, "Cannot find any RAW file") ;
        cpl_error_set(__func__, CPL_ERROR_DATA_NOT_FOUND) ;
        return -1 ;
    }

    /* Labelise the raw frames with the different settings */
    if ((labels = cpl_frameset_labelise(rawframes, cr2res_cal_dark_compare,
                &nlabels)) == NULL) {
        cpl_msg_error(__func__, "Cannot labelise input frames") ;
        cpl_frameset_delete(rawframes) ;
        hdrl_parameter_destroy(collapse_params) ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Identify if there are several DIT/NDIT/setting */
    single_dit_ndit_setting = -1 ;
    original_setting = NULL ;
    for (l=0 ; l<(int)nlabels ; l++) {
        /* Get the frames for the current setting */
        raw_one = cpl_frameset_extract(rawframes, labels, (cpl_size)l) ;

        /* Get the current setting */
        plist = cpl_propertylist_load(cpl_frame_get_filename(
                    cpl_frameset_get_position(raw_one, 0)), 0) ;
        dit = cr2res_pfits_get_dit(plist) ;
        ndit = cr2res_pfits_get_ndit(plist) ;
        setting_id = cpl_strdup(cr2res_pfits_get_wlen_id(plist)) ;
        cr2res_format_setting(setting_id) ;
        cpl_propertylist_delete(plist) ;

        /* Update single_dit_ndit_setting information */
        if (single_dit_ndit_setting < 0) {
            single_dit_ndit_setting = 1; 
            original_dit = dit ;
            original_ndit = ndit ;
            original_setting = cpl_strdup(setting_id) ;
        } else {
            if (fabs(original_dit-dit)>1e-3 || original_ndit != ndit ||
                    strcmp(original_setting, setting_id)) {
                single_dit_ndit_setting = 0 ;
            }
        }
        cpl_free(setting_id);
        cpl_frameset_delete(raw_one) ;
    }
    cpl_free(original_setting);

    /* Loop on the settings */
    for (l = 0; l < (int)nlabels; l++) {
        cpl_vector *qc_ron1, *qc_ron2, *qc_nb_bad, *qc_mean, *qc_med, *qc_sigma;
        const char *first_fname;
        int nb_frames;
        /* Get the frames for the current setting */
        raw_one = cpl_frameset_extract(rawframes, labels, (cpl_size)l) ;
        nb_frames = cpl_frameset_get_size(raw_one) ;
        first_fname = 
            cpl_frame_get_filename(cpl_frameset_get_position(raw_one, 0)) ;

        /* Get the current setting */
        plist = cpl_propertylist_load(first_fname, 0) ;
        dit = cr2res_pfits_get_dit(plist) ;
        ndit = cr2res_pfits_get_ndit(plist) ;
        setting_id = cpl_strdup(cr2res_pfits_get_wlen_id(plist)) ;
        cr2res_format_setting(setting_id) ;
        cpl_propertylist_delete(plist) ;

        cpl_msg_info(__func__, "Process SETTING %s / DIT %g / %d",
                setting_id, dit, ndit) ;
        cpl_msg_indent_more() ;

        /* Allocate qc values vectors */
        qc_ron1 = cpl_vector_new(CR2RES_NB_DETECTORS); 
        qc_ron2 = cpl_vector_new(CR2RES_NB_DETECTORS); 
        qc_nb_bad = cpl_vector_new(CR2RES_NB_DETECTORS); 
        qc_mean = cpl_vector_new(CR2RES_NB_DETECTORS); 
        qc_med = cpl_vector_new(CR2RES_NB_DETECTORS); 
        qc_sigma = cpl_vector_new(CR2RES_NB_DETECTORS); 

        /* Loop on the detectors */
        for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
            /* Initialise */
            master_darks[det_nr-1] = NULL ;
            bpms[det_nr-1] = NULL ;
            cpl_vector_set(qc_ron1, det_nr-1, -1.0) ;
            cpl_vector_set(qc_ron2, det_nr-1, -1.0) ;
            cpl_vector_set(qc_nb_bad, det_nr-1, -1.0) ;
            cpl_vector_set(qc_mean, det_nr-1, -1.0) ;
            cpl_vector_set(qc_med, det_nr-1, -1.0) ;
            cpl_vector_set(qc_sigma, det_nr-1, -1.0) ;

            /* Store the extension header for product saving */
            ext_plist[det_nr-1] = cpl_propertylist_load(first_fname,
                    cr2res_io_get_ext_idx(first_fname, det_nr, 1)) ;

            /* Compute only one detector */
            if (reduce_det != 0 && det_nr != reduce_det) continue ;

            cpl_msg_info(__func__, "Process Detector nb %i", det_nr) ;
            cpl_msg_indent_more() ;

            /* Loop on the frames */
            dark_cube = hdrl_imagelist_new();
            for (i=0; i<nb_frames ; i++) {
                /* Identify current file */
                fname=cpl_frame_get_filename(
                        cpl_frameset_get_position(raw_one, i)) ; 
                cpl_msg_info(__func__, "Load Image from File %s / Detector %i", 
                        cr2res_get_base_name(fname), det_nr) ;

                /* Load the image */
                if ((ima_data = cr2res_io_load_image(fname, det_nr)) == NULL) {
                    cpl_msg_error(__func__, 
                            "Cannot load image from File %s / Detector %d", 
                            fname, det_nr) ;
                    cpl_error_set(__func__, CPL_ERROR_DATA_NOT_FOUND) ;
                    cpl_msg_indent_less() ;
                    cpl_msg_indent_less() ;
                    cpl_frameset_delete(rawframes) ;
                    hdrl_parameter_destroy(collapse_params) ;
                    cpl_free(labels);
                    cpl_free(setting_id) ;
                    cpl_frameset_delete(raw_one) ;
                    hdrl_imagelist_delete(dark_cube) ;
                    cpl_vector_delete(qc_ron1) ;
                    cpl_vector_delete(qc_ron2) ;
                    cpl_vector_delete(qc_nb_bad) ;
                    cpl_vector_delete(qc_mean) ;
                    cpl_vector_delete(qc_med) ;
                    cpl_vector_delete(qc_sigma) ;
                    return -1 ;
                }

                /* Create the noise image */
                cpl_msg_info(__func__, "Create the associated Noise image");
                if (cr2res_detector_shotnoise_model(
                            hdrl_image_get_image(ima_data), gain, 0.0,
                            &ima_err) != CPL_ERROR_NONE) {
                    cpl_free(labels);
                    cpl_msg_error(__func__, "Cannot create the Noise image") ;
                    cpl_error_set(__func__, CPL_ERROR_DATA_NOT_FOUND) ;
                    cpl_msg_indent_less() ;
                    cpl_msg_indent_less() ;
                    cpl_frameset_delete(rawframes) ;
                    hdrl_parameter_destroy(collapse_params) ;
                    cpl_free(setting_id) ;
                    cpl_free(labels);
                    cpl_frameset_delete(raw_one) ;
                    hdrl_imagelist_delete(dark_cube) ;
                    hdrl_image_delete(ima_data); 
                    cpl_vector_delete(qc_ron1) ;
                    cpl_vector_delete(qc_ron2) ;
                    cpl_vector_delete(qc_nb_bad) ;
                    cpl_vector_delete(qc_mean) ;
                    cpl_vector_delete(qc_med) ;
                    cpl_vector_delete(qc_sigma) ;
                    return -1 ;
                }

                /* Set the new error image */
                ima_data_err =
                    hdrl_image_create(hdrl_image_get_image(ima_data), ima_err);
                cpl_image_delete(ima_err) ;
                hdrl_image_delete(ima_data) ;
                
                /* Store the hdrl image in the dark_cube */
                hdrl_imagelist_set(dark_cube, ima_data_err, i);
            }

            /* Get the proper collapsing function and do frames combination */
            if (hdrl_imagelist_collapse(dark_cube, collapse_params,
                    &(master_darks[det_nr-1]), &contrib_map) != CPL_ERROR_NONE){
                cpl_msg_warning(__func__, "Cannot collapse Detector %d",det_nr);
                master_darks[det_nr-1] = NULL ;
                contrib_map = NULL ;
            }
            cpl_image_delete(contrib_map);
       
            /* Compute BPM from the MASTER dark */
            if (master_darks[det_nr-1] != NULL) {

                /* Handle default method */
                if (bpm_method == CR2RES_BPM_UNSPECIFIED) {
                    if (cr2res_is_short_wavelength(setting_id)) {
                        my_bpm_method = CR2RES_BPM_GLOBAL_STATS ;
                        if (bpm_kappa < 0.0) 
                            my_bpm_kappa = bpm_kappa_global_default ;
                        else 
                            my_bpm_kappa = bpm_kappa ;
                    } else {
                        my_bpm_method = CR2RES_BPM_LOCAL_STATS ;
                        if (bpm_kappa < 0.0) 
                            my_bpm_kappa = bpm_kappa_local_default ;
                        else 
                            my_bpm_kappa = bpm_kappa ;
                    }
                } else {
                    my_bpm_method = bpm_method ;
                    my_bpm_kappa = bpm_kappa ;
                }
                if ((my_bpm = cr2res_bpm_compute(
                                hdrl_image_get_image(master_darks[det_nr-1]),
                                my_bpm_method, my_bpm_kappa, bpm_lines_ratio, 
                                0))==NULL) {
                    cpl_msg_warning(__func__, "Cannot create BPM") ;
                } else {
                    /* Convert mask to BPM */
                    bpms[det_nr-1] = cr2res_bpm_from_mask(my_bpm, 
                            CR2RES_BPM_DARK);
                    cpl_mask_delete(my_bpm) ;
                }
            }
                        
            /* Set the BPM in the master dark and the RAW */
            if (bpms[det_nr-1] != NULL) {
                /* Get Mask */
                bpm = cpl_mask_threshold_image_create(bpms[det_nr-1],-0.5,0.5) ;
                cpl_mask_not(bpm) ;

                /* In dark_cube */
                for (i=0; i<hdrl_imagelist_get_size(dark_cube) ; i++) {
                    hdrl_image_reject_from_mask(hdrl_imagelist_get(dark_cube, 
                                i), bpm) ;
                }

                /* In Master Dark */
                hdrl_image_reject_from_mask(master_darks[det_nr-1], bpm) ;

                cpl_mask_delete(bpm) ;
            }

            /* QCs */

            /* QCs from RAW */
            if (hdrl_imagelist_get_size(dark_cube) >= 3) {
                ron1 = cr2res_dark_qc_ron(
                        hdrl_image_get_image(hdrl_imagelist_get(dark_cube,0)), 
                        hdrl_image_get_image(hdrl_imagelist_get(dark_cube,1)), 
                        ron_hsize, ron_nsamples, ndit) ;
                ron2 = cr2res_dark_qc_ron(
                        hdrl_image_get_image(hdrl_imagelist_get(dark_cube,1)), 
                        hdrl_image_get_image(hdrl_imagelist_get(dark_cube,2)), 
                        ron_hsize, ron_nsamples, ndit) ;
                if (cpl_error_get_code()) { 
                    cpl_error_reset() ;
                } else {
                    cpl_propertylist_append_double(ext_plist[det_nr-1], 
                            CR2RES_HEADER_QC_DARK_RON1, ron1) ;
                    cpl_propertylist_append_double(ext_plist[det_nr-1], 
                            CR2RES_HEADER_QC_DARK_RON2, ron2) ;
                    cpl_vector_set(qc_ron1, det_nr-1, ron1) ;
                    cpl_vector_set(qc_ron2, det_nr-1, ron2) ;
                }

            }
            hdrl_imagelist_delete(dark_cube);

            /* QCs from MASTER DARK */
            if (master_darks[det_nr-1] != NULL) {
                 /* Compute Thresholds */
                med = cpl_image_get_median_dev(
                        hdrl_image_get_image(master_darks[det_nr-1]), &sigma) ;
                mean = cpl_image_get_mean(
                        hdrl_image_get_image(master_darks[det_nr-1])) ;

                cpl_propertylist_append_double(ext_plist[det_nr-1], 
                        CR2RES_HEADER_QC_DARK_MEAN, mean) ;
                cpl_propertylist_append_double(ext_plist[det_nr-1], 
                        CR2RES_HEADER_QC_DARK_MEDIAN, med) ;
                cpl_propertylist_append_double(ext_plist[det_nr-1], 
                        CR2RES_HEADER_QC_DARK_STDEV, sigma) ;
                cpl_vector_set(qc_mean, det_nr-1, mean) ;
                cpl_vector_set(qc_med, det_nr-1, med) ;
                cpl_vector_set(qc_sigma, det_nr-1, sigma) ;
            }
            /* QCs from BPM */
            if (bpms[det_nr-1] != NULL) {
                nb_bad = cr2res_bpm_count(bpms[det_nr-1], CR2RES_BPM_DARK) ;
                cpl_propertylist_append_int(ext_plist[det_nr-1], 
                        CR2RES_HEADER_QC_DARK_NBAD, nb_bad) ;
                cpl_vector_set(qc_nb_bad, det_nr-1, (double)nb_bad) ;
            }
            cpl_msg_indent_less() ;
        }

        /* Compute Global QCs (primary header) */
        qc_main = NULL ;
        if (reduce_det == 0) {
            qc_main = cpl_propertylist_new() ;
            cpl_propertylist_append_double(qc_main, 
                    CR2RES_HEADER_QC_DARK_RON1_AVG,
                    cpl_vector_get_mean(qc_ron1)) ;
            cpl_propertylist_append_double(qc_main, 
                    CR2RES_HEADER_QC_DARK_RON1_RMS,
                    cpl_vector_get_stdev(qc_ron1)) ;

            cpl_propertylist_append_double(qc_main, 
                    CR2RES_HEADER_QC_DARK_RON2_AVG,
                    cpl_vector_get_mean(qc_ron2)) ;
            cpl_propertylist_append_double(qc_main, 
                    CR2RES_HEADER_QC_DARK_RON2_RMS,
                    cpl_vector_get_stdev(qc_ron2)) ;

            cpl_propertylist_append_double(qc_main, 
                    CR2RES_HEADER_QC_DARK_NBAD_AVG,
                    cpl_vector_get_mean(qc_nb_bad)) ;
            cpl_propertylist_append_double(qc_main, 
                    CR2RES_HEADER_QC_DARK_NBAD_RMS,
                    cpl_vector_get_stdev(qc_nb_bad)) ;

            cpl_propertylist_append_double(qc_main, 
                    CR2RES_HEADER_QC_DARK_MEAN_AVG,
                    cpl_vector_get_mean(qc_mean)) ;
            cpl_propertylist_append_double(qc_main, 
                    CR2RES_HEADER_QC_DARK_MEAN_RMS,
                    cpl_vector_get_stdev(qc_mean)) ;

            cpl_propertylist_append_double(qc_main, 
                    CR2RES_HEADER_QC_DARK_MEDIAN_AVG,
                    cpl_vector_get_mean(qc_med)) ;
            cpl_propertylist_append_double(qc_main, 
                    CR2RES_HEADER_QC_DARK_MEDIAN_RMS,
                    cpl_vector_get_stdev(qc_med)) ;

            cpl_propertylist_append_double(qc_main, 
                    CR2RES_HEADER_QC_DARK_STDEV_AVG,
                    cpl_vector_get_mean(qc_sigma)) ;
            cpl_propertylist_append_double(qc_main, 
                    CR2RES_HEADER_QC_DARK_STDEV_RMS,
                    cpl_vector_get_stdev(qc_sigma)) ;

        }
        cpl_vector_delete(qc_ron1) ;
        cpl_vector_delete(qc_ron2) ;
        cpl_vector_delete(qc_nb_bad) ;
        cpl_vector_delete(qc_mean) ;
        cpl_vector_delete(qc_med) ;
        cpl_vector_delete(qc_sigma) ;
        
        /* Save the results */

        /* Add the Calibrations (if any) to the frameset */
        raw_one_calib = cpl_frameset_duplicate(raw_one) ;
        cpl_frameset_delete(raw_one) ;
        for (i=0 ; i<cpl_frameset_get_size(frameset) ; i++) {
            frame = cpl_frameset_get_position(frameset, i);
            if (cpl_frame_get_group(frame) == CPL_FRAME_GROUP_CALIB) {
                cpl_frameset_insert(raw_one_calib,
                        cpl_frame_duplicate(frame)) ;
            }
        }

        /* MASTER DARK */
        if (single_dit_ndit_setting) {
            filename = cpl_sprintf("%s_master.fits", RECIPE_STRING);
        } else {
            filename = cpl_sprintf("%s_%s_%gx%d_master.fits", 
                    RECIPE_STRING, setting_id, dit, ndit); 
        }

        if (cr2res_io_save_MASTER_DARK(filename, frameset, raw_one_calib, 
                    parlist, master_darks, qc_main, ext_plist, 
                    CR2RES_CAL_DARK_MASTER_PROCATG, RECIPE_STRING) != 0) {
            if (qc_main != NULL) cpl_propertylist_delete(qc_main) ;
            cpl_frameset_delete(rawframes) ;
            cpl_frameset_delete(raw_one_calib) ;
            for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
                if (bpms[det_nr-1] != NULL) 
                    cpl_image_delete(bpms[det_nr-1]);
                if (master_darks[det_nr-1] != NULL) 
                    hdrl_image_delete(master_darks[det_nr-1]);
                if (ext_plist[det_nr-1] != NULL) 
                    cpl_propertylist_delete(ext_plist[det_nr-1]);
            }
            cpl_free(labels);
            cpl_free(setting_id);
            hdrl_parameter_destroy(collapse_params) ;
            cpl_free(filename) ;
            cpl_msg_error(__func__, "Cannot save the MASTER DARK") ;
            cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
            cpl_msg_indent_less() ;
            return -1 ;

        }
        if (qc_main != NULL) cpl_propertylist_delete(qc_main) ;
        cpl_free(filename) ;

        /* BPM */
        if (single_dit_ndit_setting) {
            filename = cpl_sprintf("%s_bpm.fits", RECIPE_STRING);
        } else {
            filename = cpl_sprintf("%s_%s_%gx%d_bpm.fits", 
                        RECIPE_STRING, setting_id, dit, ndit); 
        }

        if (cr2res_io_save_BPM(filename, frameset, raw_one_calib,
                    parlist, bpms, NULL, ext_plist, 
                    CR2RES_CAL_DARK_BPM_PROCATG, RECIPE_STRING) != 0) {
            cpl_frameset_delete(rawframes) ;
            cpl_frameset_delete(raw_one_calib) ;
            for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
                if (bpms[det_nr-1] != NULL) 
                    cpl_image_delete(bpms[det_nr-1]);
                if (master_darks[det_nr-1] != NULL) 
                    hdrl_image_delete(master_darks[det_nr-1]);
                if (ext_plist[det_nr-1] != NULL) 
                    cpl_propertylist_delete(ext_plist[det_nr-1]);
            }
            cpl_free(labels);
            cpl_free(setting_id);
            hdrl_parameter_destroy(collapse_params) ;
            cpl_free(filename) ;
            cpl_msg_error(__func__, "Cannot save the BPM") ;
            cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
            cpl_msg_indent_less() ;
            return -1 ;

        }
        cpl_free(filename) ;

        /* Free */
        cpl_free(setting_id);
        cpl_frameset_delete(raw_one_calib) ;
        for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
            if (bpms[det_nr-1] != NULL) 
                cpl_image_delete(bpms[det_nr-1]);
            if (master_darks[det_nr-1] != NULL) 
                hdrl_image_delete(master_darks[det_nr-1]);
            if (ext_plist[det_nr-1] != NULL) 
                cpl_propertylist_delete(ext_plist[det_nr-1]);
        }
        cpl_msg_indent_less() ;
    }
    cpl_free(labels);
    hdrl_parameter_delete(collapse_params);
    cpl_frameset_delete(rawframes) ;

    return (int)cpl_error_get_code();
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Comparison function to identify different settings
  @param    frame1  first frame
  @param    frame2  second frame
  @return   0 if frame1!=frame2, 1 if frame1==frame2, -1 in error case
 */
/*----------------------------------------------------------------------------*/
static int cr2res_cal_dark_compare(
        const cpl_frame   *   frame1,
        const cpl_frame   *   frame2)
{
    int                     comparison ;
    cpl_propertylist    *   plist1 ;
    cpl_propertylist    *   plist2 ;
    double                  dval1, dval2 ;
    int                     ival1, ival2 ;
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

    /* Compare the DIT used */
    dval1 = cr2res_pfits_get_dit(plist1) ;
    dval2 = cr2res_pfits_get_dit(plist2) ;
    if (cpl_error_get_code()) {
        cpl_msg_error(__func__, "Cannot get the DIT");
        cpl_propertylist_delete(plist1) ;
        cpl_propertylist_delete(plist2) ;
        return -1 ;
    }
    if (fabs(dval1-dval2) > 1e-3) comparison = 0 ;

    /* Compare the NDIT used */
    ival1 = cr2res_pfits_get_ndit(plist1) ;
    ival2 = cr2res_pfits_get_ndit(plist2) ;
    if (cpl_error_get_code()) {
        cpl_msg_error(__func__, "Cannot get the NDIT");
        cpl_propertylist_delete(plist1) ;
        cpl_propertylist_delete(plist2) ;
        return -1 ;
    }
    if (ival1 != ival2) comparison = 0 ;

    cpl_propertylist_delete(plist1) ;
    cpl_propertylist_delete(plist2) ;
    return comparison ;
}



