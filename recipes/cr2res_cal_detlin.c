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
        cpl_image           *   current_bpm,
        cpl_image           *   new_bpm,
        cpl_imagelist       *   current_coeffs,
        cpl_imagelist       *   new_coeffs) ;
static int cr2res_cal_detlin_compare(
        const cpl_frame   *   frame1,
        const cpl_frame   *   frame2) ;
static int cr2res_cal_detlin_reduce(
        const cpl_frameset  *   rawframes,
        double                  bpm_kappa,
        int                     trace_degree,
        int                     trace_min_cluster,
        double                  trace_smooth,
        int                     trace_opening,
        int                     trace_collapse,
        int                     reduce_det,
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

static char cr2res_cal_detlin_description[] =
"CRIRES+ detector linearity recipe\n"
"The files listed in the Set Of Frames (sof-file) must be tagged:\n"
"raw-file.fits " CR2RES_DETLIN_RAW "\n"
" The recipe produces the following products:\n"
"cr2res_cal_detlin_coeffs.fits " CR2RES_DETLIN_COEFFS_PROCATG "\n"
"cr2res_cal_detlin_bpm.fits " CR2RES_DETLIN_BPM_PROCATG "\n"
"\n";

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
                    "cr2res_cal_detlin",
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
            "cr2res.cr2res_cal_detlin", 3.0);
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

    p = cpl_parameter_new_value("cr2res.cr2res_cal_detlin.trace_smooth",
            CPL_TYPE_DOUBLE, "Length of the smoothing kernel",
            "cr2res.cr2res_cal_detlin", 5.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "trace_smooth");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_cal_detlin.trace_opening",
            CPL_TYPE_BOOL, "Use a morphological opening to rejoin clusters",
            "cr2res.cr2res_cal_detlin", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "trace_opening");
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
                            trace_opening, reduce_det ;
    double                  bpm_kappa, trace_smooth ;
    cpl_frameset        *   rawframes ;
    cpl_size            *   labels ;
    cpl_size                nlabels ;
    cpl_frameset        *   raw_one ;
    char                *   setting_id ;
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
            "cr2res.cr2res_cal_detlin.trace_smooth");
    trace_smooth = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_detlin.trace_opening");
    trace_opening = cpl_parameter_get_bool(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_detlin.trace_collapse");
    trace_collapse = cpl_parameter_get_bool(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_cal_detlin.detector");
    reduce_det = cpl_parameter_get_int(param);

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

    /* Loop on the settings */
    for (l=0 ; l<(int)nlabels ; l++) {
        /* Get the frames for the current setting */
        raw_one = cpl_frameset_extract(rawframes, labels, (cpl_size)l) ;

        /* Get the current setting */
        plist = cpl_propertylist_load(cpl_frame_get_filename(
                    cpl_frameset_get_position(raw_one, 0)), 0) ;
        setting_id = cpl_strdup(cr2res_pfits_get_wlen_id(plist)) ;
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
                        trace_degree, trace_min_cluster, trace_smooth, 
                        trace_opening, trace_collapse, det_nr,
                        &(coeffs[det_nr-1]),
                        &(bpm[det_nr-1]),
                        &(ext_plist[det_nr-1])) == -1) {
                cpl_msg_warning(__func__, 
                        "Failed to reduce SETTING %s / det %d", 
                        setting_id, det_nr);
            } 
            cpl_msg_indent_less() ;
        }

        /* Save the products */

        /* BPM */
        out_file = cpl_sprintf("%s_%c_bpm.fits", RECIPE_STRING,
                setting_id[0]) ;
        cr2res_io_save_BPM(out_file, frameset, raw_one, parlist, bpm, NULL, 
                ext_plist, CR2RES_DETLIN_BPM_PROCATG, RECIPE_STRING) ;
        cpl_free(out_file);

        /* COEFFS */
        out_file = cpl_sprintf("%s_%c_coeffs.fits", RECIPE_STRING,
                setting_id[0]) ;
        cr2res_io_save_DETLIN_COEFFS(out_file, frameset, raw_one, parlist, 
                coeffs, NULL, ext_plist, CR2RES_DETLIN_COEFFS_PROCATG, 
                RECIPE_STRING) ;
        cpl_free(out_file);

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
    cpl_frameset_delete(rawframes) ;
    cpl_free(labels) ;

    return (int)cpl_error_get_code();
}

/*----------------------------------------------------------------------------*/
/**
  @brief  
  @param 
  @return  
 */
/*----------------------------------------------------------------------------*/
static int cr2res_cal_detlin_reduce(
        const cpl_frameset  *   rawframes,
        double                  bpm_kappa,
        int                     trace_degree,
        int                     trace_min_cluster,
        double                  trace_smooth,
        int                     trace_opening,
        int                     trace_collapse,
        int                     reduce_det,
        hdrl_imagelist      **  coeffs,
        cpl_image           **  bpm,
        cpl_propertylist    **  ext_plist)
{
    const char          *   first_file ;
    cpl_imagelist       *   imlist ;
    cpl_image           *   cur_im ;
    float               *   pcur_im ;
    cpl_vector          *   dits ;
    hdrl_image          *   collapsed ;
    cpl_image           *   collapsed_ima ;
    hdrl_image          *   master_flat_loc ;
    cpl_table           *   traces ;
    cpl_imagelist       *   coeffs_loc ;
    cpl_image           *   cur_coeffs ;
    float               *   pcur_coeffs ;
    cpl_imagelist       *   errors_loc ;
    cpl_image           *   cur_errors ;
    float               *   pcur_errors ;
    cpl_propertylist    *   plist ;
    cpl_image           *   bpm_loc ;
    int                 *   pbpm_loc ;
    cpl_image           *   trace_image ;
    int                 *   pti ;
    cpl_polynomial      *   fit1d ;
	cpl_matrix			*	samppos ;
    cpl_vector          *   fitvals ;
    cpl_vector          *   fit_residuals ;
    cpl_boolean             sampsym ;
    cpl_mask            *   bpm_mask ;
    int                     i, j, k, idx, ext_nr, order, trace_id, nx, ny ;
    cpl_size                max_degree, l ;
    double                  low_thresh, high_thresh, median, sigma ;
    int                     qc_nb_bad, qc_nbfailed, qc_nbsuccess ;
    double                  qc_median, qc_fitquality, qc_gain,
                            qc_min_level, qc_max_level ;
    
    /* Check Inputs */
    if (rawframes == NULL) return -1 ;

    /* Initialise */
    max_degree = 2 ;
    sampsym = CPL_TRUE ;

    /* Get the Extension number */
    first_file = cpl_frame_get_filename(
            cpl_frameset_get_position_const(rawframes, 0)) ;
    ext_nr = cr2res_io_get_ext_idx(first_file, reduce_det, 1) ;

    /* Load the image list */
    imlist = cpl_imagelist_load_frameset(rawframes, CPL_TYPE_FLOAT, 1, ext_nr) ;
    if (imlist == NULL) {
        cpl_msg_error(__func__, "Failed to Load the images") ;
        return -1 ;
    }

    /* Load the DITs */
    dits = cpl_vector_new(cpl_frameset_get_size(rawframes)) ;
    for (i=0 ; i< cpl_vector_get_size(dits) ; i++) {
        plist = cpl_propertylist_load(
                cpl_frame_get_filename(
                    cpl_frameset_get_position_const(rawframes, i)), 0) ;
        cpl_vector_set(dits, i, cr2res_pfits_get_dit(plist)) ;
        cpl_propertylist_delete(plist) ;
    }
            
    if (cpl_msg_get_level() == CPL_MSG_DEBUG) cpl_vector_dump(dits, stdout) ;

    /* Load the extension header for saving */
    plist = cpl_propertylist_load(first_file, ext_nr) ;
    if (plist == NULL) return -1 ;

    if (trace_collapse) {
        /* Collapse all input images */
        cpl_msg_info(__func__, "Collapse the input images") ;
        cpl_msg_indent_more() ;
        if ((collapsed_ima = cpl_imagelist_collapse_create(imlist)) == NULL) {
            cpl_msg_error(__func__, "Failed to Calibrate and collapse") ;
            cpl_imagelist_delete(imlist) ;
            cpl_vector_delete(dits); 
            cpl_propertylist_delete(plist);
            cpl_msg_indent_less() ;
            return -1 ;
        }
        cpl_msg_indent_less() ;
    } else {
        /* Only use the first image */
        collapsed_ima = cpl_image_duplicate(cpl_imagelist_get(imlist, 0)) ;
    }

    nx = cpl_image_get_size_x(collapsed_ima) ;
    ny = cpl_image_get_size_y(collapsed_ima) ;

    /* Create the HDRL collapsed */
    collapsed = hdrl_image_create(collapsed_ima, NULL) ;
    cpl_image_delete(collapsed_ima) ;

    /* Compute traces */
    cpl_msg_info(__func__, "Compute the traces") ;
    cpl_msg_indent_more() ;
    if ((traces = cr2res_trace(hdrl_image_get_image(collapsed), 
                    trace_smooth, trace_opening, trace_degree, 
                    trace_min_cluster)) == NULL) {
        cpl_msg_error(__func__, "Failed compute the traces") ;
        cpl_imagelist_delete(imlist) ;
        cpl_vector_delete(dits); 
        cpl_propertylist_delete(plist);
        hdrl_image_delete(collapsed) ;
        cpl_msg_indent_less() ;
        return -1 ;
    }
    hdrl_image_delete(collapsed) ;
    cpl_msg_indent_less() ;

    /* Allocate */
    bpm_loc = cpl_image_new(nx, ny, CPL_TYPE_INT) ;
    pbpm_loc = cpl_image_get_data_int(bpm_loc) ;
    coeffs_loc = cpl_imagelist_new() ;
    errors_loc = cpl_imagelist_new() ;

    /* Initialise the coeffs cube */
    for (l=0 ; l<=max_degree ; l++) {
        cpl_imagelist_set(coeffs_loc, cpl_image_new(nx, ny, CPL_TYPE_FLOAT),l) ;
        cpl_imagelist_set(errors_loc, cpl_image_new(nx, ny, CPL_TYPE_FLOAT),l) ;
    }

    /* Initialise the BPM as Out of Order */
    cpl_image_add_scalar(bpm_loc, CR2RES_BPM_OUTOFORDER); 

    /* Create the trace image */
    trace_image = cr2res_trace_gen_image(traces, nx, ny) ;
    pti = cpl_image_get_data_int(trace_image) ;

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
                /* We are in a trace, let's compute the linearity */

                /* Store the DITS */
				samppos = cpl_matrix_wrap(1,
                            cpl_vector_get_size(dits),
                            cpl_vector_get_data(dits)) ;

                /* Prepare values to fit */
                fitvals = cpl_vector_new(cpl_vector_get_size(dits)) ;
                for (k=0 ; k<cpl_vector_get_size(dits) ; k++) {
                    cur_im = cpl_imagelist_get(imlist,  k) ;
                    pcur_im = cpl_image_get_data_float(cur_im) ;
                    cpl_vector_set(fitvals, k, pcur_im[idx]) ;
                }

                /* Fit  */
                fit1d = cpl_polynomial_new(1);
                if (cpl_polynomial_fit(fit1d, samppos, &sampsym, fitvals, NULL,
                        CPL_FALSE, NULL, &max_degree) != CPL_ERROR_NONE) {
                    /* Failed Fit - Fill the coefficientѕ */
                    pbpm_loc[idx] = CR2RES_BPM_DETLIN ;
                    for (l=0 ; l<=max_degree ; l++) {
                        cur_coeffs = cpl_imagelist_get(coeffs_loc, l) ;
                        pcur_coeffs = cpl_image_get_data_float(cur_coeffs) ;
                        pcur_coeffs[idx] = 0.0 ;
                        cur_errors = cpl_imagelist_get(errors_loc, l) ;
                        pcur_errors = cpl_image_get_data_float(cur_errors) ;
                        pcur_errors[idx] = 0.0 ;
                    } 
                    qc_nbfailed++ ;
                    cpl_error_reset() ;
                    continue ;
                }
                qc_nbsuccess++;

                /* Compute the residuals */
                fit_residuals = cpl_vector_new(cpl_vector_get_size(dits)) ;
                cpl_vector_fill_polynomial_fit_residual(fit_residuals,
                        fitvals, NULL, fit1d, samppos, NULL) ;
                cpl_matrix_unwrap(samppos) ;
                cpl_vector_delete(fitvals) ;
                cpl_vector_delete(fit_residuals) ;

                /* Store the Coefficients in the output image list */
                pbpm_loc[idx] = 0 ;
                for (l=0 ; l<=max_degree ; l++) {
                    cur_coeffs = cpl_imagelist_get(coeffs_loc, l) ;
                    pcur_coeffs = cpl_image_get_data_float(cur_coeffs) ;
                    pcur_coeffs[idx] = cpl_polynomial_get_coeff(fit1d, &l) ;
                    cur_errors = cpl_imagelist_get(errors_loc, l) ;
                    pcur_errors = cpl_image_get_data_float(cur_errors) ;
                    /* TODO : Store error */
                    pcur_errors[idx] = 0.0 ;
                }
                cpl_polynomial_delete(fit1d) ;
            }
        }
    }
    cpl_msg_indent_less() ;

    /* Use the second coefficient stats for the BPM detection */
    cpl_msg_info(__func__, "BPM detection") ;
    cur_coeffs = cpl_imagelist_get(coeffs_loc, 1) ;
    pcur_coeffs = cpl_image_get_data_float(cur_coeffs) ;
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
    /* TODO : pass the proper inputs */
    qc_fitquality = 0.0;
    qc_median = cr2res_qc_detlin_median(coeffs_loc) ;
    qc_gain = cr2res_qc_detlin_gain(coeffs_loc) ;
    cr2res_qc_detlin_min_max_level(NULL, &qc_min_level, &qc_max_level) ;

    /* Store the QC parameters in the plist */
    cpl_propertylist_append_int(plist, "ESO QC DETLIN NBBAD", qc_nb_bad) ;
    cpl_propertylist_append_int(plist, "ESO QC DETLIN NBFAILED", qc_nbfailed) ;
    cpl_propertylist_append_int(plist, "ESO QC DETLIN NBSUCCESS",qc_nbsuccess) ;
    cpl_propertylist_append_double(plist, "ESO QC DETLIN FIT_QUALITY", 
            qc_fitquality) ;
    cpl_propertylist_append_double(plist, "ESO QC DETLIN MEDIAN", qc_median) ;
    cpl_propertylist_append_double(plist, "ESO QC DETLIN GAIN", qc_gain) ;
    cpl_propertylist_append_double(plist, "ESO QC DETLIN MIN_LEVEL",
            qc_min_level) ;
    cpl_propertylist_append_double(plist, "ESO QC DETLIN MAX_LEVEL",
            qc_max_level) ;

    /* Free */
    cpl_image_delete(trace_image) ;
    cpl_table_delete(traces) ;
    cpl_imagelist_delete(imlist) ;
    cpl_vector_delete(dits); 

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



/* TODO: Move/Use this function in cr2res_util_merge_detlin */
/*----------------------------------------------------------------------------*/
/**
  @brief Only pixels not yet computed (CR2RES_BPM_OUTOFORDER) are updated
  @param 
  @return   0 if ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_cal_detlin_update(
        cpl_image           *   current_bpm,
        cpl_image           *   new_bpm,
        cpl_imagelist       *   current_coeffs,
        cpl_imagelist       *   new_coeffs)
{
    cpl_size        i, j, k, nx, ny, ni, idx ;
    int         *   pcurrent_bpm ;
    int         *   pnew_bpm ;
    cpl_image   *   current_coeffs_ima ;
    cpl_image   *   new_coeffs_ima ;
    float       *   pcurrent_coeffs_ima ;
    float       *   pnew_coeffs_ima ;
    
    /* Check Inputs */
    if (current_bpm == NULL || new_bpm == NULL || current_coeffs == NULL
            || new_coeffs == NULL) return 0 ;

    /* Initialise */
    nx = cpl_image_get_size_x(current_bpm) ;
    ny = cpl_image_get_size_y(current_bpm) ;
    ni = cpl_imagelist_get_size(current_coeffs) ;

    if (cpl_image_get_size_x(new_bpm) != nx ||
            cpl_image_get_size_y(new_bpm) != ny ||
            cpl_imagelist_get_size(new_coeffs) != ni) {
        return -1 ;
    }

    pcurrent_bpm = cpl_image_get_data_int(current_bpm) ;
    pnew_bpm = cpl_image_get_data_int(new_bpm) ;

    /* Loop on the pixels */
    for (j=0 ; j<ny ; j++) {
        for (i=0 ; i<nx ; i++) {
            idx = i+j*nx ;
            if (pcurrent_bpm[idx] == CR2RES_BPM_OUTOFORDER &&
                    pnew_bpm[idx] != CR2RES_BPM_OUTOFORDER) {
                /* Put the new value */
                pcurrent_bpm[idx] = pnew_bpm[idx] ;
                for (k=0 ; k<ni ; k++) {
                    current_coeffs_ima = cpl_imagelist_get(current_coeffs, k) ;
                    new_coeffs_ima = cpl_imagelist_get(new_coeffs, k) ;
                    pcurrent_coeffs_ima = 
                        cpl_image_get_data_float(current_coeffs_ima) ;
                    pnew_coeffs_ima = cpl_image_get_data_float(new_coeffs_ima) ;
                    pcurrent_coeffs_ima[idx] = pnew_coeffs_ima[idx] ;
                }
            }
        }
    }
    return 0 ;
}

