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

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

#define RECIPE_STRING "cr2res_obs_1d"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_obs_1d_reduce(
        const cpl_frameset  *   rawframes,
        const cpl_frame     *   detlin_frame,
        const cpl_frame     *   master_dark_frame,
        const cpl_frame     *   master_flat_frame,
        const cpl_frame     *   bpm_frame,
        int                     calib_cosmics_corr,
        int                     extract_oversample,
        int                     extract_swath_width,
        int                     extract_height,
        double                  extract_smooth,
        int                     extract_sum_only,
        int                     reduce_det,
        hdrl_image          **  combineda,
        hdrl_image          **  combinedb,
        cpl_table           **  extracta,
        cpl_table           **  extractb,
        cpl_propertylist    **  ext_plist) ;
static int cr2res_obs_1d_create(cpl_plugin *);
static int cr2res_obs_1d_exec(cpl_plugin *);
static int cr2res_obs_1d_destroy(cpl_plugin *);
static int cr2res_obs_1d(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_obs_1d_description[] =
"CRIRES+ 1d Observation recipe\n"
"The files listed in the Set Of Frames (sof-file) must be tagged:\n"
"raw-file.fits " CR2RES_OBS_1D_RAW"\n"
"detlin.fits " CR2RES_DETLIN_COEFFS_PROTYPE "\n"
"master_dark.fits " CR2RES_MASTER_DARK_PROTYPE "\n"
"master_flat.fits " CR2RES_MASTER_FLAT_PROTYPE "\n"
"bpm.fits " CR2RES_BPM_PROTYPE "\n"
"trace_wave.fits " CR2RES_TRACE_WAVE_PROTYPE "\n"
" The recipe produces the following products:\n"
"cr2res_obs_1d_extractA.fits " CR2RES_OBS_1D_EXTRACTA_PROCATG "\n"
"cr2res_obs_1d_extractB.fits " CR2RES_OBS_1D_EXTRACTB_PROCATG "\n"
"cr2res_obs_1d_combinedA.fits " CR2RES_OBS_1D_COMBINEDA_PROCATG "\n"
"cr2res_obs_1d_combinedB.fits " CR2RES_OBS_1D_COMBINEDB_PROCATG "\n"
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
                    "cr2res_obs_1d",
                    "1D Observation recipe",
                    cr2res_obs_1d_description,
                    "Thomas Marquart, Yves Jung",
                    PACKAGE_BUGREPORT,
                    cr2res_get_license(),
                    cr2res_obs_1d_create,
                    cr2res_obs_1d_exec,
                    cr2res_obs_1d_destroy)) {    
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
static int cr2res_obs_1d_create(cpl_plugin * plugin)
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
    p = cpl_parameter_new_value("cr2res.cr2res_obs_1d.extract_oversample",
            CPL_TYPE_INT, "factor by which to oversample the extraction",
            "cr2res.cr2res_obs_1d", 2);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_oversample");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_1d.extract_swath_width",
            CPL_TYPE_INT, "The swath width", "cr2res.cr2res_obs_1d", 24);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_swath_width");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_1d.extract_height",
            CPL_TYPE_INT, "Extraction height",
            "cr2res.cr2res_obs_1d", -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_height");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_1d.extract_smooth",
            CPL_TYPE_DOUBLE,
            "Smoothing along the slit (1 for high S/N, 5 for low)",
            "cr2res.cr2res_obs_1d", 1.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "extract_smooth");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_1d.detector",
            CPL_TYPE_INT, "Only reduce the specified detector",
            "cr2res.cr2res_obs_1d", 0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "detector");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_1d.order",
            CPL_TYPE_INT, "Only reduce the specified order",
            "cr2res.cr2res_obs_1d", -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "order");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_obs_1d.trace_nb",
            CPL_TYPE_INT, "Only reduce the specified trace number",
            "cr2res.cr2res_obs_1d", -1);
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
static int cr2res_obs_1d_exec(cpl_plugin * plugin)
{
    cpl_recipe  *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1;

    return cr2res_obs_1d(recipe->frames, recipe->parameters);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_obs_1d_destroy(cpl_plugin * plugin)
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
static int cr2res_obs_1d(
        cpl_frameset            *   frameset,
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter *   param ;
    int                     extract_oversample, extract_swath_width,
                            extract_height, reduce_det, reduce_order, 
                            reduce_trace ;
    double                  extract_smooth ;
    cpl_frameset        *   rawframes ;
    const cpl_frame     *   detlin_frame ;
    const cpl_frame     *   master_dark_frame ;
    const cpl_frame     *   master_flat_frame ;
    const cpl_frame     *   bpm_frame ;
    const cpl_frame     *   trace_wave_frame ;
    hdrl_image          *   combineda[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   combinedb[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   ext_plist[CR2RES_NB_DETECTORS] ;
    cpl_table           *   extracta[CR2RES_NB_DETECTORS] ;
    cpl_table           *   extractb[CR2RES_NB_DETECTORS] ;
    char                *   out_file;
    int                     i, det_nr; 

    /* Initialise */

    /* RETRIEVE INPUT PARAMETERS */
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_1d.extract_oversample");
    extract_oversample = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_1d.extract_swath_width");
    extract_swath_width = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_1d.extract_height");
    extract_height = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_1d.extract_smooth");
    extract_smooth = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_1d.detector");
    reduce_det = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_1d.order");
    reduce_order = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_obs_1d.trace_nb");
    reduce_trace = cpl_parameter_get_int(param);

    /* Identify the RAW and CALIB frames in the input frameset */
    if (cr2res_dfs_set_groups(frameset)) {
        cpl_msg_error(__func__, "Cannot identify RAW and CALIB frames") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }
	
    /* Get Calibration frames */
    detlin_frame = cpl_frameset_find_const(frameset,
            CR2RES_DETLIN_COEFFS_PROCATG);
    master_dark_frame = cpl_frameset_find_const(frameset,
            CR2RES_MASTER_DARK_PROCATG) ; 
    master_flat_frame = cpl_frameset_find_const(frameset,
            CR2RES_FLAT_MASTER_FLAT_PROCATG) ; 
    bpm_frame = cpl_frameset_find_const(frameset,
            CR2RES_FLAT_BPM_PROCATG) ;
    trace_wave_frame = cpl_frameset_find_const(frameset,
            CR2RES_FLAT_TRACE_WAVE_PROCATG) ;

    /* Get the RAW Frames */
    rawframes = cr2res_extract_frameset(frameset, CR2RES_OBS_1D_RAW) ;
    if (rawframes == NULL) {
        cpl_msg_error(__func__, "Could not find RAW frames") ;
        return -1 ;
    }
    cpl_msg_indent_more() ;

    /* Loop on the detectors */
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
        /* Initialise */
        combineda[det_nr-1] = NULL ;
        combinedb[det_nr-1] = NULL ;
        extracta[det_nr-1] = NULL ;
        extractb[det_nr-1] = NULL ;
        ext_plist[det_nr-1] = NULL ;

        /* Compute only one detector */
        if (reduce_det != 0 && det_nr != reduce_det) continue ;
    
        cpl_msg_info(__func__, "Process Detector %d", det_nr) ;
        cpl_msg_indent_more() ;

        /* Call the reduction function */
        if (cr2res_obs_1d_reduce(rawframes, detlin_frame, 
                    master_dark_frame, master_flat_frame, bpm_frame, 0,
                    extract_oversample, extract_swath_width, extract_height, 
                    extract_smooth, 0, det_nr,
                    &(combineda[det_nr-1]),
                    &(combinedb[det_nr-1]),
                    &(extracta[det_nr-1]),
                    &(extractb[det_nr-1]),
                    &(ext_plist[det_nr-1])) == -1) {
            cpl_msg_warning(__func__, "Failed to reduce detector %d", det_nr);
        }
        cpl_msg_indent_less() ;
    }

    /* Ѕave Products */

    /* Extracted A */
    out_file = cpl_sprintf("%s_extractedA.fits", RECIPE_STRING) ;
    cr2res_io_save_EXTRACT_1D(out_file, frameset, rawframes, parlist, extracta,
            NULL, ext_plist, CR2RES_OBS_1D_EXTRACTA_PROCATG, RECIPE_STRING) ;
    cpl_free(out_file);

    /* Free */
    cpl_frameset_delete(rawframes) ;
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
        if (combineda[det_nr-1] != NULL)
            hdrl_image_delete(combineda[det_nr-1]) ;
        if (combinedb[det_nr-1] != NULL)
            hdrl_image_delete(combinedb[det_nr-1]) ;
        if (extracta[det_nr-1] != NULL) 
            cpl_table_delete(extracta[det_nr-1]) ;
        if (extractb[det_nr-1] != NULL) 
            cpl_table_delete(extractb[det_nr-1]) ;
        if (ext_plist[det_nr-1] != NULL) 
            cpl_propertylist_delete(ext_plist[det_nr-1]) ;
    }
    cpl_msg_indent_less() ;

    return (int)cpl_error_get_code();
}
 
/*----------------------------------------------------------------------------*/
/**
  @brief  
  @param 
  @return  
 */
/*----------------------------------------------------------------------------*/
static int cr2res_obs_1d_reduce(
        const cpl_frameset  *   rawframes,
        const cpl_frame     *   detlin_frame,
        const cpl_frame     *   master_dark_frame,
        const cpl_frame     *   master_flat_frame,
        const cpl_frame     *   bpm_frame,
        int                     calib_cosmics_corr,
        int                     extract_oversample,
        int                     extract_swath_width,
        int                     extract_height,
        double                  extract_smooth,
        int                     extract_sum_only,
        int                     reduce_det,
        hdrl_image          **  combineda,
        hdrl_image          **  combinedb,
        cpl_table           **  extracta,
        cpl_table           **  extractb,
        cpl_propertylist    **  ext_plist)
{
    hdrl_imagelist      *   in ;
    hdrl_imagelist      *   in_calib ;
    hdrl_imagelist      *   in_noda ;
    hdrl_imagelist      *   in_nodb ;
    cr2res_nodding_pos  *   nod_positions ;
    cpl_vector          *   dits ;
    cpl_size                nframes, i ;
    int                     det_nr ;
	
    /* Check Inputs */
    if (combineda == NULL || combinedb == NULL || extracta == NULL ||
            extractb == NULL || ext_plist == NULL || rawframes == NULL)
        return -1 ;

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
    if (master_dark_frame != NULL)  dits = cr2res_read_dits(rawframes) ;
    else                            dits = NULL ;
    if (cpl_msg_get_level() == CPL_MSG_DEBUG && dits != NULL) 
        cpl_vector_dump(dits, stdout) ;

    /* Load image list */
    in = cr2res_io_load_image_list_from_set(rawframes, reduce_det) ;

    /* Calibrate the images */
    in_calib = cr2res_calib_imagelist(in, reduce_det, 0, master_flat_frame,
            master_dark_frame, bpm_frame, detlin_frame, dits) ;

    /* Split the image lists */
    cr2res_combine_nodding_split(in_calib, nod_positions, &in_noda, 
            &in_nodb) ; 

    if (dits != NULL) cpl_vector_delete(dits) ;

    *combineda = *combinedb = *extracta = *extractb = NULL ;

    return 0 ;
}





