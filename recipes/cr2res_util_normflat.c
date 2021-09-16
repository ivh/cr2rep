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
#include "cr2res_flat.h"
#include "cr2res_bpm.h"
#include "cr2res_trace.h"
#include "cr2res_extract.h"
#include "cr2res_io.h"

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

#define RECIPE_STRING "cr2res_util_normflat"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static cpl_frameset * cr2res_util_normflat_find_RAW(const cpl_frameset * in) ;
static int cr2res_util_normflat_reduce(
        const cpl_frameset  *   rawframes,
        const cpl_frame     *   slitmodel_frame,
        double                  bpm_low,
        double                  bpm_high,
        double                  bpm_linemax,
        int                     reduce_det,
        hdrl_image          **  master_flat,
        cpl_image           **  bpm,
        cpl_propertylist    **  ext_plist) ;
static int cr2res_util_normflat_compare(
        const cpl_frame   *   frame1,
        const cpl_frame   *   frame2) ;
static int cr2res_util_normflat_create(cpl_plugin *);
static int cr2res_util_normflat_exec(cpl_plugin *);
static int cr2res_util_normflat_destroy(cpl_plugin *);
static int cr2res_util_normflat(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_util_normflat_description[] = "\
Flat normalization                                                      \n\
  The input RAW files are grouped by setting/decker and each group is   \n\
  reduced separately. For each group, a slit model file with matching   \n\
  setting/decker is expected.                                           \n\
                                                                        \n\
  Inputs                                                                \n\
    raw.fits " CR2RES_FLAT_RAW " [1 to n]                               \n\
          or " CR2RES_UTIL_CALIB_PROCATG "                              \n\
    slit_model.fits " CR2RES_CAL_FLAT_SLIT_MODEL_PROCATG " [1 to m]     \n\
                 or " CR2RES_UTIL_SLIT_MODEL_PROCATG "                  \n\
                 or " CR2RES_OBS_NODDING_SLITMODELA_PROCATG "           \n\
                 or " CR2RES_OBS_NODDING_SLITMODELB_PROCATG "           \n\
                                                                        \n\
  Outputs                                                               \n\
    cr2res_util_normflat_[setting]_[Decker]_bpm.fits "
    CR2RES_UTIL_NORM_BPM_PROCATG "\n\
    cr2res_util_normflat_[setting]_[Decker]_master.fits "
    CR2RES_UTIL_MASTER_FLAT_PROCATG "\n\
                                                                        \n\
  Algorithm                                                             \n\
    group the input frames by different settings                        \n\
    loop on groups g:                                                   \n\
      group the input frames by different decker positions              \n\
      loop on decker positions p:                                       \n\
        loop on detectors d:                                            \n\
          cr2res_util_normflat_reduce() computes (master_flat,bpm)(g,p,d)\n\
      Save master_flat(g,p)                                             \n\
      Save bpm(g,p)                                                     \n\
                                                                        \n\
    cr2res_util_normflat_reduce()                                       \n\
      Load the images list                                              \n\
      Average the images to avg                                         \n\
      Load the input slit_model with the proper setting/decker          \n\
      Compute the master flat with cr2res_master_flat(avg,              \n\
               slit_model, --bpm_low, --bpm_high, --bpm_lines_ratio)    \n\
        -> master_flat, bpm                                             \n\
                                                                        \n\
Library functions used:                                                 \n\
    cr2res_extract_frameset()                                           \n\
    cr2res_io_extract_decker_frameset()                                 \n\
    cr2res_io_find_SLIT_MODEL()                                         \n\
    cr2res_io_load_image_list_from_set()                                \n\
    cr2res_io_load_SLIT_MODEL()                                         \n\
    cr2res_master_flat()                                                \n\
    cr2res_io_save_MASTER_FLAT()                                        \n\
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
                    "Flat Normalization utility",
                    cr2res_util_normflat_description,
                    CR2RES_PIPELINE_AUTHORS,
                    PACKAGE_BUGREPORT,
                    cr2res_get_license(),
                    cr2res_util_normflat_create,
                    cr2res_util_normflat_exec,
                    cr2res_util_normflat_destroy)) {    
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
static int cr2res_util_normflat_create(cpl_plugin * plugin)
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
    p = cpl_parameter_new_value("cr2res.cr2res_util_normflat.bpm_low",
            CPL_TYPE_DOUBLE, "Low threshold for BPM detection",
            "cr2res.cr2res_util_normflat", 0.5);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "bpm_low");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_normflat.bpm_high",
            CPL_TYPE_DOUBLE, "High threshold for BPM detection",
            "cr2res.cr2res_util_normflat", 2.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "bpm_high");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_normflat.bpm_lines_ratio",
            CPL_TYPE_DOUBLE, "Maximum ratio of bad pixels per line",
            "cr2res.cr2res_util_normflat", 0.5);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "bpm_lines_ratio");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);


    p = cpl_parameter_new_value("cr2res.cr2res_util_normflat.detector",
            CPL_TYPE_INT, "Only reduce the specified detector",
            "cr2res.cr2res_util_normflat", 0);
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
static int cr2res_util_normflat_exec(cpl_plugin * plugin)
{
    cpl_recipe  *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1;

    return cr2res_util_normflat(recipe->frames, recipe->parameters);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_util_normflat_destroy(cpl_plugin * plugin)
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
static int cr2res_util_normflat(
        cpl_frameset            *   frameset,
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter *   param ;
    int                     reduce_det ;
    double                  bpm_low, bpm_high, bpm_lines_ratio ;
    cpl_frameset        *   rawframes ;
    const char          *   used_tag ;
    cpl_frameset        *   raw_one_setting ;
    cpl_frameset        *   raw_one_setting_decker ;
    cpl_size            *   labels ;
    cpl_size                nlabels ;
    cpl_propertylist    *   plist ;
    char                *   setting_id ;
    cpl_frame           *   slitmodel_frame ;
    hdrl_image          *   master_flat[CR2RES_NB_DETECTORS] ;
    cpl_image           *   bpm[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   ext_plist[CR2RES_NB_DETECTORS] ;
    char                *   out_file;
    int                     i, det_nr, l ; 

    /* Initialise */
    cr2res_decker decker_values[CR2RES_NB_DECKER_POSITIONS] = 
        {CR2RES_DECKER_NONE, CR2RES_DECKER_1_3, CR2RES_DECKER_2_4} ; 
    char * decker_desc[CR2RES_NB_DECKER_POSITIONS] =
        {"Open", "Decker1", "Decker2"} ;

    /* RETRIEVE INPUT PARAMETERS */
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_normflat.bpm_low");
    bpm_low = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_normflat.bpm_high");
    bpm_high = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_normflat.bpm_lines_ratio");
    bpm_lines_ratio = cpl_parameter_get_double(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_normflat.detector");
    reduce_det = cpl_parameter_get_int(param);

    /* Identify the RAW and CALIB frames in the input frameset */
    if (cr2res_dfs_set_groups(frameset)) {
        cpl_msg_error(__func__, "Cannot identify RAW and CALIB frames") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Extract RAW frames */
    rawframes = cr2res_util_normflat_find_RAW(frameset) ;
    if (rawframes==NULL || cpl_frameset_get_size(rawframes) <= 0) {
        cpl_msg_error(__func__, "Cannot find any RAW file") ;
        cpl_error_set(__func__, CPL_ERROR_DATA_NOT_FOUND) ;
        return -1 ;
    }
    /* Keep track of the first file tag (should all be the same) */
    used_tag = cpl_frame_get_tag(cpl_frameset_get_position_const(rawframes,0));

    /* Labelise the raw frames with the different settings */
    if ((labels = cpl_frameset_labelise(rawframes, cr2res_util_normflat_compare,
                &nlabels)) == NULL) {
        cpl_msg_error(__func__, "Cannot labelise input frames") ;
        cpl_frameset_delete(rawframes) ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Loop on the settings */
    for (l=0 ; l<(int)nlabels ; l++) {
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

            /* Get the slit model for the setting / decker */
            slitmodel_frame = cr2res_io_find_SLIT_MODEL(frameset,
                    setting_id, decker_values[i]);

            /* Loop on the detectors */
            for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
                /* Initialise */
                master_flat[det_nr-1] = NULL ;
                bpm[det_nr-1] = NULL ;
                ext_plist[det_nr-1] = NULL ;

                /* Compute only one detector */
                if (reduce_det != 0 && det_nr != reduce_det) continue ;
            
                cpl_msg_info(__func__, "Process Detector %d", det_nr) ;
                cpl_msg_indent_more() ;

                /* Call the reduction function */
                if (cr2res_util_normflat_reduce(raw_one_setting_decker, 
                            slitmodel_frame, bpm_low, bpm_high, 
                            bpm_lines_ratio, det_nr,
                            &(master_flat[det_nr-1]),
                            &(bpm[det_nr-1]),
                            &(ext_plist[det_nr-1])) == -1) {
                    cpl_msg_warning(__func__, 
                            "Failed to reduce detector %d of %s Frames", 
                            det_nr, decker_desc[i]);
                }
                cpl_msg_indent_less() ;
            }
            cpl_msg_indent_less() ;
            cpl_frame_delete(slitmodel_frame) ;

            /* Save Products */

            /* MASTER_FLAT */
            if (nlabels == 1) {
                out_file = cpl_sprintf("%s_%s_master_flat.fits", 
                        RECIPE_STRING, decker_desc[i]) ;
            } else {
                out_file = cpl_sprintf("%s_%s_%s_master_flat.fits", 
                        RECIPE_STRING, setting_id, decker_desc[i]) ;
            }
            cr2res_io_save_MASTER_FLAT(out_file, frameset,
                    raw_one_setting_decker, parlist, master_flat, NULL,
                    ext_plist, CR2RES_UTIL_MASTER_FLAT_PROCATG,
                    RECIPE_STRING);
            cpl_free(out_file);

            /* BPM */
            if (nlabels == 1) {
                out_file = cpl_sprintf("%s_%s_master_bpm.fits", 
                        RECIPE_STRING, decker_desc[i]) ;
            } else {
                out_file = cpl_sprintf("%s_%s_%s_master_bpm.fits", 
                        RECIPE_STRING, setting_id, decker_desc[i]) ;
            }
            cr2res_io_save_BPM(out_file, frameset,
                    raw_one_setting_decker, parlist, bpm, NULL,ext_plist,
                    CR2RES_UTIL_NORM_BPM_PROCATG, RECIPE_STRING) ;
            cpl_free(out_file);

            /* Free */
            cpl_frameset_delete(raw_one_setting_decker) ;
            for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
                if (master_flat[det_nr-1] != NULL)
                    hdrl_image_delete(master_flat[det_nr-1]) ;
                if (bpm[det_nr-1] != NULL)
                    cpl_image_delete(bpm[det_nr-1]) ;
                if (ext_plist[det_nr-1] != NULL)
                    cpl_propertylist_delete(ext_plist[det_nr-1]) ;
            }
        }
        cpl_msg_indent_less() ;
        cpl_frameset_delete(raw_one_setting) ;
        cpl_free(setting_id) ;
    }
    cpl_free(labels);
    cpl_frameset_delete(rawframes) ;
    return (int)cpl_error_get_code();
}

/*----------------------------------------------------------------------------*/
/**
  @brief  Normalise the flat for 1 setting, 1 decker position, 1 detector
  @param rawframes          Input raw frames (same setting, same decker)
  @param slitmodel_frame    Slit model frame for this setting/decker
  @param bpm_low            Threshold for BPM detection
  @param bpm_high           Threshold for BPM detection
  @param bpm_linemax        Max fraction of BPM per line
  @param reduce_det            The detector to compute
  @param master_flat        [out] Master flat
  @param bpm                [out] the BPM
  @param ext_plist          [out] the header for saving the products
 */
/*----------------------------------------------------------------------------*/
static int cr2res_util_normflat_reduce(
        const cpl_frameset  *   rawframes,
        const cpl_frame     *   slitmodel_frame,
        double                  bpm_low,
        double                  bpm_high,
        double                  bpm_linemax,
        int                     reduce_det,
        hdrl_image          **  master_flat,
        cpl_image           **  bpm,
        cpl_propertylist    **  ext_plist)
{
    const char          *   first_file ;
    hdrl_imagelist      *   imlist ;
    hdrl_image          *   collapsed ;
    hdrl_image          *   slit_model ;
    cpl_image           *   contrib ;
    cpl_image           *   bpm_flat ;
    cpl_propertylist    *   plist ;
    hdrl_image          *   master_flat_loc ;
    int                     i, ext_nr ;
    
    /* Check Inputs */
    if (rawframes == NULL) return -1 ;

    /* Get the Extension number */
    first_file = cpl_frame_get_filename(
            cpl_frameset_get_position_const(rawframes, 0)) ;
    ext_nr = cr2res_io_get_ext_idx(first_file, reduce_det, 1) ;

    /* Load the extension header for saving */
    plist = cpl_propertylist_load(first_file, ext_nr) ;
    if (plist == NULL) return -1 ;

    /* Load the image list */
    imlist = cr2res_io_load_image_list_from_set(rawframes, reduce_det) ;
    if (imlist == NULL) {
        cpl_msg_error(__func__, "Failed to Load the images") ;
        cpl_propertylist_delete(plist);
        return -1 ;
    }

    /* Collapse */
    cpl_msg_info(__func__, "Collapse the input images") ;
    cpl_msg_indent_more() ;
    if (hdrl_imagelist_collapse_mean(imlist, &collapsed, &contrib) !=
            CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Failed to Collapse") ;
        cpl_propertylist_delete(plist);
        hdrl_imagelist_delete(imlist) ;
        cpl_msg_indent_less() ;
        return -1 ;
    }
    hdrl_imagelist_delete(imlist) ;
    cpl_image_delete(contrib) ;
    cpl_msg_indent_less() ;

    /* Load the Model master */
    if ((slit_model = cr2res_io_load_SLIT_MODEL(
                    cpl_frame_get_filename(slitmodel_frame),
                    reduce_det)) == NULL) {
        cpl_msg_error(__func__, "Cannot load the slit model") ;
        cpl_propertylist_delete(plist);
        hdrl_image_delete(collapsed) ;
        return -1 ;
    }

    /* Compute the Master flat */
    cpl_msg_info(__func__, "Compute the master flat") ;
    cpl_msg_indent_more() ;
    if ((master_flat_loc = cr2res_master_flat(collapsed, slit_model, bpm_low, 
                    bpm_high, bpm_linemax, &bpm_flat)) == NULL) {
        cpl_msg_error(__func__, "Failed compute the Master Flat") ;
        cpl_propertylist_delete(plist);
        hdrl_image_delete(collapsed) ;
        hdrl_image_delete(slit_model) ;
        cpl_msg_indent_less() ;
        return -1 ;
    }
    cpl_msg_indent_less() ;
    hdrl_image_delete(slit_model) ;
    hdrl_image_delete(collapsed) ;

    /* Return the results */
    *master_flat = master_flat_loc ;
    *bpm = bpm_flat ;
    *ext_plist = plist ;
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
static int cr2res_util_normflat_compare(
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
  @brief    Get the RAW frames from a frameset
  @param    set     Input frame set
  @return   the RAW frameset or NULL in error case or if it is missing
    Allowed RAW types : CR2RES_FLAT_RAW
                        CR2RES_UTIL_CALIB_PROCATG
 */
/*----------------------------------------------------------------------------*/
static cpl_frameset * cr2res_util_normflat_find_RAW(const cpl_frameset * in)
{
    cpl_frameset    *   out ;

    /* Check entries */
    if (in == NULL) return NULL ;

    out = cr2res_extract_frameset(in, CR2RES_FLAT_RAW) ;
    if (out == NULL)
        out = cr2res_extract_frameset(in, CR2RES_UTIL_CALIB_PROCATG) ;
    return out ;
}



