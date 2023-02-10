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
#include "hdrl.h"

#include "cr2res_utils.h"
#include "cr2res_pfits.h"
#include "cr2res_dfs.h"
#include "cr2res_io.h"

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

#define RECIPE_STRING "cr2res_util_genlines"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_util_genlines_create(cpl_plugin *);
static int cr2res_util_genlines_exec(cpl_plugin *);
static int cr2res_util_genlines_destroy(cpl_plugin *);
static int cr2res_util_genlines(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_util_genlines_description[] = " \
Generate Lines calibration tables                                       \n\
                                                                        \n\
  Inputs                                                                \n\
    raw1.txt " CR2RES_EMISSION_LINES_TXT_RAW" [1]                       \n\
    raw2.txt " CR2RES_LINES_SELECTION_TXT_RAW" [0 to 1]                 \n\
    The ASCII file raw1.txt must contain two columns:                   \n\
      1st: Wavelengths in increasing order (the unit is corrected by    \n\
               the factor option to obtain nanometers).                 \n\
      2nd: The atmospheric emission.                                    \n\
      The file is in the catalogs/ directory of the CR2RES distribution.\n\
    The optional ASCII files raw2.txt contain a list of wavelength      \n\
      ranges (1 per line) of type: 1632.25,1632.70                      \n\
      The files are in the catalogs/selection/ directory of the CR2RES  \n\
          distribution. They are called XXX.dat, where XXX              \n\
          identifies the setting [L3244, L3340, ...].                   \n\
                                                                        \n\
  Output                                                                \n\
    cr2res_util_genlines.fits "CR2RES_EMISSION_LINES_PROCATG"           \n\
    cr2res_util_genlines_XXX.fits "CR2RES_EMISSION_LINES_PROCATG"       \n\
                                                                        \n\
  Algorithm                                                             \n\
    Parse and load the 2 columns raw1.txt file                          \n\
    Apply the --wl_factor correction                                    \n\
    if (--display) plot it                                              \n\
    Create the CPL table                                                \n\
    Save the table with all lines                                       \n\
    Loop on the selection files                                         \n\
        Only keep the lines that fall in the selection ranges (if any)  \n\
        Save the table with the selected lines                          \n\
                                                                        \n\
  Library functions used                                                \n\
    cr2res_io_save_EMISSION_LINES()                                     \n\
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
                    "Generate spectrum calibration FITS tables",
                    cr2res_util_genlines_description,
                    CR2RES_PIPELINE_AUTHORS,
                    PACKAGE_BUGREPORT,
                    cr2res_get_license(),
                    cr2res_util_genlines_create,
                    cr2res_util_genlines_exec,
                    cr2res_util_genlines_destroy)) {    
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
static int cr2res_util_genlines_create(cpl_plugin * plugin)
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
    p = cpl_parameter_new_value("cr2res_util_genlines.wl_factor", 
            CPL_TYPE_DOUBLE, "The factor used to multiply the wl",
            RECIPE_STRING, 1.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "wl_factor");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res_util_genlines.display", 
            CPL_TYPE_BOOL, "Flag to plot", RECIPE_STRING, FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "display");
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
static int cr2res_util_genlines_exec(cpl_plugin * plugin)
{
    cpl_recipe  *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1;

    return cr2res_util_genlines(recipe->frames, recipe->parameters);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_util_genlines_destroy(cpl_plugin * plugin)
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
static int cr2res_util_genlines(
        cpl_frameset            *   frameset,
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter *   par ;
    cpl_frameset        *   lines_list_frames ;
    cpl_frame           *   lines_list_frame ;
    const char          *   lines_list_fname ;
    cpl_frameset        *   selection_frames ;
    cpl_frame           *   selection_frame ;
    const char          *   selection_fname ;
    char                *   out_file;
    double                  wl_fac ;
    int                     display ;
    cpl_bivector        *   bivec ;
    cpl_bivector        *   bivec_sorted ;
    cpl_bivector        *   bivec_selec ;
    cpl_bivector        *   bivec_selected ;
    double              *   pbivec_x ;
    double              *   pbivec_y ;
    double              *   pbivec_selec_x ;
    double              *   pbivec_selec_y ;
    double              *   pbivec_selected_x ;
    double              *   pbivec_selected_y ;
    int                     nvals ;
    cpl_table           *   tab ;
    char                *   setting ;
    int                     i, j, k ;

    /* Retrieve input parameters */
    par=cpl_parameterlist_find_const(parlist, "cr2res_util_genlines.display");
    display = cpl_parameter_get_bool(par);
    par=cpl_parameterlist_find_const(parlist, "cr2res_util_genlines.wl_factor");
    wl_fac = cpl_parameter_get_double(par);

    /* Identify the RAW and CALIB frames in the input frameset */
    if (cr2res_dfs_set_groups(frameset)) {
        cpl_msg_error(__func__, "Cannot identify RAW and CALIB frames") ;
        return -1 ;
    }

    /* Get Calibration frames */
    selection_frames = cr2res_extract_frameset(frameset,
            CR2RES_LINES_SELECTION_TXT_RAW) ;
        
    /* Get the rawframes */
    lines_list_frames = cr2res_extract_frameset(frameset, 
            CR2RES_EMISSION_LINES_TXT_RAW) ;
    if (lines_list_frames==NULL || 
            cpl_frameset_get_size(lines_list_frames) != 1) {
        cpl_msg_error(__func__, "Please provide 1 and only 1 Lines list file") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        if (selection_frames != NULL) cpl_frameset_delete(selection_frames) ;
        if (lines_list_frames != NULL) cpl_frameset_delete(lines_list_frames) ;
        return -1 ;
    }
    lines_list_frame = cpl_frameset_get_position(lines_list_frames, 0);
    lines_list_fname = cpl_frame_get_filename(lines_list_frame) ;

    cpl_msg_info(__func__, "Process Lines List %s", lines_list_fname) ;
    cpl_msg_indent_more() ;

    /* Load the file */
    if ((bivec=cpl_bivector_read(lines_list_fname))==NULL){
        cpl_msg_error(__func__, "Cannot load the lines in the bivector") ;
        cpl_msg_indent_less() ;
        if (selection_frames != NULL) cpl_frameset_delete(selection_frames) ;
        cpl_frameset_delete(lines_list_frames) ;
        return -1 ;
    }

    /* Use wl_factor */
    cpl_vector_multiply_scalar(cpl_bivector_get_x(bivec), wl_fac) ;

    /* Sort if needed */
    int sort = 1 ;
    if (sort) {
        bivec_sorted = cpl_bivector_duplicate(bivec) ;
        cpl_bivector_sort(bivec_sorted, bivec, CPL_SORT_ASCENDING,
                CPL_SORT_BY_X) ;
        cpl_bivector_delete(bivec) ;
        bivec = bivec_sorted ;
        bivec_sorted = NULL;
    }

    /* Display if requested */
    if (display) {
        cpl_plot_bivector(
            "set grid;set xlabel 'Wavelength (nm)';set ylabel 'Emission';",
            "t 'Catalog lines' w lines", "", bivec);
    }

    /* Allocate the data container */
    nvals = cpl_bivector_get_size(bivec) ;
    tab = cpl_table_new(nvals) ;
    cpl_table_wrap_double(tab, cpl_bivector_get_x_data(bivec),
            CR2RES_COL_WAVELENGTH) ;
    cpl_table_wrap_double(tab, cpl_bivector_get_y_data(bivec),
            CR2RES_COL_EMISSION) ;

    /* Save the table */
    cpl_msg_info(__func__, "Saving the table with %d rows", nvals) ;
    out_file = cpl_sprintf("%s.fits",
            cr2res_get_base_name(cr2res_get_root_name(lines_list_fname)));
    if (cr2res_io_save_EMISSION_LINES(out_file, tab, parlist, frameset,
                RECIPE_STRING, NULL) == -1) {
        cpl_msg_error(__func__, "Cannot write the table") ;
        cpl_table_unwrap(tab, CR2RES_COL_WAVELENGTH) ;
        cpl_table_unwrap(tab, CR2RES_COL_EMISSION) ;
        cpl_table_delete(tab) ;
        cpl_free(out_file);
        if (selection_frames != NULL) cpl_frameset_delete(selection_frames) ;
        cpl_frameset_delete(lines_list_frames) ;
        cpl_bivector_delete(bivec) ;
        cpl_msg_indent_less() ;
        return -1 ;
    }
    cpl_free(out_file);
    cpl_table_unwrap(tab, CR2RES_COL_WAVELENGTH) ;
    cpl_table_unwrap(tab, CR2RES_COL_EMISSION) ;
    cpl_table_delete(tab) ;

    /* Loop on the selection frames */
    for (i=0 ; i<cpl_frameset_get_size(selection_frames) ; i++) {
        /* Get the Current Frame */
        selection_frame = cpl_frameset_get_position(selection_frames, i) ;
        selection_fname = cpl_frame_get_filename(selection_frame) ;
        bivec_selected = NULL ;

        /* Apply the selection */
        bivec_selec = cpl_bivector_read(selection_fname) ;
        if (bivec_selec != NULL) {
            pbivec_x = cpl_bivector_get_x_data(bivec) ;
            pbivec_y = cpl_bivector_get_y_data(bivec) ;
            pbivec_selec_x = cpl_bivector_get_x_data(bivec_selec) ;
            pbivec_selec_y = cpl_bivector_get_y_data(bivec_selec) ;
            /* Count the selected lines */
            nvals = 0 ;
            for (j=0 ; j<cpl_bivector_get_size(bivec) ; j++) {
                /* Loop on the selected ranges */
                for (k=0 ; k<cpl_bivector_get_size(bivec_selec) ; k++) {
                    /* Check if the line is in one of the ranges */
                    if (pbivec_x[j] >= pbivec_selec_x[k] &&
                            pbivec_x[j] <= pbivec_selec_y[k]) {
                        nvals++ ;
                        /* Next line */
                        break ;
                    }
                }
            }

            if (nvals > 0) {
                /* Apply selection */
                bivec_selected = cpl_bivector_new(nvals) ;
                pbivec_selected_x = cpl_bivector_get_x_data(bivec_selected);
                pbivec_selected_y = cpl_bivector_get_y_data(bivec_selected);
                nvals = 0 ;
                for (j=0 ; j<cpl_bivector_get_size(bivec) ; j++) {
                    /* Loop on the selected ranges */
                    for (k=0 ; k<cpl_bivector_get_size(bivec_selec) ; k++) {
                        /* Check if the line is in one of the ranges */
                        if (pbivec_x[j] >= pbivec_selec_x[k] &&
                                pbivec_x[j] <= pbivec_selec_y[k]) {
                            pbivec_selected_x[nvals] = pbivec_x[j] ;
                            pbivec_selected_y[nvals] = pbivec_y[j] ;
                            nvals++ ;
                            /* Next line */
                            break ;
                        }
                    }
                }
            } else {
                bivec_selected = NULL ;
            }
            cpl_bivector_delete(bivec_selec) ;
        } 

        /* Save the selection */
        if (bivec_selected != NULL) {
            /* Allocate the data container */
            nvals = cpl_bivector_get_size(bivec_selected) ;
            tab = cpl_table_new(nvals) ;
            cpl_table_wrap_double(tab, cpl_bivector_get_x_data(bivec_selected),
                    CR2RES_COL_WAVELENGTH) ;
            cpl_table_wrap_double(tab, cpl_bivector_get_y_data(bivec_selected),
                    CR2RES_COL_EMISSION) ;

            /* Save the table */
            cpl_msg_info(__func__, "Saving the table with %d rows", nvals) ;

            /* Format the setting */
            setting = cpl_strdup(
                cr2res_get_base_name(cr2res_get_root_name(selection_fname))) ;
            out_file = cpl_sprintf("%s_%s.fits",
                cr2res_get_base_name(cr2res_get_root_name(lines_list_fname)),
                setting);
            cr2res_format_setting2(setting) ;
            if (cr2res_io_save_EMISSION_LINES(out_file, tab, parlist, frameset,
                        RECIPE_STRING, setting) == -1) {
                cpl_msg_error(__func__, "Cannot write the table") ;
                cpl_table_unwrap(tab, CR2RES_COL_WAVELENGTH) ;
                cpl_table_unwrap(tab, CR2RES_COL_EMISSION) ;
                cpl_table_delete(tab) ;
                cpl_free(out_file);
                if (selection_frames != NULL) 
                    cpl_frameset_delete(selection_frames) ;
                cpl_frameset_delete(lines_list_frames) ;
                cpl_bivector_delete(bivec_selected) ;
                cpl_bivector_delete(bivec) ;
                cpl_free(setting) ;
                cpl_msg_indent_less() ;
                return -1 ;
            }
            cpl_free(setting) ;
            cpl_free(out_file);
            cpl_table_unwrap(tab, CR2RES_COL_WAVELENGTH) ;
            cpl_table_unwrap(tab, CR2RES_COL_EMISSION) ;
            cpl_table_delete(tab) ;
            cpl_bivector_delete(bivec_selected) ;
        }
    }
    cpl_msg_indent_less() ;
    if (selection_frames != NULL) cpl_frameset_delete(selection_frames) ;
    cpl_frameset_delete(lines_list_frames) ;
    cpl_bivector_delete(bivec) ;
    return 0 ;
}


