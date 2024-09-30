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
#include <cpl.h>
#include "hdrl.h"

#include "cr2res_utils.h"
#include "cr2res_dfs.h"
#include "cr2res_io.h"

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

#define RECIPE_STRING "cr2res_util_genstd"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_util_genstd_create(cpl_plugin *);
static int cr2res_util_genstd_exec(cpl_plugin *);
static int cr2res_util_genstd_destroy(cpl_plugin *);
static int cr2res_util_genstd(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                          	Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_util_genstd_description[] = "\
Photospheric flux table generation                                      \n\
                                                                        \n\
  This utility is used to generate the photospheric flux table          \n\
                                                                        \n\
  Inputs                                                                \n\
    raw.txt " CR2RES_PHOTO_FLUX_TXT_RAW " [1 to n]                      \n\
    The specified files are named after the standard star they represent\n\
    (e.g. HIP61007.txt).                                                \n\
    The first line of the file must contain the RA and DEC (hh mm ss).  \n\
    (e.g. # 13 20 35.818        -36 42 44.26).                          \n\
    The rest of the file must contain two columns:                      \n\
    1st: Wavelengths in increasing order                                \n\
    2nd: The atmospheric emission                                       \n\
                                                                        \n\
  Outputs                                                               \n\
    cr2res_util_genstd.fits " CR2RES_PHOTO_FLUX_PROCATG"                \n\
                                                                        \n\
  Algorithm                                                             \n\
                                                                        \n\
  Library functions used                                                \n\
" ;

/*-----------------------------------------------------------------------------
                                Functions code
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
                    "Generate standard star FITS tables",
                    cr2res_util_genstd_description,
                    CR2RES_PIPELINE_AUTHORS,
                    PACKAGE_BUGREPORT,
                    cr2res_get_license(),
                    cr2res_util_genstd_create,
                    cr2res_util_genstd_exec,
                    cr2res_util_genstd_destroy)) {
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
static int cr2res_util_genstd_create(cpl_plugin * plugin)
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
    p = cpl_parameter_new_value("cr2res_util_genstd.display",
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
static int cr2res_util_genstd_exec(cpl_plugin * plugin)
{
    cpl_recipe  *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1;

    return cr2res_util_genstd(recipe->frames, recipe->parameters);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_util_genstd_destroy(cpl_plugin * plugin)
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
  @brief    The FITS file creation occurs here 
  @param    framelist   the frames list
  @return   0 iff everything is ok

  The recipe expects a text file with two columns, and will create a FITS file
  out of it. 
 */
/*----------------------------------------------------------------------------*/
static int cr2res_util_genstd(
        cpl_frameset            *   framelist,
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter *   par ;
    char                line[1024] ;
    cpl_bivector    *   bivec_ref ;
    int                 ra1, ra2, dec1, dec2;
    double              ra3, dec3;
    char                isign;
    cpl_frame       *   cur_frame ;
    const char      *   cur_fname ;
    cpl_table       *   tab ;
    int                 nframes, nvals_ref, display ;
    double          *   pwave_ref ;
    cpl_array       *   array ;
    int                 i, j ;

    /* Needed for sscanf() */
    setlocale(LC_NUMERIC, "C");

    /* Retrieve input parameters */
    par=cpl_parameterlist_find_const(parlist, "cr2res_util_genstd.display");
    display = cpl_parameter_get_bool(par);

    /* Identify the RAW and CALIB frames in the input frameset */
    if (cr2res_dfs_set_groups(framelist)) {
        cpl_msg_error(__func__, "Cannot identify RAW and CALIB frames") ;
        return -1 ;
    }
    nframes = cpl_frameset_get_size(framelist) ;
  
    /* Get the first star */
    cur_frame = cpl_frameset_get_position(framelist, 0) ;
    cur_fname = cpl_frame_get_filename(cur_frame) ;
    if ((bivec_ref = cpl_bivector_read(cur_fname))==NULL){
        cpl_msg_error(__func__, "Cannot load the file in the bivector") ;
        return -1 ;
    }
    pwave_ref = cpl_bivector_get_x_data(bivec_ref) ;
    nvals_ref = cpl_bivector_get_size(bivec_ref) ;

    /* Create the table */
    tab = cpl_table_new(nframes+1) ;
    cpl_table_new_column(tab, CR2RES_COL_STDNAME, CPL_TYPE_STRING) ;
    cpl_table_new_column(tab, CR2RES_COL_RA, CPL_TYPE_DOUBLE) ;
    cpl_table_new_column(tab, CR2RES_COL_DEC, CPL_TYPE_DOUBLE) ;
    cpl_table_new_column_array(tab, CR2RES_COL_PHOTOFLUX, CPL_TYPE_DOUBLE, 
            nvals_ref) ;
    
    /* Write the wavelength */
    cpl_table_set_string(tab, CR2RES_COL_STDNAME, 0, "WAVE") ;
    cpl_table_set_double(tab, CR2RES_COL_RA, 0, -1.0) ;
    cpl_table_set_double(tab, CR2RES_COL_DEC, 0, -1.0) ;
    array = cpl_array_wrap_double(pwave_ref, nvals_ref) ;
    cpl_table_set_array(tab, CR2RES_COL_PHOTOFLUX, 0, array) ;
    cpl_array_unwrap(array) ;

    /* Loop on the input frames */
    for (i = 0; i < nframes; i++) {
        FILE *in;
        cpl_bivector *bivec;
        double ra, dec;
        double *pwave;
        double *pemiss;
        int nvals;

        /* Get the frame */
        cur_frame = cpl_frameset_get_position(framelist, i) ;
        cur_fname = cpl_frame_get_filename(cur_frame) ;

        /* Get the RA / DEC */
        if ((in = fopen(cur_fname, "r")) == NULL) {
            cpl_msg_error(__func__, "Could not open %s", cur_fname) ;
            cpl_table_delete(tab) ;
            cpl_bivector_delete(bivec_ref) ;
            return -1 ;
        }
        if (fgets(line, 1024, in) == NULL) {
            fclose(in) ;
            cpl_table_delete(tab) ;
            cpl_bivector_delete(bivec_ref) ;
            return -1 ;
        }
        if (sscanf(line, "#%d %d %lg %c%d %d %lg ", &ra1, &ra2, &ra3, &isign,
                    &dec1,  &dec2, &dec3) != 7) {
            cpl_msg_error(__func__, "Invalid first line in file %s", cur_fname);
            fclose(in) ;
            cpl_table_delete(tab) ;
            cpl_bivector_delete(bivec_ref) ;
            return -1 ;
        }
        fclose(in) ;
        ra = cr2res_ra_hms2deg(ra1, ra2, ra3) ;
        dec = cr2res_dec_hms2deg(dec1, dec2, dec3) ;
        if (isign == '-') dec *= -1.0 ;
     
        /* Load the file */
        if ((bivec = cpl_bivector_read(cur_fname))==NULL){
            cpl_msg_error(__func__, "Cannot load the file in the bivector") ;
            cpl_bivector_delete(bivec_ref) ;
            cpl_table_delete(tab) ;
            return -1 ;
        }
        pwave = cpl_bivector_get_x_data(bivec) ;
        pemiss = cpl_bivector_get_y_data(bivec) ;
        nvals = cpl_bivector_get_size(bivec) ;

        /* Check the file size */
        if (nvals != nvals_ref) {
            cpl_msg_error(__func__, "Invalid file size: %s", cur_fname) ;
            cpl_bivector_delete(bivec_ref) ;
            cpl_bivector_delete(bivec) ;
            cpl_table_delete(tab) ;
            return -1 ;
        }

        /* Check that the wavelength bins are the same */
        for (j=0 ; j<nvals ; j++) {
            if (pwave[j] != pwave_ref[j]) {
                cpl_msg_error(__func__, "Invalid bins in %s", cur_fname) ;
                cpl_bivector_delete(bivec_ref) ;
                cpl_bivector_delete(bivec) ;
                cpl_table_delete(tab) ;
                return -1 ;
            }
        }
            
        /* Display if requested */
        if (display) {
            cpl_plot_bivector(
                "set grid;set xlabel 'Wavelength (nm)';set ylabel 'Flux (jy)';",
                "t 'Photospheric flux' w lines", "", bivec);
        }
        
        /* Write the star name */
        cpl_table_set_string(tab, CR2RES_COL_STDNAME, i+1,
                cr2res_get_root_name(cr2res_get_base_name(cur_fname))) ;
        
        /* Write the RA/DEC  */
        cpl_table_set_double(tab, CR2RES_COL_RA, i+1, ra) ;
        cpl_table_set_double(tab, CR2RES_COL_DEC, i+1, dec) ;

        /* Write the signal */
        array = cpl_array_wrap_double(pemiss, nvals) ;
        cpl_table_set_array(tab, CR2RES_COL_PHOTOFLUX, i+1, array) ;
        cpl_array_unwrap(array) ;
    
        cpl_bivector_delete(bivec) ;
    }
    cpl_bivector_delete(bivec_ref) ;
    
    /* Save the table */
    cpl_msg_info(__func__, "Save the table") ;

    if (cr2res_io_save_PHOTO_FLUX("cr2res_util_genstd.fits", tab,
                parlist, framelist, RECIPE_STRING) == -1) {
        cpl_msg_error(__func__, "Cannot write the table") ;
        cpl_table_delete(tab) ;
        return -1 ;
    }
    cpl_table_delete(tab) ;

    return 0 ;
}
