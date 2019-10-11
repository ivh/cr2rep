/* $Id: crires_util_genstd.c,v 1.12 2011-11-24 08:27:46 yjung Exp $
 *
 * This file is part of the CRIRES Pipeline
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
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * $Author: yjung $
 * $Date: 2011-11-24 08:27:46 $
 * $Revision: 1.12 $
 * $Name: not supported by cvs2svn $
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/*-----------------------------------------------------------------------------
                                Includes
 -----------------------------------------------------------------------------*/

#include <locale.h>
#include "crires_recipe.h"

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

#define RECIPE_STRING "crires_util_genstd"

/*-----------------------------------------------------------------------------
                            Functions prototypes
 -----------------------------------------------------------------------------*/

static int crires_util_genstd_save(cpl_table *, const cpl_parameterlist *, 
        cpl_frameset *);

static char crires_util_genstd_description[] =
"This recipe is used to generate standard star photospheric flux tables.\n"
"The sof consists of file names tagged with "CRIRES_UTIL_GENSTD_RAW".\n"
"The specified files are named after the standard star they represent\n"
"(e.g. HIP61007.txt).\n"
"The first line of the file must contain the RA and DEC (hh mm ss).\n"
"(e.g. # 13 20 35.818        -36 42 44.26).\n"
"The rest of the file must contain two columns:\n"
"1st: Wavelengths in increasing order (the unit is corrected by\n"
"     the factor option to obtain nanometers).\n"
"2nd: The atmospheric emission.\n"
"The file is generated using the ASCII files in the catalogs/stdstar\n"
"directory of the CRIRES source-code distribution."
"\n"
"This recipe produces 1 file for each input file:\n"
"First product:     the table with the photospheric flux of the std.\n"
"                   (PRO TYPE = "CRIRES_PROTYPE_PHO_FLUX")\n" ;

CRIRES_RECIPE_DEFINE(crires_util_genstd,
        CRIRES_PARAM_PLOT,
        "Generate standard star FITS tables",
        crires_util_genstd_description) ;

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static struct {
    /* Inputs */
    int             display ;
    /* Outputs */
} crires_util_genstd_config ;

/*-----------------------------------------------------------------------------
                                Functions code
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @brief    The FITS file creation occurs here 
  @param    framelist   the frames list
  @return   0 iff everything is ok

  The recipe expects a text file with two columns, and will create a FITS file
  out of it. 
 */
/*----------------------------------------------------------------------------*/
static int crires_util_genstd(
        cpl_frameset            *   framelist,
        const cpl_parameterlist *   parlist)
{
    FILE            *   in ;
    char                line[1024] ;
    cpl_bivector    *   bivec_ref ;
    cpl_bivector    *   bivec ;
    int                 ra1, ra2, dec1, dec2;
    double              ra3, dec3;
    double              ra, dec ;
    char                isign;
    cpl_frame       *   cur_frame ;
    const char      *   cur_fname ;
    cpl_table       *   tab ;
    int                 nvals, nframes, nvals_ref ;
    double          *   pwave ;
    double          *   pemiss ;
    double          *   pwave_ref ;
    cpl_array       *   array ;
    int                 i, j ;

    /* Needed for sscanf() */
    setlocale(LC_NUMERIC, "C");

    /* Retrieve input parameters */
    crires_util_genstd_config.display = crires_parameterlist_get_bool(
            parlist, RECIPE_STRING, CRIRES_PARAM_PLOT) ;
 
    /* Identify the RAW and CALIB frames in the input frameset */
    if (crires_dfs_set_groups(framelist, NULL)) {
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
    cpl_table_new_column(tab, CRIRES_COL_STDNAME, CPL_TYPE_STRING) ;
    cpl_table_new_column(tab, CRIRES_COL_RA, CPL_TYPE_DOUBLE) ;
    cpl_table_new_column(tab, CRIRES_COL_DEC, CPL_TYPE_DOUBLE) ;
    cpl_table_new_column_array(tab, CRIRES_COL_PHOTOFLUX, CPL_TYPE_DOUBLE, 
            nvals_ref) ;
    
    /* Write the wavelength */
    cpl_table_set_string(tab, CRIRES_COL_STDNAME, 0, "WAVE") ;
    cpl_table_set_double(tab, CRIRES_COL_RA, 0, -1.0) ;
    cpl_table_set_double(tab, CRIRES_COL_DEC, 0, -1.0) ;
    array = cpl_array_wrap_double(pwave_ref, nvals_ref) ;
    cpl_table_set_array(tab, CRIRES_COL_PHOTOFLUX, 0, array) ;
    cpl_array_unwrap(array) ;

    /* Loop on the input frames */
    for (i=0 ; i<nframes ; i++) {
        
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
        ra = crires_ra_hms2deg(ra1, ra2, ra3) ;
        dec = crires_dec_hms2deg(dec1, dec2, dec3) ;
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
        if (crires_util_genstd_config.display) {
            cpl_plot_bivector(
                "set grid;set xlabel 'Wavelength (nm)';set ylabel 'Flux (jy)';",
                "t 'Photospheric flux' w lines", "", bivec);
        }
        
        /* Write the star name */
        cpl_table_set_string(tab, CRIRES_COL_STDNAME, i+1,
                crires_get_root_name(crires_get_base_name(cur_fname))) ;
        
        /* Write the RA/DEC  */
        cpl_table_set_double(tab, CRIRES_COL_RA, i+1, ra) ;
        cpl_table_set_double(tab, CRIRES_COL_DEC, i+1, dec) ;

        /* Write the signal */
        array = cpl_array_wrap_double(pemiss, nvals) ;
        cpl_table_set_array(tab, CRIRES_COL_PHOTOFLUX, i+1, array) ;
        cpl_array_unwrap(array) ;
    
        cpl_bivector_delete(bivec) ;
    }
    cpl_bivector_delete(bivec_ref) ;
    
    /* Save the table */
    cpl_msg_info(__func__, "Save the table") ;
    if (crires_util_genstd_save(tab, parlist, framelist) == -1) {
        cpl_msg_error(__func__, "Cannot write the table") ;
        cpl_table_delete(tab) ;
        return -1 ;
    }
    cpl_table_delete(tab) ;

    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save the product of the recipe
  @param    out_table   the table 
  @param    parlist     the input list of parameters
  @param    set         the input frame set
  @return   0 if everything is ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
static int crires_util_genstd_save(
        cpl_table               *   out_table,
        const cpl_parameterlist *   parlist,
        cpl_frameset            *   set)
{
    cpl_propertylist    *   plist ;

    plist = cpl_propertylist_new();
    cpl_propertylist_append_string(plist, "INSTRUME", "CRIRES") ;
    cpl_propertylist_append_string(plist, CPL_DFS_PRO_CATG,
            CRIRES_CALPRO_STD_PHOTOFLUX) ;
    cpl_propertylist_append_string(plist, CPL_DFS_PRO_TYPE,
            CRIRES_PROTYPE_PHO_FLUX) ;

    if (cpl_dfs_save_table(set, NULL, parlist, set, NULL, out_table,
                NULL, "crires_util_genstd", plist, NULL, 
                PACKAGE "/" PACKAGE_VERSION,
                "crires_util_genstd.fits") != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot save the table") ;
        return -1 ;
    }
    cpl_propertylist_delete(plist) ;

    /* Return */
    return 0 ;
}
