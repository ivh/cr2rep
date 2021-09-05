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
#include <math.h>

#include <cpl.h>

#include "cr2res_idp.h"
#include "cr2res_pfits.h"
#include "cr2res_utils.h"
#include "cr2res_dfs.h"

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_idp_copy(
    cpl_table       *   out,
    const cpl_table *   in,
    cpl_size            out_start_idx,
    int                 det_nr,
    int                 order,
    int                 tracenb) ;

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_idp     IDP related functions
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Save an IDP file
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    rawframes   The recipe used input frames
  @param    parlist     The recipe input parameters
  @param    tables      The tables to save (1 per detector)
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_idp_save(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   rawframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const char              *   recipe)
{
    cpl_table           *   idp_tab ;
    char                *   idp_filename ;
    cpl_frame           *   out_frame ;
    const cpl_frame     *   ref_frame ;
    const char          *   ref_fname ;
    cpl_propertylist    *   pri_head ;
    cpl_propertylist    *   ext_head ;
    double                  dit, exptime, texptime, mjd_start, mjd_end,
                            wmin, wmax ;
	const char			*	progid ;
    int                     err, i, nexp, nraw, obid, nrows ;

    cpl_msg_info(__func__, "Create IDPs for %s", filename) ;
    
    /* Output file name */
    idp_filename = cpl_sprintf("idp_%s", filename) ;

    /* Create the big table */
    idp_tab = cr2res_idp_create_table(tables) ;

    /* Set the units */
    cpl_table_set_column_unit(idp_tab, CR2RES_COL_WAVELENGTH, "nm") ;
    cpl_table_set_column_unit(idp_tab, CR2RES_COL_SPECTRUM, "ADU/sec") ;
    cpl_table_set_column_unit(idp_tab, CR2RES_COL_ERROR, "") ;
    cpl_table_set_column_unit(idp_tab, CR2RES_COL_QUALITY, "") ;
    cpl_table_set_column_unit(idp_tab, CR2RES_COL_ORDER, "Nb Order") ;
    cpl_table_set_column_unit(idp_tab, CR2RES_COL_TRACENB, "Nb Trace") ;
    cpl_table_set_column_unit(idp_tab, CR2RES_COL_DETECTOR, "Nb Detector") ;
    cpl_table_set_column_unit(idp_tab, CR2RES_COL_XPOS, "pixels") ;

    /* Get wmin / wmax */
    nrows = cpl_table_get_nrow(idp_tab) ;
    wmin = cpl_table_get(idp_tab, CR2RES_COL_WAVELENGTH, 0, NULL) ;
    wmax = cpl_table_get(idp_tab, CR2RES_COL_WAVELENGTH, nrows-1, NULL) ;

	/* Prepare frame */
	out_frame = cpl_frame_new();
	cpl_frame_set_filename(out_frame, idp_filename);
	cpl_frame_set_tag(out_frame, CR2RES_OBS_NODDING_IDP_PROCATG);
	cpl_frame_set_type(out_frame, CPL_FRAME_TYPE_ANY);
	cpl_frame_set_group(out_frame, CPL_FRAME_GROUP_PRODUCT);
	cpl_frame_set_level(out_frame, CPL_FRAME_LEVEL_FINAL);

    /* Create the Primary header */
    ref_frame = cpl_frameset_get_position_const(rawframes, 0) ;
    ref_fname = cpl_frame_get_filename(ref_frame) ;

    pri_head = cpl_propertylist_load(ref_fname, 0); 
    cpl_propertylist_append_string(pri_head, CPL_DFS_PRO_CATG,
            CR2RES_OBS_NODDING_IDP_PROCATG) ;
    cpl_propertylist_append_string(pri_head, CPL_DFS_PRO_TYPE,
            CR2RES_EXTRACT_1D_IDP_PROTYPE) ;
    cpl_dfs_setup_product_header(pri_head, out_frame, allframes,
            parlist, recipe, VERSION, "PRO-1.16", ref_frame);

    /* Add Keywords to primary header */
    // ESO INS FILT1 NAME      --->    FILTER
    if (cpl_propertylist_has(pri_head, "ESO INS FILT1 NAME")) {
        cpl_propertylist_update_string(pri_head, "FILTER",
                cpl_propertylist_get_string(pri_head, "ESO INS FILT1 ID")) ;
        cpl_propertylist_set_comment(pri_head, "FILTER",
                "");
    }

    /* Collect data from main header */
    dit = cr2res_pfits_get_dit(pri_head) * cr2res_pfits_get_ndit(pri_head) ; 
    cpl_propertylist_update_double(pri_head, "DIT", dit);

    nexp = cr2res_pfits_get_nexp(pri_head) ;
    exptime = dit * nexp ;
    cpl_propertylist_update_double(pri_head, "EXPTIME", exptime);
    cpl_propertylist_set_comment(pri_head, "EXPTIME", 
            "[s] Total integration time per pixel");

    texptime = exptime ;
    cpl_propertylist_update_double(pri_head, "TEXPTIME", texptime);
    cpl_propertylist_set_comment(pri_head, "TEXPTIME", 
            "[s] Total integration time of exposures");

    cr2res_idp_compute_mjd(rawframes, &mjd_start, &mjd_end) ;
    if (mjd_start > 0.0) {
        cpl_propertylist_update_double(pri_head, "MJD-OBS", mjd_start) ;
        cpl_propertylist_set_comment(pri_head, "MJD-OBS", 
                "[d] Start of observations (days)");
    }
    if (mjd_end > 0.0) {
        cpl_propertylist_update_double(pri_head, "MJD-END", mjd_end) ;
        cpl_propertylist_set_comment(pri_head, "MJD-END", 
                "[d] End of observations (days)");
    }

    progid = cr2res_pfits_get_progid(pri_head) ;
    cpl_propertylist_update_string(pri_head, "PROG_ID", progid) ;
    cpl_propertylist_set_comment(pri_head, "PROG_ID", 
            "ESO programme identification");

	obid = cr2res_pfits_get_obs_id(pri_head) ;
    cpl_propertylist_update_int(pri_head, "OBID1", obid);
    cpl_propertylist_set_comment(pri_head, "OBID1", "Observation block ID");

    cpl_propertylist_update_bool(pri_head, "M_EPOCH", 1) ;
    cpl_propertylist_set_comment(pri_head, "M_EPOCH",
            "TRUE if resulting from multiple epochs") ;

// SINGLEEXP ?

    nraw = cpl_frameset_get_size(rawframes) ;
    cpl_propertylist_update_int(pri_head, "NCOMBINE", nraw);

    cpl_propertylist_update_string(pri_head, "OBSTECH", "NODDING") ;
	cpl_propertylist_set_comment(pri_head, "OBSTECH",
			"Technique of observation") ;

	cpl_propertylist_update_string(pri_head, "FLUXCAL", "ABSOLUTE") ;
	cpl_propertylist_set_comment(pri_head, "FLUXCAL", 
            "Type of flux calibration");

    cpl_propertylist_update_string(pri_head, "PROCSOFT",
            PACKAGE "/" PACKAGE_VERSION) ;
    cpl_propertylist_set_comment(pri_head, "PROCSOFT", "ESO pipeline version");

    cpl_propertylist_update_string(pri_head, "REFERENC", "") ;
    cpl_propertylist_set_comment(pri_head, "REFERENC", "Reference publication");

    cpl_propertylist_update_string(pri_head, "PRODCATG", "SCIENCE.SPECTRUM") ;
    cpl_propertylist_set_comment(pri_head, "PRODCATG", "Data product category");

// DATE-OBS ?
 
    /* Remove the ASSON keywords */
    //cpl_propertylist_erase_regexp(pri_head, "ASSO*", 0);

	/* Save the main header */
	cpl_propertylist_save(pri_head, idp_filename, CPL_IO_CREATE);

    /* Create the first extension header */
    ext_head = cpl_propertylist_new() ;

    /* Add Keywords to extension header */
    cpl_propertylist_update_string(ext_head, "EXTNAME", "IDP_SPECTRUM") ;

    cpl_propertylist_update_string(ext_head, "VOPUB", "ESO/SAF") ;
    cpl_propertylist_set_comment(ext_head, "VOPUB", "VO Publishing Authority") ;
    cpl_propertylist_update_string(ext_head, "VOCLASS", "SPECTRUM V1.0") ;
    cpl_propertylist_set_comment(ext_head, "VOCLASS",
            "Data Model name and version") ;
    cpl_propertylist_update_int(ext_head, "NELEM", nrows) ;
    cpl_propertylist_set_comment(ext_head, "NELEM", "Length of the data array");
    cpl_propertylist_update_double(ext_head, "APERTURE", 0.0000555) ;
    cpl_propertylist_set_comment(ext_head, "APERTURE", 
            "Slit width in deg") ;

    if (mjd_end > 0 && mjd_start > 0) {
        cpl_propertylist_update_double(ext_head, "TELAPSE", mjd_end-mjd_start) ;
        cpl_propertylist_set_comment(ext_head, "TELAPSE", 
                "Total elapsed time in seconds [s]") ;

        cpl_propertylist_update_double(ext_head, "TMID",
                (mjd_end+mjd_start)/2.0) ;
        cpl_propertylist_set_comment(ext_head, "TMID", 
                "Exposure midpoint [MJD]") ;
    }

    cpl_propertylist_update_double(ext_head, "SPEC_VAL", (wmax+wmin)/2.0) ;
    cpl_propertylist_set_comment(ext_head, "SPEC_VAL", 
            "Characteristic spectral coordinate value [nm]") ;

    cpl_propertylist_update_double(ext_head, "SPEC_BW", wmax-wmin) ;
    cpl_propertylist_set_comment(ext_head, "SPEC_BW", 
            "Width of the spectrum [nm]") ;

    /* Remove keywords */
    cpl_propertylist_erase(ext_head, "CRDER3");
    cpl_propertylist_erase(ext_head, "CUNIT3");
    cpl_propertylist_erase(ext_head, "BUNIT");
    cpl_propertylist_erase(ext_head, "CTYPE1");
    cpl_propertylist_erase(ext_head, "CUNIT2");
    cpl_propertylist_erase(ext_head, "CSYER1");
    cpl_propertylist_erase(ext_head, "CSYER2");
    cpl_propertylist_erase(ext_head, "CDELT1");
    cpl_propertylist_erase(ext_head, "CRPIX1");
    cpl_propertylist_erase(ext_head, "CRVAL1");
    cpl_propertylist_erase(ext_head, "CUNIT1");

    /* Save the table */
    cpl_table_save(idp_tab, NULL, ext_head, idp_filename, CPL_IO_EXTEND) ;

    /* Insert it in the Frameset */
    cpl_frameset_insert(allframes, out_frame);

    /* Clean */
    cpl_free(idp_filename) ;
    cpl_table_delete(idp_tab) ;
    cpl_propertylist_delete(pri_head) ;
    cpl_propertylist_delete(ext_head) ;

    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief     Convert EXTRACT_1D tables into IDP table
  @param      tables    CR2RES_NB_DETECTORS EXTRACT_1D tables
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_idp_create_table(
        cpl_table               **  tables)
{
    cpl_table           *   idp_tab ;
    cpl_array           *   col_names ;
    const char          *   col_name ;
    char                *   col_type ;
    cpl_propertylist    *   sort_list ;
    int                     order, trace_nb ;
    double                  cur_val, pre_val ;
    cpl_size                i, j, ntot, ncols, nb_selected ;

    /* Check Inputs */
    if (tables == NULL) return NULL ;

    /* Count the Rows */
    ntot = 0 ;
    for (i=0 ; i<CR2RES_NB_DETECTORS ; i++) {
        if (tables[i] != NULL) {
            col_names = cpl_table_get_column_names(tables[i]);
            ncols = cpl_table_get_ncol(tables[i]) ;
            for (j=0 ; j<ncols ; j++) {
                col_name = cpl_array_get_string(col_names, j);
                col_type = cr2res_dfs_SPEC_colname_parse(col_name, 
                        &order, &trace_nb) ;
                if (col_type != NULL && 
                        !strcmp(col_type, CR2RES_COL_SPEC_SUFFIX)) {
                    /* Handle this extracted spectrum */
                    ntot += cpl_table_get_nrow(tables[i]) ;
                }
                if (col_type != NULL) cpl_free(col_type) ;
            }
            cpl_array_delete(col_names) ;
        }
    }
 
    /* Create the table */
    idp_tab = cpl_table_new(ntot) ;
    cpl_table_new_column(idp_tab, CR2RES_COL_WAVELENGTH, CPL_TYPE_DOUBLE);
    cpl_table_new_column(idp_tab, CR2RES_COL_SPECTRUM, CPL_TYPE_DOUBLE);
    cpl_table_new_column(idp_tab, CR2RES_COL_ERROR, CPL_TYPE_DOUBLE);
    cpl_table_new_column(idp_tab, CR2RES_COL_QUALITY, CPL_TYPE_DOUBLE);
    cpl_table_new_column(idp_tab, CR2RES_COL_ORDER, CPL_TYPE_INT);
    cpl_table_new_column(idp_tab, CR2RES_COL_TRACENB, CPL_TYPE_INT);
    cpl_table_new_column(idp_tab, CR2RES_COL_DETECTOR, CPL_TYPE_INT);
    cpl_table_new_column(idp_tab, CR2RES_COL_XPOS, CPL_TYPE_INT);

    /* Fill the table */
    ntot = 0 ;
    for (i=0 ; i<CR2RES_NB_DETECTORS ; i++) {
        if (tables[i] != NULL) {
            col_names = cpl_table_get_column_names(tables[i]);
            ncols = cpl_table_get_ncol(tables[i]) ;
            for (j=0 ; j<ncols ; j++) {
                col_name = cpl_array_get_string(col_names, j);
                col_type = cr2res_dfs_SPEC_colname_parse(col_name, 
                        &order, &trace_nb) ;
                if (col_type != NULL && 
                        !strcmp(col_type, CR2RES_COL_SPEC_SUFFIX)) {
                    /* Handle this extracted spectrum */
                    cr2res_idp_copy(idp_tab, tables[i], ntot, i+1, order, 
                            trace_nb) ;
                    ntot += cpl_table_get_nrow(tables[i]) ;
                }
                if (col_type != NULL) cpl_free(col_type) ;
            }
            cpl_array_delete(col_names) ;
        }
    }

    /* Sort by the wavelengths */
    sort_list = cpl_propertylist_new() ;
    cpl_propertylist_append_bool(sort_list, CR2RES_COL_WAVELENGTH, 0) ;
    cpl_propertylist_append_bool(sort_list, CR2RES_COL_XPOS, 1) ;
    cpl_table_sort(idp_tab, sort_list) ;
    cpl_propertylist_delete(sort_list) ;

    /* Identify close values */
    cpl_table_unselect_all(idp_tab) ;
    for (i=1 ; i<ntot ; i++) {
        pre_val = cpl_table_get_double(idp_tab,CR2RES_COL_WAVELENGTH,i-1,NULL);
        cur_val = cpl_table_get_double(idp_tab,CR2RES_COL_WAVELENGTH,i, NULL) ;
        if (fabs(pre_val-cur_val) < 1e-3) {
            cpl_table_select_row(idp_tab, i) ;
        }
    }

    nb_selected = cpl_table_count_selected(idp_tab) ;
    if (nb_selected > 0) {
        cpl_msg_info(__func__, 
                "IDP Creation - Double defined WLs : %"CPL_SIZE_FORMAT, 
                nb_selected) ;
    }
    cpl_table_erase_selected(idp_tab) ;

    return idp_tab ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Compute mjd start and end
  @param    fset        Raw frames
  @return
*/
/*----------------------------------------------------------------------------*/
int cr2res_idp_compute_mjd(
        cpl_frameset        *   fset,
        double              *   mjd_start,
        double              *   mjd_end)
{
    cpl_frameset        *   rawframes ;
    cpl_frame           *   cur_frame ;
    const char          *   cur_fname ;
    cpl_propertylist    *   plist ;
    double                  mjd_end_cur, mjd_end_max, mjd_obs, dit,
                            mjd_start_cur, mjd_start_min ;
    int                     i, nframes, ndit ;

    /* Initialise */
    mjd_start_cur = mjd_end_cur = -2.0 ;
    mjd_end_max = -1.0 ;
    mjd_start_min = 1e10 ;
    nframes = cpl_frameset_get_size(fset) ;

    /* Loop */
    for (i=0 ; i<nframes ; i++) {
        cur_frame = cpl_frameset_get_position(fset, i) ;
        cur_fname = cpl_frame_get_filename(cur_frame);

        /* Get header infos */
        plist = cpl_propertylist_load(cur_fname, 0) ;
        dit = cr2res_pfits_get_dit(plist) ;
        ndit = cr2res_pfits_get_ndit(plist) ;
        mjd_obs = cr2res_pfits_get_mjd_obs(plist) ;
        if (cpl_error_get_code() == CPL_ERROR_NONE) {
            /* Compute current value */
            mjd_start_cur = mjd_obs ;
            mjd_end_cur = mjd_obs + (dit*ndit/(24*60*60)) ;
        } else {
            mjd_end_cur = -2.0 ;
            mjd_start_cur = 1e11 ;
            cpl_error_reset() ;
        }
        /* Update max */
        if (mjd_end_cur > mjd_end_max)  mjd_end_max = mjd_end_cur ;
        if (mjd_start_cur < mjd_start_min) {
            mjd_start_min = mjd_start_cur ;
        }
        cpl_propertylist_delete(plist) ;
    }

    if (mjd_start_min < 1e9) {
        *mjd_start = mjd_start_min ;
    } else {
        *mjd_start = -1.0 ;
    }
    if (mjd_end_max > 0.0)      *mjd_end = mjd_end_max ;
    else                        *mjd_end = -1.0 ;
    return 0 ;
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Copy a spectrum (spec, err, wave) in a bigger table 
                - with additional entries
  @param    out             Destination table
  @param    in              Source table
  @param    out_start_idx   Start index in the destination table
  @param    det_nr          Detector nb
  @param    order           Order nb
  @param    tracenb         Trace Number
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
static int cr2res_idp_copy(
    cpl_table       *   out,
    const cpl_table *   in,
    cpl_size            out_start_idx,
    int                 det_nr,
    int                 order,
    int                 tracenb) 
{
    char            *   spec_name ;
    char            *   wave_name ;
    char            *   err_name ;
    const double    *   pspec ;
    const double    *   pwave ;
    const double    *   perr ;
    cpl_size            in_size, out_size, i ;

    /* Check entries */
    if (out == NULL  || in == NULL) return -1 ;
    in_size = cpl_table_get_nrow(in) ;
    out_size = cpl_table_get_nrow(out) ;

    /* Get the Spectrum */
    spec_name = cr2res_dfs_SPEC_colname(order, tracenb) ;
    if ((pspec = cpl_table_get_data_double_const(in, spec_name)) == NULL) {
        cpl_msg_error(__func__, "Cannot find the spectrum") ;
        cpl_free(spec_name) ;
        return -1 ;
    }
    cpl_free(spec_name) ;

    /* Get the Wavelength */
    wave_name = cr2res_dfs_WAVELENGTH_colname(order, tracenb) ;
    if ((pwave = cpl_table_get_data_double_const(in, wave_name)) == NULL) {
        cpl_msg_error(__func__, "Cannot find the wavelength") ;
        cpl_free(wave_name) ;
        return -1 ;
    }
    cpl_free(wave_name) ;

    /* Get the Spectrum Error */
    err_name = cr2res_dfs_SPEC_ERR_colname(order, tracenb) ;
    if ((perr = cpl_table_get_data_double_const(in, err_name)) == NULL) {
        cpl_msg_error(__func__, "Cannot find the spectrum error") ;
        cpl_free(err_name) ;
        return -1 ;
    }
    cpl_free(err_name) ;
    
    /* copy line by line */
    for (i=0 ; i<in_size ; i++) {
        if (out_start_idx + i <out_size) {
            if (!isnan(pwave[i])) 
                cpl_table_set_double(out, CR2RES_COL_WAVELENGTH, 
                        out_start_idx + i, pwave[i]) ;
            if (!isnan(pspec[i])) 
                cpl_table_set_double(out, CR2RES_COL_SPECTRUM, 
                        out_start_idx + i, pspec[i]) ;
            if (!isnan(perr[i])) 
                cpl_table_set_double(out, CR2RES_COL_ERROR, 
                        out_start_idx + i, perr[i]) ;
            cpl_table_set_double(out, CR2RES_COL_QUALITY, out_start_idx + i, 
                    0);
            cpl_table_set_int(out, CR2RES_COL_XPOS, out_start_idx + i, 
                    i+1) ;
            cpl_table_set_int(out, CR2RES_COL_DETECTOR, out_start_idx + i, 
                    det_nr) ;
            cpl_table_set_int(out, CR2RES_COL_ORDER, out_start_idx + i, 
                    order) ;
            cpl_table_set_int(out, CR2RES_COL_TRACENB, out_start_idx + i, 
                    tracenb) ;
        }
    }
    return 0 ;
}

