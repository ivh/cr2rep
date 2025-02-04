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
#include "cr2res_io.h"
#include "cr2res_qc.h"

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_idp_copy_spec(
    cpl_table       *   out,
    const cpl_table *   in,
    cpl_size            out_start_idx,
    int                 det_nr,
    int                 order,
    int                 tracenb,
    const char      *   setting) ;
static int cr2res_idp_copy_pol(
    cpl_table       *   out,
    const cpl_table *   in,
    cpl_size            out_start_idx,
    int                 det_nr,
    int                 order,
    const char          *   setting) ;

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_idp     IDP related functions
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Save an IDP file
  @param    filename        The FITS file name
  @param    allframes       The recipe input frames
  @param    rawframes       The recipe used input frames
  @param    parlist         The recipe input parameters
  @param    tables          The tables to save (1 per detector)
  @param    main_qc_plist   The main header QCs
  @param    ext_plist       CR2RES_NB_DETECTORS extension headers
  @param    recipe          The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_idp_save(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   rawframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        cpl_propertylist        *   main_qc_plist,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    cpl_table           *   idp_tab ;
    const cpl_array     *   wlen_arr;
    cpl_array           *   tmp_arr;
    char                *   idp_filename ;
    cpl_frame           *   out_frame ;
    const cpl_frame     *   ref_frame ;
    const char          *   ref_fname ;
    cpl_propertylist    *   pri_head ;
    cpl_propertylist    *   ext_head ;
    double                  dit, exptime, texptime, mjd_start, mjd_end,
                            wmin, wmax, resol, spec_bin ;
	const char			*	progid ;
	const char			*	slitname ;
	const char			*	setting ;
	const char			*	poltype ;
    int                     err, i, ndit, nframes, nraw, obid, nrows, ord ;
    char                *   keyname;
    char                *   tmp_string;

    poltype = NULL;

    cpl_msg_info(__func__, "Create IDPs for %s", filename) ;
    /* Output file name */
    idp_filename = cpl_sprintf("idp_%s", filename) ;

    /* Create the Primary header */
    ref_frame = cpl_frameset_get_position_const(rawframes, 0) ;
    ref_fname = cpl_frame_get_filename(ref_frame) ;
    pri_head = cpl_propertylist_load(ref_fname, 0); 


    cr2res_qc_dup_mtrlgy_key(rawframes, pri_head);
    
    /* Create the first EXTENSION header */
    ext_head = cpl_propertylist_new() ;

    /* Read spectral setting because we need it to make the table*/
    setting = cpl_propertylist_get_string(pri_head, "ESO INS WLEN ID");

    /* Create the big table */
    idp_tab = cr2res_idp_create_table(tables, recipe, setting) ;

    /* Get wmin / wmax */
    wlen_arr = cpl_table_get_array(idp_tab, CR2RES_IDP_COL_WAVE, 0);
    nrows = cpl_array_get_size(wlen_arr);
    wmin = cpl_array_get_double(wlen_arr, 0, &err) ;
    wmax = cpl_array_get_double(wlen_arr, nrows-1, &err) ;

	/* Prepare frame */
	out_frame = cpl_frame_new();
	cpl_frame_set_filename(out_frame, idp_filename);
	cpl_frame_set_type(out_frame, CPL_FRAME_TYPE_ANY);
	cpl_frame_set_group(out_frame, CPL_FRAME_GROUP_PRODUCT);
	cpl_frame_set_level(out_frame, CPL_FRAME_LEVEL_FINAL);

    // Set the PRO.CATG
    cpl_frame_set_tag(out_frame, procatg);
    cpl_propertylist_append_string(pri_head, CPL_DFS_PRO_CATG, procatg);

    // TODO: split this type by recipe
    cpl_propertylist_append_string(pri_head, CR2RES_HEADER_DRS_TYPE,
            CR2RES_EXTRACT_1D_IDP_DRSTYPE) ;

    if (!strcmp(recipe, "cr2res_obs_nodding")) {
        cpl_propertylist_update_string(pri_head, "OBSTECH", "NODDING") ;
    }
    else if (!strcmp(recipe, "cr2res_obs_staring")) {
        cpl_propertylist_update_string(pri_head, "OBSTECH", "STARING") ;
    }
    else if (!strcmp(recipe, "cr2res_obs_2d")) {
        cpl_propertylist_update_string(pri_head, "OBSTECH", "GENERIC_OFFSET") ;
    }
    else if (!strcmp(recipe, "cr2res_obs_pol")) {
        cpl_propertylist_update_string(pri_head, "OBSTECH", "POLARIMETRY") ;
    }
    cpl_propertylist_set_comment(pri_head, "OBSTECH",
			"Technique of observation") ;

    /* RADECSYS renamed to RADESYS */
    if (cpl_propertylist_has(pri_head, "RADECSYS")) {
        if (!cpl_propertylist_has(pri_head, "RADESYS")) {
            const cpl_property *_property =
                cpl_propertylist_get_property_const(pri_head, "RADECSYS");
            cpl_property *property = cpl_property_duplicate(_property);
            cpl_property_set_name(property, "RADESYS");
            cpl_propertylist_append_property(pri_head, property);
            cpl_property_delete(property);
        }
        cpl_propertylist_erase(pri_head, "RADECSYS");
    }

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
    dit = cr2res_pfits_get_dit(pri_head);
    ndit = cr2res_pfits_get_ndit(pri_head) ; 
    
    nframes = cpl_frameset_get_size(rawframes);
    nraw = 0;

    cpl_errorstate tempes = cpl_errorstate_get();

    for (i = 1; i <= nframes; i++) {
        const char *fname;
        keyname = cpl_sprintf("ESO PRO REC1 RAW%d NAME", i);
        fname = cpl_propertylist_get_string(pri_head, keyname);
        cpl_free(keyname);
        if (fname != NULL){
            keyname = cpl_sprintf("PROV%d", i);
            cpl_propertylist_update_string(pri_head, keyname, fname);
            cpl_free(keyname);
            nraw++;
        }
    }

    cpl_errorstate_set(tempes);

    /* Add QC from the first Extension  ext_plist[0] */
    if (ext_plist[0] != NULL)
        cpl_propertylist_copy_property_regexp(ext_head, ext_plist[0], "QC*", 0);

    cpl_propertylist_update_int(pri_head, "NCOMBINE", nraw);
    exptime = dit * ndit * nraw;
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

    if (mjd_end > 0 && mjd_start > 0) {
        cpl_propertylist_update_double(ext_head, "TELAPSE",
                (mjd_end-mjd_start)*24*3600) ;
        cpl_propertylist_set_comment(ext_head, "TELAPSE", 
                "Total elapsed time in seconds [s]") ;

        cpl_propertylist_update_double(ext_head, "TMID",
                (mjd_end+mjd_start)/2.0) ;
        cpl_propertylist_set_comment(ext_head, "TMID", 
                "Exposure midpoint [MJD]") ;
    }

    progid = cr2res_pfits_get_progid(pri_head) ;
    cpl_propertylist_update_string(pri_head, "PROG_ID", progid) ;
    cpl_propertylist_set_comment(pri_head, "PROG_ID", 
            "ESO programme identification");

	obid = cr2res_pfits_get_obs_id(pri_head) ;
    cpl_propertylist_update_int(pri_head, "OBID1", obid);
    cpl_propertylist_set_comment(pri_head, "OBID1", "Observation block ID");

    if (!strcmp(recipe, "cr2res_obs_2d"))
        cpl_propertylist_update_bool(pri_head, "EXT_OBJ", 1) ;
    else
        cpl_propertylist_update_bool(pri_head, "EXT_OBJ", 0) ;
    cpl_propertylist_set_comment(pri_head, "EXT_OBJ",
			"True for extended objects, cr2res_obs_2d was used") ;

	cpl_propertylist_update_string(pri_head, "FLUXCAL", "UNCALIBRATED") ;
	cpl_propertylist_set_comment(pri_head, "FLUXCAL", 
            "Type of flux calibration");

	cpl_propertylist_update_int(pri_head, "FLUXERR", -1) ;
	cpl_propertylist_set_comment(pri_head, "FLUXERR", 
            "Fractional uncertainty [%%] on the flux scale");

	cpl_propertylist_update_bool(pri_head, "TOT_FLUX", 0) ;
	cpl_propertylist_set_comment(pri_head, "TOT_FLUX", 
            "True if flux data represent the total source flux.");

    /* Find blaze to decide if continuum normalized or not */
    if ( cpl_frameset_find_const(allframes, 
            CR2RES_CAL_FLAT_EXTRACT_1D_PROCATG) != NULL) {
        cpl_propertylist_update_bool(pri_head, "CONTNORM", 1) ;
    } else {
        cpl_propertylist_update_bool(pri_head, "CONTNORM", 0) ;
    }
    cpl_propertylist_set_comment(pri_head, "CONTNORM",
            "TRUE if the spectrum has been divided by the blaze") ;


    cpl_propertylist_update_string(pri_head, "SPECSYS", "TOPOCENT") ;
	cpl_propertylist_set_comment(pri_head, "SPECSYS", 
            "Frame of reference for spectral coordinates");    

    cpl_propertylist_update_string(pri_head, "PROCSOFT",
            PACKAGE "/" PACKAGE_VERSION) ;
    cpl_propertylist_set_comment(pri_head, "PROCSOFT", "ESO pipeline version");

    cpl_propertylist_update_string(pri_head, "REFERENC", "") ;
    cpl_propertylist_set_comment(pri_head, "REFERENC", "Reference publication");

    cpl_propertylist_update_string(pri_head, "PRODCATG", "SCIENCE.SPECTRUM") ;
    cpl_propertylist_set_comment(pri_head, "PRODCATG", "Data product category");


    cpl_propertylist_update_double(pri_head, "WAVELMIN", wmin) ;
    cpl_propertylist_set_comment(pri_head, "WAVELMIN", 
            "Minimum Wavelength [nm]") ;

    cpl_propertylist_update_double(pri_head, "WAVELMAX", wmax) ;
    cpl_propertylist_set_comment(pri_head, "WAVELMAX", 
            "Maximum Wavelength [nm]") ;

    cpl_propertylist_update_double(pri_head, "SPEC_VAL", (wmax+wmin)/2.0) ;
    cpl_propertylist_set_comment(pri_head, "SPEC_VAL", 
            "Characteristic spectral coordinate value [nm]") ;

    cpl_propertylist_update_double(pri_head, "SPEC_BW", wmax-wmin) ;
    cpl_propertylist_set_comment(pri_head, "SPEC_BW", 
            "Width of the spectrum [nm]") ;

    cpl_propertylist_update_string(pri_head, "SPECSYS", "TOPOCENT") ;
    cpl_propertylist_set_comment(pri_head, "SPECSYS", 
            "") ;

    tmp_arr = cpl_array_new(cpl_array_get_size(wlen_arr)-1, CPL_TYPE_DOUBLE);
    for (i=1; i<cpl_array_get_size(wlen_arr);i++){
        cpl_array_set_double(tmp_arr, i-1, 
            cpl_array_get_double(wlen_arr,i, NULL)- 
                cpl_array_get_double(wlen_arr,i-1, NULL)
        );
    }
    spec_bin = cpl_array_get_median(tmp_arr);
    cpl_array_delete(tmp_arr);
    cpl_propertylist_update_double(pri_head, "SPEC_BIN", spec_bin) ;
    cpl_propertylist_set_comment(pri_head, "SPEC_BIN", 
            "Median spectral bin width [nm]") ;

    slitname = cpl_propertylist_get_string(pri_head,"ESO INS SLIT1 ID");
    if (!strcmp(slitname,"w_0.2")) {
        cpl_propertylist_update_double(ext_head, "APERTURE", 5.5555555E-5) ;
        cpl_propertylist_update_double(pri_head, "SPEC_RES", SPEC_RESOL_SLIT02);
    }
    else if (!strcmp(slitname,"w_0.4")) {

        cpl_propertylist_update_double(ext_head, "APERTURE", 1.1111111E-4) ;
        cpl_propertylist_update_double(pri_head, "SPEC_RES", SPEC_RESOL_SLIT04);
    }
    cpl_propertylist_set_comment(ext_head, "APERTURE", 
            "Slit width in deg") ;
    cpl_propertylist_set_comment(pri_head, "SPEC_RES",
                "Nominal resolving power for the given slit"); 

    /* Get some keys from the extension headers*/
    tmp_arr = cpl_array_new(12*CR2RES_NB_DETECTORS, CPL_TYPE_DOUBLE);
    cpl_array * tmp_arr_sig = cpl_array_new(CR2RES_NB_DETECTORS, CPL_TYPE_DOUBLE);
    for (i=0; i<CR2RES_NB_DETECTORS; i++){
        resol = cpl_propertylist_get_double(ext_plist[i], CR2RES_HEADER_QC_SIGNAL);
        if (resol==0){
            cpl_error_reset();
        } else {
            cpl_array_set_double(tmp_arr_sig, i, resol);
        }
        for (ord=0; ord < 12; ord++){
            keyname = cpl_sprintf(CR2RES_HEADER_QC_SNR, ord+1);
            resol = cpl_propertylist_get_double(ext_plist[i], keyname);
            if (resol==0){
                cpl_error_reset();
            } else {
            /* FIXME: The error propagation is not correct. PIPE-11912
             * Correction_factor = <Error_propagated_corrected> / <Error_propagated> 
             * where
             * <Error_propagated_corrected> = [((a*<Error_propagated>)**2 + b**2)**0.5]
             * where a=2.12 and b=26.1, as above.
             */
                double perr = cpl_array_get_double(tmp_arr_sig, i, NULL) / resol;
                double corr = sqrt(pow(2.12*perr,2) + 26.1*26.1) / perr;
                cpl_array_set_double(tmp_arr, ord + (i*12),resol/corr);
                cpl_propertylist_update_double(ext_head, keyname, resol/corr) ;
            }
            cpl_free(keyname);
        }
    }
    cpl_propertylist_update_double(pri_head, "SNR",
                                cpl_array_get_median(tmp_arr)) ;
    cpl_propertylist_set_comment(pri_head, "SNR",
                "Median signal-to-noise in all detector-orders"); 
    cpl_array_delete(tmp_arr);

    if (!strcmp(recipe,"cr2res_obs_pol")) {
        /* Find Pol.TYPE */
        poltype = cr2res_pfits_get_poltype(pri_head);
        tmp_string = cpl_sprintf("/I/%s/",poltype);
        cpl_propertylist_update_string(pri_head, "STOKES", tmp_string);
        cpl_free(tmp_string);
    }

    /* Remove the ASSON keywords */
    cpl_propertylist_erase_regexp(pri_head, "ASSO*", 0);

    /* Add QC from the main Extension */
    if (main_qc_plist != NULL) {
        cpl_propertylist_copy_property_regexp(pri_head, main_qc_plist, "QC*",
                                              0);
    }

    /* Save the main header */
	cpl_propertylist_save(pri_head, idp_filename, CPL_IO_CREATE);


    /* Add Keywords to extension header */
    cpl_propertylist_update_string(ext_head, "EXTNAME", "SPECTRUM") ;

    cpl_propertylist_update_double(ext_head, "RA",
                    cpl_propertylist_get_double(pri_head, "RA"));
    cpl_propertylist_set_comment(ext_head, "RA",
                    cpl_propertylist_get_comment(pri_head, "RA"));
    cpl_propertylist_update_double(ext_head, "DEC",
                    cpl_propertylist_get_double(pri_head, "DEC"));
    cpl_propertylist_set_comment(ext_head, "DEC",
                    cpl_propertylist_get_comment(pri_head, "DEC"));
    cpl_propertylist_update_string(ext_head, "OBJECT",
                    cpl_propertylist_get_string(pri_head, "OBJECT"));
    cpl_propertylist_set_comment(ext_head, "OBJECT",
                    cpl_propertylist_get_comment(pri_head, "OBJECT"));
    tmp_string = cpl_sprintf("%s - %f", 
        cpl_propertylist_get_string(pri_head, "OBJECT"), mjd_start);

    cpl_propertylist_update_string(ext_head, "TITLE", tmp_string);
    cpl_free(tmp_string);
    cpl_propertylist_set_comment(ext_head, "TITLE",
                    "Title is OBJECT and MJD at start");

    cpl_propertylist_update_double(ext_head, "SPEC_VAL", (wmax+wmin)/2.0) ;
    cpl_propertylist_set_comment(ext_head, "SPEC_VAL", 
            "Characteristic spectral coordinate value [nm]") ;
    cpl_propertylist_update_double(ext_head, "SPEC_BW", wmax-wmin) ;
    cpl_propertylist_set_comment(ext_head, "SPEC_BW", 
            "Width of the spectrum [nm]") ;
    
    cpl_propertylist_update_string(ext_head, "VOPUB", "ESO/SAF") ;
    cpl_propertylist_set_comment(ext_head, "VOPUB", "VO Publishing Authority") ;
    cpl_propertylist_update_string(ext_head, "VOCLASS", "SPECTRUM V1.0") ;
    cpl_propertylist_set_comment(ext_head, "VOCLASS",
            "Data Model name and version") ;
    cpl_propertylist_update_int(ext_head, "NELEM", nrows) ;
    cpl_propertylist_set_comment(ext_head, "NELEM", "Length of the data array");
    
    cpl_propertylist_update_string(ext_head, "TUTYP1", 
                        "spec:Data.SpectralAxis.Value");
    cpl_propertylist_update_string(ext_head, "TUCD1", "em.wl");
    cpl_propertylist_update_double(ext_head, "TDMIN1", wmin);
    cpl_propertylist_update_double(ext_head, "TDMAX1", wmax);

    cpl_propertylist_update_string(ext_head, "TUTYP2",
                                "spec:Data.FluxAxis.Value");
    cpl_propertylist_update_string(ext_head, "TUCD2", 
            "phot.flux.density;em.wl;stat.uncalib;arith.ratio");

    cpl_propertylist_update_string(ext_head, "TUTYP3",
                    "spec:Data.FluxAxis.Accuracy.StatError");
    cpl_propertylist_update_string(ext_head, "TUCD3",
                                   "stat.error;phot.flux.density;em.ql;stat.uncalib;arith.ratio");

    cpl_propertylist_update_string(ext_head, "TUTYP4",
                    "spec:Data.FluxAxis.Accuracy.QualityStatus");
    cpl_propertylist_update_string(ext_head, "TUCD4",
                        "meta.code.qual");

    cpl_propertylist_update_string(ext_head, "TUTYP5",
                    "");
    cpl_propertylist_update_string(ext_head, "TUCD5",
                        "instr.order");

    cpl_propertylist_update_string(ext_head, "TUTYP6",
                    "");
    cpl_propertylist_update_string(ext_head, "TUCD6",
                        "meta.number;instr.det"); 

    if (!strcmp(recipe,"cr2res_obs_pol")) { // POL
        cpl_propertylist_update_string(ext_head, "TUTYP7", "");
        tmp_string = cpl_sprintf("phot.count;phys.polarization.stokes.%s",
                                    poltype);
        cpl_propertylist_update_string(ext_head, "TUCD7", tmp_string);
        cpl_free(tmp_string);

        cpl_propertylist_update_string(ext_head, "TUTYP8", "");
        tmp_string = cpl_sprintf("stat.error;phys.polarization.stokes.%s",
                                 poltype);
        cpl_propertylist_update_string(ext_head, "TUCD8", tmp_string);
        cpl_free(tmp_string);
    } else { //  Nodding and staring
        cpl_propertylist_update_string(ext_head, "TUTYP7",
                        "");
        cpl_propertylist_update_string(ext_head, "TUCD7",
                            "pos.cartesian.x;instr.det");

        cpl_propertylist_update_string(ext_head, "TUTYP8",
                        "");
        cpl_propertylist_update_string(ext_head, "TUCD8",
                            "meta.number"); // TraceNb

    }
    /* For Y pixel coordinate in OBS_2D
    TUCDi = pos.cartesian.y;instr.det

    for pol
    phot.count;phys.polarization.stokes.I
    */

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
        cpl_table               **  tables,
        const char              *   recipe,
        const char              *   setting)
{
    cpl_table           *   idp_tab ;
    cpl_table           *   tmp_tab ;
    cpl_array           *   col_names ;
    const char          *   col_name ;
    char                *   col_kind ;
    cpl_type                col_type;
    cpl_propertylist    *   sort_list ;
    int                     order, trace_nb ;
    cpl_size                i, j, ntot, ncols, nrows, nb_selected ;

    /* Check Inputs */
    if (tables == NULL) return NULL ;

    /* Count the Rows */
    ntot = 0 ;
    for (i=0 ; i<CR2RES_NB_DETECTORS ; i++) {
        if (tables[i] != NULL) {
            col_names = cpl_table_get_column_names(tables[i]);
            ncols = cpl_table_get_ncol(tables[i]) ;
            if (strcmp(recipe, "cr2res_obs_pol")) { // nodding and staring
                for (j=0 ; j<ncols ; j++) {
                    col_name = cpl_array_get_string(col_names, j);
                    col_kind = cr2res_dfs_SPEC_colname_parse(col_name, 
                            &order, &trace_nb) ;
                    if (col_kind != NULL && 
                            !strcmp(col_kind, CR2RES_COL_SPEC_SUFFIX)) {
                        /* Handle this extracted spectrum */
                        ntot += cpl_table_get_nrow(tables[i]) ;
                    }
                    if (col_kind != NULL) cpl_free(col_kind) ;
                }
            } else { // POL
                for (j=0 ; j<ncols ; j++) {
                    col_name = cpl_array_get_string(col_names, j);
                    col_kind = cr2res_dfs_POL_colname_parse(col_name, 
                            &order) ;
                    if (col_kind != NULL && 
                            !strcmp(col_kind, CR2RES_COL_POL_INTENS_SUFFIX)) {
                        /* Handle this extracted spectrum */
                        ntot += cpl_table_get_nrow(tables[i]) ;
                    }
                    if (col_kind != NULL) cpl_free(col_kind) ;
                }
            }
            cpl_array_delete(col_names) ;
        }
    }
    cpl_msg_info(__func__, "Found %lld total rows.", ntot);

    /* Create the table */
    tmp_tab = cpl_table_new(ntot) ;
    cpl_table_new_column(tmp_tab, CR2RES_IDP_COL_WAVE, CPL_TYPE_DOUBLE);
    cpl_table_new_column(tmp_tab, CR2RES_IDP_COL_FLUX, CPL_TYPE_DOUBLE);
    cpl_table_new_column(tmp_tab, CR2RES_IDP_COL_ERR, CPL_TYPE_DOUBLE);
    cpl_table_new_column(tmp_tab, CR2RES_IDP_COL_QUAL, CPL_TYPE_INT);
    cpl_table_new_column(tmp_tab, CR2RES_IDP_COL_ORDER, CPL_TYPE_INT);
    cpl_table_new_column(tmp_tab, CR2RES_IDP_COL_DETEC, CPL_TYPE_INT);
    if (strcmp(recipe, "cr2res_obs_pol")) { // nodding and staring
        cpl_table_new_column(tmp_tab, CR2RES_IDP_COL_XPOS, CPL_TYPE_INT);
        cpl_table_new_column(tmp_tab, CR2RES_IDP_COL_TRACE, CPL_TYPE_INT);
    } else {
        cpl_table_new_column(tmp_tab, CR2RES_IDP_COL_STOKES, CPL_TYPE_DOUBLE);
        cpl_table_new_column(tmp_tab,
                                CR2RES_IDP_COL_STOKESERR, CPL_TYPE_DOUBLE);

    }

    /* Fill the table */
    ntot = 0 ;
    for (i=0 ; i<CR2RES_NB_DETECTORS ; i++) {
        if (tables[i] != NULL) {
            col_names = cpl_table_get_column_names(tables[i]);
            ncols = cpl_table_get_ncol(tables[i]) ;
            if (strcmp(recipe, "cr2res_obs_pol")) { // nodding and staring
                for (j=0 ; j<ncols ; j++) {
                    col_name = cpl_array_get_string(col_names, j);
                    col_kind = cr2res_dfs_SPEC_colname_parse(col_name, 
                            &order, &trace_nb) ;
                    if (col_kind != NULL && 
                            !strcmp(col_kind, CR2RES_COL_SPEC_SUFFIX)) {
                        /* Handle this extracted spectrum */
                        cr2res_idp_copy_spec(tmp_tab, tables[i], ntot, i+1, 
                                order, trace_nb, setting) ;
                        ntot += cpl_table_get_nrow(tables[i]) ;
                    }
                    if (col_kind != NULL) cpl_free(col_kind) ;
                }
            } else { 
                for (j=0 ; j<ncols ; j++) {
                    col_name = cpl_array_get_string(col_names, j);
                    col_kind = cr2res_dfs_POL_colname_parse(col_name, &order);
                    if (col_kind != NULL && 
                            !strcmp(col_kind, CR2RES_COL_POL_INTENS_SUFFIX)) {
                        cr2res_idp_copy_pol(tmp_tab, tables[i], ntot, i + 1,
                                        order, setting);
                        ntot += cpl_table_get_nrow(tables[i]) ;
                    }
                    if (col_kind != NULL) cpl_free(col_kind) ;
                }
            }
            cpl_array_delete(col_names) ;
        }
    }

    /* Sort by the wavelengths */
    sort_list = cpl_propertylist_new() ;
    cpl_propertylist_append_bool(sort_list, CR2RES_IDP_COL_WAVE, 0) ;
    //cpl_propertylist_append_bool(sort_list, CR2RES_IDP_COL_XPOS, 1) ;
    cpl_table_sort(tmp_tab, sort_list) ;
    cpl_propertylist_delete(sort_list) ;

    cpl_table_unselect_all(tmp_tab) ;
    /*Find first valid row*/
    int ii;
    int flag = 0;
    for (ii = 0; ii < ntot; ii++) {
        cpl_table_get_double(tmp_tab,CR2RES_IDP_COL_WAVE,ii, &flag) ;
        if(!flag) break;
        cpl_table_select_row(tmp_tab, ii);
    }

    /* Identify close values */
    for (i = ii+1; i < ntot; i++) {
        double cur_val, pre_val;
        flag = 0;
        pre_val = cpl_table_get_double(tmp_tab,CR2RES_IDP_COL_WAVE,i-1, NULL);
        cur_val = cpl_table_get_double(tmp_tab,CR2RES_IDP_COL_WAVE,i, &flag) ;
        if (fabs(pre_val-cur_val) < 1e-3 || flag) {
            cpl_table_select_row(tmp_tab, i) ;
        }
    }

    nb_selected = cpl_table_count_selected(tmp_tab) ;
    if (nb_selected > 0) {
        cpl_msg_info(__func__, 
                "IDP Creation - Double defined WLs : %"CPL_SIZE_FORMAT, 
                nb_selected) ;
    }
    cpl_table_erase_selected(tmp_tab) ;

    /* Transform table to single row with arrays */
    col_names = cpl_table_get_column_names(tmp_tab);
    ncols = cpl_table_get_ncol(tmp_tab);
    nrows = cpl_table_get_nrow(tmp_tab);
    idp_tab = cpl_table_new(1);
    for (j = 0; j < ncols; j++) {
        cpl_array *tmp_arr;
        col_name = cpl_array_get_string(col_names, j);
        cpl_msg_debug(__func__,"colname: %s", col_name);
        col_type = cpl_table_get_column_type(tmp_tab, col_name);
        tmp_arr = cpl_array_new(nrows, col_type);
        if (col_type==CPL_TYPE_DOUBLE){
            cpl_array_copy_data_double(tmp_arr, 
                cpl_table_get_data_double_const(tmp_tab, col_name));
        } else if (col_type==CPL_TYPE_INT){
            cpl_array_copy_data_int(tmp_arr, 
                cpl_table_get_data_int_const(tmp_tab, col_name));
        } else {continue;}
        cpl_table_new_column_array(idp_tab, col_name, col_type, nrows);
        cpl_table_set_array(idp_tab, col_name, 0, tmp_arr);
        cpl_array_delete(tmp_arr);
    }
    cpl_array_delete(col_names) ;
    cpl_table_delete(tmp_tab);


    /* Set the units */
    cpl_table_set_column_unit(idp_tab, CR2RES_IDP_COL_WAVE, "nm") ;
    cpl_table_set_column_unit(idp_tab, CR2RES_IDP_COL_DETEC, " ") ;
    cpl_table_set_column_unit(idp_tab, CR2RES_IDP_COL_ORDER, " ") ;
    cpl_table_set_column_unit(idp_tab, CR2RES_IDP_COL_QUAL, " ") ;
    cpl_table_set_column_unit(idp_tab, CR2RES_IDP_COL_FLUX, " ") ;
    cpl_table_set_column_unit(idp_tab, CR2RES_IDP_COL_ERR, " ") ;
    if (strcmp(recipe, "cr2res_obs_pol")) { // nodding and staring
        cpl_table_set_column_unit(idp_tab, CR2RES_IDP_COL_TRACE, " ") ;
        cpl_table_set_column_unit(idp_tab, CR2RES_IDP_COL_XPOS, "pixel") ;
    } else {
        cpl_table_set_column_unit(idp_tab, CR2RES_IDP_COL_STOKES, " ") ;
        cpl_table_set_column_unit(idp_tab, CR2RES_IDP_COL_STOKESERR, " ") ;
    }

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
    cpl_frame           *   cur_frame ;
    const char          *   cur_fname ;
    cpl_propertylist    *   plist ;
    double                  mjd_end_cur, mjd_end_max, mjd_obs, dit,
                            mjd_start_cur, mjd_start_min ;
    int                     i, nframes, ndit ;

    /* Initialise */
    mjd_end_max = -1.0 ;
    mjd_start_min = 1e10 ;
    nframes = cpl_frameset_get_size(fset) ;

    cpl_errorstate tempes = cpl_errorstate_get();
    cpl_error_reset();

    /* Loop */
    for (int ifr = 0; ifr <= 1; ifr++) {
        for (i = 0; i < nframes; i++) {
            cur_frame = cpl_frameset_get_position(fset, i);
            if (ifr == 0 &&
                cpl_frame_get_group(cur_frame) != CPL_FRAME_GROUP_RAW)
                continue;
            cur_fname = cpl_frame_get_filename(cur_frame);

            /* Get header infos */
            plist = cpl_propertylist_load(cur_fname, 0);
            dit = cr2res_pfits_get_dit(plist);
            ndit = cr2res_pfits_get_ndit(plist);
            mjd_obs = cr2res_pfits_get_mjd_obs(plist);
            if (cpl_error_get_code() == CPL_ERROR_NONE) {
                /* Compute current value */
                mjd_start_cur = mjd_obs;
                mjd_end_cur = mjd_obs + (dit * ndit / (24 * 60 * 60));
            }
            else {
                mjd_end_cur = -2.0;
                mjd_start_cur = 1e11;
                cpl_error_reset();
            }
            /* Update max */
            if (mjd_end_cur > mjd_end_max)
                mjd_end_max = mjd_end_cur;
            if (mjd_start_cur < mjd_start_min) {
                mjd_start_min = mjd_start_cur;
            }
            cpl_propertylist_delete(plist);
        }
        if (mjd_start_min != 1e10)
            break;
    }

    cpl_errorstate_set(tempes);

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
static int cr2res_idp_copy_spec(
    cpl_table       *   out,
    const cpl_table *   in,
    cpl_size            out_start_idx,
    int                 det_nr,
    int                 order,
    int                 tracenb,
    const char          *   setting) 
{
    char            *   spec_name ;
    char            *   wave_name ;
    char            *   err_name ;
    const double    *   pspec ;
    const double    *   pwave ;
    const double    *   perr ;
    cpl_size            in_size, out_size, i, edgepix ;

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
            //if (!isnan(pwave[i]))
                cpl_table_set_double(out, CR2RES_IDP_COL_WAVE, 
                        out_start_idx + i, pwave[i]) ;
            //if (!isnan(pspec[i])) 
                cpl_table_set_double(out, CR2RES_IDP_COL_FLUX, 
                        out_start_idx + i, pspec[i]) ;
            //if (!isnan(perr[i])){ 
            /* FIXME: PIPE-11912 - Error propagation is not correct.
             * Scale up the propagated errors (column name ERR) according to
             * the following empirically derived equation:
             * Error_propagated_corrected = ((a*Error_propagated)^2 + b^2)^0.5
             * where a=2.12 and b=26.1. 
             * For the record, these factors were measured with uncertainties
             * of 0.02 and 1.2, respectively.
             */
                double err_corrected = sqrt(pow(2.12*perr[i],2) + pow(26.1,2));
                cpl_table_set_double(out, CR2RES_IDP_COL_ERR, 
                        out_start_idx + i, err_corrected) ;
            
            if (i<5 || (i>=CR2RES_DETECTOR_SIZE-5)) edgepix = 2;
            else edgepix=0;
            cpl_table_set_int(out, CR2RES_IDP_COL_QUAL, out_start_idx + i, 
                    edgepix + cr2res_wl_is_ghost(setting, pwave[i]));
            cpl_table_set_int(out, CR2RES_IDP_COL_XPOS, out_start_idx + i, 
                    i+1) ;
            cpl_table_set_int(out, CR2RES_IDP_COL_DETEC, out_start_idx + i, 
                    det_nr) ;
            cpl_table_set_int(out, CR2RES_IDP_COL_ORDER, out_start_idx + i, 
                    order) ;
            cpl_table_set_int(out, CR2RES_IDP_COL_TRACE, out_start_idx + i, 
                    tracenb) ;
        }
    }
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Copy polarimetry (stokes, intensity, err, wave) in a bigger table 
                - with additional entries
  @param    out             Destination table
  @param    in              Source table
  @param    out_start_idx   Start index in the destination table
  @param    det_nr          Detector nb
  @param    order           Order nb
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
static int cr2res_idp_copy_pol(
    cpl_table       *   out,
    const cpl_table *   in,
    cpl_size            out_start_idx,
    int                 det_nr,
    int                 order,
    const char          *   setting) 
{
    char            *   spec_name ;
    char            *   wave_name ;
    char            *   err_name ;
    char            *   stokes_name ;
    char            *   stokeserr_name ;
    const double    *   pspec ;
    const double    *   pwave ;
    const double    *   perr ;
    const double    *   pstokes ;
    const double    *   pstokeserr ;
    cpl_size            in_size, out_size, i, edgepix ;

    /* Check entries */
    if (out == NULL  || in == NULL) return -1 ;
    in_size = cpl_table_get_nrow(in) ;
    out_size = cpl_table_get_nrow(out) ;

    /* Get the Spectrum */
    spec_name = cr2res_dfs_POL_INTENS_colname(order) ;
    if ((pspec = cpl_table_get_data_double_const(in, spec_name)) == NULL) {
        cpl_msg_error(__func__, "Cannot find the intensity spec") ;
        cpl_free(spec_name) ;
        return -1 ;
    }
    cpl_free(spec_name) ;

    /* Get the Wavelength */
    wave_name = cr2res_dfs_POL_WAVELENGTH_colname(order) ;
    if ((pwave = cpl_table_get_data_double_const(in, wave_name)) == NULL) {
        cpl_msg_error(__func__, "Cannot find the wavelength") ;
        cpl_free(wave_name) ;
        return -1 ;
    }
    cpl_free(wave_name) ;

    /* Get the Spectrum Error */
    err_name = cr2res_dfs_POL_INTENS_ERROR_colname(order) ;
    if ((perr = cpl_table_get_data_double_const(in, err_name)) == NULL) {
        cpl_msg_error(__func__, "Cannot find the intensity error") ;
        cpl_free(err_name) ;
        return -1 ;
    }
    cpl_free(err_name) ;

    /* Get the Stokes */
    stokes_name = cr2res_dfs_POL_STOKES_colname(order) ;
    if ((pstokes = cpl_table_get_data_double_const(in, stokes_name)) == NULL) {
        cpl_msg_error(__func__, "Cannot find the Stokes spectrum") ;
        cpl_free(stokes_name) ;
        return -1 ;
    }
    cpl_free(stokes_name) ;

    /* Get the stokes err */
    stokeserr_name = cr2res_dfs_POL_STOKES_ERROR_colname(order) ;
    if ((pstokeserr = cpl_table_get_data_double_const(in, stokeserr_name)) == NULL) {
        cpl_msg_error(__func__, "Cannot find the Stokes error") ;
        cpl_free(stokeserr_name) ;
        return -1 ;
    }
    cpl_free(stokeserr_name) ;

    /* copy line by line */
    for (i=0 ; i<in_size ; i++) {
        if (out_start_idx + i <out_size) {
            if (!isnan(pwave[i])) 
                cpl_table_set_double(out, CR2RES_IDP_COL_WAVE, 
                        out_start_idx + i, pwave[i]) ;
            if (!isnan(pspec[i])) 
                cpl_table_set_double(out, CR2RES_IDP_COL_FLUX, 
                        out_start_idx + i, pspec[i]) ;
            if (!isnan(perr[i])) 
                cpl_table_set_double(out, CR2RES_IDP_COL_ERR, 
                        out_start_idx + i, perr[i]) ;
            if (!isnan(pstokes[i])) 
                cpl_table_set_double(out, CR2RES_IDP_COL_STOKES, 
                        out_start_idx + i, pstokes[i]) ;
            if (!isnan(pstokeserr[i])) 
                cpl_table_set_double(out, CR2RES_IDP_COL_STOKESERR, 
                        out_start_idx + i, pstokeserr[i]) ;

            if (i<5 || (i>=CR2RES_DETECTOR_SIZE-5)) edgepix = 2;
            else edgepix=0;
            cpl_table_set_int(out, CR2RES_IDP_COL_QUAL, out_start_idx + i, 
                    edgepix + cr2res_wl_is_ghost(setting, pwave[i]));
            cpl_table_set_int(out, CR2RES_IDP_COL_XPOS, out_start_idx + i, 
                    i+1) ;
            cpl_table_set_int(out, CR2RES_IDP_COL_DETEC, out_start_idx + i, 
                    det_nr) ;
            cpl_table_set_int(out, CR2RES_IDP_COL_ORDER, out_start_idx + i, 
                    order) ;
        }
    }
    return 0 ;
}

/*---------------------------------------------------------------------------*/
/**
  @brief    check whether a wl is affected by ghost
  @param    setting         The setting in question
  @param    wl              The wavelength
  @return   0 or 1, for unaffected and affected, respectively.
 */
/*---------------------------------------------------------------------------*/
int cr2res_wl_is_ghost(const char * setting, double wl){
    cpl_vector * start;
    cpl_vector * end;
    cpl_bivector * ghosts;
    cpl_size nghosts;
    int i;

    ghosts = cr2res_get_ghosts(setting);
    if (ghosts == NULL)
        return 0;
    start = cpl_bivector_get_x(ghosts);
    end = cpl_bivector_get_y(ghosts);
    nghosts = cpl_vector_get_size(start);

    for (i=0; i<nghosts; i++){
        if ( wl > cpl_vector_get(start,i) && wl < cpl_vector_get(end,i)) {
            cpl_bivector_delete(ghosts);
            return 1;
        }
    }
    cpl_bivector_delete(ghosts);
    return 0;
}
/*---------------------------------------------------------------------------*/
/**
  @brief    Create a bivector from static values for ghosts, start and end wl
  @param    setting         The setting in question
  @return   bivector, needs to be deallocated by caller
 */
/*---------------------------------------------------------------------------*/
cpl_bivector * cr2res_get_ghosts(const char * setting){
    cpl_vector * start;
    cpl_vector * end;
    int nghost;

    if (!strcmp(setting,"Y1029")) {
        nghost = 9;
        start = cpl_vector_new(nghost);
        end = cpl_vector_new(nghost);
        cpl_vector_set(start, 0,  955.51); cpl_vector_set(end, 0, 956.59);
        cpl_vector_set(start, 1,  972.36); cpl_vector_set(end, 1, 973.46);
        cpl_vector_set(start, 2,  990.18); cpl_vector_set(end, 2, 990.81);
        cpl_vector_set(start, 3, 1008.18); cpl_vector_set(end, 3, 1008.82);
        cpl_vector_set(start, 4, 1026.81); cpl_vector_set(end, 4, 1027.47);
        cpl_vector_set(start, 5, 1046.15); cpl_vector_set(end, 5, 1046.90);
        cpl_vector_set(start, 6, 1066.29); cpl_vector_set(end, 6, 1067.05);
        cpl_vector_set(start, 7, 1087.15); cpl_vector_set(end, 7, 1087.93);
        cpl_vector_set(start, 8, 1108.83); cpl_vector_set(end, 8, 1109.63);
    } else if (!strcmp(setting,"Y1028")) {
        nghost = 9;
        start = cpl_vector_new(nghost);
        end = cpl_vector_new(nghost);
        cpl_vector_set(start, 0,  954.15); cpl_vector_set(end, 0, 955.20);
        cpl_vector_set(start, 1,  971.25); cpl_vector_set(end, 1, 971.87);
        cpl_vector_set(start, 2,  988.65); cpl_vector_set(end, 2, 989.23);
        cpl_vector_set(start, 3, 1006.61); cpl_vector_set(end, 3, 1007.24);
        cpl_vector_set(start, 4, 1025.23); cpl_vector_set(end, 4, 1025.90);
        cpl_vector_set(start, 5, 1044.60); cpl_vector_set(end, 5, 1045.25);
        cpl_vector_set(start, 6, 1064.67); cpl_vector_set(end, 6, 1065.33);
        cpl_vector_set(start, 7, 1085.46); cpl_vector_set(end, 7, 1086.18);
        cpl_vector_set(start, 8, 1107.18); cpl_vector_set(end, 8, 1107.91);
    } else if (!strcmp(setting,"J1226")) {
        nghost = 9;
        start = cpl_vector_new(nghost);
        end = cpl_vector_new(nghost);
        cpl_vector_set(start, 0, 1119.44); cpl_vector_set(end, 0, 1020.33);
        cpl_vector_set(start, 1, 1142.78); cpl_vector_set(end, 1, 1143.76);
        cpl_vector_set(start, 2, 1167.12); cpl_vector_set(end, 2, 1168.08);
        cpl_vector_set(start, 3, 1192.52); cpl_vector_set(end, 3, 1193.49);
        cpl_vector_set(start, 4, 1219.01); cpl_vector_set(end, 4, 1220.04);
        cpl_vector_set(start, 5, 1246.71); cpl_vector_set(end, 5, 1247.76);
        cpl_vector_set(start, 6, 1275.70); cpl_vector_set(end, 6, 1276.80);
        cpl_vector_set(start, 7, 1306.05); cpl_vector_set(end, 7, 1307.15);
        cpl_vector_set(start, 8, 1337.98); cpl_vector_set(end, 8, 1338.94);
    } else if (!strcmp(setting,"J1228")) {
        nghost = 9;
        start = cpl_vector_new(nghost);
        end = cpl_vector_new(nghost);
        cpl_vector_set(start, 0, 1121.21); cpl_vector_set(end, 0, 1122.12);
        cpl_vector_set(start, 1, 1144.57); cpl_vector_set(end, 1, 1145.53);
        cpl_vector_set(start, 2, 1168.98); cpl_vector_set(end, 2, 1169.91);
        cpl_vector_set(start, 3, 1194.37); cpl_vector_set(end, 3, 1195.38);
        cpl_vector_set(start, 4, 1220.92); cpl_vector_set(end, 4, 1221.95);
        cpl_vector_set(start, 5, 1248.68); cpl_vector_set(end, 5, 1249.72);
        cpl_vector_set(start, 6, 1277.71); cpl_vector_set(end, 6, 1278.83);
        cpl_vector_set(start, 7, 1308.12); cpl_vector_set(end, 7, 1309.17);
        cpl_vector_set(start, 8, 1339.95); cpl_vector_set(end, 8, 1341.03);
    } else if (!strcmp(setting,"J1232")) {
        nghost = 9;
        start = cpl_vector_new(nghost);
        end = cpl_vector_new(nghost);
        cpl_vector_set(start, 0, 1124.71); cpl_vector_set(end, 0, 1125.60);
        cpl_vector_set(start, 1, 1148.15); cpl_vector_set(end, 1, 1149.46);
        cpl_vector_set(start, 2, 1172.62); cpl_vector_set(end, 2, 1173.55);
        cpl_vector_set(start, 3, 1198.13); cpl_vector_set(end, 3, 1199.08);
        cpl_vector_set(start, 4, 1224.77); cpl_vector_set(end, 4, 1225.74);
        cpl_vector_set(start, 5, 1252.56); cpl_vector_set(end, 5, 1253.59);
        cpl_vector_set(start, 6, 1281.74); cpl_vector_set(end, 6, 1282.74);
        cpl_vector_set(start, 7, 1312.19); cpl_vector_set(end, 7, 1313.26);
        cpl_vector_set(start, 8, 1344.15); cpl_vector_set(end, 8, 1345.23);
    } else if (!strcmp(setting,"H1559")) {
        nghost = 2;
        start = cpl_vector_new(nghost);
        end = cpl_vector_new(nghost);
        cpl_vector_set(start, 0, 1682.63); cpl_vector_set(end, 0, 1682.97);
        cpl_vector_set(start, 1, 1735.50); cpl_vector_set(end, 1, 1737.42);
    } else if (!strcmp(setting,"H1567")) {
        nghost = 1;
        start = cpl_vector_new(nghost);
        end = cpl_vector_new(nghost);
        cpl_vector_set(start, 0, 1744.20); cpl_vector_set(end, 0, 1745.84);
    } else if (!strcmp(setting,"H1575")) {
        nghost = 1;
        start = cpl_vector_new(nghost);
        end = cpl_vector_new(nghost);
        cpl_vector_set(start, 0, 1753.20); cpl_vector_set(end, 0, 1754.28);
    } else {
        return NULL;
    }

    return cpl_bivector_wrap_vectors(start,end);
}
