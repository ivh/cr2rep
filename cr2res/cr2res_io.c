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

#include "irplib_wlxcorr.h"

#include "cr2res_io.h"
#include "cr2res_dfs.h"

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_io_save_imagelist(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_imagelist           **  data,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        cpl_type                    type,
        const char              *   recipe,
        const char              *   procatg,
        const char              *   protype) ;

static int cr2res_io_save_image(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  data,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        cpl_type                    type,
        const char              *   recipe,
        const char              *   procatg,
        const char              *   protype) ;
static int cr2res_io_save_table(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  slit_func,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   recipe,
        const char              *   procatg,
        const char              *   protype) ;

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_io  IO related functions
 *
 * TBD
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Create Extname
  @param    detector    The wished detector (1 to CR2RES_NB_DETECTORS)
  @param    data        1 for the data image, 0 for the error
  @return   the newly allocated string with the EXTNAME
 */
/*----------------------------------------------------------------------------*/
char * cr2res_io_create_extname(
        int             detector,
        int             data)
{
    char                *   wished_extname ;

    /* Check entries */
    if (detector < 1 || detector > CR2RES_NB_DETECTORS) return NULL ;

    /* Create wished EXTNAME */
    
    if (data)   wished_extname = cpl_sprintf("CHIP%d.INT1", detector) ;
    else        wished_extname = cpl_sprintf("CHIP%dERR", detector) ;

    return wished_extname ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the wished extension number for a detector
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to CR2RES_NB_DETECTORS)
  @param    data        1 for the data image, 0 for the error
  @return   the Extension number or -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_get_ext_idx(
        const char  *   filename,
        int             detector,
        int             data)
{
    const char          *   extname ;
    char                *   wished_extname ;
    int                     wished_ext_nb = -1 ;
    cpl_propertylist    *   pl ;
    int                     nb_ext, i ;

    /* Check entries */
    if (filename == NULL) return -1 ;
    if (detector < 1 || detector > CR2RES_NB_DETECTORS) return -1 ;

    /* Create wished EXTNAME */
    wished_extname = cr2res_io_create_extname(detector, data) ;

    /* Get the number of extensions */
    nb_ext = cpl_fits_count_extensions(filename) ;

    /* Loop on the extensions */
    for (i=1 ; i<=nb_ext ; i++) {
        /* Get the header */
        pl = cpl_propertylist_load(filename, i) ;
        /* Read the EXTNAME */
        extname = cpl_propertylist_get_string(pl, "EXTNAME");

        /* Compare to the wished one */
        if (strcmp(extname, wished_extname)==0) wished_ext_nb = i ;
        cpl_propertylist_delete(pl) ;
    }
    cpl_free(wished_extname) ;

    return wished_ext_nb ;
}

/*----------------------------------------------------------------------------*/
/*--------------------       LOADING FUNCTIONS       -------------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Load an EMISSION_LINES bivector
  @param    filename    The FITS file name
  @return   The returned object needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_bivector * cr2res_io_load_EMISSION_LINES(
        const char  *   filename)
{
    cpl_bivector    *   lines ;
    double          *   lines_x ;
    double          *   lines_y ;
    cpl_table       *   lines_tab ;
    double              val ;
    int                 i, tab_size, log_flag ;

    /* Check Entries */
    if (filename == NULL) return NULL ;

    /* Initialise */
    log_flag = 0 ;

    /* Load the file in a table */
    lines_tab = cpl_table_load(filename, 1, 1) ;
    tab_size = cpl_table_get_nrow(lines_tab) ;

    /* Create the bivector */
    lines = cpl_bivector_new(tab_size) ;
    lines_x = cpl_bivector_get_x_data(lines) ;
    lines_y = cpl_bivector_get_y_data(lines) ;
    for (i=0 ; i<tab_size ; i++) {
        lines_x[i] = cpl_table_get(lines_tab, CR2RES_COL_WAVELENGTH, i, NULL) ;
        val = cpl_table_get(lines_tab, CR2RES_COL_EMISSION, i, NULL) ;
        if (log_flag && val > 0)    lines_y[i] = log10(val) ;
        else                        lines_y[i] = val ;
    }

    /* Free and return */
    cpl_table_delete(lines_tab) ;
    return lines ;
}
/*----------------------------------------------------------------------------*/
/**
  @brief    Load an image from a MASTER_DARK
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to CR2RES_NB_DETECTORS)
  @param    data        1 for the data image, 0 for the error
  @return   A float type image or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_image * cr2res_io_load_MASTER_DARK(
        const char  *   filename,
        int             detector,
        int             data)
{
    int     wished_ext_nb ;

    /* Check entries */
    if (filename == NULL) return NULL ;
    if (detector < 1 || detector > CR2RES_NB_DETECTORS) return NULL ;

    /* Get the extension number for this detector */
    wished_ext_nb = cr2res_io_get_ext_idx(filename, detector, data) ;

    /* The wished extension was not found */
    if (wished_ext_nb < 0) return NULL ;

    return cpl_image_load(filename, CPL_TYPE_FLOAT, 0, wished_ext_nb);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load an image from a BPM
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to CR2RES_NB_DETECTORS)
  @return   An integer type image or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_image * cr2res_io_load_BPM(
        const char  *   filename,
        int             detector)
{
    int     wished_ext_nb ;

    /* Check entries */
    if (filename == NULL) return NULL ;
    if (detector < 1 || detector > CR2RES_NB_DETECTORS) return NULL ;

    /* Get the extension number for this detector */
    wished_ext_nb = cr2res_io_get_ext_idx(filename, detector, 1) ;

    /* The wished extension was not found */
    if (wished_ext_nb < 0) return NULL ;

    return cpl_image_load(filename, CPL_TYPE_INT, 0, wished_ext_nb);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load the detlin coefficients
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to CR2RES_NB_DETECTORS)
  @return   A float image list with the polynomial coefficients for each
              pixel of the wished detector. The returned object list
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_imagelist * cr2res_io_load_DETLIN_COEFFS(
        const char  *   filename,
        int             detector)
{
    int     wished_ext_nb ;

    /* Check entries */
    if (filename == NULL) return NULL ;
    if (detector < 1 || detector > CR2RES_NB_DETECTORS) return NULL ;

    /* Get the extension number for this detector */
    wished_ext_nb = cr2res_io_get_ext_idx(filename, detector, 1) ;

    /* The wished extension was not found */
    if (wished_ext_nb < 0) return NULL ;

    return cpl_imagelist_load(filename, CPL_TYPE_FLOAT, wished_ext_nb);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load an image from a MASTER_FLAT
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to CR2RES_NB_DETECTORS)
  @param    data        1 for the data image, 0 for the error
  @return   A float image or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_image * cr2res_io_load_MASTER_FLAT(
        const char  *   filename,
        int             detector,
        int             data)
{
    int     wished_ext_nb ;

    /* Check entries */
    if (filename == NULL) return NULL ;
    if (detector < 1 || detector > CR2RES_NB_DETECTORS) return NULL ;

    /* Get the extension number for this detector */
    wished_ext_nb = cr2res_io_get_ext_idx(filename, detector, data) ;

    /* The wished extension was not found */
    if (wished_ext_nb < 0) return NULL ;

    return cpl_image_load(filename, CPL_TYPE_FLOAT, 0, wished_ext_nb);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load a table from a TRACE_WAVE
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to CR2RES_NB_DETECTORS)
  @return   A table or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_io_load_TRACE_WAVE(
        const char  *   filename,
        int             detector)
{
    int                     wished_ext_nb ;

    /* Check entries */
    if (filename == NULL) return NULL ;
    if (detector < 1 || detector > CR2RES_NB_DETECTORS) return NULL ;

    /* Get the extension number for this detector */
    wished_ext_nb = cr2res_io_get_ext_idx(filename, detector, 1) ;

    /* The wished extension was not found */
    if (wished_ext_nb < 0) return NULL ;

    return cpl_table_load(filename, wished_ext_nb, 1);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load an image from a SLIT MODEL
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to CR2RES_NB_DETECTORS)
  @param    data        1 for the data image, 0 for the error
  @return   A float image or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_image * cr2res_io_load_SLIT_MODEL(
        const char  *   filename,
        int             detector,
        int             data)
{
    return NULL ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load an image from a WAVE_MAP
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to CR2RES_NB_DETECTORS)
  @param    data        1 for the data image, 0 for the error
  @return   A float image or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_image * cr2res_io_load_WAVE_MAP(
        const char  *   filename,
        int             detector,
        int             data)
{
    return NULL ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load an image from a SLIT_CURV_MAP
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to CR2RES_NB_DETECTORS)
  @param    data        1 for the data image, 0 for the error
  @return   A float image or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_image * cr2res_io_load_SLIT_CURV_MAP(
        const char  *   filename,
        int             detector,
        int             data)
{
    return NULL ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load a table from a SLIT_CURV
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to CR2RES_NB_DETECTORS)
  @return   A table or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_io_load_SLIT_CURV(
        const char  *   filename,
        int             detector)
{
    return NULL ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load a table from a EXTRACT_1D
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to CR2RES_NB_DETECTORS)
  @return   A table or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_io_load_EXTRACT_1D(
        const char  *   filename,
        int             detector)
{
    int                     wished_ext_nb ;

    /* Check entries */
    if (filename == NULL) return NULL ;
    if (detector < 1 || detector > CR2RES_NB_DETECTORS) return NULL ;

    /* Get the extension number for this detector */
    wished_ext_nb = cr2res_io_get_ext_idx(filename, detector, 1) ;

    /* The wished extension was not found */
    if (wished_ext_nb < 0) return NULL ;

    return cpl_table_load(filename, wished_ext_nb, 1);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load a table from a SPLICED_1D
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to CR2RES_NB_DETECTORS)
  @return   A table or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_io_load_SPLICED_1D(
        const char  *   filename,
        int             detector)
{
    return NULL ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load a table from a EXTRACT_2D
  @param    filename    The FITS file name
  @return   A table or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_io_load_EXTRACT_2D(
        const char  *   filename)
{
    return NULL ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load a table from a EXTRACT_POL
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to CR2RES_NB_DETECTORS)
  @return   A table or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_io_load_EXTRACT_POL(
        const char  *   filename,
        int             detector)
{
    return NULL ;
}

/*----------------------------------------------------------------------------*/
/*---------------------       SAVING FUNCTIONS       -------------------------*/
/*----------------------------------------------------------------------------*/
/* TODO ? set frame levels to mark temp, intermediate and final frames?       */

/*----------------------------------------------------------------------------*/
/**
  @brief    Save EMISSION_LINES file
  @param   	table		The table to save
  @param    parlist     The recipe input parameters
  @param    set         The recipe input frames
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_EMISSION_LINES(
        cpl_table               *   out_table,
        const cpl_parameterlist *   parlist,
        cpl_frameset            *   set,
        const char              *   recipe)
{
    cpl_propertylist    *   plist ;
    char                *   fname ;

    plist = cpl_propertylist_new();
    cpl_propertylist_append_string(plist, "INSTRUME", "CR2RES") ;
    cpl_propertylist_append_string(plist, CPL_DFS_PRO_CATG,
            CR2RES_EMISSION_LINES_PROCATG) ;
    cpl_propertylist_append_string(plist, CPL_DFS_PRO_TYPE,
            CR2RES_PROTYPE_CATALOG) ;

    fname = cpl_sprintf("%s.fits", recipe) ;
    if (cpl_dfs_save_table(set, NULL, parlist, set, NULL, out_table,
                NULL, recipe, plist, NULL,
                PACKAGE "/" PACKAGE_VERSION, fname) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot save the table") ;
        cpl_free(fname);
        return -1 ;
    }
    cpl_free(fname);
    cpl_propertylist_delete(plist) ;

    /* Return */
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a MASTER_DARK
  @param    allframes   The recipe input frames
  @param    filename    The FITS file name
  @param    parlist     The recipe input parameters
  @param    master_darks  The data/error master darks (1 per detector)
  @param    qc_list     The QC parameters
  @param    ext_plist   The extensions property lists
  @param    procatg     The PRO CATG value
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_MASTER_DARK(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  master_darks,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    return cr2res_io_save_image(filename, allframes, parlist,
            master_darks, qc_list, ext_plist, CPL_TYPE_FLOAT, recipe,
            procatg, CR2RES_MASTER_DARK_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a DETLIN COEFFS
  @param    allframes   The recipe input frames
  @param    filename    The FITS file name
  @param    parlist     The recipe input parameters
  @param    coeffs      The detlin coefficients (1 per detector)
  @param    qc_list     The QC parameters
  @param    ext_plist   The extensions property lists
  @param    procatg     The PRO CATG value
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_DETLIN_COEFFS(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_imagelist           **  coeffs,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    return cr2res_io_save_imagelist(filename, allframes, parlist,
            coeffs, qc_list, ext_plist, CPL_TYPE_FLOAT, recipe,
            procatg, CR2RES_DETLIN_COEFFS_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a BPM
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    bpms        The BPMs (1 per detector)
  @param    qc_list     The QC parameters
  @param    ext_plist   The extensions property lists
  @param    procatg     The PRO CATG value
  @param    recipe      The recipe name

  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_BPM(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_image               **  bpms,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    hdrl_image      *   hdrl_bpms[CR2RES_NB_DETECTORS] ;
    int                 det_nr, ret ;
            
    /* Convert to HDRL images */
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
        if (bpms[det_nr-1] == NULL) 
            hdrl_bpms[det_nr-1] = NULL ;
        else
            hdrl_bpms[det_nr-1] = hdrl_image_create(bpms[det_nr-1], NULL) ;
    }

    /* Save */
    ret = cr2res_io_save_image(filename, allframes, parlist,
            hdrl_bpms, qc_list, ext_plist, CPL_TYPE_INT, recipe,
            procatg, CR2RES_BPM_PROTYPE) ;

    /* Free and return */
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
        hdrl_image_delete(hdrl_bpms[det_nr-1]) ;
    }
    return ret ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a CALIBRATED frame
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    calib_collapsed The data/error (1 per detector)
  @param    qc_list     The QC parameters
  @param    ext_plist   The extensions property lists
  @param    procatg     The PRO CATG value
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_CALIBRATED(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  calib_collapsed,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    return cr2res_io_save_image(filename, allframes, parlist,
            calib_collapsed, qc_list, ext_plist, CPL_TYPE_FLOAT, recipe,
            procatg, CR2RES_CALIBRATED_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a MASTER_FLAT
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    master_flats The data/error FLATs (1 per detector)
  @param    qc_list     The QC parameters
  @param    ext_plist   The extensions property lists
  @param    procatg     The PRO CATG value
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_MASTER_FLAT(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  master_flats,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    return cr2res_io_save_image(filename, allframes, parlist,
            master_flats, qc_list, ext_plist, CPL_TYPE_FLOAT, recipe,
            procatg, CR2RES_MASTER_FLAT_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a TRACE_WAVE
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    tables      The tables to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    ext_plist   The extensions property lists
  @param    procatg     The PRO CATG value
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_TRACE_WAVE(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    return cr2res_io_save_table(filename, allframes, parlist, tables,
            qc_list, ext_plist, recipe, procatg, CR2RES_TRACE_WAVE_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a 1D extracted spectrum   
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    tables      The tables to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    ext_plist   The extensions property lists
  @param    procatg     The PRO CATG value
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_EXTRACT_1D(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    return cr2res_io_save_table(filename, allframes, parlist, tables,
            qc_list, ext_plist, recipe, procatg, CR2RES_EXTRACT_1D_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a SLIT_FUNC
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    data        The tables to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    ext_plist   The extensions property lists
  @param    procatg     The PRO CATG value
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_SLIT_FUNC(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  slit_func,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    return cr2res_io_save_table(filename, allframes, parlist, slit_func,
            qc_list, ext_plist, recipe, procatg, CR2RES_SLIT_FUNC_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a SLIT_MODEL
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    data        The data images to save (DATA and ERROR per detector)
  @param    qc_list     The QC parameters
  @param    ext_plist   The extensions property lists
  @param    procatg     The PRO CATG value
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_SLIT_MODEL(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  data,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    return cr2res_io_save_image(filename, allframes, parlist,
            data, qc_list, ext_plist, CPL_TYPE_FLOAT, recipe, procatg, 
            CR2RES_SLIT_MODEL_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a WAVE_MAP
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    data        The data images to save (DATA and ERROR per detector)
  @param    qc_list     The QC parameters
  @param    ext_plist   The extensions property lists
  @param    procatg     The PRO CATG value
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_WAVE_MAP(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  data,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    return cr2res_io_save_image(filename, allframes, parlist,
            data, qc_list, ext_plist, CPL_TYPE_FLOAT, recipe,
            procatg, CR2RES_WAVE_MAP_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a SLIT_CURV_MAP
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    data        The data images to save (DATA and ERROR per detector)
  @param    qc_list     The QC parameters
  @param    ext_plist   The extensions property lists
  @param    procatg     The PRO CATG value
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_SLIT_CURV_MAP(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  data,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    return cr2res_io_save_image(filename, allframes, parlist,
            data, qc_list, ext_plist, CPL_TYPE_FLOAT, recipe,
            procatg, CR2RES_SLIT_CURV_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a SLIT_CURV
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    tables      The tables to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    ext_plist   The extensions property lists
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_SLIT_CURV(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    return cr2res_io_save_table(filename, allframes, parlist, tables,
            qc_list, ext_plist, recipe, procatg, CR2RES_SLIT_CURV_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a SPLICED_1D
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    tables      The tables to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    ext_plist   The extensions property lists
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_SPLICED_1D(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   recipe)
{
    return -1 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a EXTRACT_2D
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    table       The table to save
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_EXTRACT_2D(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               *   table,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe)
{
    return -1 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a EXTRACT_POL
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    tables      The tables to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    ext_plist   The extensions property lists
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_EXTRACT_POL(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   recipe)
{
    return -1 ;
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a multi extension table
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    tab         The tables to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    ext_plist   The extensions property lists
  @param    recipe      The recipe name
  @param    procatg     PRO.CATG
  @param    protype     PRO.TYPE
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
static int cr2res_io_save_table(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tab,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   recipe,
        const char              *   procatg,
        const char              *   protype)
{
    cpl_propertylist    *   pro_list ;
    cpl_propertylist    *   ext_head ;
    char                *   wished_extname ;
    int                     i ;

    /* Test entries */
    if (allframes == NULL || filename == NULL || ext_plist == NULL) return -1 ;

    /* Add the PRO keys */
    if (qc_list != NULL) pro_list = cpl_propertylist_duplicate(qc_list) ;
    else pro_list = cpl_propertylist_new() ;

    /* Add PRO Keys */
    cpl_propertylist_append_string(pro_list, CPL_DFS_PRO_CATG, procatg) ;
    cpl_propertylist_append_string(pro_list, CPL_DFS_PRO_TYPE, protype) ;

    /* Create the first extension header */
    if (ext_plist[0]==NULL) {
        ext_head = cpl_propertylist_new() ;
    } else {
        ext_head = cpl_propertylist_duplicate(ext_plist[0]);
        cpl_propertylist_erase(ext_head, "EXTNAME");
    }
    wished_extname = cr2res_io_create_extname(1, 1) ;
    cpl_propertylist_update_string(ext_head, "EXTNAME", wished_extname) ;
    cpl_free(wished_extname) ;

    /* Save the first extension */
    if (tab[0] != NULL) {
        if (cpl_dfs_save_table(allframes, NULL, parlist, allframes, NULL,
                    tab[0], ext_head, recipe, pro_list, NULL,
                    PACKAGE "/" PACKAGE_VERSION, filename) != CPL_ERROR_NONE) {
            cpl_msg_error(__func__, "Cannot save the first extension table") ;
            cpl_propertylist_delete(ext_head) ;
            cpl_propertylist_delete(pro_list) ;
            return -1 ;
        }
    } else {
        if (cpl_dfs_save_propertylist(allframes, NULL, parlist,
                    allframes, NULL, recipe, pro_list, NULL, PACKAGE "/"
                    PACKAGE_VERSION, filename) != CPL_ERROR_NONE) {
            cpl_msg_error(__func__, "Cannot save the empty HDU") ;
            cpl_propertylist_delete(ext_head) ;
            cpl_propertylist_delete(pro_list) ;
            return -1 ;
        }
        cpl_propertylist_save(ext_head, filename, CPL_IO_EXTEND) ;
    }
    cpl_propertylist_delete(ext_head) ;
    cpl_propertylist_delete(pro_list) ;

    /* Save the extensions */
    for (i=1 ; i<CR2RES_NB_DETECTORS ; i++) {
        /* Create the first extension header */
        if (ext_plist[i] == NULL) {
            ext_head = cpl_propertylist_new() ;
        } else {
            ext_head = cpl_propertylist_duplicate(ext_plist[i]);
            cpl_propertylist_erase(ext_head, "EXTNAME");
        }
        wished_extname = cr2res_io_create_extname(i+1, 1) ;
        cpl_propertylist_update_string(ext_head, "EXTNAME", wished_extname) ;
        cpl_free(wished_extname) ;

        if (tab[i] != NULL) {
            cpl_table_save(tab[i], NULL, ext_head, filename, CPL_IO_EXTEND) ;
        } else {
            cpl_propertylist_save(ext_head, filename, CPL_IO_EXTEND) ;
        }
        cpl_propertylist_delete(ext_head) ;
    }
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a multi extension images
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    data        The images to save (data and error per detector)
  @param    qc_list     The QC parameters
  @param    ext_plist   The extensions property lists
  @param    type        CPL_TYPE_FLOAT, CPL_TYPE_INT,...
  @param    recipe      The recipe name
  @param    procatg     PRO.CATG
  @param    protype     PRO.TYPE
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
static int cr2res_io_save_image(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  data,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        cpl_type                    type,
        const char              *   recipe,
        const char              *   procatg,
        const char              *   protype)
{
    cpl_propertylist    *   qclist_loc ;
    cpl_image           *   to_save ;
    char          		*   wished_extname ;
    int                     det_nr ;

    /* Create a local QC list and add the PRO.CATG */
    if (qc_list == NULL) {
        qclist_loc = cpl_propertylist_new();
    } else {
        qclist_loc = cpl_propertylist_duplicate(qc_list) ;
    }
    cpl_propertylist_update_string(qclist_loc, CPL_DFS_PRO_CATG, procatg);
    cpl_propertylist_update_string(qclist_loc, CPL_DFS_PRO_TYPE, protype);

    /* Create the Primary Data Unit without data */
    if (cpl_dfs_save_image(allframes, NULL, parlist, allframes, NULL, NULL,
                type, recipe, qclist_loc, NULL,
                PACKAGE "/" PACKAGE_VERSION, filename) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot save the empty primary HDU") ;
        cpl_propertylist_delete(qclist_loc) ;
        return -1 ;
    }
    /* Delete PRO LIST */
    cpl_propertylist_delete(qclist_loc) ;

    /* Save the extensions */
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
        if (ext_plist[det_nr-1] == NULL) {
            qclist_loc = cpl_propertylist_new();
        } else {
            qclist_loc = cpl_propertylist_duplicate(ext_plist[det_nr-1]) ;
            cpl_propertylist_erase(qclist_loc, "EXTNAME");
        }

        /* Save the DATA */
        wished_extname = cr2res_io_create_extname(det_nr, 1) ;
        cpl_propertylist_prepend_string(qclist_loc, "EXTNAME", wished_extname) ;

        if (data[det_nr-1] == NULL)
            to_save = NULL ;
        else                        
            to_save = hdrl_image_get_image(data[det_nr-1]);
        cpl_image_save(to_save, filename, type, qclist_loc, CPL_IO_EXTEND) ;
        cpl_free(wished_extname) ;

        /* Save the NOISE */
        wished_extname = cr2res_io_create_extname(det_nr, 0) ;
        cpl_propertylist_update_string(qclist_loc, "EXTNAME", wished_extname) ;
        if (data[det_nr-1] == NULL)    
            to_save = NULL ;
        else                        
            to_save = hdrl_image_get_error(data[det_nr-1]);
        cpl_image_save(to_save, filename, type, qclist_loc, CPL_IO_EXTEND) ;
        cpl_propertylist_delete(qclist_loc) ;
        cpl_free(wished_extname) ;
    }

	return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a multi extension imagelists
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    data        The imagelists to save (only data for the moment)
  @param    qc_list     The QC parameters
  @param    ext_plist   The extensions property lists
  @param    type        CPL_TYPE_FLOAT, CPL_TYPE_INT,...
  @param    recipe      The recipe name
  @param    procatg     PRO.CATG
  @param    protype     PRO.TYPE
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
static int cr2res_io_save_imagelist(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_imagelist           **  data,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        cpl_type                    type,
        const char              *   recipe,
        const char              *   procatg,
        const char              *   protype)
{
    cpl_propertylist    *   qclist_loc ;
    char          		*   wished_extname ;
    int                     det_nr ;

    /* Create a local QC list and add the PRO.CATG */
    if (qc_list == NULL) {
        qclist_loc = cpl_propertylist_new();
    } else {
        qclist_loc = cpl_propertylist_duplicate(qc_list) ;
    }
    cpl_propertylist_update_string(qclist_loc, CPL_DFS_PRO_CATG, procatg);
    cpl_propertylist_update_string(qclist_loc, CPL_DFS_PRO_TYPE, protype);

    /* Create the Primary Data Unit without data */
    if (cpl_dfs_save_image(allframes, NULL, parlist, allframes, NULL, NULL,
                type, recipe, qclist_loc, NULL,
                PACKAGE "/" PACKAGE_VERSION, filename) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot save the empty primary HDU") ;
        cpl_propertylist_delete(qclist_loc) ;
        return -1 ;
    }
    /* Delete PRO LIST */
    cpl_propertylist_delete(qclist_loc) ;

    /* Save the extensions */
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {
        if (ext_plist[det_nr-1] == NULL) {
            qclist_loc = cpl_propertylist_new();
        } else {
            qclist_loc = cpl_propertylist_duplicate(ext_plist[det_nr-1]) ;
            cpl_propertylist_erase(qclist_loc, "EXTNAME");
        }

        /* Save the DATA */
        wished_extname = cr2res_io_create_extname(det_nr, 1) ;
        cpl_propertylist_prepend_string(qclist_loc, "EXTNAME", wished_extname) ;

        if (data[det_nr-1] == NULL) {
            cpl_propertylist_save(qclist_loc, filename, CPL_IO_EXTEND) ;
        } else {                        
            cpl_imagelist_save(data[det_nr-1], filename, type, qclist_loc, 
                    CPL_IO_EXTEND) ;
        }
        cpl_free(wished_extname) ;

        /* Save the NOISE */
        wished_extname = cr2res_io_create_extname(det_nr, 0) ;
        cpl_propertylist_update_string(qclist_loc, "EXTNAME", wished_extname) ;
        if (data[det_nr-1] == NULL) {
            cpl_propertylist_save(qclist_loc, filename, CPL_IO_EXTEND) ;
        } else {                        
            cpl_imagelist_save(data[det_nr-1], filename, type, qclist_loc, 
                    CPL_IO_EXTEND) ;
        }
        cpl_propertylist_delete(qclist_loc) ;
        cpl_free(wished_extname) ;
    }

	return 0 ;
}
