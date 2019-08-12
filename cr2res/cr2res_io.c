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
#include "cr2res_pfits.h"

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_io_check_pro_type(
        const char  *   filename,
        const char  *   expected_protype) ;
static cpl_imagelist * cr2res_hdrl_to_cpl_list(
        const hdrl_imagelist    *   in, 
        int                         data) ;
static int cr2res_table_check_column(
        const cpl_table     *   tab,
        const char          *   col) ;
static int cr2res_io_save_imagelist(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        hdrl_imagelist          **  data,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        cpl_type                    type,
        const char              *   recipe,
        const char              *   procatg,
        const char              *   protype) ;

static int cr2res_io_save_image(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
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
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  slit_func,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   recipe,
        const char              *   procatg,
        const char              *   protype) ;
static int cr2res_io_save_one_table(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        cpl_table               *   tab,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        *   ext_plist,
        const char              *   extname,
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
  @brief    Get the first SLIT_MODEL frame from a frameset 
  @param    set     Input frame set
  @param    setting The setting to match 
  @param    decker  The decker position to match
  @return   The new allocated frame or NULL in error case or if it is missing
 */
/*----------------------------------------------------------------------------*/
cpl_frame * cr2res_io_find_SLIT_MODEL(
        const cpl_frameset  *   in,
        const char          *   setting_id,
        cr2res_decker           decker)
{
    cpl_frameset        *   fset ;
    cpl_frame           *   out ;
    const cpl_frame     *   cur_frame ;
    const char          *   cur_fname ;
    cpl_propertylist    *   plist ;
    cr2res_decker           cur_decker ;
    char                *   cur_setting ;
    int                     i ;

    /* Check entries */
    if (in == NULL) return NULL ;

    /* Initialise */
    out = NULL ;

    /* Get the slit model frames */
    fset=cr2res_extract_frameset(in, CR2RES_CAL_FLAT_SLIT_MODEL_PROCATG) ;
    if (fset == NULL) 
        fset=cr2res_extract_frameset(in, CR2RES_UTIL_SLIT_MODEL_PROCATG) ;
    if (fset == NULL) 
        fset=cr2res_extract_frameset(in,CR2RES_OBS_NODDING_SLITMODELA_PROCATG) ;
    if (fset == NULL) 
        fset=cr2res_extract_frameset(in,CR2RES_OBS_NODDING_SLITMODELB_PROCATG) ;
    if (fset == NULL) return NULL ;

    /* Find out if there is a matching one */
    for (i=0 ; i<cpl_frameset_get_size(fset) ; i++) {
        if (out == NULL) {
			/* Get the Current Frame */
			cur_frame = cpl_frameset_get_position(fset, i) ;
			cur_fname = cpl_frame_get_filename(cur_frame) ;
            plist = cpl_propertylist_load(cur_fname, 0);
            cur_setting = cpl_strdup(cr2res_pfits_get_wlen_id(plist)) ;
            cr2res_format_setting(cur_setting) ;
            cur_decker = cr2res_pfits_get_decker_position(plist) ;
            cpl_propertylist_delete(plist) ;
            if (!strcmp(cur_setting, setting_id) && cur_decker == decker) 
                out = cpl_frame_duplicate(cur_frame) ;
            cpl_free(cur_setting) ;
        }
    }
    cpl_frameset_delete(fset) ;
    return out ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the first TRACE_WAVE frame from a frameset
  @param    set     Input frame set
  @return   the frame reference or NULL in error case or if it is missing
 */
/*----------------------------------------------------------------------------*/
const cpl_frame * cr2res_io_find_TRACE_WAVE(const cpl_frameset * in)
{
    const cpl_frame *   out ;

    /* Check entries */
    if (in == NULL) return NULL ;

    out=cpl_frameset_find_const(in, CR2RES_CAL_FLAT_TW_PROCATG) ;
    if (out == NULL) 
        out=cpl_frameset_find_const(in, CR2RES_CAL_FLAT_TW_MERGED_PROCATG) ;
    if (out == NULL) 
        out=cpl_frameset_find_const(in, CR2RES_UTIL_TRACE_TW_PROCATG) ;
    if (out == NULL) 
        out=cpl_frameset_find_const(in, CR2RES_UTIL_WAVE_TW_PROCATG) ;
    if (out == NULL) 
        out=cpl_frameset_find_const(in, CR2RES_CAL_WAVE_TW_PROCATG) ;
    if (out == NULL) 
       out=cpl_frameset_find_const(in,CR2RES_UTIL_SLIT_CURV_TW_PROCATG);
    return out ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the first BPM frame from a frameset
  @param    set     Input frame set
  @return   the frame reference or NULL in error case or if it is missing
 */
/*----------------------------------------------------------------------------*/
const cpl_frame * cr2res_io_find_BPM(const cpl_frameset * in)
{
    const cpl_frame *   out ;

    /* Check entries */
    if (in == NULL) return NULL ;

    out=cpl_frameset_find_const(in, CR2RES_CAL_DARK_BPM_PROCATG) ;
    if (out == NULL) 
        out = cpl_frameset_find_const(in, CR2RES_CAL_FLAT_BPM_PROCATG) ;
    if (out == NULL) 
        out = cpl_frameset_find_const(in, CR2RES_CAL_DETLIN_BPM_PROCATG) ;
    if (out == NULL) 
        out = cpl_frameset_find_const(in, CR2RES_UTIL_BPM_SPLIT_PROCATG) ;
    if (out == NULL) 
        out = cpl_frameset_find_const(in, CR2RES_UTIL_NORM_BPM_PROCATG) ;
    return out ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the DITS from a frame set
  @param    set     Input frame set
  @return   the DITS or NULL in error case
 */
/*----------------------------------------------------------------------------*/
cpl_vector * cr2res_io_read_dits(const cpl_frameset * in)
{
    cpl_vector          *   dits ;
    cpl_propertylist    *   plist ;
    cpl_size                i ;

    /* Check entries */
    if (in == NULL) return NULL ;

    /* Allocate the vector */
    dits = cpl_vector_new(cpl_frameset_get_size(in)) ;

    /* Loop on the frames */
    for (i=0 ; i< cpl_vector_get_size(dits) ; i++) {
        plist = cpl_propertylist_load(cpl_frame_get_filename(
                    cpl_frameset_get_position_const(in, i)), 0) ;
        cpl_vector_set(dits, i, cr2res_pfits_get_dit(plist)) ;
        cpl_propertylist_delete(plist) ;
    }

    return dits ;
}
/*----------------------------------------------------------------------------*/
/**
  @brief    Get the decker positions from a frame set
  @param    set     Input frame set
  @return   the DECKER positions or NULL in error case
 */
/*----------------------------------------------------------------------------*/
cr2res_decker * cr2res_io_read_decker_positions(const cpl_frameset * in)
{
    cr2res_decker       *   out ;
    cpl_propertylist    *   plist ;
    const char          *   fname ;
    cpl_size                nframes, i ;

    /* Check entries */
    if (in == NULL) return NULL ;

    /* Initialise */
    nframes = cpl_frameset_get_size(in) ;

    /* Allocate the vector */
    out = cpl_malloc(nframes * sizeof(cr2res_decker)) ;

    /* Loop on the frames */
    for (i=0 ; i< nframes ; i++) {
        plist = cpl_propertylist_load(cpl_frame_get_filename(
                    cpl_frameset_get_position_const(in, i)), 0) ;
        out[i] = cr2res_pfits_get_decker_position(plist) ;
        cpl_propertylist_delete(plist) ;
    }
    return out ;
}

/*----------------------------------------------------------------------------*/
/**
   @brief   Extract the frames with the given tag and Decker position
   @param   in      A non-empty frameset
   @param   tag     The tag of the requested frames
   @param   decker  CR2RES_DECKER_NONE,CR2RES_DECKER_1_3 or CR2RES_DECKER_2_4
   @return  The newly created frameset or NULL on error

   The returned frameset must be de allocated with cpl_frameset_delete()
 */
/*----------------------------------------------------------------------------*/
cpl_frameset * cr2res_io_extract_decker_frameset(
        const cpl_frameset  *   in,
        const char          *   tag,
        cr2res_decker           decker)
{
    cpl_frameset        *   out ;
    const cpl_frame     *   cur_frame ;
    cpl_frame           *   loc_frame ;
    int                     nbframes;
    cpl_propertylist    *   plist ;
    int                     i ;

    /* Test entries */
    if (in == NULL) return NULL ;
    if (tag == NULL) return NULL ;
    if (decker != CR2RES_DECKER_NONE && decker != CR2RES_DECKER_1_3 &&
            decker != CR2RES_DECKER_2_4)
        return NULL ;

    /* Initialise */
    nbframes = cpl_frameset_get_size(in) ;

    /* Count the frames with the tag */
    if ((cpl_frameset_count_tags(in, tag)) == 0) return NULL ;

    /* Create the output frameset */
    out = cpl_frameset_new() ;

    /* Loop on the requested frames and store them in out */
    for (i=0 ; i<nbframes ; i++) {
        cur_frame = cpl_frameset_get_position_const(in, i) ;
        if (!strcmp(cpl_frame_get_tag(cur_frame), tag)) {
            /* Get the propertylist */
            plist = cpl_propertylist_load(cpl_frame_get_filename(cur_frame), 0);
            if (cr2res_pfits_get_decker_position(plist) == decker) {
                loc_frame = cpl_frame_duplicate(cur_frame) ;
                cpl_frameset_insert(out, loc_frame) ;
            }
            cpl_propertylist_delete(plist) ;
        }
    }
    /* No matching frame found */
    if (cpl_frameset_get_size(out) == 0){
        cpl_frameset_delete(out) ;
        return NULL ;
    }
    return out ;
}
/*----------------------------------------------------------------------------*/
/**
  @brief    Convert the order to the keyword index
  @param    order   Order (-49 to 50)
  @return   the order index or a negative value in error case
            (00 to 99)
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_convert_order_to_idx(int order)
{
    /* Check entries */
    if (order < -49 || order > 50) return -1 ;

    /* Conversion order <-> keyword Index */
    if (order < 0)  return order + 100 ;
    else            return order ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Convert the keyword index to the order
  @param    order_idx   the order index (00 to 99)
  @return   Order (-50 to 50)
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_convert_idx_to_order(int order_idx)
{
    /* Check entries */
    if (order_idx < 0 || order_idx > 99) return -1 ;

    /* Conversion order <-> keyword Index */
    if (order_idx > 50) return order_idx - 100 ;
    else                return order_idx ;
}

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
    /* If this changes, update warning message in cr2res_io_get_ext_idx() */
    if (data)   wished_extname = cpl_sprintf("CHIP%d.INT1", detector) ;
    else        wished_extname = cpl_sprintf("CHIP%dERR.INT1", detector) ;

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

    /* EXTNAME expectation */
    if (wished_ext_nb < 0) {
        /* Warning only for data - error are sometimeѕ optional */
        if (data == 1)
            cpl_msg_warning(__func__,
        "EXTNAME is supposed to match CHIPn.INT1 or CHIPnERR.INT1 (n=1/2/3)") ;
    } 

    return wished_ext_nb ;
}

/*----------------------------------------------------------------------------*/
/*--------------------       LOADING FUNCTIONS       -------------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Load an hdrl image from a image file
  @param    fname       The input file name
  @param    detector    The wished detector (1 to CR2RES_NB_DETECTORS)
  @return   A hdrl image or NULL in error case. 
            The returned object needs to be deallocated
  This function load imageѕ files (also where the error is missing)
  with the proper EXTNAME convention
 */
/*----------------------------------------------------------------------------*/
hdrl_image * cr2res_io_load_image(
        const char  *   in,
        int             detector)
{
    hdrl_image      *   out ;
    cpl_image       *   data ;
    cpl_image       *   err ;
    int                 ext_nr_data, ext_nr_err ;

    /* Check entries */
    if (in == NULL) return NULL ;
    if (detector < 1 || detector > CR2RES_NB_DETECTORS) return NULL ;

    /* Get the extension numbers for this detector */
    ext_nr_data = cr2res_io_get_ext_idx(in, detector, 1) ;
    ext_nr_err = cr2res_io_get_ext_idx(in, detector, 0) ;

    /* The wished extension was not found */
    if (ext_nr_data < 0) return NULL ;
    
    /* Load the image */
    data = cpl_image_load(in, CPL_TYPE_DOUBLE, 0, ext_nr_data) ;
    if (ext_nr_err >= 0) 
        err = cpl_image_load(in, CPL_TYPE_DOUBLE, 0, ext_nr_err) ;
    else    
        err = NULL ;
    out = hdrl_image_create(data, err) ;
    if (data != NULL) cpl_image_delete(data) ;
    if (err != NULL) cpl_image_delete(err) ;

    /* Return  */
    return out ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load an hdrl image list from a cube file
  @param    fname       The input file name
  @param    detector    The wished detector (1 to CR2RES_NB_DETECTORS)
  @return   A hdrl image list or NULL in error case. 
            The returned object needs to be deallocated
  This function loads image list from a cube (also where the error is missing)
  with the proper EXTNAME convention
 */
/*----------------------------------------------------------------------------*/
hdrl_imagelist * cr2res_io_load_image_list(
        const char  *   in,
        int             detector)
{
    hdrl_imagelist  *   out ;
    cpl_imagelist   *   data ;
    cpl_imagelist   *   error ;
    int                 wished_ext_nb_data, wished_ext_nb_error ;

    /* Check entries */
    if (in == NULL) return NULL ;
    if (detector < 1 || detector > CR2RES_NB_DETECTORS) return NULL ;

    /* Get the extension number for this detector */
    wished_ext_nb_data = cr2res_io_get_ext_idx(in, detector, 1) ;
    wished_ext_nb_error = cr2res_io_get_ext_idx(in, detector, 0) ;

    /* The wished extension was not found */
    if (wished_ext_nb_data < 0) return NULL ;
    
    /* Load the image list */
    data = cpl_imagelist_load(in, CPL_TYPE_DOUBLE, wished_ext_nb_data);
    if (wished_ext_nb_error >= 0) 
        error = cpl_imagelist_load(in, CPL_TYPE_DOUBLE, wished_ext_nb_error);
    else 
        error = NULL ;
    out = hdrl_imagelist_create(data, error) ;
    if (data != NULL) cpl_imagelist_delete(data) ;
    if (error != NULL) cpl_imagelist_delete(error) ;

    /* Return  */
    return out ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load an hdrl image list from an images frameset
  @param    fset        The input frame set
  @param    detector    The wished detector (1 to CR2RES_NB_DETECTORS)
  @return   A hdrl imagelist or NULL in error case. The returned object
              needs to be deallocated
  The returned hdrl image list contains the list of all data images
  for a given detector from a list of input image frames
  This function load imageѕ files (also where the error is missing)
  with the proper EXTNAME convention
 */
/*----------------------------------------------------------------------------*/
hdrl_imagelist * cr2res_io_load_image_list_from_set(
        const cpl_frameset  *   in,
        int                     detector)
{
    const char      *   first_file ;
    hdrl_imagelist  *   out ;
    cpl_imagelist   *   data ;
    cpl_imagelist   *   err ;
    int                 ext_nr_data, ext_nr_err ;

    /* Check entries */
    if (in == NULL) return NULL ;
    if (cpl_frameset_get_size(in) < 1) return NULL ;
    if (detector < 1 || detector > CR2RES_NB_DETECTORS) return NULL ;

    /* Get the extension number for this detector */
    first_file = cpl_frame_get_filename(cpl_frameset_get_position_const(in,0)) ;
    ext_nr_data = cr2res_io_get_ext_idx(first_file, detector, 1) ;
    ext_nr_err = cr2res_io_get_ext_idx(first_file, detector, 0) ;

    /* The wished extension was not found */
    if (ext_nr_data < 0) return NULL ;
    
    /* Load the image list */
    data = cpl_imagelist_load_frameset(in, CPL_TYPE_DOUBLE, 1, ext_nr_data) ;
    if (ext_nr_err >= 0) 
        err = cpl_imagelist_load_frameset(in, CPL_TYPE_DOUBLE, 1, ext_nr_err) ;
    else 
        err = NULL ;
    out = hdrl_imagelist_create(data, err) ;
    if (data != NULL) cpl_imagelist_delete(data) ;
    if (err != NULL) cpl_imagelist_delete(err) ;

    /* Return  */
    return out ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load the table accordingly
  @param    in          the input file
  @param    det_nr      the detector number
  @param    pmin        the first pixel to load (-1 if all)
  @param    pmax        the last pixel to load (-1 if all)
  @return   the loaded table or NULL in error case
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_load_table(
        const char  *   in,
        int             det_nr,
        int             pmin,
        int             pmax)
{
    cpl_table   *   out ;
    cpl_table   *   out_tmp ;
    int             ext_nr ;

    /* Check entries */
    if (in == NULL) return NULL ;
    if (det_nr < 1 || det_nr > CR2RES_NB_DETECTORS) return NULL ;

    /* Get the extension number for this detector */
    ext_nr = cr2res_io_get_ext_idx(in, det_nr, 1) ;

    /* Load the table */
    if ((out = cpl_table_load(in, det_nr, 0)) == NULL) {
        cpl_msg_error(__func__, "Cannot load %s as a table", in) ;
        return NULL ;
    }

    /* Select only between pmin and pmax */
    if (pmin>0 && pmax>0 && pmax >= pmin && pmax <= cpl_table_get_nrow(out)) {
        out_tmp = cpl_table_extract(out, pmin, pmax-pmin+1) ;
        if (out_tmp != NULL) {
            cpl_table_delete(out) ;
            out = out_tmp ;
        }
    }
    return out ;
}

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

    /* Check PRO.TYPE */
    if (cr2res_io_check_pro_type(filename, CR2RES_PROTYPE_CATALOG) != 1) {
        cpl_msg_info(__func__, "File check failed for %s", filename);
        return NULL ;
    }

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
  @brief    Load an image from a BPM
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to CR2RES_NB_DETECTORS)
  @param    data        1 for the data image, 0 for the error
  @return   An integer type image or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_image * cr2res_io_load_BPM(
        const char  *   filename,
        int             detector,
        int             data)
{
    int                     wished_ext_nb ;

    /* Check entries */
    if (filename == NULL) return NULL ;
    if (detector < 1 || detector > CR2RES_NB_DETECTORS) return NULL ;

    /* Check PRO.TYPE */
    if (cr2res_io_check_pro_type(filename, CR2RES_BPM_PROTYPE) != 1)
        return NULL ;

    /* Get the extension number for this detector */
    wished_ext_nb = cr2res_io_get_ext_idx(filename, detector, data) ;

    /* The wished extension was not found */
    if (wished_ext_nb < 0) return NULL ;

    return cpl_image_load(filename, CPL_TYPE_INT, 0, wished_ext_nb);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load an image from a MASTER_DARK
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to CR2RES_NB_DETECTORS)
  @return   A hdrl image or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
hdrl_image * cr2res_io_load_MASTER_DARK(
        const char  *   filename,
        int             detector)
{
    hdrl_image          *   master_dark ;

    /* Check entries */
    if (filename == NULL) return NULL ;
    if (detector < 1 || detector > CR2RES_NB_DETECTORS) return NULL ;

    /* Check PRO.TYPE */
    if (cr2res_io_check_pro_type(filename, CR2RES_MASTER_DARK_PROTYPE) != 1)
        return NULL ;

    /* Load */
    master_dark = cr2res_io_load_image(filename, detector) ;

    /* Error must exist */
    if (hdrl_image_get_error(master_dark) == NULL) {
        cpl_msg_error(__func__, "The error is missing") ;
        hdrl_image_delete(master_dark) ;
        return NULL ;
    }

    /* Return  */
    return master_dark ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load the detlin coefficients
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to CR2RES_NB_DETECTORS)
  @return   An image list with the polynomial coefficients for each
              pixel of the wished detector. The returned object list
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
hdrl_imagelist * cr2res_io_load_DETLIN_COEFFS(
        const char  *   filename,
        int             detector)
{
    hdrl_imagelist      *   detlin_coeffs ;

    /* Check entries */
    if (filename == NULL) return NULL ;
    if (detector < 1 || detector > CR2RES_NB_DETECTORS) return NULL ;

    /* Check PRO.TYPE */
    if (cr2res_io_check_pro_type(filename, CR2RES_DETLIN_COEFFS_PROTYPE) != 1)
        return NULL ;

    /* Load */
    detlin_coeffs = cr2res_io_load_image_list(filename, detector) ;

    /* TODO - Error must exist */

    /* Return  */
    return detlin_coeffs ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load an hdrl image from a MASTER_FLAT
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to CR2RES_NB_DETECTORS)
  @return   A hdrl image or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
hdrl_image * cr2res_io_load_MASTER_FLAT(
        const char  *   filename,
        int             detector)
{
    hdrl_image          *   master_flat ;

    /* Check entries */
    if (filename == NULL) return NULL ;
    if (detector < 1 || detector > CR2RES_NB_DETECTORS) return NULL ;

    /* Check PRO.TYPE */
    if (cr2res_io_check_pro_type(filename, CR2RES_MASTER_FLAT_PROTYPE) != 1)
        return NULL ;

    /* Load */
    master_flat = cr2res_io_load_image(filename, detector) ;

    /* Error must exist */
    if (hdrl_image_get_error(master_flat) == NULL) {
        cpl_msg_error(__func__, "The error is missing") ;
        hdrl_image_delete(master_flat) ;
        return NULL ;
    }

    /* Return  */
    return master_flat ;
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
    cpl_table           *   trace_wave_tab ;

     /* Check entries */
    if (filename == NULL) return NULL ;
    if (detector < 1 || detector > CR2RES_NB_DETECTORS) return NULL ;

    /* Check PRO.TYPE */
    if (cr2res_io_check_pro_type(filename, CR2RES_TW_PROTYPE) != 1)
        return NULL ;

    /* Load the table */
    trace_wave_tab = cr2res_load_table(filename, detector, -1, -1) ;

    /* Return  */
    return trace_wave_tab ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load an hdrl image from a SLIT MODEL
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to CR2RES_NB_DETECTORS)
  @return   A hdrl image or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
hdrl_image * cr2res_io_load_SLIT_MODEL(
        const char  *   filename,
        int             detector) 
{
    hdrl_image          *   slit_model ;

     /* Check entries */
    if (filename == NULL) return NULL ;
    if (detector < 1 || detector > CR2RES_NB_DETECTORS) return NULL ;

    /* Check PRO.TYPE */
    if (cr2res_io_check_pro_type(filename, CR2RES_SLIT_MODEL_PROTYPE) != 1)
        return NULL ;

    /* Load */
    slit_model = cr2res_io_load_image(filename, detector) ;

    /* Error must exist */
    if (hdrl_image_get_error(slit_model) == NULL) {
        cpl_msg_error(__func__, "The error is missing") ;
        hdrl_image_delete(slit_model) ;
        return NULL ;
    }

    /* Return  */
    return slit_model ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load an hdrl image from a TRACE_MAP
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to CR2RES_NB_DETECTORS)
  @return   A hdrl image or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
hdrl_image * cr2res_io_load_TRACE_MAP(
        const char  *   filename,
        int             detector)
{
    hdrl_image          *   trace_map ;

     /* Check entries */
    if (filename == NULL) return NULL ;
    if (detector < 1 || detector > CR2RES_NB_DETECTORS) return NULL ;

    /* Check PRO.TYPE */
    if (cr2res_io_check_pro_type(filename, CR2RES_TRACE_MAP_PROTYPE) != 1)
        return NULL ;

    /* Load */
    trace_map = cr2res_io_load_image(filename, detector) ;

    /* Return  */
    return trace_map ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load an hdrl image from a WAVE_MAP
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to CR2RES_NB_DETECTORS)
  @return   A hdrl image or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
hdrl_image * cr2res_io_load_WAVE_MAP(
        const char  *   filename,
        int             detector)
{
    hdrl_image          *   wave_map ;

    /* Check entries */
    if (filename == NULL) return NULL ;
    if (detector < 1 || detector > CR2RES_NB_DETECTORS) return NULL ;

    /* Check PRO.TYPE */
    if (cr2res_io_check_pro_type(filename, CR2RES_WAVE_MAP_PROTYPE) != 1)
        return NULL ;

    /* Load */
    wave_map = cr2res_io_load_image(filename, detector) ;

    /* Return  */
    return wave_map ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load an hdrl image from a SLIT_CURV_MAP
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to CR2RES_NB_DETECTORS)
  @return   A hdrl image or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
hdrl_image * cr2res_io_load_SLIT_CURV_MAP(
        const char  *   filename,
        int             detector)
{
    hdrl_image          *   slit_curv_map ;

    /* Check entries */
    if (filename == NULL) return NULL ;
    if (detector < 1 || detector > CR2RES_NB_DETECTORS) return NULL ;

    /* Check PRO.TYPE */
    if (cr2res_io_check_pro_type(filename, CR2RES_SLIT_CURV_MAP_PROTYPE) != 1)
        return NULL ;

    /* Load */
    slit_curv_map = cr2res_io_load_image(filename, detector) ;

    /* Return  */
    return slit_curv_map ;
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
    cpl_table           *   slit_curv_tab ;

    /* Check entries */
    if (filename == NULL) return NULL ;
    if (detector < 1 || detector > CR2RES_NB_DETECTORS) return NULL ;

    /* Check PRO.TYPE */
    if (cr2res_io_check_pro_type(filename, CR2RES_SLIT_CURV_PROTYPE) != 1)
        return NULL ;

    /* Load the table */
    slit_curv_tab = cr2res_load_table(filename, detector, -1, -1) ;

    /* Return  */
    return slit_curv_tab ;
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
    cpl_table           *   extract_1D_tab ;

    /* Check entries */
    if (filename == NULL) return NULL ;
    if (detector < 1 || detector > CR2RES_NB_DETECTORS) return NULL ;

    /* Check PRO.TYPE */
    if (cr2res_io_check_pro_type(filename, CR2RES_EXTRACT_1D_PROTYPE) != 1)
        return NULL ;

    /* Load the table */
    extract_1D_tab = cr2res_load_table(filename, detector, -1, -1) ;

    /* Return  */
    return extract_1D_tab ;
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
    cpl_table           *   spliced_1D_tab ;

    /* Check entries */
    if (filename == NULL) return NULL ;
    if (detector < 1 || detector > CR2RES_NB_DETECTORS) return NULL ;

    /* Check PRO.TYPE */
    if (cr2res_io_check_pro_type(filename, CR2RES_SPLICED_1D_PROTYPE) != 1)
        return NULL ;

    /* Load the table */
    spliced_1D_tab = cr2res_load_table(filename, detector, -1, -1) ;

    /* Return  */
    return spliced_1D_tab ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load a table from a EXTRACT_2D
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to CR2RES_NB_DETECTORS)
  @return   A table or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_io_load_EXTRACT_2D(
        const char  *   filename,
        int             detector)
{
    cpl_table           *   extract_2D_tab ;

    /* Check PRO.TYPE */
    if (cr2res_io_check_pro_type(filename, CR2RES_EXTRACT_2D_PROTYPE) != 1)
        return NULL ;

    /* Load the table */
    extract_2D_tab = cr2res_load_table(filename, detector, -1, -1) ;

    /* Return  */
    return extract_2D_tab ;
}

/*----------------------------------------------------------------------------*/
/*---------------------       SAVING FUNCTIONS       -------------------------*/
/*----------------------------------------------------------------------------*/
/* TODO ? set frame levels to mark temp, intermediate and final frames?       */

/*----------------------------------------------------------------------------*/
/**
  @brief    Save EMISSION_LINES file
  @param    filename    The file name
  @param   	table		The table to save
  @param    parlist     The recipe input parameters
  @param    set         The recipe input frames
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_EMISSION_LINES(
        const char              *   filename,
        cpl_table               *   out_table,
        const cpl_parameterlist *   parlist,
        cpl_frameset            *   set,
        const char              *   recipe)
{
    cpl_propertylist    *   plist ;

    plist = cpl_propertylist_new();
    cpl_propertylist_append_string(plist, CR2RES_HEADER_INSTRUMENT, "CR2RES") ;
    cpl_propertylist_append_string(plist, CPL_DFS_PRO_CATG,
            CR2RES_EMISSION_LINES_PROCATG) ;
    cpl_propertylist_append_string(plist, CPL_DFS_PRO_TYPE,
            CR2RES_PROTYPE_CATALOG) ;

    if (cpl_dfs_save_table(set, NULL, parlist, set, NULL, out_table,
                NULL, recipe, plist, NULL,
                PACKAGE "/" PACKAGE_VERSION, filename) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot save the table") ;
        cpl_propertylist_delete(plist) ;
        return -1 ;
    }
    cpl_propertylist_delete(plist) ;

    /* Return */
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a MASTER_DARK
  @param    filename    The file name
  @param    allframes   The recipe input frames
  @param    inframes    The recipe used input frames
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
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  master_darks,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    return cr2res_io_save_image(filename, allframes, inframes, parlist,
            master_darks, qc_list, ext_plist, CPL_TYPE_DOUBLE, recipe,
            procatg, CR2RES_MASTER_DARK_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a DETLIN COEFFS
  @param    filename    The file name
  @param    allframes   The recipe input frames
  @param    inframes    The recipe used input frames
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
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        hdrl_imagelist          **  coeffs,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    return cr2res_io_save_imagelist(filename, allframes, inframes, parlist,
            coeffs, qc_list, ext_plist, CPL_TYPE_DOUBLE, recipe,
            procatg, CR2RES_DETLIN_COEFFS_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a BPM
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    inframes    The recipe used input frames
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
        cpl_frameset            *   inframes,
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
    ret = cr2res_io_save_image(filename, allframes, inframes, parlist,
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
  @param    inframes    The recipe used input frames
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
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  calib_collapsed,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    return cr2res_io_save_image(filename, allframes, inframes, parlist,
            calib_collapsed, qc_list, ext_plist, CPL_TYPE_DOUBLE, recipe,
            procatg, CR2RES_CALIBRATED_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a MASTER_FLAT
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    inframes    The recipe used input frames
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
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  master_flats,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    return cr2res_io_save_image(filename, allframes, inframes, parlist,
            master_flats, qc_list, ext_plist, CPL_TYPE_DOUBLE, recipe,
            procatg, CR2RES_MASTER_FLAT_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a TRACE_WAVE
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    inframes    The recipe used input frames
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
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    return cr2res_io_save_table(filename, allframes, inframes, parlist, tables,
            qc_list, ext_plist, recipe, procatg, CR2RES_TW_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a LINES_DIAGNOSICS
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    inframes    The recipe used input frames
  @param    parlist     The recipe input parameters
  @param    tables      The tables to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    ext_plist   The extensions property lists
  @param    procatg     The PRO CATG value
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_LINES_DIAGNOSTICS(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    return cr2res_io_save_table(filename, allframes, inframes, parlist, tables,
            qc_list, ext_plist, recipe, procatg,
            CR2RES_LINES_DIAGNOSTICS_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a 1D extracted spectrum   
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    inframes    The recipe used input frames
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
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    return cr2res_io_save_table(filename, allframes, inframes, parlist, tables,
            qc_list, ext_plist, recipe, procatg, CR2RES_EXTRACT_1D_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a SLIT_FUNC
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    inframes    The recipe used input frames
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
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  slit_func,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    return cr2res_io_save_table(filename, allframes, inframes, parlist, 
            slit_func, qc_list, ext_plist, recipe, procatg, 
            CR2RES_SLIT_FUNC_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a SLIT_MODEL
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    inframes    The recipe used input frames
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
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  data,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    return cr2res_io_save_image(filename, allframes, inframes, parlist,
            data, qc_list, ext_plist, CPL_TYPE_DOUBLE, recipe, procatg, 
            CR2RES_SLIT_MODEL_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a COMBINED
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    inframes    The recipe used input frames
  @param    parlist     The recipe input parameters
  @param    data        The data images to save (DATA and ERROR per detector)
  @param    qc_list     The QC parameters
  @param    ext_plist   The extensions property lists
  @param    procatg     The PRO CATG value
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_COMBINED(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  data,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    return cr2res_io_save_image(filename, allframes, inframes, parlist,
            data, qc_list, ext_plist, CPL_TYPE_DOUBLE, recipe, procatg, 
            CR2RES_COMBINED_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a TRACE_MAP
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    inframes    The recipe used input frames
  @param    parlist     The recipe input parameters
  @param    data        The data images to save (DATA and ERROR per detector)
  @param    qc_list     The QC parameters
  @param    ext_plist   The extensions property lists
  @param    procatg     The PRO CATG value
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_TRACE_MAP(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  data,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    return cr2res_io_save_image(filename, allframes, inframes, parlist,
            data, qc_list, ext_plist, CPL_TYPE_DOUBLE, recipe,
            procatg, CR2RES_TRACE_MAP_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a WAVE_MAP
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    inframes    The recipe used input frames
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
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  data,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    return cr2res_io_save_image(filename, allframes, inframes, parlist,
            data, qc_list, ext_plist, CPL_TYPE_DOUBLE, recipe,
            procatg, CR2RES_WAVE_MAP_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a SLIT_CURV_MAP
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    inframes    The recipe used input frames
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
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  data,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    return cr2res_io_save_image(filename, allframes, inframes, parlist,
            data, qc_list, ext_plist, CPL_TYPE_DOUBLE, recipe,
            procatg, CR2RES_SLIT_CURV_MAP_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a SLIT_CURV
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    inframes    The recipe used input frames
  @param    parlist     The recipe input parameters
  @param    tables      The tables to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    ext_plist   The extensions property lists
  @param    procatg     The PRO CATG value
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_SLIT_CURV(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    return cr2res_io_save_table(filename, allframes, inframes, parlist, tables,
            qc_list, ext_plist, recipe, procatg, CR2RES_SLIT_CURV_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a SPLICED_1D
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    inframes    The recipe used input frames
  @param    parlist     The recipe input parameters
  @param    spliced_1d  The table to save
  @param    qc_list     The QC parameters
  @param    ext_plist   The extension property list
  @param    procatg     The PRO CATG value
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_SPLICED_1D(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        cpl_table               *   spliced_1d,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        *   ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    return cr2res_io_save_one_table(filename, allframes, inframes, parlist, 
            spliced_1d, qc_list, ext_plist, "SPLICED", recipe, procatg, 
            CR2RES_SPLICED_1D_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a 2D extracted spectrum   
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    inframes    The recipe used input frames
  @param    parlist     The recipe input parameters
  @param    tables      The tables to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    ext_plist   The extensions property lists
  @param    procatg     The PRO CATG value
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_EXTRACT_2D(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    return cr2res_io_save_table(filename, allframes, inframes, parlist, tables,
            qc_list, ext_plist, recipe, procatg, CR2RES_EXTRACT_2D_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a Polarymetry spectrum   
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    inframes    The recipe used input frames
  @param    parlist     The recipe input parameters
  @param    tables      The tables to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    ext_plist   The extensions property lists
  @param    procatg     The PRO CATG value
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_POL_SPEC(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe)
{
    return cr2res_io_save_table(filename, allframes, inframes, parlist, tables,
            qc_list, ext_plist, recipe, procatg, CR2RES_POL_SPEC_PROTYPE) ;
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a multi extension table
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    inframes    The recipe used input frames
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
        cpl_frameset            *   inframes,
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
        if (cpl_dfs_save_table(allframes, NULL, parlist, inframes, NULL,
                    tab[0], ext_head, recipe, pro_list, NULL,
                    PACKAGE "/" PACKAGE_VERSION, filename) != CPL_ERROR_NONE) {
            cpl_msg_error(__func__, "Cannot save the first extension table") ;
            cpl_propertylist_delete(ext_head) ;
            cpl_propertylist_delete(pro_list) ;
            return -1 ;
        }
    } else {
        if (cpl_dfs_save_propertylist(allframes, NULL, parlist,
                    inframes, NULL, recipe, pro_list, NULL, PACKAGE "/"
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
  @brief    Save a single extension table
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    inframes    The recipe used input frames
  @param    parlist     The recipe input parameters
  @param    tab         The table to save
  @param    qc_list     The QC parameters
  @param    ext_plist   The extensions property list
  @param    recipe      The recipe name
  @param    extname     The extension name
  @param    procatg     PRO.CATG
  @param    protype     PRO.TYPE
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
static int cr2res_io_save_one_table(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        cpl_table               *   tab,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        *   ext_plist,
        const char              *   extname,
        const char              *   recipe,
        const char              *   procatg,
        const char              *   protype)
{
    cpl_propertylist    *   pro_list ;
    cpl_propertylist    *   ext_head ;
    int                     i ;

    /* Test entries */
    if (allframes == NULL || filename == NULL || ext_plist == NULL ||
            extname == NULL) return -1 ;

    /* Add the PRO keys */
    if (qc_list != NULL) pro_list = cpl_propertylist_duplicate(qc_list) ;
    else pro_list = cpl_propertylist_new() ;

    /* Add PRO Keys */
    cpl_propertylist_append_string(pro_list, CPL_DFS_PRO_CATG, procatg) ;
    cpl_propertylist_append_string(pro_list, CPL_DFS_PRO_TYPE, protype) ;

    /* Create the first extension header */
	ext_head = cpl_propertylist_duplicate(ext_plist);
	cpl_propertylist_erase(ext_head, "EXTNAME");
    cpl_propertylist_update_string(ext_head, "EXTNAME", extname) ;

    /* Save the first extension */
	if (cpl_dfs_save_table(allframes, NULL, parlist, inframes, NULL,
				tab, ext_head, recipe, pro_list, NULL,
				PACKAGE "/" PACKAGE_VERSION, filename) != CPL_ERROR_NONE) {
		cpl_msg_error(__func__, "Cannot save the first extension table") ;
		cpl_propertylist_delete(ext_head) ;
		cpl_propertylist_delete(pro_list) ;
		return -1 ;
	}
    cpl_propertylist_delete(ext_head) ;
    cpl_propertylist_delete(pro_list) ;

    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a multi extension images
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    inframes    The recipe used input frames
  @param    parlist     The recipe input parameters
  @param    data        The images to save (data and error per detector)
  @param    qc_list     The QC parameters
  @param    ext_plist   The extensions property lists
  @param    type        CPL_TYPE_DOUBLE, CPL_TYPE_INT,...
  @param    recipe      The recipe name
  @param    procatg     PRO.CATG
  @param    protype     PRO.TYPE
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
static int cr2res_io_save_image(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
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
    if (cpl_dfs_save_propertylist(allframes, NULL, parlist, inframes, NULL, 
                recipe, qclist_loc, NULL,
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
  @param    inframes    The recipe used input frames
  @param    parlist     The recipe input parameters
  @param    data        The imagelists to save (only data for the moment)
  @param    qc_list     The QC parameters
  @param    ext_plist   The extensions property lists
  @param    type        CPL_TYPE_DOUBLE, CPL_TYPE_INT,...
  @param    recipe      The recipe name
  @param    procatg     PRO.CATG
  @param    protype     PRO.TYPE
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
static int cr2res_io_save_imagelist(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        hdrl_imagelist          **  data,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        cpl_type                    type,
        const char              *   recipe,
        const char              *   procatg,
        const char              *   protype)
{
    cpl_propertylist    *   qclist_loc ;
    cpl_imagelist       *   list_data ;
    cpl_imagelist       *   list_noise ;
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
    if (cpl_dfs_save_propertylist(allframes, NULL, parlist, inframes, NULL,
                recipe, qclist_loc, NULL,
                PACKAGE "/" PACKAGE_VERSION, filename) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot save the empty primary HDU") ;
        cpl_propertylist_delete(qclist_loc) ;
        return -1 ;
    }
    /* Delete PRO LIST */
    cpl_propertylist_delete(qclist_loc) ;

    /* Save the extensions */
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {

        /* Extract separately data / noise */
        if (data[det_nr-1] == NULL) {
            list_data = list_noise = NULL ;
        } else {
            list_data = cr2res_hdrl_to_cpl_list(data[det_nr-1], 1) ;
            list_noise = cr2res_hdrl_to_cpl_list(data[det_nr-1], 0) ;
        }

        /* Prepare header */
        if (ext_plist[det_nr-1] == NULL) {
            qclist_loc = cpl_propertylist_new();
        } else {
            qclist_loc = cpl_propertylist_duplicate(ext_plist[det_nr-1]) ;
            cpl_propertylist_erase(qclist_loc, "EXTNAME");
        }

        /* Save the DATA */
        wished_extname = cr2res_io_create_extname(det_nr, 1) ;
        cpl_propertylist_prepend_string(qclist_loc, "EXTNAME", wished_extname) ;
        if (list_data == NULL) {
            cpl_propertylist_save(qclist_loc, filename, CPL_IO_EXTEND) ;
        } else {                        
            cpl_imagelist_save(list_data, filename, type, qclist_loc, 
                    CPL_IO_EXTEND) ;
            cpl_imagelist_delete(list_data) ;
        }
        cpl_free(wished_extname) ;

        /* Save the NOISE */
        wished_extname = cr2res_io_create_extname(det_nr, 0) ;
        cpl_propertylist_update_string(qclist_loc, "EXTNAME", wished_extname) ;
        if (list_noise == NULL) {
            cpl_propertylist_save(qclist_loc, filename, CPL_IO_EXTEND) ;
        } else {                        
            cpl_imagelist_save(list_noise, filename, type, qclist_loc, 
                    CPL_IO_EXTEND) ;
            cpl_imagelist_delete(list_noise) ;
        }
        cpl_propertylist_delete(qclist_loc) ;
        cpl_free(wished_extname) ;
    }

	return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief Create a cpl_imagelist from the data/error of a hdrl_imagelist
  @param    in      the hdrl_imageḻist
  @param    data    1 -> data, 0 -> error
  @return  the newly allocated cpl_imagelist
 */
/*----------------------------------------------------------------------------*/
static cpl_imagelist * cr2res_hdrl_to_cpl_list(
        const hdrl_imagelist    *   in, 
        int                         data) 
{
    cpl_imagelist       *   out ;
    const hdrl_image    *   cur_im ;
    const cpl_image     *   cur_cpl_im ;
    cpl_size            i ;

    /* Check Entries */
    if (in == NULL) return NULL ;
    if (data != 0 && data != 1) return NULL ;

    /* Create output list */
    out = cpl_imagelist_new() ;
    for (i=0 ; i< hdrl_imagelist_get_size(in) ; i++) {
        cur_im = hdrl_imagelist_get(in, i) ;
        if (data == 1) cur_cpl_im = hdrl_image_get_image_const(cur_im) ;
        if (data == 0) cur_cpl_im = hdrl_image_get_error_const(cur_im) ;

        cpl_imagelist_set(out, cpl_image_duplicate(cur_cpl_im),
                cpl_imagelist_get_size(out)) ;
    }

    /* Check output */
    if (hdrl_imagelist_get_size(in) != cpl_imagelist_get_size(out)) {
        cpl_imagelist_delete(out) ;
        return NULL ;
    }
    /* Return  */
    return out ;
}

static int cr2res_table_check_column(
        const cpl_table     *   tab,
        const char          *   col)
{
    int         ret ;
    if (!cpl_table_has_column(tab, col)) {
        cpl_msg_error(__func__, "Column %s is missing", col) ;
        ret = 1 ;
    } else ret = 0 ;
    return ret ;
}

static int cr2res_io_check_pro_type(
        const char  *   filename,
        const char  *   expected_protype)
{
    cpl_propertylist    *   plist ;
    const char          *   protype ;

    /* Check entries */
    if (filename == NULL) return -1 ;

    /* Check the PRO.TYPE */
    plist = cpl_propertylist_load(filename, 0) ;
    if (plist == NULL) return -1;
    protype = cr2res_pfits_get_protype(plist) ;
    if (strcmp(protype, expected_protype)) {
        cpl_msg_error(__func__, "Unexpected PRO.TYPE: %s != %s",
                protype, expected_protype) ;
        cpl_propertylist_delete(plist) ;
        return 0 ;
    }
    cpl_propertylist_delete(plist) ;
    return 1 ;
}

