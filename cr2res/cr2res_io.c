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

#include "cr2res_io.h"
#include "cr2res_dfs.h"

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_io_save_image(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  data,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   procatg,
        const char              *   protype) ;
static int cr2res_io_save_table(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  slit_func,
        const cpl_propertylist  *   qc_list,
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
/*--------------------       LOADING FUNCTIONS       -------------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Load an image from a MASTER_DARK
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to 3)
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
    return NULL ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load an image from a MASTER_BPM
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to 3)
  @return   An integer type image or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_image * cr2res_io_load_MASTER_BPM(
        const char  *   filename,
        int             detector)
{
    return NULL ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load the detlin coefficients
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to 3)
  @return   A float image list with the polynomial coefficients for each
              pixel of the wished detector. The returned object list
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_imagelist * cr2res_io_load_DETLIN_COEFFS(
        const char  *   filename,
        int             detector)
{
        return NULL ;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Load an image from a MASTER_FLAT
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to 3)
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
        return NULL ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load a table from a TRACE_OPEN
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to 3)
  @return   A table or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_io_load_TRACE_OPEN(
        const char  *   filename,
        int             detector)
{
    const char          *   extname ;
    char                *   wished_extname ;
    int                     wished_ext_nb ;
    cpl_propertylist    *   pl ;
    int                     nb_ext, i ;

    /* Check entries */
    if (filename == NULL) return NULL ;
    if (detector < 1 || detector > CR2RES_NB_DETECTORS) return NULL ;

    /* Initialise */
    wished_ext_nb = -1 ;

    /* Create wished EXTNAME */
    wished_extname = cpl_sprintf("CHIP%d", detector) ;

    /* Get the number of extensions */
    nb_ext = cpl_fits_count_extensions(filename) ;

    /* Loop on the extensions */
    for (i=1 ; i<=nb_ext ; i++) {
        /* Get the header */
        pl = cpl_propertylist_load(filename, i) ;
        /* Read the ext EXTNAME */
        extname = cpl_propertylist_get_string(pl, "EXTNAME");

        /* Compare to the wished one */
        if (!strcmp(extname, wished_extname)) wished_ext_nb = i ;
        cpl_propertylist_delete(pl) ;
    }
    cpl_free(wished_extname) ;

    /* The wished extension was not found */
    if (wished_ext_nb < 0) return NULL ;

    return cpl_table_load(filename, wished_ext_nb, 1);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load a table from a TRACE_DECKER
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to 3)
  @param    decker_type [out] CR2RES_DECKER_1_3 or CR2RES_DECKER_2_4
  @return   A table or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_io_load_TRACE_DECKER(
        const char  *   filename,
        int             detector,
        cr2res_decker * decker_type)
{
        return NULL ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Load a table from a BLAZE
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to 3)
  @return   A table or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_io_load_BLAZE(
        const char  *   filename,
        int             detector)
{
        return NULL ;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Load an image from a BLAZE_IMAGE
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to 3)
  @return   A float image or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_image * cr2res_io_load_BLAZE_IMAGE(
        const char  *   filename,
        int             detector)
{
        return NULL ;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Load an image from a SLIT_MODEL
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to 3)
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
  @param    detector    The wished detector (1 to 3)
  @return   A float image or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_image * cr2res_io_load_WAVE_MAP(
        const char  *   filename,
        int             detector)
{
        return NULL ;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Load a table from a WAVE_SUB_ORDER
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to 3)
  @return   A table or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_io_load_WAVE_SUB_ORDER(
        const char  *   filename,
        int             detector)
{
        return NULL ;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Load an image from a SLITPOS_MAP
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to 3)
  @return   A float image or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_image * cr2res_io_load_SLITPOS_MAP(
        const char  *   filename,
        int             detector)
{
        return NULL ;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Load an image from a TILT_MAP
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to 3)
  @return   A float image or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_image * cr2res_io_load_TILT_MAP(
        const char  *   filename,
        int             detector)
{
        return NULL ;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Load a table from a TILT_POLY
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to 3)
  @return   A table or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_io_load_TILT_POLY(
        const char  *   filename,
        int             detector)
{
        return NULL ;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Load a table from a EXTRACT_1D
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to 3)
  @return   A table or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_io_load_EXTRACT_1D(
        const char  *   filename,
        int             detector)
{
        return NULL ;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Load a table from a SPLICED_1D
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to 3)
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
  @param    detector    The wished detector (1 to 3)
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

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a MASTER_DARK
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    master_darks  The data/error master darks (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_MASTER_DARK(
        cpl_frameset            *   allframes,
        const char              *   filename,
        cpl_frameset            *   used_frames,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  master_darks,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe)
{
    return cr2res_io_save_image(filename, allframes, parlist,
            master_darks, qc_list, recipe, CR2RES_MASTER_DARK_PROCATG, "") ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a MASTER_BPM
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    data        The data images to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_MASTER_BPM(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_imagelist           *   data,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe)
{
            return -1 ;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Save a DETLIN_COEFFS
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    data        The data imagelists to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_DETLIN_COEFFS(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_imagelist           **  data,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe)
{
            return -1 ;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Save a MASTER_FLAT
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    data        The data images to save (1 per detector)
  @param    errors      The error images to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_MASTER_FLAT(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_imagelist           *   data,
        cpl_imagelist           *   errors,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe)
{
            return -1 ;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Save a TRACE_OPEN
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    tables      The tables to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_TRACE_OPEN(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe)
{
            return -1 ;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Save a TRACE_DECKER
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    tables      The tables to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_TRACE_DECKER(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe)
{
            return -1 ;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Save a BLAZE
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    tables      The tables to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_BLAZE(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe)
{
            return -1 ;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Save a BLAZE_IMAGE
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    data        The data images to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_BLAZE_IMAGE(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_imagelist           *   data,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe)
{
            return -1 ;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Save a SLIT_FUNC
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    data        The tables to save (1 per detector)
  @param    qc_list     The QC parameters
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
        const char              *   recipe)
{
    return cr2res_io_save_table(filename, allframes, parlist, slit_func, 
            qc_list, recipe, CR2RES_SLIT_FUNC_PROCATG,
            CR2RES_SLIT_FUNC_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a SLIT_MODEL
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    data        The data images to save (DATA and ERROR per detector)
  @param    qc_list     The QC parameters
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
        const char              *   recipe)
{
    return cr2res_io_save_image(filename, allframes, parlist,
            data, qc_list, recipe, CR2RES_SLIT_MODEL_PROCATG,
            CR2RES_SLIT_MODEL_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a WAVE_MAP
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    data        The data images to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_WAVE_MAP(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_imagelist           *   data,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe)
{
            return -1 ;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Save a WAVE_SUB_ORDER
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    tables      The tables to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_WAVE_SUB_ORDER(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe)
{
            return -1 ;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Save a SLITPOS_MAP
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    data        The data images to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_SLITPOS_MAP(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_imagelist           *   data,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe)
{
            return -1 ;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Save a TILT_MAP
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    data        The data images to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_TILT_MAP(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_imagelist           *   data,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe)
{
            return -1 ;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Save a TILT_POLY
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    tables      The tables to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_TILT_POLY(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe)
{
            return -1 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a EXTRACT_1D
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    tables      The tables to save (1 per detector)
  @param    qc_list     The QC parameters
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
        const char              *   recipe)
{
    return cr2res_io_save_table(filename, allframes, parlist, tables, 
            qc_list, recipe, CR2RES_EXTRACT_1D_PROCATG,
            CR2RES_EXTRACT_1D_PROTYPE) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a SPLICED_1D
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    tables      The tables to save (1 per detector)
  @param    qc_list     The QC parameters
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
        const char              *   recipe,
        const char              *   procatg,
        const char              *   protype)
{
    cpl_propertylist    *   pro_list ;
    cpl_propertylist    *   ext_head ;
    char                    sval[16] ;
    int                     i ;

    /* Test entries */
    if (allframes == NULL || filename == NULL) return -1 ;

    /* Add the PRO keys */
    if (qc_list != NULL) pro_list = cpl_propertylist_duplicate(qc_list) ;
    else pro_list = cpl_propertylist_new() ;
        
    /* Add PRO Keys */
    cpl_propertylist_append_string(pro_list, CPL_DFS_PRO_CATG, procatg) ;
    cpl_propertylist_append_string(pro_list, CPL_DFS_PRO_TYPE, protype) ;

    /* Create the first extension header */
    ext_head = cpl_propertylist_new() ;
    cpl_propertylist_append_string(ext_head, "EXTNAME", "CHIP1") ;

    /* Save the first extension */
    if (cpl_dfs_save_table(allframes, NULL, parlist, allframes, NULL,
                tab[0], ext_head, recipe, pro_list, NULL,
                PACKAGE "/" PACKAGE_VERSION, filename) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot save the first extension table") ;
        cpl_propertylist_delete(ext_head) ;
        cpl_propertylist_delete(pro_list) ;
        return -1 ;
    }
    cpl_propertylist_delete(ext_head) ;
    cpl_propertylist_delete(pro_list) ;

    /* Save the extensions */
    for (i=1 ; i<CR2RES_NB_DETECTORS ; i++) {
        ext_head = cpl_propertylist_new() ;
        sprintf(sval, "CHIP%d", i+1) ;
        cpl_propertylist_prepend_string(ext_head, "EXTNAME", sval) ;
        cpl_table_save(tab[i], NULL, ext_head, filename, CPL_IO_EXTEND) ;
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
        const char              *   recipe,
        const char              *   procatg,
        const char              *   protype)
{
    cpl_propertylist    *   qclist_loc ;
    char          		*   sval ;
    int                     ext ;

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
                CPL_BPP_IEEE_FLOAT, recipe, qclist_loc, NULL,
                PACKAGE "/" PACKAGE_VERSION, filename) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot save the empty primary HDU") ;
        cpl_propertylist_delete(qclist_loc) ;
        return -1 ;
    }
    /* Delete PRO LIST */
    cpl_propertylist_delete(qclist_loc) ;

    /* Save the extensions */
    for (ext=1 ; ext<=CR2RES_NB_DETECTORS ; ext++) {
        /* Save the DATA */
        qclist_loc = cpl_propertylist_new() ;
        sval = cpl_sprintf("CHIP%d", ext) ;
        cpl_propertylist_prepend_string(qclist_loc, "EXTNAME", sval) ;
        cpl_image_save(hdrl_image_get_image(data[ext-1]),
                filename, CPL_BPP_IEEE_FLOAT, qclist_loc, CPL_IO_EXTEND) ;
        cpl_propertylist_delete(qclist_loc) ;
        cpl_free(sval) ;

        /* Save the NOISE */
        qclist_loc = cpl_propertylist_new() ;
        sval = cpl_sprintf("CHIP%dERR", ext) ;
        cpl_propertylist_prepend_string(qclist_loc, "EXTNAME", sval) ;
        cpl_image_save(hdrl_image_get_error(data[ext-1]),
                filename, CPL_BPP_IEEE_FLOAT, qclist_loc, CPL_IO_EXTEND) ;
        cpl_propertylist_delete(qclist_loc) ;
        cpl_free(sval) ;
    }

	return 0 ;
}
