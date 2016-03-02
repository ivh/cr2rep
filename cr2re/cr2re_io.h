

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

/*----------------------------------------------------------------------------*/
/**
  @brief    Load an image from a SLIT_ILLUM
  @param    filename    The FITS file name
  @param    detector    The wished detector (1 to 3)
  @param    data        1 for the data image, 0 for the error
  @return   A float image or NULL in error case. The returned object
              needs to be deallocated
 */
/*----------------------------------------------------------------------------*/
cpl_image * cr2res_io_load_SLIT_ILLUM(
        const char  *   filename,
        int             detector,
        int             data)

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

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a MASTER_DARK
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    data        The data images to save (1 per detector)
  @param    errors      The error images to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @param    pipe_id     PACKAGE "/" PACKAGE_VERSION
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_MASTER_DARK(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_imagelist           *   data,
        cpl_imagelist           *   errors,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id)

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a MASTER_BPM
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    data        The data images to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @param    pipe_id     PACKAGE "/" PACKAGE_VERSION
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_MASTER_BPM(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_imagelist           *   data,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id)

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a DETLIN_COEFFS
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    data        The data imagelists to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @param    pipe_id     PACKAGE "/" PACKAGE_VERSION
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_DETLIN_COEFFS(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_imagelist           **  data,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id)

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
  @param    pipe_id     PACKAGE "/" PACKAGE_VERSION
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
        const char              *   recipe,
        const char              *   pipe_id)

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a TRACE_OPEN
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    tables      The tables to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @param    pipe_id     PACKAGE "/" PACKAGE_VERSION
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_TRACE_OPEN(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id)

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a TRACE_DECKER
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    tables      The tables to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @param    pipe_id     PACKAGE "/" PACKAGE_VERSION
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_TRACE_DECKER(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id)

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a BLAZE
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    tables      The tables to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @param    pipe_id     PACKAGE "/" PACKAGE_VERSION
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_BLAZE(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id)

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a BLAZE_IMAGE
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    data        The data images to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @param    pipe_id     PACKAGE "/" PACKAGE_VERSION
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_BLAZE_IMAGE(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_imagelist           *   data,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id)

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a SLIT_ILLUM
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    data        The data images to save (1 per detector)
  @param    errors      The error images to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @param    pipe_id     PACKAGE "/" PACKAGE_VERSION
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_SLIT_ILLUM(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_imagelist           *   data,
        cpl_imagelist           *   errors,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id)

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a WAVE_MAP
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    data        The data images to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @param    pipe_id     PACKAGE "/" PACKAGE_VERSION
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_WAVE_MAP(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_imagelist           *   data,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id)

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a WAVE_SUB_ORDER
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    tables      The tables to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @param    pipe_id     PACKAGE "/" PACKAGE_VERSION
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_WAVE_SUB_ORDER(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id)

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a SLITPOS_MAP
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    data        The data images to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @param    pipe_id     PACKAGE "/" PACKAGE_VERSION
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_SLITPOS_MAP(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_imagelist           *   data,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id)

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a TILT_MAP
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    data        The data images to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @param    pipe_id     PACKAGE "/" PACKAGE_VERSION
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_TILT_MAP(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_imagelist           *   data,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id)

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a TILT_POLY
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    tables      The tables to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @param    pipe_id     PACKAGE "/" PACKAGE_VERSION
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_TILT_POLY(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id)

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a EXTRACT_1D
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    tables      The tables to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @param    pipe_id     PACKAGE "/" PACKAGE_VERSION
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_EXTRACT_1D(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id)

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a SPLICED_1D
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    tables      The tables to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @param    pipe_id     PACKAGE "/" PACKAGE_VERSION
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_SPLICED_1D(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id)

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a EXTRACT_2D
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    table       The table to save
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @param    pipe_id     PACKAGE "/" PACKAGE_VERSION
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_EXTRACT_2D(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               *   table,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id)

/*----------------------------------------------------------------------------*/
/**
  @brief    Save a EXTRACT_POL
  @param    filename    The FITS file name
  @param    allframes   The recipe input frames
  @param    parlist     The recipe input parameters
  @param    tables      The tables to save (1 per detector)
  @param    qc_list     The QC parameters
  @param    recipe      The recipe name
  @param    pipe_id     PACKAGE "/" PACKAGE_VERSION
  @return   0 if ok, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_io_save_EXTRACT_POL(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id)





