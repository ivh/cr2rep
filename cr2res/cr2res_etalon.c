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
#include <math.h>
#include "cr2res_dfs.h"
#include "cr2res_etalon.h"
#include "cr2res_wave.h"
#include "cr2res_pfits.h"

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

#define SPEED_OF_LIGHT 299792.458
#define modulo(x, y) ((x) - ((y) * trunc((x)/(y))))
#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

#define CR2RES_ETALON_GAP_H 9.9950e6
#define CR2RES_ETALON_GAP_J 9.9900e6
#define CR2RES_ETALON_GAP_K 9.9999e6
#define CR2RES_ETALON_GAP_Y 9.9913e6

#define CR2RES_ETALON_DEGREE_H 2
#define CR2RES_ETALON_DEGREE_J 1
#define CR2RES_ETALON_DEGREE_K 1
#define CR2RES_ETALON_DEGREE_Y 2

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

#ifdef CR2RES_UNUSED
static cpl_mask * cr2res_etalon_binary_image(const cpl_image *) ;
#endif
static cpl_vector * cr2res_etalon_get_peaks_gaussian(
    cpl_bivector        *  spectra,
    cpl_bivector        *  spectra_err,
    cpl_polynomial      *  wavesol_init,
    cpl_array           *  wavesol_init_err,
    cpl_vector          *  peaks,
    int                    display,
    cpl_vector         **  sigma,
    cpl_vector         **  heights,
    cpl_vector         **  fit_error);

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_etalon       Etalon related
 */
/*----------------------------------------------------------------------------*/

/**@{*/
#ifdef CR2RES_UNUSED
/*----------------------------------------------------------------------------*/
/**
  @brief    Get the etalon Image, and detects the blobs
  @param    in      The Etalon image
  @return   The labels int image or NULL in error case
 */
/*----------------------------------------------------------------------------*/
cpl_image * cr2res_etalon_computation(const cpl_image * in)
{
    cpl_mask        *   mask ;
    cpl_size            nlabels ;
    cpl_image       *   labels ;
    cpl_apertures   *   aperts ;

    /* Compute Binary image */
    if ((mask = cr2res_etalon_binary_image(in)) == NULL) {
        cpl_msg_error(__func__, "Cannot Binarise the image") ;
        return NULL ;
    }
        
    /* Debugging */
    cpl_mask_save(mask, "mask.fits", NULL, CPL_IO_CREATE) ;

    /* Labelise the different detected apertures */
    if ((labels = cpl_image_labelise_mask_create(mask, &nlabels))==NULL) {
        cpl_msg_error(cpl_func, "Cannot Labelise") ;
        cpl_mask_delete(mask) ;
        return NULL ;
    }
    cpl_mask_delete(mask) ;

    cpl_msg_debug(__func__, "Number of Apertures: %"CPL_SIZE_FORMAT, nlabels) ;

    /* Create the detected apertures list */
    if ((aperts = cpl_apertures_new_from_image(in, labels)) == NULL) {
        cpl_msg_error(cpl_func, "Cannot Compute the apertures") ;
        cpl_image_delete(labels) ;
        return NULL ;
    }

    /* Display Apertures */
    cpl_apertures_dump(aperts, stdout) ;

    /* Free and return */
    cpl_apertures_delete(aperts) ;
    return labels ;
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Create a mask identifying the blogs (all separated)
  @param    in      The Etalon image
  @return   The mask
 */
/*----------------------------------------------------------------------------*/
static cpl_mask * cr2res_etalon_binary_image(const cpl_image * in)
{
    cpl_mask        *   mask ;
    cpl_vector      *   collapsed_vec ;
    cpl_image       *   collapsed_im ;
    cpl_vector      *   local_maxima ;
    double          *   pmax ;
    int i;

    /* Check Entries */
    if (in == NULL) return NULL; 

    /*************************/
    /* Get Local Maxima list */
    /*************************/
    /* Collapse in y */
    collapsed_im = cpl_image_collapse_create(in, 0) ;
    
    /* Smooth */
    /* TODO */

    /* Convert as vector */
    collapsed_vec = cpl_vector_new_from_image_row(collapsed_im, 1) ;
    cpl_image_delete(collapsed_im) ;

    /* Debugging */
    cpl_plot_vector(NULL, "w lines", NULL, collapsed_vec) ;

    /* Store local maxima positions of the plot */
    local_maxima = cr2res_etalon_get_maxpos(collapsed_vec) ;
    pmax = cpl_vector_get_data(local_maxima) ;
    cpl_vector_delete(collapsed_vec) ;

    /* Create the output mask  */
    mask = cpl_mask_new(cpl_image_get_size_x(in), cpl_image_get_size_y(in)) ;
   
    /* Loop on the Maxima positions and isolate the blob */
    for (i = 0; i < cpl_vector_get_size(local_maxima); i++) {
        cpl_image *blob_image;
        cpl_mask *blob_mask;
        double threshold;
        int j, start_x, end_x;

        /* Get start_x end_x */
        start_x = 1;
        end_x = cpl_image_get_size_x(in) ;
        if (i>0) 
            start_x = (int)(pmax[i-1] + (pmax[i]-pmax[i-1])/2.) ;
        if (i<cpl_vector_get_size(local_maxima)-1)
            end_x = (int)(pmax[i] + (pmax[i+1]-pmax[i])/2.) ;

        cpl_msg_debug(__func__, "Extract %d -> %d", start_x, end_x) ;

        /* Extract the current blob */
        blob_image = cpl_image_extract(in, start_x, 1, end_x,
                cpl_image_get_size_y(in)) ;

        /* Compute the threshold */
        /* TODO */
        threshold = cpl_image_get_mean(blob_image) ;
        
        cpl_msg_debug(__func__, "Threshold: %g", threshold) ;

        /* Binarise the image */
        blob_mask = cpl_mask_threshold_image_create(blob_image, threshold, 
                DBL_MAX) ;
        cpl_image_delete(blob_image); 

        /* Set the left column to 0 to separate from the neighbor */
        for (j=0 ; j<cpl_mask_get_size_y(blob_mask) ; j++) 
            cpl_mask_set(blob_mask, 1, j+1, CPL_BINARY_0) ;

        /* Fill the Binary with the current blob */
        cpl_mask_copy(mask, blob_mask, start_x, 1); 

        cpl_mask_delete(blob_mask) ;
    }
    /* Free and return */
    cpl_vector_delete(local_maxima) ;

    return mask ;
}

/**@}*/
/*----------------------------------------------------------------------------*/
/**
  @brief    Detect maxima from a 1d periodic signal and store their positions
  @param    in      The 1d signal as a vector
  @return   The vector with the maxima positions
 */
/*----------------------------------------------------------------------------*/
cpl_vector * cr2res_etalon_get_maxpos(const cpl_vector * in)
{
    const double    *   pin ;
    cpl_vector      *   maxima_pos ;
    double          *   pmax ;
    int                 i, nb_max ;

    /* Check Entries */
    if (in==NULL) return NULL ;

    /* Initialise */
    pin = cpl_vector_get_data_const(in) ;

    /* Count the number of Max positions */
    nb_max = 0 ;
    for (i=1 ; i<cpl_vector_get_size(in)-1 ; i++) {
        if (pin[i] > pin[i+1] && pin[i] > pin[i-1]) nb_max++ ;
    }

    /* Create the output vector */
    maxima_pos = cpl_vector_new(nb_max) ;
    pmax = cpl_vector_get_data(maxima_pos) ;
    nb_max = 0 ;
    for (i=1 ; i<cpl_vector_get_size(in)-1 ; i++) {
        if (pin[i] > pin[i+1] && pin[i] > pin[i-1]) {
            /* Store the Max X position */
            pmax[nb_max] = i+1 ;
            nb_max++ ;
        }
    }
    return maxima_pos ;
}

#endif
/*----------------------------------------------------------------------------*/
/**
  @brief    Find local maxima in a 1D array. This function finds all local 
            maxima in a 1D array and returns the indices for their edges and 
            midpoints (rounded down for even plateau sizes).
  @param    in      The 1d signal as a vector
  @param    left_edges [out] left edge positions of the peaks
  @param    right_edges [out] right edge positions of the peaks
  @return   The vector with the maxima positions

  The left_edge and right_edge vectors need to be deleted afterwards
  Note that the edges here refer to points with exactly the same value, not
  the width of the gaussian peak

 */
/*----------------------------------------------------------------------------*/
cpl_vector * cr2res_etalon_get_local_maxima(
    const cpl_vector * in, 
    cpl_vector ** left_edges, 
    cpl_vector ** right_edges){
    /*
    Find local maxima in a 1D array.
    This function finds all local maxima in a 1D array and returns the indices
    for their edges and midpoints (rounded down for even plateau sizes).
    Parameters
    ----------
    x : ndarray
        The array to search for local maxima.
    Returns
    -------
    midpoints : ndarray
        Indices of midpoints of local maxima in `x`.
    left_edges : ndarray
        Indices of edges to the left of local maxima in `x`.
    right_edges : ndarray
        Indices of edges to the right of local maxima in `x`.
    Notes
    -----
    - Compared to `argrelmax` this function is significantly faster and can
      detect maxima that are more than one sample wide. However this comes at
      the cost of being only applicable to 1D arrays.
    - A maxima is defined as one or more samples of equal value that are
      surrounded on both sides by at least one smaller sample.
    .. versionadded:: 1.1.0
    */
    
    cpl_vector * midpoints;
    int m, i, i_ahead, i_max;

    int width = cpl_vector_get_size(in);

    // Preallocate, there can't be more maxima than half the size of `x`
    midpoints = cpl_vector_new(width / 2);
    *left_edges = cpl_vector_new(width / 2);
    *right_edges = cpl_vector_new(width / 2);
    m = 0;  // Pointer to the end of valid area in allocated arrays

    i = 1;  // Pointer to current sample, first one can't be maxima
    i_max = width - 1;  // Last sample can't be maxima
    while (i < i_max){
        // Test if previous sample is smaller
        if (cpl_vector_get(in, i - 1) < cpl_vector_get(in, i)){
            i_ahead = i + 1;  // Index to look ahead of current sample

            // Find next sample that is unequal to x[i]
            while (i_ahead < i_max && cpl_vector_get(in, i_ahead) == cpl_vector_get(in, i)){
                i_ahead += 1;
            }

            // Maxima is found if next unequal sample is smaller than x[i]
            if (cpl_vector_get(in, i_ahead) < cpl_vector_get(in, i)){
                cpl_vector_set(*left_edges, m, i);
                cpl_vector_set(*right_edges, m, i_ahead - 1);
                cpl_vector_set(midpoints, m, (i + i_ahead - 1) / 2);
                m += 1;
                // Skip samples that can't be maximum
                i = i_ahead;
            }
        }
        i++;
    }

    if (m <1) {
        cpl_msg_warning(__func__, "Cannot find maxima");
        cpl_vector_delete(midpoints);
        cpl_vector_delete(*left_edges);
        cpl_vector_delete(*right_edges);
        return NULL;
    }
    // Keep only valid part of array memory.
    cpl_vector_set_size(midpoints, m);

    cpl_vector_set_size(*left_edges, m);
    cpl_vector_set_size(*right_edges, m);

    return midpoints;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Selects the highest peaks that are atleast distance pixel apart
  @param    peaks          The 1d signal as a vector
  @param    peak_heights   The heights (priority) of each peak
  @param    distance       The minimum distance between peaks to keep
  @return   A boolean vector where 1 means to keep the peak, or 0 to discard it
 */
/*----------------------------------------------------------------------------*/
cpl_vector * cr2res_etalon_select_by_peak_distance(const cpl_vector * peaks,
                             const cpl_vector * peak_heights,
                             float distance)
{
    /*
    Evaluate which peaks fulfill the distance condition.
    Parameters
    ----------
    peaks : ndarray
        Indices of peaks in `vector`.
    priority : ndarray
        An array matching `peaks` used to determine priority of each peak. A
        peak with a higher priority value is kept over one with a lower one.
    distance : np.float64
        Minimal distance that peaks must be spaced.
    Returns
    -------
    keep : ndarray[bool]
        A boolean mask evaluating to true where `peaks` fulfill the distance
        condition.
    Notes
    -----
    Declaring the input arrays as C-contiguous doesn't seem to have performance
    advantages.
    .. versionadded:: 1.1.0
    */
    

    cpl_vector * keep;
    cpl_vector * priority_to_position;
    cpl_bivector * bivector_tmp;
    int i, k, peaks_size, distance_;

    peaks_size = cpl_vector_get_size(peaks);
    // Round up because actual peak distance can only be natural number
    distance_ = distance;
    // Prepare array of flags
    keep = cpl_vector_new(peaks_size);
    cpl_vector_fill(keep, 1);

    // Create map from `i` (index for `peaks` sorted by `priority`) to `j` (index
    // for `peaks` sorted by position). This allows to iterate `peaks` and `keep`
    // with `j` by order of `priority` while still maintaining the ability to
    // step to neighbouring peaks with (`j` + 1) or (`j` - 1).
    
    bivector_tmp = cpl_bivector_new(peaks_size);
    priority_to_position = cpl_bivector_get_x(bivector_tmp);

    for (i = 0; i < peaks_size; i++){
        cpl_vector_set(priority_to_position, i, i);
    }
    cpl_vector_copy(cpl_bivector_get_y(bivector_tmp), peak_heights);
    
    cpl_bivector_sort(bivector_tmp, bivector_tmp, CPL_SORT_ASCENDING, CPL_SORT_BY_Y);


    // Highest priority first -> iterate in reverse order (decreasing)
    for (i = peaks_size - 1; i> -1; i--){
        int j;
        // "Translate" `i` to `j` which points to current peak whose
        // neighbours are to be evaluated
        j = cpl_vector_get(priority_to_position, i);
        if (cpl_vector_get(keep, j) == 0){
            // Skip evaluation for peak already marked as "don't keep"
            continue;
        }

        k = j - 1;
        // Flag "earlier" peaks for removal until minimal distance is exceeded
        while (0 <= k &&  cpl_vector_get(peaks, j) - cpl_vector_get(peaks,k) < distance_){
            cpl_vector_set(keep, k, 0);
            k--;
        }

        k = j + 1;
        // Flag "later" peaks for removal until minimal distance is exceeded
        while (k < peaks_size && cpl_vector_get(peaks, k) - cpl_vector_get(peaks, j) < distance_){
            cpl_vector_set(keep, k, 0);
            k++;
        }
    }
    cpl_bivector_delete(bivector_tmp);

    return keep;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Detect peaks from a 1d periodic signal and store their positions
  @param    in       The 1d signal as a vector
  @param    height   The minimum height of each peak
  @param    distance The minimum distance between peaks
  @return   The vector with the maxima positions

  Loosely based on scipy.signal.find_peaks
 */
/*----------------------------------------------------------------------------*/
cpl_vector * cr2res_etalon_find_peaks(
    const cpl_vector * in, 
    double height, 
    double distance){

    cpl_vector * right_edges;
    cpl_vector * left_edges;
    cpl_vector * peaks;
    cpl_vector * peak_heights;
    cpl_vector * peaks_out;
    cpl_vector * peak_distance;
    
    double peak, peak_height ;
    int npeaks;
    int k;

    peaks = cr2res_etalon_get_local_maxima(in, &left_edges, &right_edges);
    if (peaks == NULL) return NULL;
    npeaks = cpl_vector_get_size(peaks);

    peak_heights = cpl_vector_new(npeaks);
    for (cpl_size i = 0; i < npeaks; i++)
    {
        peak = cpl_vector_get(peaks, i);
        peak_height = cpl_vector_get(in, peak);
        cpl_vector_set(peak_heights, i, peak_height);
    }

    peak_distance = cr2res_etalon_select_by_peak_distance(peaks, 
                        peak_heights, distance);

    // Evaluate height condition
    peaks_out = cpl_vector_new(npeaks);
    k = 0;
    for (cpl_size i = 0; i < npeaks; i++)
    {
        peak = cpl_vector_get(peaks, i);
        peak_height = cpl_vector_get(in, peak);

        if ((peak_height > height) && (cpl_vector_get(peak_distance, i))){
            cpl_vector_set(peaks_out, k, peak);
            k++;
        }
    }
    cpl_vector_set_size(peaks_out, k);
    cpl_vector_delete(peaks);
    cpl_vector_delete(left_edges);
    cpl_vector_delete(right_edges);
    cpl_vector_delete(peak_heights);
    cpl_vector_delete(peak_distance);

    return peaks_out;
}

static cpl_vector * cr2res_etalon_get_peaks_gaussian(
    cpl_bivector        *  spectra,
    cpl_bivector        *  spectra_err,
    cpl_polynomial      *  wavesol_init,
    cpl_array           *  wavesol_init_err,
    cpl_vector          *  peaks,
    int                    display,
    cpl_vector         **  sigma,
    cpl_vector         **  heights,
    cpl_vector         **  fit_error)
{
    cpl_bivector    *  linelist;
    cpl_matrix      *  px;
    cpl_vector      *  py;
    cpl_vector      *  sigma_loc;
    cpl_vector      *  heights_loc;
    cpl_vector      *  fit_error_loc;
    cpl_vector      *  wavelength;
    cpl_vector      *  peak_height;
    cpl_vector      *  new_peaks;
    const cpl_vector      *  in;
    cpl_size j, npeaks;

    in = cpl_bivector_get_y_const(spectra);
    npeaks = cpl_vector_get_size(peaks);
    linelist = cpl_bivector_new(npeaks);

    wavelength = cpl_bivector_get_x(linelist);
    peak_height = cpl_bivector_get_y(linelist);
    for (j = 0; j < npeaks; j++)
    {
        double wave, height;
        wave = cpl_polynomial_eval_1d(wavesol_init, 
                    cpl_vector_get(peaks, j), NULL);
        height = cpl_vector_get(in, (cpl_size) cpl_vector_get(peaks, j));
        cpl_vector_set(wavelength, j, wave);
        cpl_vector_set(peak_height, j, height);
    }

    if (cr2res_wave_extract_lines(spectra, spectra_err, wavesol_init, 
            wavesol_init_err, linelist, 40, 4, display, &px, &py, &sigma_loc, 
            &heights_loc, &fit_error_loc) == -1)
        {
        // Abort
        cpl_msg_warning(__func__, "No lines could be found");
        cpl_bivector_delete(linelist);
        return NULL;
        };

    // Make px the new peaks
    npeaks = cpl_matrix_get_nrow(px);
    new_peaks = cpl_vector_wrap(npeaks, cpl_matrix_get_data(px));
    cpl_matrix_unwrap(px);
    // Delete all the other stuff we don't need
    cpl_bivector_delete(linelist);
    cpl_vector_delete(py);

    if (sigma != NULL) *sigma = sigma_loc;
    else cpl_vector_delete(sigma_loc);
    if (heights != NULL) *heights = heights_loc;
    else cpl_vector_delete(heights_loc);
    if (fit_error_loc != NULL) *fit_error = fit_error_loc;
    else cpl_vector_delete(fit_error_loc);

    return new_peaks;
}

#define signum(x) ((x > 0) ? 1 : ((x < 0) ? -1 : 0))

double rofunc(
    double * x, 
    double * y, 
    int ndata, 
    double * a, 
    double * abdev, 
    const double b) 
{
    // Evaluates the right-hand side of equation (15.7.16) for a given value of b.
    int j;
    double sum=0.0;
    double arr[ndata];
    cpl_vector * tmp_vec;

    for (j=0;j<ndata;j++) arr[j]=y[j]-b*x[j];
    tmp_vec = cpl_vector_wrap(ndata, arr);
    *a = cpl_vector_get_median(tmp_vec);
    cpl_vector_unwrap(tmp_vec);

    *abdev=0.0;
    for (j=0;j<ndata;j++) 
    {
        double d;
        d=y[j]-(b*x[j]+*a);
        *abdev += *abdev + fabs(d);
        if (y[j] != 0.0) 
            d /= fabs(y[j]);
        if (fabs(d) > DBL_EPSILON) 
            sum += (d >= 0.0 ? x[j] : -x[j]);
    }
    return sum;
}

// This is taken from Numerical Recipes Chapter 15.7
// with some minor changes to make it work in C instead of C++
cpl_polynomial * cr2res_robust_polynomial_fit(cpl_matrix *px, cpl_vector *py)
{
    // Object for fitting a straight line y = a+bx to a set of points .xi ; yi /, by the criterion of least
    // absolute deviations. Call the constructor to calculate the fit. The answers are then available as
    // the variables a, b, and abdev (the mean absolute deviation of the points from the line).
    // Constructor. Given a set of data points xx[0..ndata-1], yy[0..ndata-1], sets a, b, and abdev.

    double a, b, abdev;
    int j;
    double b1, del, f1, sigb;
    double sx=0.0, sy=0.0, sxy=0.0, sxx=0.0, chisq=0.0;
    double *x, *y;
    int ndata;
    cpl_polynomial * poly;
    cpl_size deg;

    x = cpl_matrix_get_data(px);
    y = cpl_vector_get_data(py);
    ndata = cpl_vector_get_size(py);

    if(ndata < 1)
        return NULL;

    // As a first guess for a and b, we will find the
    // least-squares fitting line.
    for (j=0;j<ndata;j++) 
    { 
        sx += x[j];
        sy += y[j];
        sxy += x[j] * y[j];
        sxx += x[j] * x[j];
    }
    del = ndata * sxx - sx * sx;
    a = (sxx * sy - sx * sxy) / del; // Least-squares solutions.
    b = (ndata * sxy - sx * sy) / del;

    for (j = 0; j < ndata; j++){
        double temp;
        temp = y[j] - (a + b * x[j]);
        chisq += temp * temp;
    }
    sigb = sqrt(chisq / del); 
    // The standard deviation will give some idea of
    // how big an iteration step to take.
    b1 = b;
    f1 = rofunc(x, y, ndata, &a, &abdev, b1);
    if (sigb > 0.0) 
    {
        double f2, b2;
        //Guess bracket as 3-sigma away, in the downhill direction known from f1.
        b2 = b + 3.0 * sigb * signum(f1); 
        f2 = rofunc(x, y, ndata, &a, &abdev, b2);
        if (b2 == b1) 
        {
            abdev /= ndata;
            poly = cpl_polynomial_new(1);
            deg = 0;
            cpl_polynomial_set_coeff(poly, &deg, a);
            deg = 1;
            cpl_polynomial_set_coeff(poly, &deg, b);
            return poly;
        }
        //Bracketing
        while (f1 * f2 > 0.0) 
        { 
            b = b2 + 1.6 * (b2-b1);
            b1 = b2;
            f1 = f2;
            b2 = b;
            f2 = rofunc(x, y, ndata, &a, &abdev, b2);
        }
        sigb = 0.01 * sigb;
        while (fabs(b2 - b1) > sigb) 
        {
            double f;
            b = b1 + 0.5 * (b2 - b1); //Bisection.
            if (b == b1 || b == b2) 
                break;
            f=rofunc(x, y, ndata, &a, &abdev, b);
            if (f * f1 >= 0.0) 
            {
                f1 = f;
                b1 = b;
            } else {
                //f2 = f; Never used
                b2 = b;
            }
        }
    }
    abdev /= ndata;

    poly = cpl_polynomial_new(1);
    deg = 0;
    cpl_polynomial_set_coeff(poly, &deg, a);
    deg = 1;
    cpl_polynomial_set_coeff(poly, &deg, b);

    return poly;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Create the 2d wavecal fit using etalon peaks
  @param    spectra         List of extracted spectra
  @param    spectra_err     List of extracted spectra errors
  @param    wavesol_init    List of Initial wavelength solutions
  @param    wavesol_init_err List of Initial wavelength error (can be NULL)
  @param    orders          List of orders of the various spectra
  @param    ninputs         Number of entries in the previous parameters
  @param    degree_x        The polynomial degree in x
  @param    degree_y        The polynomial degree in y
  @return   Wavelength solution, i.e. polynomial that translates pixel
            values to wavelength.

 */
/*----------------------------------------------------------------------------*/
cpl_polynomial * cr2res_etalon_wave_2d(
    cpl_bivector        **  spectra,
    cpl_bivector        **  spectra_err,
    cpl_polynomial      **  wavesol_init,
    cpl_array           **  wavesol_init_err,
    int                 *   orders,
    int                 *   traces_nb,
    int                     ninputs,
    cpl_size                degree_x,
    cpl_size                degree_y,
    int                     zp_order,
    int                     display,
    cpl_array           **  wavelength_error,
    cpl_table           **  line_diagnostics)
{
    const cpl_vector * in;
    cpl_vector * peaks;
    cpl_vector * peaks_new;
    cpl_vector * freq;
    cpl_vector * wcen;
    cpl_vector * mpos;
    cpl_polynomial * poly;
    cpl_matrix * px;
    cpl_vector * pxa;
    cpl_vector * pxb;
    cpl_vector * py;
    cpl_vector ** heights;
    cpl_vector ** sigmas;
    cpl_vector ** fit_errors;
    cpl_vector ** fpe_xobs;
    cpl_vector ** fpe_wobs;
    cpl_vector ** fpe_freq;
    cpl_vector ** fpe_mord;
    cpl_vector ** fpe_cord;
    cpl_polynomial * result;
    cpl_vector * pos;
    cpl_vector * diff;

    cpl_vector * tmp_vec;
    cpl_vector * corr;
    cpl_size degree_2d[2];
    cpl_size i, j, k, deg, npeaks, npeaks_total, npoints, setting_deg;
    double wave, gap, tmp;
    double f0, fr, m;
    char * path;
    cpl_table * lines_diagnostics_loc;
    double offset;
    int pad;

    /* Check Inputs */
    if (spectra==NULL || spectra_err==NULL || wavesol_init==NULL || 
            orders==NULL)
        return NULL ;
    for (i=0 ; i<ninputs ; i++) {
        if (spectra[i]==NULL || spectra_err[i]==NULL || wavesol_init[i]==NULL)
            return NULL ;
    }
    *wavelength_error = NULL;

    // Determine the setting based on the initial wavelength guess
    deg = 0;
    tmp = cpl_polynomial_get_coeff(wavesol_init[0], &deg);
    switch ((int)(tmp/100))
    {
    case 9:
    case 10:
    case 11:
        gap = CR2RES_ETALON_GAP_Y;
        setting_deg = CR2RES_ETALON_DEGREE_Y;
        break;
    case 12:
    case 13:
    case 14:
        gap = CR2RES_ETALON_GAP_J;
        setting_deg = CR2RES_ETALON_DEGREE_J;
        break;
    case 16:
    case 17:
    case 18:
        gap = CR2RES_ETALON_GAP_H;
        setting_deg = CR2RES_ETALON_DEGREE_H;
        break;
    case 23:
    case 24:
    case 25:
        gap = CR2RES_ETALON_GAP_K;
        setting_deg = CR2RES_ETALON_DEGREE_K;
        break;
    default:
        gap = -1;
        setting_deg = 1;
        break;
    }


    fpe_xobs = cpl_malloc(ninputs * sizeof(cpl_vector *));
    fpe_wobs = cpl_malloc(ninputs * sizeof(cpl_vector *));
    fpe_freq = cpl_malloc(ninputs * sizeof(cpl_vector *));
    fpe_mord = cpl_malloc(ninputs * sizeof(cpl_vector *));
    fpe_cord = cpl_malloc(ninputs * sizeof(cpl_vector *));

    heights = cpl_malloc(ninputs * sizeof(cpl_vector *));
    sigmas = cpl_malloc(ninputs * sizeof(cpl_vector *));
    fit_errors = cpl_malloc(ninputs * sizeof(cpl_vector *));


    npeaks_total = 0;

    for (i = 0; i < ninputs; i++){
        if (spectra[i] == NULL){
            // Skip this spectrum
            fpe_xobs[i] = NULL;
            fpe_wobs[i] = NULL;
            fpe_freq[i] = NULL;
            fpe_mord[i] = NULL;
            fpe_cord[i] = NULL;
            sigmas[i] = NULL;
            heights[i] = NULL;
            fit_errors[i] = NULL;
            continue;
        }
        // Find peaks in the etalon spectra
        in = cpl_bivector_get_y_const(spectra[i]);
        if ((peaks = cr2res_etalon_find_peaks(in, cpl_vector_get_mean(in), 3))
                 == NULL){
            fpe_xobs[i] = NULL;
            fpe_wobs[i] = NULL;
            fpe_freq[i] = NULL;
            fpe_mord[i] = NULL;
            fpe_cord[i] = NULL;
            sigmas[i] = NULL;
            heights[i] = NULL;
            fit_errors[i] = NULL;
            continue;
        }

        // get the peak using the gaussian fit
        if ((peaks_new = cr2res_etalon_get_peaks_gaussian(spectra[i], spectra_err[i],
                         wavesol_init[i], wavesol_init_err[i], peaks, display,
                         &sigmas[i], &heights[i], &fit_errors[i])) == NULL){
            cpl_vector_delete(peaks);
            fpe_xobs[i] = NULL;
            fpe_wobs[i] = NULL;
            fpe_freq[i] = NULL;
            fpe_mord[i] = NULL;
            fpe_cord[i] = NULL;
            sigmas[i] = NULL;
            heights[i] = NULL;
            fit_errors[i] = NULL;
            continue;
        }
        cpl_vector_delete(peaks);
        npeaks = cpl_vector_get_size(peaks_new);
        npeaks_total += npeaks;

        /*
        Adjust frequency sequence so that freq_i = i * const
        Then determine FPE order number m for each peak in a given CRIRES+ order i.
        This is done with the following smart algebra based on diffraction equation:

             m * lambda_m     = const (1)
         (m+1) * lambda_(m+1) = const

        Subtracting:

          m * [lambda_m - lambda_(m+1)] = lambda_(m+1)
          m = lambda_(m+1) / [lambda_m - lambda_(m+1)] = lambda_(m+1) / delta lambda_m

        To appreciate how smart this is note that const is not really constant but a slow
        function of the wavelength. For the two adjacent lines it is very much the same
        helping us avoiding 9th order polynomial fit to the const as e.g. in Cersullo et al.,
        2019, A&A 624, 122.

        We can reformulate the formula to frequencies by using
          freq = c_light / lambda
          freq_i = f0 + n_i * fr 
        where we can fit f0 and fr to the detected peaks, to:
          m = abs(f0 / fr + n_i + 1)
        */
        
        // Create the frequencies
        freq = cpl_vector_new(npeaks);
        wcen = cpl_vector_new(npeaks);
        for ( j = 0; j < npeaks; j++)
        {
            wave = cpl_polynomial_eval_1d(wavesol_init[i], 
                        cpl_vector_get(peaks_new, j), NULL);
            cpl_vector_set(freq, j, SPEED_OF_LIGHT / wave);
            cpl_vector_set(wcen, j, wave);
        }

        // get first guesses for the polynomial parameters fr
        tmp_vec = cpl_vector_new(npeaks - 1);
        for ( j = 0; j < npeaks - 1; j++)
        {
            cpl_vector_set(tmp_vec, j, 
                cpl_vector_get(freq, j + 1) 
                - cpl_vector_get(freq, j));
        }
        fr = fabs(cpl_vector_get_median(tmp_vec));
        cpl_vector_delete(tmp_vec);

        // and a first guess for fd
        tmp_vec = cpl_vector_new(npeaks);
        for (j = 0; j < npeaks; j++)
        {
            cpl_vector_set(tmp_vec, j,
                modulo(cpl_vector_get(freq, j), fr)); 
        }
        f0 = cpl_vector_get_median(tmp_vec);
        cpl_vector_delete(tmp_vec);

        // Fit a 1d polynomial to the frequencies
        // determine the peak numbers
        px = cpl_matrix_new(1, npeaks);
        for (j = 0; j < npeaks; j++)
        {
            m = (cpl_vector_get(freq, j) - f0) / fr;
            cpl_matrix_set(px, 0, j, round(m));
        }

        deg = 1;
        poly = cpl_polynomial_new(1);
        cpl_polynomial_fit(poly, px, NULL, freq, NULL, CPL_FALSE, NULL, &deg);
        // and evaluate it at each peak
        for ( j = 0; j < npeaks; j++)
        {
            m = cpl_matrix_get(px, 0, j);
            // actually the frequency, just reusing the parameter name
            wave = cpl_polynomial_eval_1d(poly, m, NULL);
            cpl_vector_set(freq, j, wave);
            cpl_vector_set(wcen, j, SPEED_OF_LIGHT / wave);
        }
        deg = 0;
        //f0 = cpl_polynomial_get_coeff(poly, &deg); Never used
        deg = 1;
        //fr = cpl_polynomial_get_coeff(poly, &deg); Never used
        cpl_polynomial_delete(poly);
        // Determine M
        mpos = cpl_vector_new(npeaks);
        for ( j = 0; j < npeaks; j++)
        {
            // This is really sensitive to the wavelength solution...
            // TODO: m is awkwardly close to 0.5 before rounding...
            m = cpl_matrix_get(px, 0, j);
            cpl_vector_set(mpos, j, round(m));
            // m = cpl_vector_get(wcen, j-1)/ 
            //         fabs(cpl_vector_get(wcen, j) - cpl_vector_get(wcen, j-1));
            // cpl_vector_set(mpos, j, m);
        }
        // cpl_vector_set(mpos, 0, 1 + cpl_vector_get(mpos, 1));
        cpl_matrix_delete(px);

        // 
        tmp_vec = cpl_vector_new(npeaks);
        for (j = 1; j < npeaks; j++){
            m = cpl_vector_get(wcen, j-1)/ 
                    fabs(cpl_vector_get(wcen, j) - cpl_vector_get(wcen, j-1));
            cpl_vector_set(tmp_vec, j, m - cpl_vector_get(mpos, j));
        }
        offset = cpl_vector_get_median(tmp_vec);
        cpl_vector_add_scalar(mpos, round(offset));
        cpl_vector_delete(tmp_vec);

        // Store vectors for later
        fpe_xobs[i] = peaks_new;
        fpe_wobs[i] = wcen;
        fpe_freq[i] = freq;
        fpe_mord[i] = mpos;
        // This is only used for the debug output
        fpe_cord[i] = cpl_vector_new(npeaks);
        cpl_vector_fill(fpe_cord[i], orders[i] + zp_order);
    }

    if (npeaks_total == 0){
        cpl_msg_error(__func__, "No peaks found for Etalon wavecal");
        for (i = 0; i < ninputs; i++)
        {
            if (fpe_mord[i] == NULL) continue;
            cpl_vector_delete(fpe_xobs[i]);
            cpl_vector_delete(fpe_wobs[i]);
            cpl_vector_delete(fpe_freq[i]);
            cpl_vector_delete(fpe_mord[i]);
            cpl_vector_delete(fpe_cord[i]);
            cpl_vector_delete(heights[i]);
            cpl_vector_delete(sigmas[i]);
            cpl_vector_delete(fit_errors[i]);
        }
        cpl_free(fpe_xobs);
        cpl_free(fpe_wobs);
        cpl_free(fpe_freq);
        cpl_free(fpe_mord);
        cpl_free(fpe_cord);
        cpl_free(heights);
        cpl_free(sigmas);
        cpl_free(fit_errors);
        return NULL;
    }

    if (cpl_msg_get_level() == CPL_MSG_DEBUG){
        path = cpl_sprintf("debug_etalon_mord.fits");
        cpl_vector_save(fpe_mord[0], path, CPL_TYPE_DOUBLE, NULL, CPL_IO_CREATE);
        for (i = 1; i < ninputs; i++)
        {
            cpl_vector_save(fpe_mord[i], path, CPL_TYPE_DOUBLE, NULL, CPL_IO_EXTEND);
        }
        cpl_free(path);
        path = cpl_sprintf("debug_etalon_wobs.fits");
        cpl_vector_save(fpe_wobs[0], path, CPL_TYPE_DOUBLE, NULL, CPL_IO_CREATE);
        for (i = 1; i < ninputs; i++)
        {
            cpl_vector_save(fpe_wobs[i], path, CPL_TYPE_DOUBLE, NULL, CPL_IO_EXTEND);
        }
        cpl_free(path);
    }

    /*
    If the initial wavelength solution has an offset larger than the FPE line
    spacing the m's defined above will not match across the CRIRES+ orders.
    Here with find the offsets using eq. 1 and apply them:
    */
    // gap here is the constant m * w
    if (gap < 0)
    { 
        cpl_vector * fpe_gap;
        k = 0;
        fpe_gap = cpl_vector_new(npeaks_total);
        for (i = 0; i < ninputs; i++)
        {
            if (fpe_mord[i] == NULL) continue;
            for (j = 0; j < cpl_vector_get_size(fpe_mord[i]); j++)
            {
                cpl_vector_set(fpe_gap, k,
                    cpl_vector_get(fpe_mord[i], j)
                    * cpl_vector_get(fpe_wobs[i], j));
                k++;
            }
        }
        gap = cpl_vector_get_median(fpe_gap);
        cpl_vector_delete(fpe_gap);
        cpl_msg_debug(__func__, "Median gap is: %g", gap);
    }
    
    // Calculate and apply the correction
    corr = cpl_vector_new(ninputs);
    for (i = 0; i < ninputs; i++)
    {
        if (fpe_mord[i] == NULL) continue;
        tmp_vec = cpl_vector_new(cpl_vector_get_size(fpe_mord[i]));
        for (j = 0; j < cpl_vector_get_size(fpe_mord[i]); j++)
        {
            m = cpl_vector_get(fpe_mord[i], j);
            wave = cpl_vector_get(fpe_wobs[i], j);
            tmp = gap / wave - m;
            cpl_vector_set(tmp_vec, j, tmp);
        }
        cpl_vector_set(corr, i, round(cpl_vector_get_median(tmp_vec)));
        cpl_vector_delete(tmp_vec);
        // Apply the correction
        cpl_vector_add_scalar(fpe_mord[i], cpl_vector_get(corr, i));
    }
    cpl_vector_delete(corr);

    // In a second step we assume that m*wave varies linearly across all orders
    // and then correct for that
    npoints = 0;
    pad = min(2, ninputs-2);
    for (i = pad; i < ninputs - pad; i++)
    {
        if (fpe_mord[i] == NULL) continue;
        npoints += cpl_vector_get_size(fpe_mord[i]);
    }
    
    py = cpl_vector_new(npoints);
    px = cpl_matrix_new(1, npoints);
    k = 0;
    for (i = pad; i < ninputs - pad; i++)
    {
        if (fpe_mord[i] == NULL) continue;
        for (j = 0; j < cpl_vector_get_size(fpe_mord[i]); j++)
        {
            wave = cpl_vector_get(fpe_wobs[i], j);
            m = cpl_vector_get(fpe_mord[i], j);
            cpl_matrix_set(px, 0, k, m);
            cpl_vector_set(py, k, m * wave);
            k++;
        }
    }
    poly = cpl_polynomial_new(1);
    deg = 1;
    cpl_polynomial_fit(poly, px, NULL, py, NULL, CPL_FALSE, NULL, &deg);
    cpl_matrix_delete(px);
    cpl_vector_delete(py);

    // Calculate and apply the correction
    corr = cpl_vector_new(ninputs);
    for (i = 0; i < ninputs; i++)
    {
        if (fpe_mord[i] == NULL) continue;
        npeaks = cpl_vector_get_size(fpe_mord[i]);
        tmp_vec = cpl_vector_new(npeaks);
        for (j = 0; j < npeaks; j++)
        {
            m = cpl_vector_get(fpe_mord[i], j);
            wave = cpl_vector_get(fpe_wobs[i], j);
            gap = cpl_polynomial_eval_1d(poly, m, NULL);
            tmp = gap / wave - m;
            cpl_vector_set(tmp_vec, j, tmp);
        }
        cpl_vector_set(corr, i, round(cpl_vector_get_median(tmp_vec)));
        cpl_vector_delete(tmp_vec);
        // Apply the correction
        cpl_vector_add_scalar(fpe_mord[i], cpl_vector_get(corr, i));
    }
    cpl_vector_delete(corr);
    cpl_polynomial_delete(poly);

    /*
    Fit m * wave of all orders with a linear fit, i.e. assuming it only varies
    slowly between orders. Then determine new wavelengths based on that.
    */
    py = cpl_vector_new(npeaks_total);
    px = cpl_matrix_new(1, npeaks_total);
    k = 0;
    for (i = 0; i < ninputs; i++)
    {
        if (fpe_mord[i] == NULL) continue;
        for (j = 0; j < cpl_vector_get_size(fpe_mord[i]); j++)
        {
            wave = cpl_vector_get(fpe_wobs[i], j);
            m = cpl_vector_get(fpe_mord[i], j);
            cpl_matrix_set(px, 0, k, m);
            cpl_vector_set(py, k, m * wave);
            k++;
        }
    }

    if (setting_deg == 1){

        cpl_matrix * px_tmp;
        cpl_vector * py_tmp;

        double m0, ymax, ymin;
        // Do a robust (L1 norm) fit to the data
        // use robust to avoid outliers
        // However the method is not as numerically stable
        // so we normalize the x and y data first
        px_tmp = cpl_matrix_duplicate(px);
        py_tmp = cpl_vector_duplicate(py);

        m0 = cpl_matrix_get_min(px_tmp);
        cpl_matrix_subtract_scalar(px_tmp, m0);
        ymin = cpl_vector_get_min(py_tmp);
        cpl_vector_subtract_scalar(py_tmp, ymin);
        ymax = cpl_vector_get_max(py_tmp);
        cpl_vector_divide_scalar(py_tmp, ymax);
        poly = cr2res_robust_polynomial_fit(px_tmp, py_tmp);

        cpl_matrix_delete(px_tmp);
        cpl_vector_delete(py_tmp);

        // Fix orders to be on the same linear fit
        k = 0;
        for (i = 0; i < ninputs; i++)
        {
            double tmp1, tmp2;
            if (fpe_mord[i] == NULL) continue;
            // mean(mord * wobs)
            tmp1 = 0;
            for (j = 0; j < cpl_vector_get_size(fpe_mord[i]); j++)
            {
                m = cpl_vector_get(fpe_mord[i], j);
                wave = cpl_vector_get(fpe_wobs[i], j);
                tmp1 += m * wave;
            }
            tmp1 /= cpl_vector_get_size(fpe_mord[i]);

            // mean(poly(mord))
            tmp2 = 0;
            for (j = 0; j < cpl_vector_get_size(fpe_mord[i]); j++)
            {
                m = cpl_vector_get(fpe_mord[i], j);
                // Remember to use the normalization factors created above
                tmp2 += cpl_polynomial_eval_1d(poly, m - m0, NULL) * ymax + ymin;
            }
            tmp2 /= cpl_vector_get_size(fpe_mord[i]);

            cpl_msg_debug(__func__, 
                "Difference between measured and expected in order %lli: %f", 
                i, fabs(tmp1 - tmp2));

            for (j = 0; j < cpl_vector_get_size(fpe_mord[i]); j++)
            {
                wave = cpl_vector_get(fpe_wobs[i], j);
                m = cpl_vector_get(fpe_mord[i], j);
                wave += (tmp2 - tmp1) / m;
                cpl_vector_set(py, k, m * wave);
                cpl_vector_set(fpe_wobs[i], j, wave);
                cpl_vector_set(fpe_freq[i], j, SPEED_OF_LIGHT / wave);
                k++;
            }
        }
        cpl_polynomial_delete(poly);
    }

    // Do yet another fit to the points
    // and recalculate the wavelength on now (hopefully) consistent data
    deg = setting_deg;
    poly = cpl_polynomial_new(1);
    cpl_polynomial_fit(poly, px, NULL, py, NULL, CPL_FALSE, NULL, &deg);
    
    cpl_matrix_delete(px);
    cpl_vector_delete(py);
    
    for (i = 0; i < ninputs; i++)
    {
        if (fpe_mord[i] == NULL) continue;
        for (j = 0; j < cpl_vector_get_size(fpe_mord[i]); j++)
        {
            m = cpl_vector_get(fpe_mord[i], j);
            wave = cpl_polynomial_eval_1d(poly, m, NULL) / m;
            cpl_vector_set(fpe_freq[i], j, SPEED_OF_LIGHT / wave);
            cpl_vector_set(fpe_wobs[i], j, wave);
        }
    }
    cpl_polynomial_delete(poly);


    // Do the 2d fit
    pxa = cpl_vector_new(npeaks_total);
    pxb = cpl_vector_new(npeaks_total);
    py = cpl_vector_new(npeaks_total);
    k = 0;
    for (i = 0; i < ninputs; i++)
    {
        if (fpe_mord[i] == NULL) continue;
        npeaks = cpl_vector_get_size(fpe_mord[i]);
        for (j = 0; j < npeaks; j++){
            cpl_vector_set(pxa, k, cpl_vector_get(fpe_xobs[i], j));
            cpl_vector_set(pxb, k, orders[i] + zp_order);
            cpl_vector_set(py, k, cpl_vector_get(fpe_wobs[i], j));
            k++;
        }
    }

    degree_2d[0] = degree_x ;
    degree_2d[1] = degree_y ;
    result = cr2res_polyfit_2d(pxa, pxb, py, degree_2d);


    if (result == NULL){
        cpl_msg_error(__func__, "Error in Etalon polynomial fit");
        // cpl_polynomial_delete(result);
        cpl_vector_delete(pxa);
        cpl_vector_delete(pxb);
        cpl_vector_delete(py);
        for (i = 0; i < ninputs; i++)
        {
            if (fpe_mord[i] == NULL) continue;
            cpl_vector_delete(fpe_xobs[i]);
            cpl_vector_delete(fpe_wobs[i]);
            cpl_vector_delete(fpe_freq[i]);
            cpl_vector_delete(fpe_mord[i]);
            cpl_vector_delete(fpe_cord[i]);
            cpl_vector_delete(heights[i]);
            cpl_vector_delete(sigmas[i]);
            cpl_vector_delete(fit_errors[i]);
        }
        cpl_free(fpe_xobs);
        cpl_free(fpe_wobs);
        cpl_free(fpe_freq);
        cpl_free(fpe_mord);
        cpl_free(fpe_cord);
        cpl_free(heights);
        cpl_free(sigmas);
        cpl_free(fit_errors);
        return NULL;
    }

    if (cpl_msg_get_level() == CPL_MSG_DEBUG){
        cpl_table * tmp_table;
        path = cpl_sprintf("debug_etalon_final_mord.fits");
        cpl_vector_save(fpe_mord[0], path, CPL_TYPE_DOUBLE, NULL, CPL_IO_CREATE);
        for (i = 1; i < ninputs; i++)
        {
            cpl_vector_save(fpe_mord[i], path, CPL_TYPE_DOUBLE, NULL, CPL_IO_EXTEND);
        }
        cpl_free(path);
        path = cpl_sprintf("debug_etalon_final_wobs.fits");
        cpl_vector_save(fpe_wobs[0], path, CPL_TYPE_DOUBLE, NULL, CPL_IO_CREATE);
        for (i = 1; i < ninputs; i++)
        {
            cpl_vector_save(fpe_wobs[i], path, CPL_TYPE_DOUBLE, NULL, CPL_IO_EXTEND);
        }
        cpl_free(path);
        path = cpl_sprintf("debug_etalon_final_xobs.fits");
        cpl_vector_save(fpe_xobs[0], path, CPL_TYPE_DOUBLE, NULL, CPL_IO_CREATE);
        for (i = 1; i < ninputs; i++)
        {
            cpl_vector_save(fpe_xobs[i], path, CPL_TYPE_DOUBLE, NULL, CPL_IO_EXTEND);
        }
        cpl_free(path);
        path = cpl_sprintf("debug_etalon_final_cord.fits");
        cpl_vector_save(fpe_cord[0], path, CPL_TYPE_DOUBLE, NULL, CPL_IO_CREATE);
        for (i = 1; i < ninputs; i++)
        {
            cpl_vector_save(fpe_cord[i], path, CPL_TYPE_DOUBLE, NULL, CPL_IO_EXTEND);
        }
        cpl_free(path);

        tmp_table = cpl_table_new((degree_x + 1) * (degree_y + 1));
        cpl_table_new_column(tmp_table, "DEGREE_X", CPL_TYPE_INT);
        cpl_table_new_column(tmp_table, "DEGREE_Y", CPL_TYPE_INT);
        cpl_table_new_column(tmp_table, "COEFFICIENT", CPL_TYPE_DOUBLE);
        k = 0;
        for (i = 0; i <= degree_x; i++)
        {
            for (j = 0; j <= degree_y; j++)
            {
                degree_2d[0] = i;
                degree_2d[1] = j;
                tmp = cpl_polynomial_get_coeff(result, &degree_2d[0]);
                cpl_table_set_int(tmp_table, "DEGREE_X", k, i);
                cpl_table_set_int(tmp_table, "DEGREE_Y", k, j);
                cpl_table_set_double(tmp_table, "COEFFICIENT", k, tmp);
                k++;
            }
        }
        path = cpl_sprintf("debug_etalon_final_poly.fits");
        cpl_table_save(tmp_table, NULL, NULL, path, CPL_IO_CREATE);
        cpl_table_delete(tmp_table);
        cpl_free(path);
    }

    /* Create / Fill / Merge the lines diagnostics table  */
    if (line_diagnostics != NULL) {
        *line_diagnostics = NULL;
        cpl_msg_debug(__func__, "Number of lines: %"CPL_SIZE_FORMAT, npeaks_total);
        for (i = 0; i < ninputs; i++)
        {
            cpl_polynomial * wavesol_loc;
            if (fpe_mord[i] == NULL) continue;
            npeaks = cpl_vector_get_size(fpe_mord[i]);
            /* Create */
            lines_diagnostics_loc =
                cr2res_dfs_create_lines_diagnostics_table(npeaks) ;
            wavesol_loc = cr2res_wave_poly_2d_to_1d(result, 
                    orders[i] + zp_order);
            /* Fill */
            for (j=0 ; j < npeaks ; j++) {
                double pix_pos, lambda_cat, lambda_meas, line_width, line_intens, fit_error;
                pix_pos = cpl_vector_get(fpe_xobs[i], j);
                lambda_cat = cpl_vector_get(fpe_wobs[i], j) ;
                lambda_meas = cpl_polynomial_eval_1d(wavesol_loc, pix_pos,
                        NULL) ;
                line_width = cpl_vector_get(sigmas[i], j) ;
                line_intens = cpl_vector_get(heights[i], j) ;
                fit_error = cpl_vector_get(fit_errors[i], j) ;
                m = cpl_vector_get(fpe_mord[i], j);
                cpl_table_set_int(lines_diagnostics_loc,
                        CR2RES_COL_ORDER, j, orders[i]) ;
                cpl_table_set_int(lines_diagnostics_loc,
                        CR2RES_COL_TRACENB, j, traces_nb[i]) ;
                cpl_table_set_double(lines_diagnostics_loc,
                        CR2RES_COL_MEASURED_LAMBDA, j, lambda_meas) ;
                cpl_table_set_double(lines_diagnostics_loc,
                        CR2RES_COL_CATALOG_LAMBDA, j, lambda_cat);
                cpl_table_set_double(lines_diagnostics_loc,
                        CR2RES_COL_DELTA_LAMBDA, j, lambda_cat-lambda_meas);
                cpl_table_set_double(lines_diagnostics_loc,
                        CR2RES_COL_MEASURED_PIXEL, j, pix_pos);
                cpl_table_set_double(lines_diagnostics_loc,
                        CR2RES_COL_LINE_WIDTH, j, line_width) ;
                cpl_table_set_double(lines_diagnostics_loc,
                        CR2RES_COL_FIT_QUALITY, j, fit_error) ;
                cpl_table_set_double(lines_diagnostics_loc,
                        CR2RES_COL_INTENSITY, j, line_intens) ;
                cpl_table_set_double(lines_diagnostics_loc,
                        CR2RES_COL_FPET_M, j, m) ;
            }
            cpl_polynomial_delete(wavesol_loc);
            /* Merge */
            if (*line_diagnostics == NULL) {
                *line_diagnostics = lines_diagnostics_loc ;
                lines_diagnostics_loc = NULL ;
            } else if (lines_diagnostics_loc != NULL) {
                /* Merge with previous */
                cpl_table_insert(*line_diagnostics, lines_diagnostics_loc,
                        cpl_table_get_nrow(*line_diagnostics)) ;
                cpl_table_delete(lines_diagnostics_loc) ;
            }
        }
    }
    cpl_msg_debug(__func__,"%s",cpl_error_get_where());
    npeaks = cpl_vector_get_size(py);
    // Calculate absolute difference between polynomial and
    // catalog value for each line
    // use px and py, so that only good lines are used
    diff = cpl_vector_new(npeaks);
    pos = cpl_vector_new(2);
    for (i = 0; i < npeaks; i++){
        cpl_vector_set(pos, 0, cpl_vector_get(pxa, i));
        cpl_vector_set(pos, 1, cpl_vector_get(pxb, i));
        cpl_vector_set(diff, i, fabs(
            cpl_polynomial_eval(result, pos)
            - cpl_vector_get(py, i)));
    }

    if (*wavelength_error == NULL)
        *wavelength_error = cpl_array_new(2, CPL_TYPE_DOUBLE);
    cpl_array_set_double(*wavelength_error, 0,
            cpl_vector_get_mean(diff));
    cpl_array_set_double(*wavelength_error, 1,
            cpl_vector_get_max(diff));

    cpl_vector_delete(diff);
    cpl_vector_delete(pos);

    cpl_vector_delete(pxa);
    cpl_vector_delete(pxb);
    cpl_vector_delete(py);

    for (i = 0; i < ninputs; i++)
    {
        if (fpe_mord[i] == NULL) continue;
        cpl_vector_delete(fpe_xobs[i]);
        cpl_vector_delete(fpe_wobs[i]);
        cpl_vector_delete(fpe_freq[i]);
        cpl_vector_delete(fpe_mord[i]);
        cpl_vector_delete(fpe_cord[i]);
        cpl_vector_delete(heights[i]);
        cpl_vector_delete(sigmas[i]);
        cpl_vector_delete(fit_errors[i]);
    }
    cpl_free(fpe_xobs);
    cpl_free(fpe_wobs);
    cpl_free(fpe_freq);
    cpl_free(fpe_mord);
    cpl_free(fpe_cord);
    cpl_free(heights);
    cpl_free(sigmas);
    cpl_free(fit_errors);

    return result;
    

}
