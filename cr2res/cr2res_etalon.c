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
#include "cr2res_etalon.h"

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

#define SPEED_OF_LIGHT 299792.458

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static cpl_mask * cr2res_etalon_binary_image(const cpl_image *) ;


/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_etalon       Etalon related
 */
/*----------------------------------------------------------------------------*/

/**@{*/

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
    cpl_image       *   blob_image ;
    cpl_mask        *   blob_mask ;
    double              threshold ;
    int                 i, j, start_x, end_x ;

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
    for (i=0 ; i<cpl_vector_get_size(local_maxima) ; i++) {

        /* Get start_x end_x */
        start_x = 1;
        end_x = cpl_image_get_size_x(in) ;
        if (i>0) 
            start_x = (int)(pmax[i-1] + (pmax[i]-pmax[i-1])/2.) ;
        if (i<cpl_vector_get_size(local_maxima)-1)
            end_x = (int)(pmax[i] + (pmax[i+1]-pmax[i])/2.) ;

        cpl_msg_debug(__func__, "Extract %d -> %d", start_x, end_x) ;

        /* Extact the current blob */
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
    int i, j, k, peaks_size, distance_;

    peaks_size = cpl_vector_get_size(peaks);
    // Round up because actual peak distance can only be natural number
    distance_ = distance;
    // Prepare array of flags
    keep = cpl_vector_new(peaks_size);
    for (i = 0; i < peaks_size; i++){
        cpl_vector_set(keep, i, 1);
    }

    // Create map from `i` (index for `peaks` sorted by `priority`) to `j` (index
    // for `peaks` sorted by position). This allows to iterate `peaks` and `keep`
    // with `j` by order of `priority` while still maintaining the ability to
    // step to neighbouring peaks with (`j` + 1) or (`j` - 1).
    
    cpl_bivector * bivector_tmp = cpl_bivector_new(peaks_size);
    priority_to_position = cpl_bivector_get_x(bivector_tmp);

    for (i = 0; i < peaks_size; i++){
        cpl_vector_set(priority_to_position, i, i);
    }
    cpl_vector_copy(cpl_bivector_get_y(bivector_tmp), peak_heights);
    
    cpl_bivector_sort(bivector_tmp, bivector_tmp, CPL_SORT_ASCENDING, CPL_SORT_BY_Y);


    // Highest priority first -> iterate in reverse order (decreasing)
    for (i = peaks_size - 1; i> -1; i--){
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
    
    double peak, peak_height, peak_width, left, right;
    int npeaks;
    int k;

    peaks = cr2res_etalon_get_local_maxima(in, &left_edges, &right_edges);
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
    int                     ninputs,
    cpl_size                degree_x,
    cpl_size                degree_y)
{
    const cpl_vector * in;
    cpl_vector * peaks;
    cpl_vector * freq_peaks;
    cpl_vector * n_peaks;
    cpl_vector * temp;
    cpl_vector * pf;
    cpl_vector * pn;
    cpl_polynomial * result;
    cpl_polynomial * wavesol;
    cpl_matrix * samppos;
    cpl_matrix * pxo;
    cpl_size i, j, npeaks, npeaks_old, total_peaks;
    cpl_size degree, degree_2d[2];
    cpl_error_code error;
    double freq, wave, n, offset;
    double f0, fr, fd;


    /* Check Inputs */
    if (spectra==NULL || spectra_err==NULL || wavesol_init==NULL ||
            orders==NULL)
        return NULL ;

    // Initialize data
    result = cpl_polynomial_new(2);
    wavesol = cpl_polynomial_new(1);
    // TODO: up to 1000 peaks per order?
    total_peaks = ninputs * 1000;
    pxo = cpl_matrix_new(2, total_peaks);
    pf = cpl_vector_new(total_peaks);
    pn = cpl_vector_new(total_peaks);
    npeaks_old = 0;

    for (i = 0; i < ninputs; i++){
        // Find peaks in the etalon spectra
        in = cpl_bivector_get_y_const(spectra[i]);
        peaks = cr2res_etalon_find_peaks(in, cpl_vector_get_mean(in), 3);
        npeaks = cpl_vector_get_size(peaks);
        // get the frequency of each peak as estimated from the initial guess
        freq_peaks = cpl_vector_new(npeaks);
        n_peaks = cpl_vector_new(npeaks);
        for (j = 0; j < npeaks; j++)
        {
            wave = cpl_polynomial_eval_1d(wavesol_init[i], 
                        cpl_vector_get(peaks, j), NULL);
            freq = SPEED_OF_LIGHT / wave;
            cpl_vector_set(freq_peaks, j, freq);
        }

        // get first guesses for the polynomial parameters fr
        temp = cpl_vector_new(npeaks - 1);
        for ( j = 0; j < npeaks - 1; j++)
        {
            cpl_vector_set(temp, j, 
                cpl_vector_get(freq_peaks, j + 1) 
                - cpl_vector_get(freq_peaks, j));
        }
        fr = cpl_vector_get_median(temp);
        cpl_vector_delete(temp);

        // and a first guess for fd
        temp = cpl_vector_new(npeaks);
        for (j = 0; j < npeaks; j++)
        {
            cpl_vector_set(temp, j, 
                fmod(cpl_vector_get(freq_peaks, j), fr)); 
        }
        fd = cpl_vector_get_median(temp);
        cpl_vector_delete(temp);

        // determine the peak numbers
        for (j = 0; j < npeaks; j++)
        {
            n = (cpl_vector_get(freq_peaks, j) - fd) / fr;
            cpl_vector_set(n_peaks, j, round(n)); 
        }
        cpl_vector_subtract_scalar(n_peaks, cpl_vector_get(n_peaks, 0));

        // fit a polynomial to the peaks and peak numbers
        degree = 1;
        samppos = cpl_matrix_wrap(1, npeaks, cpl_vector_get_data(n_peaks));
        error = cpl_polynomial_fit(wavesol, samppos, NULL, freq_peaks, NULL, 
                                CPL_FALSE, NULL, &degree);
        cpl_matrix_unwrap(samppos);

        if (error != CPL_ERROR_NONE){
            cpl_msg_error(__func__, "%s", cpl_error_get_message());
            cpl_error_reset();
            cpl_polynomial_delete(result);
            cpl_matrix_delete(pxo);
            cpl_vector_delete(pf);
            cpl_vector_delete(pn);
            cpl_polynomial_delete(wavesol);
            cpl_vector_delete(freq_peaks);
            cpl_vector_delete(n_peaks);
            return NULL;
        }

        degree = 0;
        fd = cpl_polynomial_get_coeff(wavesol, &degree);
        degree = 1;
        fr = cpl_polynomial_get_coeff(wavesol, &degree);

        // The first order determines the anchor frequency
        // this is arbitrary
        if (i == 0){
            f0 = fd;
        } else {
            // determine the offset to the anchor frequency
            offset = round((f0 - fd) / fr);
            cpl_vector_subtract_scalar(n_peaks, offset);
        }

        for (j = 0; j < npeaks; j++){
            cpl_vector_set(pf, npeaks_old + j, cpl_vector_get(freq_peaks, j));
            cpl_matrix_set(pxo, 0, npeaks_old + j, cpl_vector_get(peaks, j));
            cpl_matrix_set(pxo, 1, npeaks_old + j, orders[i]);
            cpl_vector_set(pn, npeaks_old + j, 
                            cpl_vector_get(n_peaks, j));
        }
        npeaks_old += npeaks;
        cpl_vector_delete(freq_peaks);
        cpl_vector_delete(n_peaks);
        cpl_vector_delete(peaks);
    }
    total_peaks = npeaks_old;

    if (total_peaks == 0){
        cpl_msg_error(__func__, "No peaks found for etalon wavecal");
        cpl_polynomial_delete(result);
        cpl_matrix_delete(pxo);
        cpl_vector_delete(pf);
        cpl_vector_delete(pn);
        cpl_polynomial_delete(wavesol);
        return NULL;
    }

    // set the vector sizes
    cpl_matrix_set_size(pxo, 2, total_peaks);
    cpl_vector_set_size(pf, total_peaks);
    cpl_vector_set_size(pn, total_peaks);

    // Fit all peaks with a single polynomial
    degree = 1;
    samppos = cpl_matrix_wrap(1, total_peaks, cpl_vector_get_data(pn));
    error = cpl_polynomial_fit(wavesol, samppos, NULL, pf, NULL, 
                CPL_FALSE, NULL, &degree);
    cpl_matrix_unwrap(samppos);

    if (error != CPL_ERROR_NONE){
        cpl_msg_error(__func__, "Polynomial fit failed in etalon wavecal. %s", 
            cpl_error_get_message());
        cpl_error_reset();
        cpl_polynomial_delete(result);
        cpl_matrix_delete(pxo);
        cpl_vector_delete(pf);
        cpl_vector_delete(pn);
        cpl_polynomial_delete(wavesol);
        return NULL;
    }

    // reevaluate the peak wavelengths based on the fit
    for (j = 0; j < total_peaks; j++)
    {
        freq = cpl_polynomial_eval_1d(wavesol, 
                    cpl_vector_get(pn, j), NULL);
        wave = SPEED_OF_LIGHT / freq;
        cpl_vector_set(pf, j, wave);
    }
    cpl_polynomial_delete(wavesol);
    cpl_vector_delete(pn);


    // Do the 2d fit
    degree_2d[0] = degree_x ;
    degree_2d[1] = degree_y ;
    error = cpl_polynomial_fit(result, pxo, NULL, pf, NULL, TRUE, NULL,
                    degree_2d);

    cpl_matrix_delete(pxo);
    cpl_vector_delete(pf);

    if (error != CPL_ERROR_NONE){
        cpl_msg_error(__func__, "Polynomial fit for etalon wavecal failed");
        cpl_error_reset();
        cpl_polynomial_delete(result);
        return NULL;
    }

    return result;
}