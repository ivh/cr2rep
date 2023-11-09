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

#include "cr2res_pol.h"
#include "cr2res_dfs.h"
#include "cr2res_calib.h"
#include "cr2res_trace.h"
#include "cr2res_pfits.h"

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

int cr2res_pol_resample(cpl_vector ** intens,
                        cpl_vector ** wl,
                        cpl_vector ** errors,
                        int           n,
                        cpl_vector ** xmin,
                        cpl_vector ** xmax);

#define CR2RES_POL_MODE_1U "1u"
#define CR2RES_POL_MODE_1D "1d"
#define CR2RES_POL_MODE_2U "2u"
#define CR2RES_POL_MODE_2D "2d"
#define CR2RES_POL_MODE_3U "3u"
#define CR2RES_POL_MODE_3D "3d"
#define CR2RES_POL_MODE_4U "4u"
#define CR2RES_POL_MODE_4D "4d"
#define CR2RES_POL_MODE_ERROR "??"

#define CR2RES_POL_MODE(n) \
    (((n) == 0) ? CR2RES_POL_MODE_1U :     \
    (((n) == 1) ? CR2RES_POL_MODE_1D :     \
    (((n) == 2) ? CR2RES_POL_MODE_2U :     \
    (((n) == 3) ? CR2RES_POL_MODE_2D :     \
    (((n) == 4) ? CR2RES_POL_MODE_3U :     \
    (((n) == 5) ? CR2RES_POL_MODE_3D :     \
    (((n) == 6) ? CR2RES_POL_MODE_4U :     \
    (((n) == 7) ? CR2RES_POL_MODE_4D :     \
                  CR2RES_POL_MODE_ERROR    \
    ))))))))

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_pol      Polarimetry
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/* TODO demod fuinctions : rebin all 8 vectors to the same wl  */
/*                                     --> wl[0] as reference */
/*----------------------------------------------------------------------------*/
/**
  @brief    Demodulate extracted spectra into Stokes parameter
  @param    intens      Array of n extracted intensities
  @param    wl          Array of n extracted wavelengths
  @param    errors      Array of n extracted errors
  @param    n           Length of intens, wl and erros [needs to be 8]
  @return   cpl_bivector with Stokes parameter spectrum (P/I) and error, needs
            to be de-allocated by caller.

  The input list of the n spectra needs to come in this order:
    1u, 1d, 2u , 2d, 3u, 3d, 4u, 4d
  i.e. first exposure upper beam, then down, then second exposure etc.

  The demodulation formula is P/I = (R^1/4 - 1) / (R^1/4 + 1) with
    R = 1u/1d * 2d/2u * 3d/3u * 4u/4d
see equation 2 in Donati et al. 1997 bibcode: 1997MNRAS.291..658D
  Important : the first of the 8 input wavelength vectors (wl[0]) is the
  reference one for which the output parameters shall be computed
 */
/*----------------------------------------------------------------------------*/
cpl_bivector * cr2res_pol_demod_stokes(
        cpl_vector  **  intens,
        cpl_vector  **  wl,
        cpl_vector  **  errors,
        int             n)
{
  cpl_bivector * result;
  cpl_vector * outspec;
  cpl_vector * outerr;
  cpl_size size;
  cpl_vector * xmin, * xmax;

  cpl_vector ** intens_local;
  cpl_vector ** errors_local;
  cpl_vector ** wl_local;
  double r, s, e, tmp;
  cpl_size i, j;

  /* Check entries */
  if (intens == NULL || wl == NULL || errors == NULL) return NULL;
  if (n != 8) {
    cpl_msg_error(__func__, 
          "Got %i spectra, but expected 8 for polarimetry!", n);
    return NULL;
  }

  if (intens[0]== NULL) {
    cpl_msg_error(__func__, 
          "Spectrum %s is missing for polarimetry", CR2RES_POL_MODE(0));
    return NULL;
  }
  size = cpl_vector_get_size(intens[0]);
  for (i = 0; i < n; i++) {
    if (intens[i] == NULL) {
      cpl_msg_error(__func__, 
            "Spectrum %s is missing for polarimetry", CR2RES_POL_MODE(i));
      return NULL;
    }
    if (wl[i] == NULL) {
      cpl_msg_error(__func__, 
            "Wavelength %s is missing for polarimetry", CR2RES_POL_MODE(i));
      return NULL;
    }
    if (errors[i] == NULL) {
      cpl_msg_error(__func__, 
            "Errors %s are missing for polarimetry", CR2RES_POL_MODE(i));
      return NULL;
    }
    if (cpl_vector_get_size(intens[i]) != size) {
      cpl_msg_error(__func__, "Spectra have different sizes");
      return NULL;
    }
    if (cpl_vector_get_size(wl[i]) != size) {
      cpl_msg_error(__func__, "Wavelengths have different sizes");
      return NULL;
    }
    if (cpl_vector_get_size(errors[i]) != size) {
      cpl_msg_error(__func__, "Errors have different sizes");
      return NULL;
    }
  }

  // Resample to common wavelength grid
  // This modifies intens and errors, so we should copy them first
  intens_local = cpl_malloc(n * sizeof(cpl_vector*));
  errors_local = cpl_malloc(n * sizeof(cpl_vector*));
  wl_local = cpl_malloc(n * sizeof(cpl_vector*));
  for (i = 0; i < n; i++){
    intens_local[i] = cpl_vector_duplicate(intens[i]);
    errors_local[i] = cpl_vector_duplicate(errors[i]);
    wl_local[i] = cpl_vector_duplicate(wl[i]);
  }
  xmin = cpl_vector_new(n);
  xmax = cpl_vector_new(n);
  if (cr2res_pol_resample(intens_local, wl_local, 
                  errors_local, n, &xmin, &xmax) == -1){
    cpl_msg_error(__func__, 
        "Could not resample polarimetry spectra to the same wavelength scale");
    for (i = 0; i < n; i++)
    {
      cpl_vector_delete(intens_local[i]);
      cpl_vector_delete(errors_local[i]);
      cpl_vector_delete(wl_local[i]);
    }
    cpl_free(intens_local);
    cpl_free(errors_local);
    cpl_free(wl_local);
    cpl_vector_delete(xmin);
    cpl_vector_delete(xmax);
    return NULL;
  }

  result = cpl_bivector_new(size);
  outspec = cpl_bivector_get_x(result);
  outerr = cpl_bivector_get_y(result);

  // Initialize to 0
  for (i = 0; i < size; i++)
  {
    cpl_vector_set(outspec, i, 0.);
    cpl_vector_set(outerr, i, 0.);
  }

  for (i = cpl_vector_get_max(xmin);
        i < cpl_vector_get_min(xmax) + 1; i++){
    // Calculate R
    r = cpl_vector_get(intens_local[0], i);  // 1u
    r /= cpl_vector_get(intens_local[2], i); // 2u
    r *= cpl_vector_get(intens_local[3], i); // 2d
    r /= cpl_vector_get(intens_local[1], i); // 1d
    r *= cpl_vector_get(intens_local[5], i); // 3d
    r /= cpl_vector_get(intens_local[7], i); // 4d
    r *= cpl_vector_get(intens_local[6], i); // 4u
    r /= cpl_vector_get(intens_local[4], i); // 3u
    r = pow(r, 0.25); // 0.25 = 2/n

    // Calculate Spectrum
    s = (r - 1) / (r + 1);

    // Calculate Error
    // sum((err / spec)**2)
    e = 0;
    for (j = 0; j < n; j++){
      tmp = cpl_vector_get(errors_local[j], i) /
          cpl_vector_get(intens_local[j], i);
      e += tmp * tmp;
    }

    // 0.5 * R / (R+1)**2 * sqrt(err)
    e = sqrt(e);
    e *= 0.5 * r / ((r + 1) * (r + 1));

    cpl_vector_set(outspec, i, s);
    cpl_vector_set(outerr, i, e);
  }

  // We need to remove pointless values outside the wavelength overlap
  for (j = 0; j < cpl_vector_get_max(xmin); j++)
  {
    cpl_vector_set(outspec, j, NAN);
    cpl_vector_set(outerr, j, NAN);
  }
  for (j = cpl_vector_get_min(xmax) + 1; j < size; j++)
  {
    cpl_vector_set(outspec, j, NAN);
    cpl_vector_set(outerr, j, NAN);
  }

  for (i = 0; i < n; i++){
    cpl_vector_delete(intens_local[i]);
    cpl_vector_delete(errors_local[i]);
    cpl_vector_delete(wl_local[i]);
  }
  cpl_free(intens_local);
  cpl_free(errors_local);
  cpl_free(wl_local);
  cpl_vector_delete(xmin);
  cpl_vector_delete(xmax);

  if (cpl_error_get_code() != CPL_ERROR_NONE) {
    cpl_msg_error(__func__, "Error message: %s", cpl_error_get_message());
    cpl_bivector_delete(result);
    cpl_error_reset();
    return NULL;
  }

  return result;
}

/*----------------------------------------------------------------------------*/
/**
 * @brief Resample all spectra to the same wavelength grid
 * @param intens vectors to be resampled
 * @param wl wavelengths of those vectors
 * @param errors uncertainties of the vectors
 * @param n number of spectra
 * @return 0 on success, != 0 on failure

  Important : the first of the 8 input wavelength vectors (wl[0]) is the
  reference one for which the output parameters shall be computed
 */
/*----------------------------------------------------------------------------*/
int cr2res_pol_resample(cpl_vector ** intens,
                        cpl_vector ** wl,
                        cpl_vector ** errors,
                        int           n,
                        cpl_vector  ** xmin,
                        cpl_vector  ** xmax)
{
  // TODO: Use some other grid as baseline?
  double wmin, wmax;
  cpl_vector * master_wl;
  cpl_bivector * tmp;
  cpl_bivector * target;
  cpl_size i, j, m;

  // Construct master_wl
  // so that we only include points from within ALL wavelength ranges
  master_wl = cpl_vector_duplicate(wl[0]);
  wmin = -INFINITY;
  wmax = INFINITY;
  for (i = 0; i < n; i++)
  {
    if (cpl_vector_get(wl[i], 0) > wmin){
      wmin = cpl_vector_get(wl[i], 0);
    }
    if (cpl_vector_get(wl[i], cpl_vector_get_size(wl[i]) - 1) < wmax){
      wmax = cpl_vector_get(wl[i], cpl_vector_get_size(wl[i]) - 1);
    }
  }

  for (i = 0; i < n; i++){
    cpl_vector_set(*xmin, i, 0);
    cpl_vector_set(*xmax, i, cpl_vector_get_size(wl[i]) - 1);

    for (j = 0; j < cpl_vector_get_size(wl[i]); j++){
      if (cpl_vector_get(wl[i], j) >= wmin){
        cpl_vector_set(*xmin, i, j);
        break;
      }
    }

    for (j = cpl_vector_get_size(wl[i]) - 1; j >= 0; j--){
      if (cpl_vector_get(wl[i], j) <= wmax){
        cpl_vector_set(*xmax, i, j);
        break;
      }
    }
  }

  // We need to limit master_wl to the valid range, since bivector interpolate
  // does not extrapolate at all, but gives an error (?)
  m = cpl_vector_get_size(wl[0]);
  for (i = 0; i < m; i++)
  {
    if (cpl_vector_get(master_wl, i) < wmin){
      cpl_vector_set(master_wl, i, wmin);
    } else if (cpl_vector_get(master_wl, i) > wmax){
      cpl_vector_set(master_wl, i, wmax);
    }
  }

  // Do we need to copy the arrays before interpolation
  // or does it work in place?
  for (i = 0; i < n; i++)
  {
    tmp = cpl_bivector_wrap_vectors(wl[i], intens[i]);
    target = cpl_bivector_wrap_vectors(master_wl, intens[i]);
    cpl_bivector_interpolate_linear(target, tmp);
    cpl_bivector_unwrap_vectors(tmp);
    cpl_bivector_unwrap_vectors(target);

    tmp = cpl_bivector_wrap_vectors(wl[i], errors[i]);
    target = cpl_bivector_wrap_vectors(master_wl, errors[i]);
    cpl_bivector_interpolate_linear(target, tmp);
    cpl_bivector_unwrap_vectors(tmp);
    cpl_bivector_unwrap_vectors(target);

    // We copy the wavelength here, even though it is currently noy used in
    // the rest of the code
    cpl_vector_copy(wl[i], master_wl);

    for (j = 0; j < cpl_vector_get(*xmin, i); j++)
    {
      cpl_vector_set(intens[i], j, 1);
      cpl_vector_set(errors[i], j, 1);
    }
    for (j = cpl_vector_get(*xmax, i) + 1;
         j < cpl_vector_get_size(intens[i]); j++)
    {
      cpl_vector_set(intens[i], j, 1);
      cpl_vector_set(errors[i], j, 1);
    }
  }

  cpl_vector_delete(master_wl);

  return 0;
}

/**
  @brief    Demodulate extracted spectra into Null spectrum
  @param    intens      Array of n extracted intensities
  @param    wl          Array of n extracted wavelengths
  @param    errors      Array of n extracted errors
  @param    n           Length of intens, wl and erros [needs to be 8]
  @return   cpl_bivector with Null spectrum and error, needs
            to be de-allocated by caller.

  The input list of spectra needs to come in this order:
    1u, 1d, 2u , 2d, 3u, 3d, 4u, 4d
  i.e. first exposure upper beam, then down, then second exposure etc.

  Demodulation formula is N = (R^1/4 - 1) / (R^1/4 + 1) , with
    R = 1u/1d * 2u/2d * 3d/3u * 4d/4u
  see equation 3 in Donati et al. 1997 bibcode: 1997MNRAS.291..658D
  This is anologous to the Stokes demodulation but the last two ratios are
  inverted which makes the true polarization signal cancel out. Thus the
  Null spectrum's deviation from zero is an inverse measure of quality.
  We simply re-use the Stokes demodulation after switching the spectra
  in the input order.

  Important : the first of the 8 input wavelength vectors (wl[0]) is the
  reference one for which the output parameters shall be computed
 */
/*----------------------------------------------------------------------------*/
cpl_bivector * cr2res_pol_demod_null(
        cpl_vector  **  intens,
        cpl_vector  **  wl,
        cpl_vector  **  errors,
        int             n)
{
  cpl_vector  **  swapintens;
  cpl_vector  **  swapwl;
  cpl_vector  **  swaperrors;
  cpl_vector  *   tmpintens;
  cpl_vector  *   tmpwl;
  cpl_vector  *   tmperrors;
  cpl_bivector  * out ;
  cpl_size        size;
  int             i;

  /* Check entries */
  if (intens == NULL || wl == NULL || errors == NULL) return NULL;
  if (n != 8) {
    cpl_msg_error(__func__, 
          "Got %i spectra, but expected 8 for polarimetry!", n);
    return NULL;
  }

  if (intens[0]== NULL) {
    cpl_msg_error(__func__, 
          "Spectrum %s is missing for polarimetry", CR2RES_POL_MODE(0));
    return NULL;
  }
  size = cpl_vector_get_size(intens[0]);
  for (i = 0; i < n; i++) {
    if (intens[i] == NULL) {
      cpl_msg_error(__func__, 
            "Spectrum %s is missing for polarimetry", CR2RES_POL_MODE(i));
      return NULL;
    }
    if (wl[i] == NULL) {
      cpl_msg_error(__func__, 
            "Wavelength %s is missing for polarimetry", CR2RES_POL_MODE(i));
      return NULL;
    }
    if (errors[i] == NULL) {
      cpl_msg_error(__func__, 
            "Errors %s are missing for polarimetry", CR2RES_POL_MODE(i));
      return NULL;
    }
    if (cpl_vector_get_size(intens[i]) != size) {
      cpl_msg_error(__func__, "Spectra have different sizes");
      return NULL;
    }
    if (cpl_vector_get_size(wl[i]) != size) {
      cpl_msg_error(__func__, "Wavelengths have different sizes");
      return NULL;
    }
    if (cpl_vector_get_size(errors[i]) != size) {
      cpl_msg_error(__func__, "Errors have different sizes");
      return NULL;
    }
  }

  // copy list to leave original unchanged
  swapintens = cpl_malloc(n * sizeof(cpl_vector *));
  swapwl = cpl_malloc(n * sizeof(cpl_vector *));
  swaperrors = cpl_malloc(n * sizeof(cpl_vector *));
  for (i=0 ; i<n ; i++) {
      swapintens[i]=intens[i];
      swapwl[i]=wl[i];
      swaperrors[i]=errors[i];
  }

  // swap index 2 and 3, i.e. 2u and 2d
  tmpintens = swapintens[2];
  tmpwl = swapwl[2];
  tmperrors = swaperrors[2];
  swapintens[2] = swapintens[3];
  swapwl[2] = swapwl[3];
  swaperrors[2] = swaperrors[3];
  swapintens[3] = tmpintens;
  swapwl[3] = tmpwl;
  swaperrors[3] = tmperrors;

  // swap index 7 and 6, i.e. 4d and 4u
  tmpintens = swapintens[7];
  tmpwl = swapwl[7];
  tmperrors = swaperrors[7];
  swapintens[7] = swapintens[6];
  swapwl[7] = swapwl[6];
  swaperrors[7] = swaperrors[6];
  swapintens[6] = tmpintens;
  swapwl[6] = tmpwl;
  swaperrors[6] = tmperrors;

  out = cr2res_pol_demod_stokes(swapintens, swapwl, swaperrors, n);

  cpl_free(swapintens);
  cpl_free(swapwl);
  cpl_free(swaperrors);

  return out ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Combine extracted spectra into Intensity spectrum
  @param    intens      Array of n extracted intensities
  @param    wl          Array of n extracted wavelengths
  @param    errors      Array of n extracted errors
  @param    n           Length of intens, wl and erros [needs to be 8]
  @return   cpl_bivector with intensity spectrum and error, needs
            to be de-allocated by caller.

  The calculation is a simple sum of input spectra, divided by half the
  number of sepctra, since two pol-beams together make up one unit intensity.

  Important : the first of the 8 input wavelength vectors (wl[0]) is the
  reference one for which the output parameters shall be computed
 */
/*----------------------------------------------------------------------------*/
cpl_bivector * cr2res_pol_demod_intens(
        cpl_vector  **  intens,
        cpl_vector  **  wl,
        cpl_vector  **  errors,
        int             n)
{
    cpl_bivector    *   result;
    cpl_vector      *   outspec;
    cpl_vector      *   outerr;
    double          *   pouterr ;
    cpl_vector      *   tmp;
    int                 ncorrections ;
    cpl_size            size, i;
    cpl_vector * xmin, * xmax;

    cpl_vector ** intens_local;
    cpl_vector ** errors_local;
    cpl_vector ** wl_local;

    /* Check entries */
    if (intens == NULL || wl == NULL || errors == NULL) return NULL;
    if (n != 8) {
      cpl_msg_error(__func__, 
            "Got %i spectra, but expected 8 for polarimetry!", n);
      return NULL;
    }

    if (intens[0]== NULL) {
      cpl_msg_error(__func__, 
            "Spectrum %s is missing for polarimetry", CR2RES_POL_MODE(0));
      return NULL;
    }
    size = cpl_vector_get_size(intens[0]);
    for (i = 0; i < n; i++) {
      if (intens[i] == NULL) {
        cpl_msg_error(__func__, 
              "Spectrum %s is missing for polarimetry", CR2RES_POL_MODE(i));
        return NULL;
      }
      if (wl[i] == NULL) {
        cpl_msg_error(__func__, 
              "Wavelength %s is missing for polarimetry", CR2RES_POL_MODE(i));
        return NULL;
      }
      if (errors[i] == NULL) {
        cpl_msg_error(__func__, 
              "Errors %s are missing for polarimetry", CR2RES_POL_MODE(i));
        return NULL;
      }
      if (cpl_vector_get_size(intens[i]) != size) {
        cpl_msg_error(__func__, "Spectra have different sizes");
        return NULL;
      }
      if (cpl_vector_get_size(wl[i]) != size) {
        cpl_msg_error(__func__, "Wavelengths have different sizes");
        return NULL;
      }
      if (cpl_vector_get_size(errors[i]) != size) {
        cpl_msg_error(__func__, "Errors have different sizes");
        return NULL;
      }
    }

    // Resample to common wavelength grid
    // This modifies intens and errors, so we should copy them first
    intens_local = cpl_malloc(n * sizeof(cpl_vector*));
    errors_local = cpl_malloc(n * sizeof(cpl_vector*));
    wl_local = cpl_malloc(n * sizeof(cpl_vector*));
    for (i = 0; i < n; i++){
      intens_local[i] = cpl_vector_duplicate(intens[i]);
      errors_local[i] = cpl_vector_duplicate(errors[i]);
      wl_local[i] = cpl_vector_duplicate(wl[i]);
    }
    xmin = cpl_vector_new(n);
    xmax = cpl_vector_new(n);

    if (cr2res_pol_resample(intens_local, wl_local, 
                    errors_local, n, &xmin, &xmax) == -1){
      cpl_msg_error(__func__, 
          "Could not resample polarimetry spectra to the same wavelength scale");
      for (i = 0; i < n; i++)
      {
        cpl_vector_delete(intens_local[i]);
        cpl_vector_delete(errors_local[i]);
        cpl_vector_delete(wl_local[i]);
      }
      cpl_free(intens_local);
      cpl_free(errors_local);
      cpl_free(wl_local);
      cpl_vector_delete(xmin);
      cpl_vector_delete(xmax);
      return NULL;
    }

    result = cpl_bivector_new(size);
    outspec = cpl_bivector_get_x(result);
    outerr = cpl_bivector_get_y(result);

    // Initialize to 0
    for (i = 0; i < size; i++) {
      cpl_vector_set(outspec, i, 0.);
      cpl_vector_set(outerr, i, 0.);
    }

    /* Allocate */
    for (i=0 ; i<n ; i++) {
        if (i==0) {
            cpl_vector_copy(outspec, intens[i]);
            cpl_vector_copy(outerr, errors[i]);
            cpl_vector_power(outerr, 2.0);
        } else {
            cpl_vector_add(outspec, intens[i]);
            tmp = cpl_vector_duplicate(errors[i]);
            cpl_vector_power(tmp, 2.0);
            cpl_vector_add(outerr, tmp);
            cpl_vector_delete(tmp);
        }
    }

    /* Clean Errors */
    ncorrections = 0 ;
    pouterr = cpl_vector_get_data(outerr) ;
    for (i=0 ; i<size ; i++) {
        if (isnan(pouterr[i]) || pouterr[i] < 0.0) {
            ncorrections++ ;
            pouterr[i] = 0.0 ;
        }
    }
    if (ncorrections > 20)
        cpl_msg_warning(__func__,
                "The Errors vector contained %d negative values",
                ncorrections) ;

    cpl_vector_power(outerr, 0.5);

    for (i = 0; i < n; i++){
      cpl_vector_delete(intens_local[i]);
      cpl_vector_delete(errors_local[i]);
      cpl_vector_delete(wl_local[i]);
    }
    cpl_free(intens_local);
    cpl_free(errors_local);
    cpl_free(wl_local);
    cpl_vector_delete(xmin);
    cpl_vector_delete(xmax);

    if (cpl_error_get_code() != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Error message: %s", cpl_error_get_message());
        cpl_bivector_delete(result);
        cpl_error_reset();
        return NULL;
    }

    return result;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Create the POL_SPEC table to be saved
  @param    orders  List of orders
  @param    wl      Wavelength for the different orders
  @param    stokes  Stokes parameters for the different orders with errors
  @param    null    Null parameters for the different orders with errors
  @param    intens  Intensity for the different orders with errors
  @param    norders Number of orders
  @return   the POL_SPEC table or NULL
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_pol_POL_SPEC_create(
        int             *   orders,
        cpl_vector      **  wl,
        cpl_bivector    **  stokes,
        cpl_bivector    **  null,
        cpl_bivector    **  intens,
        int                 norders)
{
    cpl_table       *   out ;
    char            *   col_name ;
    const double    *   px ;
    const double    *   py ;
    int                 all_null, i, nrows ;

    /* Check entries */
    if (orders==NULL || wl==NULL || stokes==NULL || null==NULL || intens==NULL)
        return NULL ;

    /* Check if all bivectors are not null */
    all_null = 1 ;
    for (i=0 ; i<norders ; i++)
        if (wl[i]!=NULL && stokes[i]!=NULL && null[i]!=NULL && intens[i]!=NULL){
            nrows = cpl_vector_get_size(wl[i]) ;
            if (cpl_bivector_get_size(stokes[i]) != nrows ||
                    cpl_bivector_get_size(null[i]) != nrows ||
                    cpl_bivector_get_size(intens[i]) != nrows) {
                cpl_msg_error(__func__, "Invalid Input Sizes") ;
                return NULL ;
            }
            all_null = 0 ;
        }
    if (all_null == 1) return NULL ;

    /* Create the table */
    out = cpl_table_new(nrows);
    for (i=0 ; i<norders ; i++) {
        if (wl[i]!=NULL && stokes[i]!=NULL && null[i]!=NULL && intens[i]!=NULL){
            /* Create POL_WAVELENGTH column */
            col_name = cr2res_dfs_POL_WAVELENGTH_colname(orders[i]) ;
            cpl_table_new_column(out, col_name, CPL_TYPE_DOUBLE);
            cpl_free(col_name) ;

            /* Create POL_STOKES column */
            col_name = cr2res_dfs_POL_STOKES_colname(orders[i]) ;
            cpl_table_new_column(out, col_name, CPL_TYPE_DOUBLE);
            cpl_free(col_name) ;

            /* Create POL_STOKES_ERROR column */
            col_name = cr2res_dfs_POL_STOKES_ERROR_colname(orders[i]) ;
            cpl_table_new_column(out, col_name, CPL_TYPE_DOUBLE);
            cpl_free(col_name) ;

            /* Create POL_NULL column */
            col_name = cr2res_dfs_POL_NULL_colname(orders[i]) ;
            cpl_table_new_column(out, col_name, CPL_TYPE_DOUBLE);
            cpl_free(col_name) ;

            /* Create POL_NULL_ERROR column */
            col_name = cr2res_dfs_POL_NULL_ERROR_colname(orders[i]) ;
            cpl_table_new_column(out, col_name, CPL_TYPE_DOUBLE);
            cpl_free(col_name) ;

            /* Create POL_INTENS column */
            col_name = cr2res_dfs_POL_INTENS_colname(orders[i]) ;
            cpl_table_new_column(out, col_name, CPL_TYPE_DOUBLE);
            cpl_free(col_name) ;

            /* Create POL_INTENS_ERROR column */
            col_name = cr2res_dfs_POL_INTENS_ERROR_colname(orders[i]) ;
            cpl_table_new_column(out, col_name, CPL_TYPE_DOUBLE);
            cpl_free(col_name) ;
        }
    }
    /* Fill the table */
    for (i=0 ; i<norders ; i++) {
        if (wl[i]!=NULL && stokes[i]!=NULL && null[i]!=NULL && intens[i]!=NULL){
            /* Fill POL_WAVELENGTH column */
            px = cpl_vector_get_data_const(wl[i]) ;
            col_name = cr2res_dfs_POL_WAVELENGTH_colname(orders[i]) ;
            cpl_table_copy_data_double(out, col_name, px) ;
            cpl_free(col_name) ;

            /* Fill POL_STOKES column */
            px = cpl_bivector_get_x_data_const(stokes[i]) ;
            col_name = cr2res_dfs_POL_STOKES_colname(orders[i]) ;
            cpl_table_copy_data_double(out, col_name, px) ;
            cpl_free(col_name) ;
            /* Fill POL_STOKES_ERROR column */
            py = cpl_bivector_get_y_data_const(stokes[i]) ;
            col_name = cr2res_dfs_POL_STOKES_ERROR_colname(orders[i]) ;
            cpl_table_copy_data_double(out, col_name, py) ;
            cpl_free(col_name) ;

            /* Fill POL_NULL column */
            px = cpl_bivector_get_x_data_const(null[i]) ;
            col_name = cr2res_dfs_POL_NULL_colname(orders[i]) ;
            cpl_table_copy_data_double(out, col_name, px) ;
            cpl_free(col_name) ;
            /* Fill POL_NULL_ERROR column */
            py = cpl_bivector_get_y_data_const(null[i]) ;
            col_name = cr2res_dfs_POL_NULL_ERROR_colname(orders[i]) ;
            cpl_table_copy_data_double(out, col_name, py) ;
            cpl_free(col_name) ;

            /* Fill POL_INTENS column */
            px = cpl_bivector_get_x_data_const(intens[i]) ;
            col_name = cr2res_dfs_POL_INTENS_colname(orders[i]) ;
            cpl_table_copy_data_double(out, col_name, px) ;
            cpl_free(col_name) ;
            /* Fill POL_INTENS_ERROR column */
            py = cpl_bivector_get_y_data_const(intens[i]) ;
            col_name = cr2res_dfs_POL_INTENS_ERROR_colname(orders[i]) ;
            cpl_table_copy_data_double(out, col_name, py) ;
            cpl_free(col_name) ;
        }
    }
    return out ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Compute the positions of the passed frames
  @param   frame1   Frame #1
  @param   frame2   Frame #2
  @param   frame3   Frame #3
  @param   frame4   Frame #4
  @return   an array of 4 integer indices or NULL in error case

  If the Out array is giving the new position : ordering[0] is the
  position of frame1, ordering[1] is the position of frame2, etc...
  If the Frames order needs to be  Frame #4, #1, #3, #2, ordering will
  contain [3, 0, 2, 1]

  The returned positions correspond to the 1,2,3,4 inputs in the demod
  fuctions.
 */
/*----------------------------------------------------------------------------*/
int * cr2res_pol_sort_frames(
        const cpl_frame *   frame1,
        const cpl_frame *   frame2,
        const cpl_frame *   frame3,
        const cpl_frame *   frame4)
{
    int                 *   idx_order ;
    const char	        *	poltype_first ;
    const char	        *	poltype_curr ;
    const char	        *	fname ;
    cpl_propertylist    *   pri_head ;
    cpl_propertylist    *   tmp_head ;
    const cpl_frame     *   frames_holder[4] ;
    int                     nframes, i ;

    /* Check Inputs */
    if (frame1 == NULL || frame2 == NULL || frame3 == NULL || frame4 == NULL)
        return NULL;

    /* Initialise */
    nframes = 4 ;
    frames_holder[0] = frame1 ;
    frames_holder[1] = frame2 ;
    frames_holder[2] = frame3 ;
    frames_holder[3] = frame4 ;

    /* Read POL.TYPE */
    fname = cpl_frame_get_filename(frames_holder[0]); 
    pri_head = cpl_propertylist_load(fname, 0);
    poltype_first = cr2res_pfits_get_poltype(pri_head);
    if (cpl_error_get_code()) {
        cpl_propertylist_delete(pri_head);
        cpl_msg_error(__func__, "Missing POL.TYPE in header") ;
        return NULL ;
    }

    /* Make sure the others are the same */
    for (i=1 ; i<nframes ; i++) {
        fname = cpl_frame_get_filename(frames_holder[i]); 
        tmp_head = cpl_propertylist_load(fname, 0);
        poltype_curr = cr2res_pfits_get_poltype(tmp_head);
        if (cpl_error_get_code()) {
            cpl_propertylist_delete(tmp_head);
            cpl_propertylist_delete(pri_head);
            cpl_msg_error(__func__, "Missing POL.TYPE in header") ;
            return NULL ;
        }
        if (strcmp(poltype_first, poltype_curr)) {
            cpl_msg_error(__func__, "Different POL.TYPE! %s != %s", 
                    poltype_first, poltype_curr) ; 
            return NULL ;  
        }
        cpl_propertylist_delete(tmp_head);
    }

    idx_order = cpl_malloc(nframes * sizeof(int)) ;

    if (!strcmp(poltype_first, "V")){
        cpl_msg_info(__func__, "Set up file order for circular pol.");
        idx_order[0] = 0 ;
        idx_order[1] = 1 ;
        idx_order[2] = 2 ;
        idx_order[3] = 3 ;
    } else {
        cpl_msg_info(__func__, "Set up file order for linear pol.");
        idx_order[0] = 0 ;
        idx_order[1] = 3 ;
        idx_order[2] = 2 ;
        idx_order[3] = 1 ;
    }
    cpl_propertylist_delete(pri_head);
    return idx_order ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Merge several POL_SPEC tables together, by averaging
  @param    pol_spec_list   The list of tables
  @param    pol_spec_nb     The number of tables in list
  @return   A merged table of the same format.
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_pol_spec_pol_merge(
        const cpl_table **  pol_spec_list,
        int                 pol_spec_nb)
{
    int i, j, k;
    int n_cols, n_rows;
    cpl_table * merged_table;
    double sum, avg;
    cpl_array * column_names_array;

    /* Check Inputs */
    if (pol_spec_list == NULL || pol_spec_nb <= 0) return NULL;
    for (i = 0; i < pol_spec_nb; i++)
        if (pol_spec_list[i] == NULL) return NULL;

    n_cols = cpl_table_get_ncol(pol_spec_list[0]);
    n_rows = cpl_table_get_nrow(pol_spec_list[0]);
    merged_table = cpl_table_new(n_rows);

    column_names_array = cpl_table_get_column_names(pol_spec_list[0]);

    for (i = 0; i < n_cols; i++) {
        const char *col_name = cpl_array_get_string(column_names_array, i);
        cpl_table_new_column(merged_table, col_name, CPL_TYPE_DOUBLE);

        /* Loop through the rows */
        for (j = 0; j < n_rows; j++) {
            sum = 0.0;
            for (k = 0; k < pol_spec_nb; k++) {
                sum += cpl_table_get_double(pol_spec_list[k], col_name, j, NULL);
            }
            avg = sum / pol_spec_nb;
            cpl_table_set_double(merged_table, col_name, j, avg);
        }
    }

    /* Clean up */
    cpl_array_delete(column_names_array);

    return merged_table;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Compute the traces for the polarimetric beams
  @param    tw_in        The input traces
  @param    decker_name   The decker name
  @param    up_down       Upper or lower beam
  @return   The newly computed trace 

Uses hardcoded correction polynomials to account for the diverging beams

    */
/*----------------------------------------------------------------------------*/

cpl_table *cr2res_pol_get_beam_trace(
    const cpl_table *tw_in,
    cr2res_decker   decker_position,
    int up_or_down)
{

    cpl_table       *   tw_out;
    cpl_polynomial  *   wl_poly;
    cpl_polynomial  *   trace_poly;
    cpl_polynomial  *   upper_poly;
    cpl_polynomial  *   lower_poly;
    cpl_polynomial  *   corr_poly;
    double              wl, y_mid, y_up, y_lo, y_corr, slope_corr, halfSF, newSF;
    cpl_array       *   tmp_arr;
    cpl_size            pow0=0;
    cpl_size            pow1=1;
    char            *   band;
    int                 nb_traces, i, o, trace_id;


    wl_poly = cr2res_get_trace_wave_poly(tw_in, CR2RES_COL_WAVELENGTH, 5, 1);
    wl = cpl_polynomial_eval_1d(wl_poly, 1024, NULL);
    cpl_polynomial_delete(wl_poly);

    /* Correction polynomial gets filled with hardcoded values. */    
    /* These were derived as linear difference from the order mid-line */
    /* for the upper/lower beams and A/B nodding position, in each band.*/
    corr_poly = cpl_polynomial_new(1);
    if (wl < 1125) {
        band=cpl_sprintf("Y");
                            /*
                            'YuA': array([  0.03988915, -22.78903764]),
                            'YuB': array([ 0.04808708, 13.12469913]),
                            'YdA': array([-0.04965059, -6.51155269]),
                            'YdB': array([-0.04019722, 28.22704312]),
                            */
        if (decker_position == CR2RES_DECKER_2_4){ // nodd A
            if (up_or_down == 1){ // UP
                cpl_polynomial_set_coeff(corr_poly, &pow0, -22.78903764);
                cpl_polynomial_set_coeff(corr_poly, &pow1, 0.03988915);
            } else { // DOWN
                cpl_polynomial_set_coeff(corr_poly, &pow0, -6.51155269);
                cpl_polynomial_set_coeff(corr_poly, &pow1, -0.04965059);
            }
        } else if (decker_position == CR2RES_DECKER_1_3){ // nodd B
            if (up_or_down == 1){ // UP
                cpl_polynomial_set_coeff(corr_poly, &pow0, 13.12469913);
                cpl_polynomial_set_coeff(corr_poly, &pow1, 0.04808708);
            } else { // DOWN
                cpl_polynomial_set_coeff(corr_poly, &pow0, 28.22704312);
                cpl_polynomial_set_coeff(corr_poly, &pow1, -0.04019722);
            }
        }
    }
    else if (wl < 1360) {
        band=cpl_sprintf("J");
                            /*
                            'JuA': array([  0.03903698, -26.55176271]),
                            'JuB': array([ 0.04496791, 10.27999294]),
                            'JdA': array([ -0.047482  , -10.87508078]),
                            'JdB': array([-0.04101687, 25.44306868]),
                            */
        if (decker_position == CR2RES_DECKER_2_4){ // nodd A
            if (up_or_down == 1){ // UP
                cpl_polynomial_set_coeff(corr_poly, &pow0, -26.55176271);
                cpl_polynomial_set_coeff(corr_poly, &pow1, 0.03903698);
            } else { // DOWN
                cpl_polynomial_set_coeff(corr_poly, &pow0, -10.87508078);
                cpl_polynomial_set_coeff(corr_poly, &pow1, -0.047482);
            }
        } else if (decker_position == CR2RES_DECKER_1_3){ // nodd B
            if (up_or_down == 1){ // UP
                cpl_polynomial_set_coeff(corr_poly, &pow0, 10.27999294);
                cpl_polynomial_set_coeff(corr_poly, &pow1, 0.04496791);
            } else { // DOWN
                cpl_polynomial_set_coeff(corr_poly, &pow0, 25.44306868);
                cpl_polynomial_set_coeff(corr_poly, &pow1, -0.04101687);
            }
        }
    
    }
    else if (wl < 1850) {
        band=cpl_sprintf("H");
                            /*
                            'HuA': array([ 1.59902724e-02, -2.17975984e+01]),
                            'HuB': array([ 0.0199041 , 17.44054878]),
                            'HdA': array([ -0.02831508, -23.23946728]),
                            'HdB': array([-0.02552662, 18.31347609]),
                            */
        if (decker_position == CR2RES_DECKER_2_4){ // nodd A
            if (up_or_down == 1){ // UP
                cpl_polynomial_set_coeff(corr_poly, &pow0, -2.17975984e+01);
                cpl_polynomial_set_coeff(corr_poly, &pow1, 1.59902724e-02);
            } else { // DOWN
                cpl_polynomial_set_coeff(corr_poly, &pow0, -23.23946728);
                cpl_polynomial_set_coeff(corr_poly, &pow1, -0.02831508);
            }
        } else if (decker_position == CR2RES_DECKER_1_3){ // nodd B
            if (up_or_down == 1){ // UP
                cpl_polynomial_set_coeff(corr_poly, &pow0, 17.44054878);
                cpl_polynomial_set_coeff(corr_poly, &pow1, 0.0199041);
            } else { // DOWN
                cpl_polynomial_set_coeff(corr_poly, &pow0, 18.31347609);
                cpl_polynomial_set_coeff(corr_poly, &pow1, -0.02552662);
            }
        }
    }
    else if (wl < 2600) {
        band=cpl_sprintf("K");
                            /*
                            'KuA': array([ 2.24337404e-02, -2.30870052e+01]),
                            'KuB': array([ 0.02467974, 16.43404529]),
                            'KdA': array([ -0.02667301, -13.55392217]),
                            'KdB': array([-2.42857387e-02,  2.58313952e+01])
                            */
        if (decker_position == CR2RES_DECKER_2_4){ // nodd A
            if (up_or_down == 1){ // UP
                cpl_polynomial_set_coeff(corr_poly, &pow0, -2.30870052e+01);
                cpl_polynomial_set_coeff(corr_poly, &pow1, 2.24337404e-02);
            } else { // DOWN
                cpl_polynomial_set_coeff(corr_poly, &pow0, -13.55392217);
                cpl_polynomial_set_coeff(corr_poly, &pow1, -0.02667301);
            }
        } else if (decker_position == CR2RES_DECKER_1_3){ // nodd B
            if (up_or_down == 1){ // UP
                cpl_polynomial_set_coeff(corr_poly, &pow0, 16.43404529);
                cpl_polynomial_set_coeff(corr_poly, &pow1, 0.02467974);
            } else { // DOWN
                cpl_polynomial_set_coeff(corr_poly, &pow0, 2.58313952e+01);
                cpl_polynomial_set_coeff(corr_poly, &pow1, -2.42857387e-02);
            }
        }
    }
    else {
        cpl_msg_error(__func__, "No proper WL found to identify band.");
        cpl_polynomial_delete(corr_poly);
        return NULL;
    }
    cpl_msg_info(__func__, "Found wl=%g nm, therefore applying beam "
                    "correction for %s-band.", wl, band);


    /* Loop through all traces in input TW */
    nb_traces = cpl_table_get_nrow(tw_in) ;
    tw_out = cpl_table_duplicate(tw_in);
    for (i=0 ; i<nb_traces ; i++) {
        o = cpl_table_get(tw_in, CR2RES_COL_ORDER, i, NULL) ;
        trace_id = cpl_table_get(tw_in, CR2RES_COL_TRACENB, i, NULL);
        if (trace_id =! 1){
            cpl_msg_error(__func__, "More than one input-trace per order"
                            "is not supported.");
            cpl_table_delete(tw_out);
            cpl_free(band);
            cpl_polynomial_delete(corr_poly);
            return NULL;
        }

        /* Evaluate at middle of detector to get average positions */
        wl_poly = cr2res_get_trace_wave_poly(tw_in, CR2RES_COL_WAVELENGTH, o,1);
        trace_poly = cr2res_get_trace_wave_poly(tw_in, CR2RES_COL_ALL, o, 1);
        upper_poly = cr2res_get_trace_wave_poly(tw_in, CR2RES_COL_UPPER, o, 1);
        lower_poly = cr2res_get_trace_wave_poly(tw_in, CR2RES_COL_LOWER, o, 1);
        y_mid = cpl_polynomial_eval_1d(trace_poly, 1024, NULL);
        y_up = cpl_polynomial_eval_1d(upper_poly, 1024, NULL);
        y_lo = cpl_polynomial_eval_1d(lower_poly, 1024, NULL);
        y_corr = cpl_polynomial_eval_1d(corr_poly, wl, NULL);

        /* Calculate new slit-fraction, important for WL-correction */
        tmp_arr = cpl_array_new(3,CPL_TYPE_DOUBLE);
        newSF = 0.5 + (y_corr) / (y_up - y_lo); // 
        halfSF = CR2RES_POLARIMETRY_DEFAULT_HEIGHT / (y_up - y_lo) / 2;
        cpl_array_set(tmp_arr, 0, newSF-halfSF);
        cpl_array_set(tmp_arr, 1, newSF);
        cpl_array_set(tmp_arr, 2, newSF+halfSF);
        cpl_table_set_array(tw_out, CR2RES_COL_SLIT_FRACTION, i, tmp_arr);
        cpl_array_delete(tmp_arr);

        /* Recalculate wavelengths at the new slit-fraction */
        tw_out = cr2res_trace_shift_wavelength(tw_out, 0.5, o, 1);

        /* Evaluate WL at detector edges, because corr-poly is P(WL)*/
        wl = cpl_polynomial_eval_1d(wl_poly, 1, NULL);
        y_corr = cpl_polynomial_eval_1d(corr_poly, wl, NULL);
        wl = cpl_polynomial_eval_1d(wl_poly, CR2RES_DETECTOR_SIZE, NULL);
        slope_corr = (cpl_polynomial_eval_1d(corr_poly, wl, NULL)
                    - y_corr) /  CR2RES_DETECTOR_SIZE;  // SLOPE!

        cpl_msg_debug(__func__, 
                "y_mid, y_up, y_lo, newSF, halfSF, slope_corr, y_corr: "
                "%g, %g, %g, %g, %g, %g, %g ",
                 y_mid, y_up, y_lo, newSF, halfSF, slope_corr, y_corr);

        /* Add the correction to the existing trace polys */
        tmp_arr = cpl_array_duplicate (
                    cpl_table_get_array(tw_out, CR2RES_COL_ALL,i));
        cpl_array_set(tmp_arr, 0, cpl_array_get(tmp_arr, 0, NULL) + y_corr);
        cpl_array_set(tmp_arr, 1, cpl_array_get(tmp_arr, 1, NULL) + slope_corr);
        cpl_table_set_array(tw_out, CR2RES_COL_ALL, i, tmp_arr);
        
        /* Edge polynomials are just shifted up/down by height/2 */
        cpl_array_set(tmp_arr, 0, cpl_array_get(tmp_arr, 0, NULL) + 
                            CR2RES_POLARIMETRY_DEFAULT_HEIGHT / 2);
        cpl_table_set_array(tw_out, CR2RES_COL_UPPER, i, tmp_arr);
        cpl_array_set(tmp_arr, 0, cpl_array_get(tmp_arr, 0, NULL) - 
                            CR2RES_POLARIMETRY_DEFAULT_HEIGHT);
        cpl_table_set_array(tw_out, CR2RES_COL_LOWER, i, tmp_arr);
        cpl_array_delete(tmp_arr);

        cpl_polynomial_delete(trace_poly);
        cpl_polynomial_delete(upper_poly);
        cpl_polynomial_delete(lower_poly);
        cpl_polynomial_delete(wl_poly);
    }


    cpl_free(band);
    cpl_polynomial_delete(corr_poly);
    return tw_out;
}

//

/**@}*/
