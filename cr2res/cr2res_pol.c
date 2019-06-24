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

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_pol      Polarimetry
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Demodulate extracted spectra into Stokes parameter
  @param    intens      Array of n extracted intenѕities
  @param    wl          Array of n extracted wavelengths
  @param    errors      Array of n extracted errors
  @param    n           Length of intens, wl and erros [needs to be 8]
  @return   cpl_bivector with Stokes parameter spectrum (P/I) and error, needs
            to be de-allocated by caller.

  The input list of the n spectra needs to come in this order:
    1u, 1d, 2u , 2d, 3u, 3d, 4u, 4d
  i.e. first exposure upper beam, then down, then second exposure etc.

  Demodulation formula is P/I = (R^1/4 - 1) / (R^1/4 + 1) with
    R = 1u/2u * 2d/1d * 3u/4u * 4d/3d
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
  cpl_vector * R;
  cpl_vector * tmp;
  cpl_size size;

  /* Check entries */
  if (intens == NULL || wl == NULL || errors == NULL) return NULL;
  if (n != 8) {
    cpl_msg_error(__func__, "Need 8 spectra!");
    return NULL;
  }
  for (cpl_size i = 0; i < n; i++) {
    if (intens[i] == NULL) return NULL;
    if (wl[i] == NULL) return NULL;
    if (errors[i] == NULL) return NULL;
  }
  /* TODO : Check that all vectors have the same size */

  size = cpl_vector_get_size(intens[0]);
  result = cpl_bivector_new(size);
  outspec = cpl_bivector_get_x(result);
  outerr = cpl_bivector_get_y(result);

  // Initialize to 0
  for (cpl_size i = 0; i < size; i++)
  {
    // cpl_vector_set(outspec, i, 0.);  
    cpl_vector_set(outerr, i, 0.);
  }

  R = cpl_vector_duplicate(intens[0]); // 1u
  cpl_vector_divide(R, intens[2]); // 2u

  tmp = cpl_vector_duplicate(intens[3]); // 2d
  cpl_vector_divide(tmp, intens[1]); // 1d
  cpl_vector_multiply(R, tmp);
  cpl_vector_delete(tmp);

  tmp = cpl_vector_duplicate(intens[4]); // 3u
  cpl_vector_divide(tmp, intens[6]); // 4u
  cpl_vector_multiply(R, tmp);
  cpl_vector_delete(tmp);
  
  tmp = cpl_vector_duplicate(intens[7]); // 3d
  cpl_vector_divide(tmp, intens[5]); // 4d
  cpl_vector_multiply(R, tmp);
  cpl_vector_delete(tmp);

  cpl_vector_power(R, 0.25); // 0.25 = 2/n
  cpl_vector_copy(outspec, R);
  cpl_vector_subtract_scalar(outspec, 1.0);
  tmp = cpl_vector_duplicate(R);
  cpl_vector_add_scalar(tmp, 1.0);
  cpl_vector_divide(outspec, tmp);
  cpl_vector_delete(tmp);

  // Calculate Error
  // sum((sigma / spec)**2)
  for (cpl_size i = 0; i < n; i++) {
    tmp = cpl_vector_duplicate(errors[i]);
    cpl_vector_divide(tmp, intens[i]);
    cpl_vector_multiply(tmp, tmp);
    cpl_vector_add(outerr, tmp);
    cpl_vector_delete(tmp);
  }
  
  // 0.5 * R / (R+1)**2 * sqrt(err)
  cpl_vector_sqrt(outerr);
  tmp = cpl_vector_duplicate(R);
  cpl_vector_add_scalar(tmp, 1);
  cpl_vector_multiply(tmp, tmp);
  cpl_vector_divide(R, tmp);
  cpl_vector_multiply(outerr, R);
  cpl_vector_multiply_scalar(outerr, 0.5);

  cpl_vector_delete(R);
  cpl_vector_delete(tmp);

  if (cpl_error_get_code() != CPL_ERROR_NONE) {
    cpl_msg_error(__func__, "Error code: %i", cpl_error_get_code());
    return NULL;
  }

  return result;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Demodulate extracted spectra into Null spectrum
  @param    intens      Array of n extracted intenѕities
  @param    wl          Array of n extracted wavelengths
  @param    errors      Array of n extracted errors
  @param    n           Length of intens, wl and erros [needs to be 8]
  @return   cpl_bivector with Null spectrum and error, needs
            to be de-allocated by caller.

  The input list of spectra needs to come in this order:
    1u, 1d, 2u , 2d, 3u, 3d, 4u, 4d
  i.e. first exposure upper beam, then down, then second exposure etc.

  Demodulation formula is N = (R^1/4 - 1) / (R^1/4 + 1) , with
    R = 1u/2u * 2d/1d * 4u/3u * 3d/4d
  This is anologous to the Stokes demodulation but the last two ratios are
  inverted which makes the true polarization signal cancel out. Thus the
  Null spectrum's deviation from zero is an inverse measure of quality.
  We simply re-use the Stokes demodulation after switching the spectra
  in the input order.
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
    int             i;

  // copy list to leave original unchanged
  swapintens = cpl_malloc(n * sizeof(cpl_vector *));
  swapwl = cpl_malloc(n * sizeof(cpl_vector *));
  swaperrors = cpl_malloc(n * sizeof(cpl_vector *));
  for (i=0 ; i<n ; i++) {
      swapintens[i]=intens[i];
      swapwl[i]=wl[i];
      swaperrors[i]=errors[i];
  }

  // swap index 6 and 4, i.e. 4u and 3u
  tmpintens = swapintens[6];
  tmpwl = swapwl[6];
  tmperrors = swaperrors[6];
  swapintens[6] = swapintens[4];
  swapwl[6] = swapwl[4];
  swaperrors[6] = swaperrors[4];
  swapintens[4] = tmpintens;
  swapwl[4] = tmpwl;
  swaperrors[4] = tmperrors;

  // swap index 7 and 5, i.e. 4d and 3d
  tmpintens = swapintens[7];
  tmpwl = swapwl[7];
  tmperrors = swaperrors[7];
  swapintens[7] = swapintens[5];
  swapwl[7] = swapwl[5];
  swaperrors[7] = swaperrors[5];
  swapintens[5] = tmpintens;
  swapwl[5] = tmpwl;
  swaperrors[5] = tmperrors;

  out = cr2res_pol_demod_stokes(swapintens, swapwl, swaperrors, n);

  cpl_free(swapintens);
  cpl_free(swapwl);
  cpl_free(swaperrors);

  return out ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Combine extracted spectra into Intensity spectrum
  @param    intens      Array of n extracted intenѕities
  @param    wl          Array of n extracted wavelengths
  @param    errors      Array of n extracted errors
  @param    n           Length of intens, wl and erros [needs to be 8]
  @return   cpl_bivector with intensity spectrum and error, needs
            to be de-allocated by caller.

  The calculation is a simple sum of input spectra, divided by half the
  number of sepctra, since two pol-beams together make up one unit intensity.
 */
/*----------------------------------------------------------------------------*/
cpl_bivector * cr2res_pol_demod_intens(
        cpl_vector  **  intens, 
        cpl_vector  **  wl, 
        cpl_vector  **  errors, 
        int             n)
{
  cpl_vector * outspec;
  cpl_vector * outerr;
  cpl_vector * tmp;
  int i;

  for (i=0;i<n;i++) {
    if (i==0) {
      outspec = cpl_vector_duplicate(intens[i]);
      outerr = cpl_vector_duplicate(errors[i]);
      cpl_vector_power(outerr, 2.0);
    } else {
      cpl_vector_add(outspec, intens[i]);
      tmp = cpl_vector_duplicate(errors[i]);
      cpl_vector_power(tmp, 2.0);
      cpl_vector_add(outerr, tmp);
    }
  }
  cpl_vector_delete(tmp);
  
  cpl_vector_divide_scalar(outspec, (double)n/2.0);
  cpl_vector_power(outerr, 0.5);

  if (cpl_error_get_code() != CPL_ERROR_NONE) {
    cpl_msg_error(__func__, "Error code: %i", cpl_error_get_code());
    return NULL;
  }

  return cpl_bivector_wrap_vectors(outspec, outerr);
}
