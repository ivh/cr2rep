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
  @param    speclist  Array of bivectors
  @param    n         Length of speclist, needs to be 8
  @return   cpl_bivector with Stokes parameter spectrum (P/I) and error, needs
            to be de-allocated by caller.

  The input list of spectra needs to come in this order:
    1u, 1d, 2u , 2d, 3u, 3d, 4u, 4d
  i.e. first exposure upper beam, then down, then second exposure etc.

  Demodulation formula is P/I = (R^1/4 - 1) / (R^1/4 + 1) with
    R = 1u/2u * 2d/1d * 3u/4u * 4d/3d
 */
/*----------------------------------------------------------------------------*/
cpl_bivector * cr2res_pol_demod_stokes(cpl_bivector ** speclist, int n)
{
  cpl_bivector * result;
  cpl_vector * outspec;
  cpl_vector * outerr;
  cpl_vector * R;
  cpl_vector * tmp;
  cpl_size size;

  if (speclist == NULL) return NULL;

  if (n != 8){
    cpl_msg_error(__func__, "Need 8 spectra!");
    return NULL;
  }

  for (cpl_size i = 0; i < n; i++)
  {
    if (speclist[i] == NULL) return NULL;
  }

  size = cpl_bivector_get_size(speclist[0]);
  result = cpl_bivector_new(size);
  outspec = cpl_bivector_get_x(result);
  outerr = cpl_bivector_get_y(result);

  // Initialize to 0
  for (cpl_size i = 0; i < size; i++)
  {
    // cpl_vector_set(outspec, i, 0.);  
    cpl_vector_set(outerr, i, 0.);
  }

  R = cpl_vector_duplicate(cpl_bivector_get_x(speclist[0])); // 1u
  cpl_vector_divide(R, cpl_bivector_get_x(speclist[2])); // 2u

  tmp = cpl_vector_duplicate(cpl_bivector_get_x(speclist[3])); // 2d
  cpl_vector_divide(tmp, cpl_bivector_get_x(speclist[1])); // 1d
  cpl_vector_multiply(R, tmp);
  cpl_vector_delete(tmp);

  tmp = cpl_vector_duplicate(cpl_bivector_get_x(speclist[4])); // 3u
  cpl_vector_divide(tmp, cpl_bivector_get_x(speclist[6])); // 4u
  cpl_vector_multiply(R, tmp);
  cpl_vector_delete(tmp);
  
  tmp = cpl_vector_duplicate(cpl_bivector_get_x(speclist[7])); // 3d
  cpl_vector_divide(tmp, cpl_bivector_get_x(speclist[5])); // 4d
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
  for (cpl_size i = 0; i < n; i++)
  {
    tmp = cpl_vector_duplicate(cpl_bivector_get_y(speclist[i]));
    cpl_vector_divide(tmp, cpl_bivector_get_x(speclist[i]));
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
  @param    speclist  Array of bivectors
  @param    n         Length of speclist, needs to be 8
  @return   cpl_bivector with Null spectrum and error, needs
            to be de-allocated by caller.

  The input list of spectra needs to come in this order:
    1u, 1d, 2u , 2d, 3u, 3d, 4u, 4d
  i.e. first exposure upper beam, then down, then second exposure etc.

  Demodulation formula is N = ...  , with
    R = 1u/2u * 2d/1d * 3u/4u * 4d/3d
 */
/*----------------------------------------------------------------------------*/
cpl_bivector * cr2res_pol_demod_null(cpl_bivector ** speclist, int n)
{
  cpl_vector * outspec;
  cpl_vector * outerr;
  cpl_vector * R;
  cpl_vector * tmp;

  /* TODO: Implement */

  return NULL;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Combine extracted spectra into Intensity spectrum
  @param    speclist  Array of bivectors
  @param    n         Length of speclist
  @return   cpl_bivector with intensity spectrum and error, needs
            to be de-allocated by caller.

  The calculation is a simple sum of input spectra, divided by half the
  number of sepctra, since two pol-beams together make up one unit intensity.
 */
/*----------------------------------------------------------------------------*/
cpl_bivector * cr2res_pol_demod_intens(cpl_bivector ** speclist, int n)
{
  cpl_vector * outspec;
  cpl_vector * outerr;
  cpl_vector * tmp;
  int i;

  for (i=0;i<n;i++){
    if(i==0){
      outspec = cpl_vector_duplicate(cpl_bivector_get_x(speclist[i]));
      outerr = cpl_vector_duplicate(cpl_bivector_get_y(speclist[i]));
      cpl_vector_power(outerr, 2.0);
    }
    else {
      cpl_vector_add(outspec, cpl_bivector_get_x(speclist[i]));
      tmp = cpl_vector_duplicate(cpl_bivector_get_y(speclist[i]));
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
