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

/* TODO demod fuinctions : rebin all 8 vectors to the same wl  */
/*                                     --> wl[0] as reference */
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
  cpl_vector * R;
  cpl_vector * tmp;
  cpl_size size;

  /* Check entries */
  if (intens == NULL || wl == NULL || errors == NULL) return NULL;
  if (n != 8) {
    cpl_msg_error(__func__, "Need 8 spectra!");
    return NULL;
  }

  if (intens[0]== NULL) return NULL;
  size = cpl_vector_get_size(intens[0]);
  for (cpl_size i = 0; i < n; i++) {
    if (intens[i] == NULL) return NULL;
    if (wl[i] == NULL) return NULL;
    if (errors[i] == NULL) return NULL;
    if (cpl_vector_get_size(intens[i]) != size) return NULL;
    if (cpl_vector_get_size(wl[i]) != size) return NULL;
    if (cpl_vector_get_size(errors[i]) != size) return NULL;
  }

  result = cpl_bivector_new(size);
  outspec = cpl_bivector_get_x(result);
  outerr = cpl_bivector_get_y(result);

  // Initialize to 0
  for (cpl_size i = 0; i < size; i++)
  {
    cpl_vector_set(outspec, i, 0.);  
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
    cpl_bivector_delete(result);
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
    cpl_msg_error(__func__, "Need 8 spectra!");
    return NULL;
  }
  if (intens[0]== NULL) return NULL;
  size = cpl_vector_get_size(intens[0]);
  for (cpl_size i = 0; i < n; i++) {
    if (intens[i] == NULL) return NULL;
    if (wl[i] == NULL) return NULL;
    if (errors[i] == NULL) return NULL;
    if (cpl_vector_get_size(intens[i]) != size) return NULL;
    if (cpl_vector_get_size(wl[i]) != size) return NULL;
    if (cpl_vector_get_size(errors[i]) != size) return NULL;
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
  cpl_bivector * result;
  cpl_vector * outspec;
  cpl_vector * outerr;
  cpl_vector * tmp;
  cpl_size size;
  int i;

  /* Check entries */
  if (intens == NULL || wl == NULL || errors == NULL) return NULL;
  if (n != 8) {
    cpl_msg_error(__func__, "Need 8 spectra!");
    return NULL;
  }
  if (intens[0]== NULL) return NULL;
  size = cpl_vector_get_size(intens[0]);
  for (cpl_size i = 0; i < n; i++) {
    if (intens[i] == NULL) return NULL;
    if (wl[i] == NULL) return NULL;
    if (errors[i] == NULL) return NULL;
    if (cpl_vector_get_size(intens[i]) != size) return NULL;
    if (cpl_vector_get_size(wl[i]) != size) return NULL;
    if (cpl_vector_get_size(errors[i]) != size) return NULL;
  }

  result = cpl_bivector_new(size);
  outspec = cpl_bivector_get_x(result);
  outerr = cpl_bivector_get_y(result);

  for (i=0;i<n;i++) {
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
  
  cpl_vector_divide_scalar(outspec, (double)n/2.0);
  cpl_vector_power(outerr, 0.5);

  if (cpl_error_get_code() != CPL_ERROR_NONE) {
    cpl_msg_error(__func__, "Error code: %i", cpl_error_get_code());
    cpl_bivector_delete(result);
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

    /* Check if all bivectorѕ are not null */
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

  If the Out array is giving the new position : ordering[0] iѕ the
  position of frame1, ordering[1] is the position of frame2, etc...
  If the Frames order needs to be  Frame #4, #1, #3, #2, orderіng will
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
    int     *   idx_order ;

    /* Check Inputs */
    if (frame1 == NULL || frame2 == NULL || frame3 == NULL || frame4 == NULL)
        return NULL ;

    /* TODO */

    idx_order = cpl_malloc(4 * sizeof(int)) ;
    idx_order[0] = 0 ;
    idx_order[1] = 1 ;
    idx_order[2] = 2 ;
    idx_order[3] = 3 ;
    return idx_order ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Merge several POL_SPEC tables together
  @param    
  @param   frame2   Frame #2
  @param   frame3   Frame #3
  @param   frame4   Frame #4
  @return   an array of 4 integer indices or NULL in error case

  If the Out array is giving the new position : ordering[0] iѕ the
  position of frame1, ordering[1] is the position of frame2, etc...
  If the Frames order needs to be  Frame #2, #4, #3, #1, orderіng will
  contain [3, 0, 2, 1]

  The returned positions correspond to the 1,2,3,4 inputs in the demod
  fuctions.
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_pol_spec_pol_merge(
        const cpl_table **  pol_spec_list,
        int                 pol_spec_nb)
{
    int         i ;

    /* Check Inputs */
    if (pol_spec_list == NULL || pol_spec_nb <= 0) return NULL ;
    for (i=0 ; i<pol_spec_nb ; i++)
        if (pol_spec_list[i] == NULL) return NULL ;

    /* TODO */
    /* Currently return the first table */
    return cpl_table_duplicate(pol_spec_list[0]) ;
}

/**@}*/
