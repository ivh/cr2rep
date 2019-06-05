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

#include <math.h>
#include <cpl.h>
#include "cr2res_detlin.h"

/*-----------------------------------------------------------------------------
                                   	Defines
 -----------------------------------------------------------------------------*/

#define pow2(x) (x)*(x)
#define pow3(x) (x)*(x)*(x)

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @defgroup cr2res_detlin
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Apply the detector linearity correction
  @param    in      the input image 
  @param    detlin  the detlin coeffs
  @return   0 if everything is ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
int cr2res_detlin_correct(
        hdrl_image              *   in,
        const hdrl_imagelist    *   detlin)
{
    const cpl_image     *   ima ;
    const cpl_image     *   erra ;
    const cpl_image     *   imb ;
    const cpl_image     *   errb ;
    const cpl_image     *   imc ;
    const cpl_image     *   errc ;
    const double        *   pima ;
    const double        *   perra ;
    const double        *   pimb ;
    const double        *   perrb ;
    const double        *   pimc ;
    const double        *   perrc ;
    cpl_image           *   cur_ima ;
    double              *   pdata ;
    double              *   perr ;
    int                     nx, ny ;
    double                  correction_factor;
    int                     i, j ;

    /* Test entries */
    if (!in || !detlin) return -1 ;

    /* Initialise */
    pdata = cpl_image_get_data_double(hdrl_image_get_image(in)) ;
    perr = cpl_image_get_data_double(hdrl_image_get_error(in)) ;


    /* Load the 3 coeffs images */
    ima = hdrl_image_get_image_const(hdrl_imagelist_get_const(detlin, 0)) ;
    erra = hdrl_image_get_error_const(hdrl_imagelist_get_const(detlin, 0)) ;
    imb = hdrl_image_get_image_const(hdrl_imagelist_get_const(detlin, 1)) ;
    errb = hdrl_image_get_error_const(hdrl_imagelist_get_const(detlin, 1)) ;
    imc = hdrl_image_get_image_const(hdrl_imagelist_get_const(detlin, 2)) ;
    errc = hdrl_image_get_error_const(hdrl_imagelist_get_const(detlin, 2)) ;
    
    if (!ima || !imb || !imc) {
        cpl_msg_error(cpl_func, "Cannot access the detlin images") ;
        return -1 ;
    }
    pima = cpl_image_get_data_double_const(ima) ;
    pimb = cpl_image_get_data_double_const(imb) ;
    pimc = cpl_image_get_data_double_const(imc) ;
    perra = cpl_image_get_data_double_const(erra);
    perrb = cpl_image_get_data_double_const(errb);
    perrc = cpl_image_get_data_double_const(errc);

    /* Test sizes */
    cur_ima = hdrl_image_get_image(in) ;
    nx = cpl_image_get_size_x(cur_ima) ;
    ny = cpl_image_get_size_y(cur_ima) ;
    if ((cpl_image_get_size_x(ima) != nx) ||
            (cpl_image_get_size_x(imb) != nx) ||
            (cpl_image_get_size_x(imc) != nx) ||
            (cpl_image_get_size_y(ima) != ny) ||
            (cpl_image_get_size_y(imb) != ny) ||
            (cpl_image_get_size_y(imc) != ny)) {
        cpl_msg_error(cpl_func, "Incompatible sizes") ;
        return -1 ;
    }

    /* Loop on pixels */
    for (i=0 ; i<nx*ny ; i++) {
        // for each pixel p' = a + b * p + c * p * p

        perr[i] = pow2(perra[i] * pdata[i]) + pow2(perrb[i] * pow2(pdata[i])) + pow2(perrc[i] * pow3(pdata[i])) 
                + pow2(perr[i] * (pima[i] + 2. * pimb[i] * pdata[i] + 3. * pimc[i] * pow2(pdata[i])));
        perr[i] = sqrt(perr[i]);

        correction_factor = pima[i] + (pimb[i] + pimc[i] * pdata[i]) * pdata[i];
        pdata[i] = pdata[i] * correction_factor;
    }
    /* return */
    return 0 ;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Fits the response of a given pixel to the illumination increase
  @param    dits        Vector with the DIT values
  @param    values      Vector with illumination values
  @param    max_degree  Maximum degree for the fit 
  @param    fitted      [out] The fitted polynomial
  @param    error       [out] the errors vector
  @return   0 if ok, -1 in error case

  The input dits and values vectors must have the same size
  The *fitted polynomial coefficients are the values stored in the
  DETLIN_COEFFS product for a given pixel.
  The *error vector size must match the *fitted polynomial number of
  coefficients. Its values are stored in the error extension of the 
  DETLIN_COEFFS product.
 */
/*----------------------------------------------------------------------------*/
int cr2res_cal_detlin_fit(
        const cpl_vector    *   dits,
        const cpl_vector    *   values,
        cpl_size                max_degree,
        cpl_polynomial      **  fitted,
        cpl_vector          **  error)
{
    cpl_matrix          *   samppos ;
    cpl_boolean             sampsym ;
    cpl_polynomial      *   fitted_local ;
    cpl_vector          *   error_local ;
    cpl_vector          *   y;
    cpl_bivector        *   first_dits, *tmp;
    double   x0, x1, y0, y1, m;
    double expected_counts;

    /* Test entries */
    if (fitted == null || dits == NULL || values == NULL) return -1 ;
    if (cpl_vector_get_size(dits) != cpl_vector_get_size(values))
        return -1 ;

    /* Initialise */
    sampsym = CPL_TRUE ;

    // Find linear coefficient from the first two dits
    first_dits = cpl_bivector_new(cpl_vector_get_size(dits));
    tmp = cpl_bivector_wrap_vectors((cpl_vector*)dits, (cpl_vector*)values);
    cpl_bivector_sort(first_dits, tmp, CPL_SORT_ASCENDING, CPL_SORT_BY_X);
    cpl_bivector_unwrap_vectors(tmp);

    // x: dits, y: values
    x0 = cpl_bivector_get_x_data(first_dits)[0];
    x1 = cpl_bivector_get_x_data(first_dits)[1];
    y0 = cpl_bivector_get_y_data(first_dits)[0];
    y1 = cpl_bivector_get_y_data(first_dits)[1];
    m = (y1-y0)/(x1-x0);

    cpl_bivector_delete(first_dits);

    /* Store the DITS */
    samppos = cpl_matrix_wrap(1,
                cpl_vector_get_size(values),
                cpl_vector_get_data((cpl_vector*)values)) ;

    y = cpl_vector_new(cpl_vector_get_size(dits));
    for(cpl_size i = 0; i < cpl_vector_get_size(dits); i++)
    {
        expected_counts = y0 + m * (cpl_vector_get(dits, i) - x0);
        cpl_vector_set(y, i, expected_counts / cpl_vector_get(values, i));
    }


    /* Fit  */
    fitted_local = cpl_polynomial_new(1);
    if (cpl_polynomial_fit(fitted_local, samppos, &sampsym, y, NULL,
            CPL_FALSE, NULL, &max_degree) != CPL_ERROR_NONE) {
        /* Failed Fit - Fill the coefficientѕ */
        cpl_matrix_unwrap(samppos) ;
        cpl_polynomial_delete(fitted_local) ;
        cpl_error_reset() ;
        return -1 ;
    }
    cpl_matrix_unwrap(samppos) ;
    cpl_vector_delete(y);

    /* Compute the error */
    /* TODO */
    error_local = cpl_vector_new(max_degree+1) ;






    /* Return */
    *fitted = fitted_local ;
    if (error != NULL) *error = error_local ;
    else cpl_vector_delete(error_local) ;
    return 0 ;
}

/**@}*/
