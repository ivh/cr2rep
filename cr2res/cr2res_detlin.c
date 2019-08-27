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
#include "cr2res_pfits.h"

/*-----------------------------------------------------------------------------
                                   	Defines
 -----------------------------------------------------------------------------*/

#define pow2(x) (x)*(x)
#define pow3(x) (x)*(x)*(x)

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_detlin_frames_dit_compare(
        const cpl_frame *   in1,
        const cpl_frame *   in2) ;
static void cr2res_matrix_fill_normal_vandermonde(cpl_matrix * self,
                                               cpl_matrix * mx,
                                               const cpl_vector * xhat,
                                               cpl_boolean is_eqdist,
                                               cpl_size mindeg,
                                               const cpl_vector * values);

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

        perr[i] = pow2(perra[i] * pdata[i]) + pow2(perrb[i] * pow2(pdata[i]))
                + pow2(perrc[i] * pow3(pdata[i])) 
                + pow2(perr[i] * (pima[i] + 2. * pimb[i] * pdata[i] 
                + 3. * pimc[i] * pow2(pdata[i])));
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
int cr2res_detlin_compute(
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
    if (fitted == NULL || dits == NULL || values == NULL) return -1 ;
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
        cpl_vector_delete(y);
        cpl_polynomial_delete(fitted_local) ;
        cpl_error_reset() ;
        return -1 ;
    }

    /* Compute the error */
    error_local = cpl_vector_new(max_degree+1) ;

    cpl_size nc = max_degree + 1;
    cpl_size ndata = cpl_vector_get_size(y);
    if (nc >= ndata){
        // No uncertainty as fit should be perfectly aligned with data points
        for (cpl_size i = 0; i < max_degree + 1; i++)
        {
            cpl_vector_set(error_local, i, 0);
        }
    } else {
        // lhs = vandermode(x, order)
        // hankel = dot(lhs.T, lhs)
        // cov = inv(hankel)
        // cov *= resids / (len(x) - order)
        // error_local = diag(cov)
        cpl_matrix * hankel = cpl_matrix_new(nc, nc);
        cpl_matrix * mx = cpl_matrix_new(nc, 1); // just a temporary matrix

        // this actually returns the hankel matrix, not vandermode
        // also directly copied from the cpl source code (cpl_polynomial.c)
        cr2res_matrix_fill_normal_vandermonde(hankel, mx, values, CPL_FALSE, 
            0, y);
        cpl_matrix * inverse = cpl_matrix_invert_create(hankel);
        cpl_vector * resids = cpl_vector_new(ndata);

        cpl_vector_fill_polynomial_fit_residual(resids, y, NULL, fitted_local,
            samppos, NULL);
        cpl_matrix_multiply_scalar(inverse, 
            cpl_vector_get_sum(resids) / (double)(ndata - nc));

        for (cpl_size i = 0; i < max_degree + 1; i++) {
            cpl_vector_set(error_local, i, cpl_matrix_get(inverse, i, i));
        }
        cpl_matrix_delete(hankel);
        cpl_matrix_delete(mx);
        cpl_matrix_delete(inverse);
        cpl_vector_delete(resids);
    }
    cpl_vector_delete(y);
    cpl_matrix_unwrap(samppos) ;    

    /* Return */
    *fitted = fitted_local ;
    if (error != NULL) *error = error_local ;
    else cpl_vector_delete(error_local) ;
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Sort the frames by increaing DIT
  @param    in  The input frameset
  @return   the newly allocated sorted frameset
 */
/*----------------------------------------------------------------------------*/
cpl_frameset * cr2res_detlin_sort_frames(
        const cpl_frameset  *   in)
{
    cpl_frameset    *   sorted ;

    /* Check Inputs */
    if (in == NULL) return NULL ;

    sorted = cpl_frameset_duplicate(in) ;
    if (cpl_frameset_sort(sorted, 
                cr2res_detlin_frames_dit_compare) != CPL_ERROR_NONE) {
        cpl_frameset_delete(sorted) ;
        return NULL ;
    }
    return sorted ;
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Frames comparison using DIT values
  @param    in1 the first frame
  @param    in2 the second frame
  @return   -1, 0, or 1 if in1 is less than, equal or greater than in2
 */
/*----------------------------------------------------------------------------*/
static int cr2res_detlin_frames_dit_compare(
        const cpl_frame *   in1,
        const cpl_frame *   in2)
{
    cpl_propertylist    *   plist1 ;
    cpl_propertylist    *   plist2 ;
    double                  dit1, dit2 ;
    
    /* Test entries */
    if (in1==NULL || in2==NULL) return 0 ;

    /* Get property lists */
    plist1=cpl_propertylist_load(cpl_frame_get_filename(in1),0) ;
    plist2=cpl_propertylist_load(cpl_frame_get_filename(in2),0) ;

    /* Get DITs */
    dit1 = cr2res_pfits_get_dit(plist1) ;
    dit2 = cr2res_pfits_get_dit(plist2) ;
    if (plist1 != NULL) cpl_propertylist_delete(plist1) ;
    if (plist2 != NULL) cpl_propertylist_delete(plist2) ;
    if (cpl_error_get_code()) return 0 ;

    if (dit1<dit2)          return -1 ;
    else if (dit1>dit2)     return 1 ;
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @internal
  @brief    Fill the Hankel Matrix H=V'*V, where V is a 1D-Vandermonde matrix
  @param    self      The matrix H
  @param    mx        A right multiplication with V', mx = V' * values
  @param    xhat      The mean-transformed x-values
  @param    is_eqdist True iff xhat contains equidistant points
  @param    mindeg    The non-negative minimum fitting degree
  @param    values    The values to be interpolated
  @return   void
  @note self must have its elements initialized to zero iff is_eqdist is true.

 */
/*----------------------------------------------------------------------------*/
static void cr2res_matrix_fill_normal_vandermonde(cpl_matrix * self,
                                               cpl_matrix * mx,
                                               const cpl_vector * xhat,
                                               cpl_boolean is_eqdist,
                                               cpl_size mindeg,
                                               const cpl_vector * values)
{


    const double * dval = cpl_vector_get_data_const(values);
    const double * xval = cpl_vector_get_data_const(xhat);
    cpl_vector   * phat = cpl_vector_duplicate(xhat); /* Powers of xhat */
    cpl_vector   * qhat = NULL;                       /* mindeg Power of xhat */
    double       * dhat = cpl_vector_get_data(phat);
    double       * ehat = NULL;
    const cpl_size nc   = cpl_matrix_get_ncol(self);
    const cpl_size np   = cpl_vector_get_size(xhat);
    cpl_size       i,j;

    /* Fill Hankel matrix from top-left to main skew diagonal
       - on and above (non-skew) main diagonal */
    /* Also compute transpose(V) * b */
    /* Peel off 1st iteration */
    if (mindeg > 0) {
        double hsum = 0.0;
        cpl_size k;

        qhat = mindeg == 1 ? cpl_vector_duplicate(xhat) : cpl_vector_new(np);
        ehat = cpl_vector_get_data(qhat);

        /* Raise xhat to the power of mindeg */
        for (k=0; k < np; k++) {
            const double x = xval[k];

            if (mindeg > 1) ehat[k] = pow(x, (int)mindeg);
            dhat[k] *= ehat[k];

            hsum += ehat[k] * ehat[k];
        }
        cpl_matrix_set(self, 0, 0, hsum);
    } else {
        cpl_matrix_set(self, 0, 0, (double)np);
    }
    /* qhat is xhat to the power of mindeg, iff mindeg > 0 */
    /* dhat is xhat to the power of 1+mindeg, iff mindeg > 0 */
    for (j=1; j < 2; j++) {
        double vsum0 = 0.0;
        double hsum = 0.0;
        double vsum = 0.0;
        cpl_size k;

        for (k=0; k < np; k++) {
            const double y = dval[k];

            hsum += mindeg > 0 ? ehat[k] * dhat[k] : dhat[k];
            vsum += y * dhat[k];
            vsum0 += mindeg > 0 ? ehat[k] * y : y;
        }
        cpl_matrix_set(mx, 0, 0, vsum0);
        cpl_matrix_set(mx, j, 0, vsum);
        if (is_eqdist) continue;
        k = j;
        for (i=0; i <= k; i++, k--) {
            cpl_matrix_set(self, i, k, hsum);
        }
    }
    for (; j < nc; j++) {
        double   hsum = 0.0;
        double   vsum = 0.0;
        cpl_size k;

        for (k=0; k < np; k++) {
            const double x = xval[k];
            const double y = dval[k];

            dhat[k] *= x;
            hsum += mindeg > 0 ? ehat[k] * dhat[k] : dhat[k];
            vsum += y * dhat[k];
        }
        cpl_matrix_set(mx, j, 0, vsum);
        if (is_eqdist && (j&1)) continue;
        k = j;
        for (i=0; i <= k; i++, k--) {
            cpl_matrix_set(self, i, k, hsum);
        }
    }
    /* Fill remaining Hankel matrix - on and above (non-skew) main diagonal */

    if (mindeg > 0) {
        cpl_vector_multiply(phat, qhat);
        cpl_vector_delete(qhat);
    }

    for (i = 1; i < nc; i++) {
        cpl_size k;
        double   hsum = 0.0;

        if (is_eqdist && ((i+nc)&1)==0) {
            cpl_vector_multiply(phat, xhat);
            continue;
        }

        for (k=0; k < np; k++) {
            const double x = xval[k];

            dhat[k] *= x;
            hsum += dhat[k];
        }
        k = i;
        for (j = nc-1; k <= j; k++, j--) {
            cpl_matrix_set(self, k, j, hsum);
        }
    }

    cpl_vector_delete(phat);

}

