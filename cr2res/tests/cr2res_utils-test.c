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

#include <stdlib.h>
#include <string.h>
#include <cpl.h>
#include <cr2res_utils.h>

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static void test_cr2res_vector_get_rest(void) ;
static void test_cr2res_vector_get_int(void);

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_utils-test    Unit test of cr2res_utils
 *
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief
  @param
  @return
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_vector_get_int(void)
{
    int i;
    double d;
    int n=1000;
    cpl_vector * in = cpl_vector_new(n);
    int * res;

    for (i=0;i<n;i++) {
        d = (double)i ;
        cpl_vector_set(in, i, d + (d/(n+1)));
    }

    cpl_test( res = cr2res_vector_get_int(in) );

    for (i=0;i<n;i++) {
        cpl_test_eq(i, res[i]);
    }


    cpl_vector_delete(in);
    cpl_free(res);

    return;
}

/*----------------------------------------------------------------------------*/
/**
  @brief
  @param
  @return
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_vector_get_rest(void)
{
    int i;
    double d;
    int n=1000;
    cpl_vector * in = cpl_vector_new(n);
    cpl_vector * out = cpl_vector_new(n);
    double * res;

    for (i=0;i<n;i++) {
        d = (double)i ;
        cpl_vector_set(in, i, d + (d/(n+1)));
        cpl_vector_set(out, i, (d/(n+1)));
    }

    cpl_test( res = cr2res_vector_get_rest(in) );
    cpl_vector_delete(in);
    in = cpl_vector_wrap(n, res);
    cpl_test_vector_abs(in, out, DBL_EPSILON * n );

    cpl_vector_delete(in);
    cpl_vector_delete(out);

    return;
}

/*----------------------------------------------------------------------------*/
/**
  @brief
  @param
  @return
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_image_cut_rectify(void)
{
    return; // TODO: implement
}

/*----------------------------------------------------------------------------*/
/**
  @brief
  @param
  @return
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_polynomial_eval_vector(void)
{
    int i;
    double p0=1.1, p1=2.2, p2=3.3;
    double d, val;
    int n=1000;
    cpl_vector * in = cpl_vector_new(n);
    cpl_vector * out = cpl_vector_new(n);
    cpl_vector * res;
    cpl_polynomial * poly = cpl_polynomial_new(1);

    i=0; cpl_polynomial_set_coeff(poly, (cpl_size *)&i, p0 );
    i=1; cpl_polynomial_set_coeff(poly, (cpl_size *)&i, p1 );
    i=2; cpl_polynomial_set_coeff(poly, (cpl_size *)&i, p2 );

    for (i=0;i<n;i++) {
        d = (double)i ;
        val = d + (d/(n+1));
        cpl_vector_set(in, i, d);
        val = (p2*d*d)+(p1*d)+p0;
        cpl_vector_set(out, i, val);
    }

    cpl_test( res= cr2res_polynomial_eval_vector(poly, in));

    cpl_test_vector_abs(res, out, DBL_EPSILON * n * n * 10);

    cpl_vector_delete(in);
    cpl_vector_delete(out);
    cpl_vector_delete(res);
    cpl_polynomial_delete(poly);

    return;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Run the Unit tests
 */
/*----------------------------------------------------------------------------*/
int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_WARNING);

    test_cr2res_vector_get_rest() ;
    test_cr2res_vector_get_int() ;
    test_cr2res_polynomial_eval_vector();

    return cpl_test_end(0);
}

/**@}*/

