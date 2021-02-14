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
#include <math.h>
#include <cpl.h>
#include <hdrl.h>
#include <cr2res_dfs.h>
#include <cr2res_detlin.h>
#include <cr2res_utils.h>

#define CR2RES_DETECTOR_SIZE            2048

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_detlin-test    Unit test of cr2res_cal_detlin
 *
 */
/*----------------------------------------------------------------------------*/
/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Test detlin correction polynomial
  @param
  @return foo

  foo
 */
/*----------------------------------------------------------------------------*/
static void test_cr2res_detlin_compute(){
 /*int cr2res_detlin_compute(
        const cpl_vector    *   dits,
        const cpl_vector    *   values,
        cpl_size                max_degree,
        cpl_polynomial      **  fitted,
        cpl_vector          **  error) */

  int n = 100;
  double i;
  cpl_size pow;
  double aduPsec = 666.6;
  double adu;
  cpl_size max_degree=5;
  cpl_vector    *   dits  = cpl_vector_new(n);
  cpl_vector    *   adus = cpl_vector_new(n);
  cpl_polynomial      *  poly_fitted;
  cpl_polynomial      *  poly_nonlin = cpl_polynomial_new(1);
  cpl_vector          *  error;
  cpl_vector          *  adus_corr;

pow=0;
cpl_polynomial_set_coeff(poly_nonlin,&pow,0.0);
pow=1;
cpl_polynomial_set_coeff(poly_nonlin,&pow,0.0);
pow=2;
cpl_polynomial_set_coeff(poly_nonlin,&pow,5E-5);

/* Fill data */
for (i=0; i<n; i++){
  cpl_vector_set(dits, i, i+1);
  adu=(i+1)*aduPsec / (1+cpl_polynomial_eval_1d(poly_nonlin,i+1,NULL));
  cpl_vector_set(adus, i, adu);
}

  cr2res_detlin_compute(dits, adus, max_degree, &poly_fitted, &error);
  adus_corr = cr2res_polynomial_eval_vector(poly_fitted,adus);
  cpl_vector_dump(adus_corr,stdout);
  cpl_vector_multiply(adus_corr,adus);
  cpl_polynomial_dump(poly_fitted, stdout);
  
    /* Check if we are close to true aduPsec*/
  cpl_vector_divide_scalar(adus_corr,aduPsec);
  cpl_vector_dump(adus_corr,stdout);
  cpl_test_vector_abs(adus_corr,dits,1);

  cpl_vector_delete(dits);
  cpl_vector_delete(adus);
  cpl_vector_delete(adus_corr);
  cpl_polynomial_delete(poly_nonlin);
  if (error != NULL) cpl_vector_delete(error);
  if (poly_fitted != NULL) cpl_polynomial_delete(poly_fitted);
  return ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Run the Unit tests
 */
/*----------------------------------------------------------------------------*/
int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_DEBUG);

    test_cr2res_detlin_compute();
    
    return cpl_test_end(0);
}

/**@}*/
