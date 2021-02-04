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
  double aduPsec = 1000.1;
  cpl_size max_degree=3;
  cpl_vector    *   dits  = cpl_vector_new(n);
  cpl_vector    *   values = cpl_vector_new(n);
  cpl_polynomial      *  fitted;
  cpl_vector          *  error;




  cr2res_detlin_compute(dits, values, max_degree, &fitted, &error);

  //cpl_vector_dump(fitted);

  cpl_vector_delete(dits);
  cpl_vector_delete(values);
  if (fitted != NULL) cpl_vector_delete(error);
  if (error != NULL) cpl_polynomial_delete(fitted);
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
