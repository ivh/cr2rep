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
#include <cr2res_pol.h>

#define CR2RES_DETECTOR_SIZE            2048

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/


static void test_cr2res_pol_demod_stokes(void);

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_wave-test    Unit test of cr2res_pol
 *
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
static void test_cr2res_pol_demod_stokes(){
    int n = 8; // as it has to be
    cpl_bivector ** speclist = cpl_malloc(n * sizeof(cpl_bivector*));
    cpl_bivector * pol;
    double value = 0;
    double error = 1. / sqrt(8);

    // spec = sigma = 1
    // -> pol = 0
    // -> sigma pol = 1 / sqrt(8)
    cpl_bivector * spec = cpl_bivector_new(CR2RES_DETECTOR_SIZE);
    for (cpl_size i = 0; i < CR2RES_DETECTOR_SIZE; i++)
    {
        cpl_vector_set(cpl_bivector_get_x(spec), i, 1);
        cpl_vector_set(cpl_bivector_get_y(spec), i, 1);
    }
    
    // just use the same spectrum 8 times
    for (int i = 0; i < n; i++)
    {
        speclist[i] = spec;
    }

    cpl_test_nonnull(pol = cr2res_pol_demod_stokes(speclist, n));

    // Check that the results are as expected    
    for (cpl_size i = 0; i < CR2RES_DETECTOR_SIZE; i++)
    {
        cpl_test_abs(value, cpl_vector_get(cpl_bivector_get_x(pol), i), DBL_EPSILON);
        cpl_test_abs(error, cpl_vector_get(cpl_bivector_get_y(pol), i), DBL_EPSILON);
    }

    // Clean memory
    cpl_bivector_delete(pol);
    cpl_bivector_delete(spec);
    cpl_free(speclist);
}


/*----------------------------------------------------------------------------*/
/**
  @brief    Run the Unit tests
 */
/*----------------------------------------------------------------------------*/
int main(void)
{
    cpl_test_init(PACKAGE_BUGREPORT, CPL_MSG_DEBUG);

    test_cr2res_pol_demod_stokes();

    return cpl_test_end(0);
}

/**@}*/
