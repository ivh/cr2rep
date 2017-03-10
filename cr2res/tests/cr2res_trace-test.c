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
#include <cr2res_trace.h>

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

static void test_cr2res_trace_cpl(void) ;

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_trace-test    Unit test of cr2res_trace
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
static void test_cr2res_trace_cpl(void)
{

    cpl_image   *   ima ;
    cpl_table   *   out ;

    /* NULL Input */
    ima = NULL ;
    out = cr2res_trace_cpl(ima, CR2RES_DECKER_NONE, 1.0, 1, 6, 500) ;
    cpl_test_null(out);

    /* Empty Input */
    ima = cpl_image_new(2048, 2048, CPL_TYPE_FLOAT) ;
    /* out = cr2res_trace_cpl(ima, CR2RES_DECKER_NONE, 1.0, 1, 6, 500) ; */
    /* cpl_test_null(out); */

    cpl_image_delete(ima) ;

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

    test_cr2res_trace_cpl() ;
    /* test_cr2res_trace_nocpl() ; */
    /* test_cr2res_trace_detect() ; */
    /* test_cr2res_trace_labelize() ; */
    /* test_cr2res_trace_fit() ; */
    /* test_cr2res_trace_compare() ; */
    /* test_cr2res_trace_combine() ; */
    /* test_cr2res_trace_gen_image() ; */
    /* test_cr2res_trace_get_order_numbers() ; */
    /* test_cr2res_trace_open_get_polynomials() ; */
    /* test_cr2res_trace_compute_middle() ; */
    /* test_cr2res_trace_compute_height() ; */
        
    return cpl_test_end(0);
}

/**@}*/

