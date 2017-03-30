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
    cpl_table   *   traces ;
	int				poly_order, norders ;
    cpl_array   *   array ;
    cpl_image   *   trace_ima ;
    cpl_image   *   ima ;
    cpl_table   *   out ;
    int             i ;

    /* Initialise */
    poly_order = 2 ;
    norders = 9 ;

    /* NULL Input */
    ima = NULL ;
    out = cr2res_trace_cpl(ima, CR2RES_DECKER_NONE, 1.0, 1, 6, 500) ;
    cpl_test_null(out);

    /* Create a label image with 9 orders */
    traces = cpl_table_new(norders) ;
   	cpl_table_new_column_array(traces, "All", CPL_TYPE_DOUBLE, poly_order) ;
    cpl_table_new_column_array(traces, "Upper", CPL_TYPE_DOUBLE, poly_order) ;
    cpl_table_new_column_array(traces, "Lower", CPL_TYPE_DOUBLE, poly_order) ;
    cpl_table_new_column(traces, "Order", CPL_TYPE_INT) ;

    /*
                 All|               Upper|               Lower|  Order
 54.3148, 0.00738623|  110.379, 0.0159501|1.63328, -8.95809e-05|      1
  226.289, 0.0169145|  311.885, 0.0167957|   139.106, 0.017396|      2
  437.881, 0.0172448|  524.126, 0.0171958|  350.398, 0.0170009|      3
  660.954, 0.0178211|  747.931, 0.0179547|  572.986, 0.0177137|      4
  897.266, 0.0185319|   985.06, 0.0186544|  808.823, 0.0185517|      5
   1148.31, 0.019388|  1237.01, 0.0193813|  1059.39, 0.0194215|      6
  1415.88, 0.0202877|  1505.34, 0.0202763|  1326.43, 0.0202534|      7
   1701.85, 0.021292|  1792.03, 0.0213662|  1611.65, 0.0212178|      8
  1982.59, 0.0111388|2047.88, 8.66878e-05|   1917.3, 0.0221835|      9

  */
    double all_1[] = {54.3148, 226.289, 437.881, 660.954, 897.266,
           1148.31, 1415.88, 1701.85, 1982.59} ;
    double all_2[] = {0.00738623, 0.0169145, 0.0172448, 0.0178211,
           0.0185319, 0.019388, 0.0202877, 0.021292, 0.0111388} ;
    double upper_1[] = {110.379, 311.885, 524.126, 747.931, 985.06,
           1237.01, 1505.34, 1792.03, 2047.88} ;
    double upper_2[] = {0.0159501, 0.0167957, 0.0171958, 0.0179547,
           0.0186544, 0.0193813, 0.0202763, 0.0213662, 8.66878e-05} ;
    double lower_1[] = {1.63328, 139.106, 350.398, 572.986, 808.823,
           1059.39, 1326.43, 1611.65, 1917.3} ;
    double lower_2[] = {-8.95809e-05, 0.017396, 0.0170009, 0.0177137,
           0.0185517, 0.0194215, 0.0202534, 0.0212178, 0.0221835} ;
    array = cpl_array_new(poly_order, CPL_TYPE_DOUBLE) ;
    for (i=0 ; i<norders ; i++) {
        cpl_array_set(array, 0, all_1[i]) ; cpl_array_set(array, 1, all_2[i]) ;
        cpl_table_set_array(traces, "All", i, array);
        cpl_array_set(array, 0, upper_1[i]); cpl_array_set(array,1,upper_2[i]);
        cpl_table_set_array(traces, "Upper", i, array);
        cpl_array_set(array, 0, lower_1[i]); cpl_array_set(array,1,lower_2[i]);
        cpl_table_set_array(traces, "Lower", i, array);
        cpl_table_set(traces, "Order", i, i+1);
    }
    cpl_array_delete(array) ;

    trace_ima = cr2res_trace_gen_image(traces, 2048, 2048) ;

    cpl_table_delete(traces) ;

    cpl_image_save(trace_ima, "TEST.fits", CPL_TYPE_INT, NULL, CPL_IO_CREATE);

    out = cr2res_trace_cpl(trace_ima, CR2RES_DECKER_NONE, 1.0, 1, 6, 500) ;
    
    cpl_table_save(out, NULL, NULL, "TEST2.fits", CPL_IO_CREATE);
    cpl_table_delete(out) ;
    /* cpl_test_null(out); */

    cpl_image_delete(trace_ima) ;

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

