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
#include <hdrl.h>
#include <cr2res_dfs.h>
#include <cr2res_extract.h>
#include <cr2res_trace.h>

static void test_cr2res_slitdec_vert(void);


static cpl_table *create_test_table()
{
    int poly_order = 2;
    int n_orders = 2;
    cpl_table * traces = cpl_table_new(n_orders);
    cpl_table_new_column_array(traces, CR2RES_COL_ALL, CPL_TYPE_DOUBLE, poly_order);
    cpl_table_new_column_array(traces, CR2RES_COL_UPPER, CPL_TYPE_DOUBLE, poly_order);
    cpl_table_new_column_array(traces, CR2RES_COL_LOWER, CPL_TYPE_DOUBLE, poly_order);
    cpl_table_new_column(traces, CR2RES_COL_ORDER, CPL_TYPE_INT);
    cpl_table_new_column(traces, CR2RES_COL_TRACENB, CPL_TYPE_INT);

    double all_1[] = {86.6279, 175.5738};
    double all_2[] = {0.01699, 0.07512};
    double upper_1[] = {108.5065, 197.3485};
    double upper_2[] = {0.016601, 0.0184364};
    double lower_1[] = {64.05477, 153.7987};
    double lower_2[] = {0.017355, 0.01659297};

    cpl_array * array = cpl_array_new(poly_order, CPL_TYPE_DOUBLE);
    for (int i = 0; i < n_orders; i++)
    {
        cpl_array_set(array, 0, all_1[i]);
        cpl_array_set(array, 1, all_2[i]);
        cpl_table_set_array(traces, CR2RES_COL_ALL, i, array);
        cpl_array_set(array, 0, upper_1[i]);
        cpl_array_set(array, 1, upper_2[i]);
        cpl_table_set_array(traces, CR2RES_COL_UPPER, i, array);
        cpl_array_set(array, 0, lower_1[i]);
        cpl_array_set(array, 1, lower_2[i]);
        cpl_table_set_array(traces, CR2RES_COL_LOWER, i, array);
        cpl_table_set(traces, CR2RES_COL_ORDER, i, 7);
        cpl_table_set(traces, CR2RES_COL_TRACENB, i, i + 1);
    }

    cpl_array_delete(array);
    return traces;
}

static cpl_image *create_test_image()
{
    cpl_image *img = cpl_image_load("./tests/cr2res_utils_test_image.fits", CPL_TYPE_INT, 0, 1);
    return img;
}



static void test_cr2res_slitdec_vert(void)
{
    cpl_image * img_in = create_test_image();
    cpl_table * trace_table = create_test_table();
    int order = 7;
    int trace = 1;
    int height = 70;
    int swath = 50;
    int oversample = 1;
    double smooth_slit = 0.1;
    int width = cpl_image_get_size_x(img_in);

    cpl_vector ** slit_func;
    cpl_vector ** spec;
    hdrl_image ** model;

    cr2res_extract_slitdec_vert(img_in, trace_table, order, trace, height, swath, oversample, smooth_slit, slit_func, spec, model);

    cpl_image_delete(img_in);
    cpl_table_delete(trace_table);
    
    for(int i = 0; i < width/swath+2; i++)
    {
        cpl_vector_delete(slit_func[i]);
        cpl_vector_delete(spec[i]);
        hdrl_image_delete(model[i]);
    }
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Run the Unit tests
 */
/*----------------------------------------------------------------------------*/
int main(void)
{
    cpl_test_init("test@bugreport.se", CPL_MSG_DEBUG);

    test_cr2res_slitdec_vert();

    return cpl_test_end(0);
}
