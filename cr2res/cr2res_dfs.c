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

#include <locale.h>
#include <string.h>
#include <math.h>

#include <cpl.h>

#include "cr2res_dfs.h"
#include "cr2res_io.h"
#include "cr2res_utils.h"

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_dfs  DFS related functions
 *
 * TBD
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Set the group as RAW or CALIB in a frameset
  @param    set     the input frameset
  @return   CPL_ERROR_NONE iff OK
 */
/*----------------------------------------------------------------------------*/
cpl_error_code cr2res_dfs_set_groups(cpl_frameset * set)
{

    cpl_frame   *   frame ;
    const char  *   tag ; 
    int             nframes, i ;
    
    /* Check entries */
    if (set == NULL) return -1 ;
    
    /* Initialize */
    nframes = cpl_frameset_get_size(set) ;
    
    /* Loop on frames */
    for (i = 0 ; i < nframes ; i++) {
        frame = cpl_frameset_get_position(set, i);
        tag = cpl_frame_get_tag(frame);

        if (tag == NULL) {
            cpl_msg_warning(cpl_func, "Frame %d has no tag", i);
        } else if (!strcmp(tag, CR2RES_COMMAND_LINE) ||
                /* DO.CATG Tags for the RAW files */
                !strcmp(tag, CR2RES_EMISSION_LINES_TXT_RAW) ||
                !strcmp(tag, CR2RES_DETLIN_RAW) ||
                !strcmp(tag, CR2RES_DARK_RAW) ||
                !strcmp(tag, CR2RES_FLAT_RAW) ||
                !strcmp(tag, CR2RES_WAVE_RAW) ||
                !strcmp(tag, CR2RES_OBS_NODDING_RAW) ||
                !strcmp(tag, CR2RES_OBS_2D_RAW) ||
                !strcmp(tag, CR2RES_OBS_POL_RAW) ||
                /* PRO.TYPE tags that can be used as input RAWs */
                /* For cr2res_util_bpm_split */
                !strcmp(tag, CR2RES_BPM_PROTYPE) ||
                /* For cr2res_util_extract */
                /* For cr2res_util_trace */
                /* For cr2res_util_normflat */
                !strcmp(tag, CR2RES_CALIBRATED_PROTYPE) ||
                /* For cr2res_util_wave */
                !strcmp(tag, CR2RES_EXTRACT_1D_PROTYPE) ||
                /* For cr2res_util_slit_curv  */
                /*     cr2res_util_trace_map  */
                !strcmp(tag, CR2RES_TW_PROTYPE)) {
            /* RAW frames */
            cpl_frame_set_group(frame, CPL_FRAME_GROUP_RAW);
        } else if (!strcmp(tag, CR2RES_CAL_DETLIN_COEFFS_PROCATG) ||
                /* Produced by cr2res_cal_detlin */
                !strcmp(tag, CR2RES_CAL_DETLIN_BPM_PROCATG) ||
                /* Produced by cr2res_cal_dark */
                !strcmp(tag, CR2RES_CAL_DARK_MASTER_PROCATG) ||
                !strcmp(tag, CR2RES_CAL_DARK_BPM_PROCATG) ||
                /* Produced by cr2res_cal_flat */
                !strcmp(tag, CR2RES_CAL_FLAT_BPM_PROCATG) ||
                !strcmp(tag, CR2RES_CAL_FLAT_EXTRACT_1D_PROCATG) ||
                !strcmp(tag, CR2RES_CAL_FLAT_SLIT_MODEL_PROCATG) ||
                !strcmp(tag, CR2RES_CAL_FLAT_SLIT_FUNC_PROCATG) ||
                !strcmp(tag, CR2RES_CAL_FLAT_TW_PROCATG) ||
                !strcmp(tag, CR2RES_CAL_FLAT_TW_MERGED_PROCATG) ||
                !strcmp(tag, CR2RES_CAL_FLAT_MASTER_PROCATG) ||
                /* Produced by cr2res_cal_wave */
                !strcmp(tag, CR2RES_CAL_WAVE_TW_PROCATG) ||
                !strcmp(tag, CR2RES_CAL_WAVE_MAP_PROCATG) ||
                !strcmp(tag, CR2RES_CAL_WAVE_LINES_DIAGNOSTICS_PROCATG) ||
                /* Produced by cr2res_obs_nodding */
                !strcmp(tag, CR2RES_OBS_NODDING_EXTRACTA_PROCATG) ||
                !strcmp(tag, CR2RES_OBS_NODDING_COMBINEDA_PROCATG) ||
                !strcmp(tag, CR2RES_OBS_NODDING_SLITFUNCA_PROCATG) ||
                !strcmp(tag, CR2RES_OBS_NODDING_SLITMODELA_PROCATG) ||
                !strcmp(tag, CR2RES_OBS_NODDING_EXTRACTB_PROCATG) ||
                !strcmp(tag, CR2RES_OBS_NODDING_COMBINEDB_PROCATG) ||
                !strcmp(tag, CR2RES_OBS_NODDING_SLITFUNCB_PROCATG) ||
                !strcmp(tag, CR2RES_OBS_NODDING_SLITMODELB_PROCATG) ||
                /* Produced by cr2res_obs_2d */
                !strcmp(tag, CR2RES_OBS_2D_EXTRACT_PROCATG) ||
                /* Produced by cr2res_obs_pol */
                !strcmp(tag, CR2RES_OBS_POL_SPECA_PROCATG) ||
                !strcmp(tag, CR2RES_OBS_POL_SPECB_PROCATG) ||
                /* Produced by cr2res_util_genlines */
                !strcmp(tag, CR2RES_EMISSION_LINES_PROCATG) ||
                /* Produced by cr2res_util_calib */
                !strcmp(tag, CR2RES_UTIL_CALIB_PROCATG) ||
                /* Produced by cr2res_util_bpm_split */
                !strcmp(tag, CR2RES_UTIL_BPM_SPLIT_PROCATG) ||
                /* Produced by cr2res_util_trace */
                !strcmp(tag, CR2RES_UTIL_TRACE_TW_PROCATG) ||
                /* Produced by cr2res_util_extract */
                !strcmp(tag, CR2RES_UTIL_SLIT_FUNC_PROCATG) ||
                !strcmp(tag, CR2RES_UTIL_SLIT_MODEL_PROCATG) ||
                !strcmp(tag, CR2RES_UTIL_EXTRACT_1D_PROCATG) ||
                /* Produced by cr2res_util_normflat */
                !strcmp(tag, CR2RES_UTIL_MASTER_FLAT_PROCATG) ||
                !strcmp(tag, CR2RES_UTIL_NORM_BPM_PROCATG) ||
                /* Produced by cr2res_util_slit_curv */
                !strcmp(tag, CR2RES_UTIL_SLIT_CURV_MAP_PROCATG) ||
                !strcmp(tag, CR2RES_UTIL_SLIT_CURV_PROCATG) ||
                !strcmp(tag, CR2RES_UTIL_SLIT_CURV_TW_PROCATG) ||
                /* Produced by cr2res_util_wave */
                !strcmp(tag, CR2RES_UTIL_WAVE_TW_PROCATG) ||
                !strcmp(tag, CR2RES_UTIL_WAVE_MAP_PROCATG) ||
                !strcmp(tag, CR2RES_UTIL_WAVE_XCORR_PROCATG) ||
                !strcmp(tag, CR2RES_UTIL_WAVE_LINES_DIAGNOSTICS_PROCATG) ||
                /* Produced by cr2res_util_trace_maps */
                !strcmp(tag, CR2RES_UTIL_TRACE_MAP_SLIT_CURVE_PROCATG) ||
                !strcmp(tag, CR2RES_UTIL_TRACE_MAP_WL_PROCATG) ||
                !strcmp(tag, CR2RES_UTIL_TRACE_MAP_TRACE_PROCATG) ||
                /* Produced by cr2res_util_splice */
                !strcmp(tag, CR2RES_UTIL_SPLICE_SPLICED_1D_PROCATG)) {
            /* CALIB frames */
            cpl_frame_set_group(frame, CPL_FRAME_GROUP_CALIB);
        }
    }
    return CPL_ERROR_NONE;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the SLIT_CURV table column name for a given order
  @param    order       The order number (1->) 
  @return   the column name or NULL in error case
  The return string needs to be deallocated with cpl_free() 
 */
/*----------------------------------------------------------------------------*/
char * cr2res_dfs_SLIT_CURV_colname(int order, int trace)
{
    int         order_loc ;
    if ((order_loc = cr2res_io_convert_order_to_idx(order)) < 0) return NULL ;
    return cpl_sprintf("%02d_%02d_%s", order_loc, trace, 
            CR2RES_COL_SLIT_CURV_SUFFIX);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the POL_WAVELENGTH column name for a given order
  @param    order       The order number (1->) 
  @return   the column name or NULL in error case
  The return string needs to be deallocated with cpl_free() 
 */
/*----------------------------------------------------------------------------*/
char * cr2res_dfs_POL_WAVELENGTH_colname(int order)
{
    int         order_loc ;
    if ((order_loc = cr2res_io_convert_order_to_idx(order)) < 0) return NULL ;
    return cpl_sprintf("%02d_%s", order_loc, CR2RES_COL_WL_SUFFIX) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the POL_STOKES column name for a given order
  @param    order       The order number (1->) 
  @return   the column name or NULL in error case
  The return string needs to be deallocated with cpl_free() 
 */
/*----------------------------------------------------------------------------*/
char * cr2res_dfs_POL_STOKES_colname(int order)
{
    int         order_loc ;
    if ((order_loc = cr2res_io_convert_order_to_idx(order)) < 0) return NULL ;
    return cpl_sprintf("%02d_%s", order_loc, CR2RES_COL_POL_STOKES_SUFFIX) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the POL_STOKES_ERROR column name for a given order
  @param    order       The order number (1->) 
  @return   the column name or NULL in error case
  The return string needs to be deallocated with cpl_free() 
 */
/*----------------------------------------------------------------------------*/
char * cr2res_dfs_POL_STOKES_ERROR_colname(int order)
{
    int         order_loc ;
    if ((order_loc = cr2res_io_convert_order_to_idx(order)) < 0) return NULL ;
    return cpl_sprintf("%02d_%s", order_loc, 
            CR2RES_COL_POL_STOKES_ERROR_SUFFIX) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the POL_NULL column name for a given order
  @param    order       The order number (1->) 
  @return   the column name or NULL in error case
  The return string needs to be deallocated with cpl_free() 
 */
/*----------------------------------------------------------------------------*/
char * cr2res_dfs_POL_NULL_colname(int order)
{
    int         order_loc ;
    if ((order_loc = cr2res_io_convert_order_to_idx(order)) < 0) return NULL ;
    return cpl_sprintf("%02d_%s", order_loc, CR2RES_COL_POL_NULL_SUFFIX) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the POL_NULL_ERROR column name for a given order
  @param    order       The order number (1->) 
  @return   the column name or NULL in error case
  The return string needs to be deallocated with cpl_free() 
 */
/*----------------------------------------------------------------------------*/
char * cr2res_dfs_POL_NULL_ERROR_colname(int order)
{
    int         order_loc ;
    if ((order_loc = cr2res_io_convert_order_to_idx(order)) < 0) return NULL ;
    return cpl_sprintf("%02d_%s", order_loc, CR2RES_COL_POL_NULL_ERROR_SUFFIX) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the POL_INTENS column name for a given order
  @param    order       The order number (1->) 
  @return   the column name or NULL in error case
  The return string needs to be deallocated with cpl_free() 
 */
/*----------------------------------------------------------------------------*/
char * cr2res_dfs_POL_INTENS_colname(int order)
{
    int         order_loc ;
    if ((order_loc = cr2res_io_convert_order_to_idx(order)) < 0) return NULL ;
    return cpl_sprintf("%02d_%s", order_loc, CR2RES_COL_POL_INTENS_SUFFIX) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the POL_INTENS_ERROR column name for a given order
  @param    order       The order number (1->) 
  @return   the column name or NULL in error case
  The return string needs to be deallocated with cpl_free() 
 */
/*----------------------------------------------------------------------------*/
char * cr2res_dfs_POL_INTENS_ERROR_colname(int order)
{
    int         order_loc ;
    if ((order_loc = cr2res_io_convert_order_to_idx(order)) < 0) return NULL ;
    return cpl_sprintf("%02d_%s", order_loc,
            CR2RES_COL_POL_INTENS_ERROR_SUFFIX) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the SPEC column name for a given order/trace
  @param    order       The order number (1->) 
  @param    trace       The trace number (1->)
  @return   the column name or NULL in error case
  The return string needs to be deallocated with cpl_free() 
 */
/*----------------------------------------------------------------------------*/
char * cr2res_dfs_SPEC_colname(int order, int trace)
{
    int         order_loc ;
    if ((order_loc = cr2res_io_convert_order_to_idx(order)) < 0) return NULL ;
    return cpl_sprintf("%02d_%02d_%s", order_loc, trace,
            CR2RES_COL_SPEC_SUFFIX);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the WAVELENGTH column name for a given order/trace
  @param    order       The order number (1->) 
  @param    trace       The trace number (1->)
  @return   the column name or NULL in error case
  The return string needs to be deallocated with cpl_free() 
 */
/*----------------------------------------------------------------------------*/
char * cr2res_dfs_WAVELENGTH_colname(int order, int trace)
{
    int         order_loc ;
    if ((order_loc = cr2res_io_convert_order_to_idx(order)) < 0) return NULL ;
    return cpl_sprintf("%02d_%02d_%s", order_loc, trace,
            CR2RES_COL_WL_SUFFIX);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the ERR column name for a given order/trace
  @param    order       The order number (1->) 
  @param    trace       The trace number (1->)
  @return   the column name or NULL in error case
  The return string needs to be deallocated with cpl_free() 
 */
/*----------------------------------------------------------------------------*/
char * cr2res_dfs_SPEC_ERR_colname(int order, int trace)
{
    int         order_loc ;
    if ((order_loc = cr2res_io_convert_order_to_idx(order)) < 0) return NULL ;
    return cpl_sprintf("%02d_%02d_%s", order_loc, trace,
            CR2RES_COL_ERROR_SUFFIX);
}
 
/*----------------------------------------------------------------------------*/
/**
  @brief    Get the SLIT_FUNC table column name for a given order/trace
  @param    order       The order number (1->)
  @param    trace       The trace number (1->)
  @return   the column name or NULL in error case
  The return string needs to be deallocated with cpl_free() 
 */
/*----------------------------------------------------------------------------*/
char * cr2res_dfs_SLIT_FUNC_colname(int order, int trace)
{
    int         order_loc ;
    if ((order_loc = cr2res_io_convert_order_to_idx(order)) < 0) return NULL ;
    return cpl_sprintf("%02d_%02d_%s", order_loc,trace,
            CR2RES_COL_SLIT_FUNC_SUFFIX);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the POSITIONX table column name for a given order/trace
  @param    order       The order number (1->)
  @param    trace       The trace number (1->)
  @return   the column name or NULL in error case
  The return string needs to be deallocated with cpl_free() 
 */
/*----------------------------------------------------------------------------*/
char * cr2res_dfs_POSITIONX_colname(int order, int trace)
{
    int         order_loc ;
    if ((order_loc = cr2res_io_convert_order_to_idx(order)) < 0) return NULL ;
    return cpl_sprintf("%02d_%02d_%s", order_loc,trace,
            CR2RES_COL_POSITIONX_SUFFIX);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the POSITIONY table column name for a given order/trace
  @param    order       The order number (1->)
  @param    trace       The trace number (1->)
  @return   the column name or NULL in error case
  The return string needs to be deallocated with cpl_free() 
 */
/*----------------------------------------------------------------------------*/
char * cr2res_dfs_POSITIONY_colname(int order, int trace)
{
    int         order_loc ;
    if ((order_loc = cr2res_io_convert_order_to_idx(order)) < 0) return NULL ;
    return cpl_sprintf("%02d_%02d_%s", order_loc,trace,
            CR2RES_COL_POSITIONY_SUFFIX);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the SLIT_FRACTION table column name for a given order/trace
  @param    order       The order number (1->)
  @param    trace       The trace number (1->)
  @return   the column name or NULL in error case
  The return string needs to be deallocated with cpl_free() 
 */
/*----------------------------------------------------------------------------*/
char * cr2res_dfs_SLIT_FRACTION_colname(int order, int trace)
{
    int         order_loc ;
    if ((order_loc = cr2res_io_convert_order_to_idx(order)) < 0) return NULL ;
    return cpl_sprintf("%02d_%02d_%s", order_loc,trace,
            CR2RES_COL_SLIT_FRACTION_SUFFIX);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Parse a column name ORDER_TRACE_TYPE format
  @param    colname     The column name to parse
  @param    order       [out] The order number (1->) 
  @param    trace       [out] The trace number (1->)
  @return   the column TYPE or NULL in error case
  The return string needs to be deallocated with cpl_free() 
 */
/*----------------------------------------------------------------------------*/
char * cr2res_dfs_SPEC_colname_parse(
        const char  * colname, 
        int         * order, 
        int         * trace)
{
    char    col_type[1024] ;
    if (colname == NULL || order == NULL || trace == NULL) return NULL ;

    /* Needed for sscanf() */
    setlocale(LC_NUMERIC, "C");

    if (sscanf(colname, "%02d_%02d_%s", order, trace, col_type) != 3)
        return NULL ;
    return cpl_strdup(col_type) ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Create an empty LINES DIAGNOSTICS table
  @param    nrows       The wished number of rows
  @return   a new table
 */
/*----------------------------------------------------------------------------*/
cpl_table * cr2res_dfs_create_lines_diagnostics_table(int nrows)
{
    cpl_table   *   out ;

    /* Check entries */
    if (nrows < 1) return NULL;

    out = cpl_table_new(nrows) ;

    cpl_table_new_column(out, CR2RES_COL_ORDER, CPL_TYPE_INT) ;
    cpl_table_new_column(out, CR2RES_COL_TRACENB, CPL_TYPE_INT);
    cpl_table_new_column(out, CR2RES_COL_MEASURED_LAMBDA, CPL_TYPE_DOUBLE);
    cpl_table_new_column(out, CR2RES_COL_CATALOG_LAMBDA, CPL_TYPE_DOUBLE);
    cpl_table_new_column(out, CR2RES_COL_DELTA_LAMBDA, CPL_TYPE_DOUBLE);
    cpl_table_new_column(out, CR2RES_COL_MEASURED_PIXEL, CPL_TYPE_DOUBLE);
    cpl_table_new_column(out, CR2RES_COL_LINE_WIDTH, CPL_TYPE_DOUBLE);
    cpl_table_new_column(out, CR2RES_COL_FIT_QUALITY, CPL_TYPE_DOUBLE);
    cpl_table_new_column(out, CR2RES_COL_INTENSITY, CPL_TYPE_DOUBLE);
    return out ;
}
/*----------------------------------------------------------------------------*/
/**
  @brief    Check completeness of trace table
  @param    trace       The trace table to check
  @return   1 if complete, 0 if not, -1 in error case
 */
/*----------------------------------------------------------------------------*/
int cr2res_dfs_check_traces_table(const cpl_table * traces)
{
    /* Check entries */
    if (traces == NULL) return -1 ;

    /* Check completeness */
    if (!cpl_table_has_column(traces, CR2RES_COL_UPPER) ||
            !cpl_table_has_column(traces, CR2RES_COL_LOWER) ||
            !cpl_table_has_column(traces, CR2RES_COL_ALL) ||
            !cpl_table_has_column(traces, CR2RES_COL_WAVELENGTH) ||
            !cpl_table_has_column(traces, CR2RES_COL_WAVELENGTH_ERROR) ||
            !cpl_table_has_column(traces, CR2RES_COL_ORDER) ||
            !cpl_table_has_column(traces, CR2RES_COL_TRACENB) ||
            !cpl_table_has_column(traces, CR2RES_COL_SLIT_CURV_A) ||
            !cpl_table_has_column(traces, CR2RES_COL_SLIT_CURV_B) ||
            !cpl_table_has_column(traces, CR2RES_COL_SLIT_CURV_C) ||
            !cpl_table_has_column(traces, CR2RES_COL_SLIT_FRACTION)) {
        return 0 ;
    }
    return 1 ;
}

/**@}*/
