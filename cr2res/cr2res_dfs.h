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

#ifndef CR2RES_DFS_H
#define CR2RES_DFS_H

/*-----------------------------------------------------------------------------
                                   Define
 -----------------------------------------------------------------------------*/

/******************************************/
/* DEFINE HERE THE DIFFERENT COLUMN NAMES */
/*** EMISSION_LINES Table ***/
#define CR2RES_COL_EMISSION         "Emission"      /* No Unit */
/*** EMISSION_LINES and TRACE_WAVE Table ***/
#define CR2RES_COL_WAVELENGTH       "Wavelength"    /* In Nanometers */
/*** TRACE_WAVE Table ***/
#define CR2RES_COL_UPPER            "Upper"         /* pixel position */
#define CR2RES_COL_LOWER            "Lower"         /* pixel position */
#define CR2RES_COL_ALL              "All"           /* pixel position */
#define CR2RES_COL_ORDER            "Order"         /* 1 .. nborders */
#define CR2RES_COL_TRACENB          "TraceNb"       /* 1 .. traces */
/*** Clusters Table ***/
#define CR2RES_COL_XS               "xs"            /* pixel position */
#define CR2RES_COL_YS               "ys"            /* pixel position */
#define CR2RES_COL_CLUSTERS         "clusterѕ"      /* cluster label */

/* SLIT_FUNC Table */
#define CR2RES_COL_SLIT_FUNC_SUFFIX "SLIT_FUNC"     /* Intensity */

/* EXTRACT_1D Table */
/* BLAZE Table */
#define CR2RES_COL_SPEC_SUFFIX      "SPEC"          /* Intensity */

/* TILT_POLY Table */
#define CR2RES_COL_TILT_SUFFIX      "TILT"          /* Polynomial */

/* SPLICED_1D Table */
/* EXTRACT_2D Table */
/* EXTRACT_POL Table */
/******************************************/

/*************************************/
/* Define here the PRO.TYPE keywords */

/* For  CR2RES_EMISSION_LINES_PROCATG */
/* Table with columns CR2RES_COL_EMISSION / CR2RES_COL_WAVELENGTH */
#define CR2RES_PROTYPE_CATALOG          "CATALOG"

/* For CR2RES_DETLIN_COEFFS_PROCATG */
#define CR2RES_DETLIN_COEFFS_PROTYPE    "DETLIN_COEFFS"

/* For  CR2RES_FLAT_BPM_PROCATG */
/*      CR2RES_DARK_BPM_PROCATG */
/* BPM Image with values CR2RES_BPM_DARK / CR2RES_BPM_FLAT... */
#define CR2RES_BPM_PROTYPE              "BPM"

/* For CR2RES_MASTER_DARK_PROCATG */
#define CR2RES_MASTER_DARK_PROTYPE      "MASTER_DARK"

/* For  CR2RES_BLAZE_PROCATG */
/*      CR2RES_EXTRACT_1D_PROCATG */
/* Table with columns cr2res_dfs_SPEC_colname() */
#define CR2RES_EXTRACT_1D_PROTYPE       "EXTRACT_1D"

/* For  CR2RES_FLAT_ЅLIT_MODEL_PROCATG */
/*      CR2RES_SLIT_MODEL_PROCATG */
#define CR2RES_SLIT_MODEL_PROTYPE       "SLIT_MODEL"

/* For  CR2RES_FLAT_SLIT_FUNC_PROCATG */
/*      CR2RES_UTIL_SLIT_FUNC_PROCATG */
/* Table with columns from cr2res_dfs_SLIT_FUNC_colname() */
#define CR2RES_SLIT_FUNC_PROTYPE        "SLIT_FUNC"

/* For  CR2RES_FLAT_TRACE_WAVE_PROCATG */
/*      CR2RES_UTIL_TRACE_WAVE_PROCATG */
/* Table with Traces polynomials, orders, trace_nb wavelengths */
/*                  1 trace per Row */
#define CR2RES_TRACE_WAVE_PROTYPE       "TRACE_WAVE"

/* For CR2RES_MASTER_FLAT_PROCATG */
/* Master Flat image with values around 1 */
#define CR2RES_MASTER_FLAT_PROTYPE      "MASTER_FLAT"

/* For CR2RES_CALIBRATED_PROCATG */
#define CR2RES_CALIBRATED_PROTYPE       "CALIBRATED"

/* For  CR2RES_TILT_COEFFS_PROCATG  */
/* Table with columns from cr2res_dfs_TILT_colname() */
#define CR2RES_TILT_COEFFS_PROTYPE      "TILT_COEFFS"

/* For  CR2RES_CALIB_COLLAPSED_PROCATG */
#define CR2RES_CALIB_COLLAPSED_PROTYPE  "CALIB_COLLAPSED"

/*************************************/
/* Define here the PRO.CATG keywords */
/*************************************/
/* Produced by cr2res_util_genlines */
#define CR2RES_EMISSION_LINES_PROCATG   "EMISSION_LINES"

/* Produced by cr2res_cal_detlin */
#define CR2RES_DETLIN_COEFFS_PROCATG    "DETLIN_COEFFS"

/* Produced by cr2res_cal_dark */
#define CR2RES_MASTER_DARK_PROCATG      "MASTER_DARK"
#define CR2RES_DARK_BPM_PROCATG         "DARK_BPM"

/* Produced by cr2res_cal_flat */
#define CR2RES_FLAT_BPM_PROCATG         "FLAT_BPM"
#define CR2RES_BLAZE_PROCATG            "BLAZE"
#define CR2RES_FLAT_SLIT_MODEL_PROCATG  "FLAT_SLIT_MODEL"
#define CR2RES_FLAT_SLIT_FUNC_PROCATG   "FLAT_SLIT_FUNC"
#define CR2RES_FLAT_TRACE_WAVE_PROCATG  "FLAT_TRACE_WAVE"
#define CR2RES_MASTER_FLAT_PROCATG      "MASTER_FLAT"

/* Produced by cr2res_util_calib */
#define CR2RES_CALIBRATED_PROCATG       "CALIBRATED"

/* Produced by cr2res_util_trace */
#define CR2RES_UTIL_TRACE_WAVE_PROCATG  "UTIL_TRACE_WAVE"

/* Produced by cr2res_util_extract */
#define CR2RES_UTIL_SLIT_FUNC_PROCATG   "UTIL_SLIT_FUNC"
#define CR2RES_UTIL_SLIT_MODEL_PROCATG  "UTIL_SLIT_MODEL"

/* TODO */
/* Produced by cr2res_util_tilt */
#define CR2RES_TILT_COEFFS_PROCATG      "TILT_COEFFS"

/* Produced by cr2res_util_normflat */
/* Produced by cr2res_util_wave */

/* TODO */
#define CR2RES_WAVE_COEFFS_PROCATG      "WAVE_COEFFS"
#define CR2RES_WAVE_MAP_PROCATG         "WAVE_MAP"
#define CR2RES_EXTRACT_1D_PROCATG       "EXTRACT_1D"

/* Define here the DO.CATG keywords */
#define CR2RES_COMMAND_LINE             "COMMAND_LINE"
#define CR2RES_DETLIN_RAW               "DETLIN"
#define CR2RES_DARK_RAW                 "DARK"
#define CR2RES_FLAT_RAW                 "FLAT"
#define CR2RES_CALIB_RAW                "CALIB"

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

cpl_error_code cr2res_dfs_set_groups(cpl_frameset *);
char * cr2res_dfs_SPEC_colname(int, int) ;
char * cr2res_dfs_SLIT_FUNC_colname(int, int) ;
char * cr2res_dfs_TILT_colname(int) ;

#endif
