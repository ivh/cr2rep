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
#define CR2RES_COL_WAVELENGTH_ERROR "Wavelength_Error"  /* In Nanometers */
/*** TRACE_WAVE Table ***/
#define CR2RES_COL_UPPER            "Upper"         /* pixel position */
#define CR2RES_COL_LOWER            "Lower"         /* pixel position */
#define CR2RES_COL_ALL              "All"           /* pixel position */
#define CR2RES_COL_ORDER            "Order"         /* 1 .. nborders */
#define CR2RES_COL_TRACENB          "TraceNb"       /* 1 .. traces */
#define CR2RES_COL_SLIT_CURV_A      "SlitPolyA"     /* 3 Coefficients Poly */
#define CR2RES_COL_SLIT_CURV_B      "SlitPolyB"     /* 3 Coefficients Poly */
#define CR2RES_COL_SLIT_CURV_C      "SlitPolyC"     /* 3 Coefficients Poly */
#define CR2RES_COL_SLIT_FRACTION    "SlitFraction"  /* Position on the Slit */

/*** LINES_DIAGNOSTICS Table ***/
#define CR2RES_COL_MEASURED_LAMBDA  "Measured WL"   /* In Nanometers */
#define CR2RES_COL_CATALOG_LAMBDA   "Catalog WL"    /* In Nanometers */
#define CR2RES_COL_DELTA_LAMBDA     "Delta WL"      /* In Nanometers */
#define CR2RES_COL_MEASURED_PIXEL   "Measured Pix"  /* In Pixels */
#define CR2RES_COL_LINE_WIDTH       "Line Width"    /* In Pixels */
#define CR2RES_COL_FIT_QUALITY      "Fit Quality"   /*  */
#define CR2RES_COL_INTENSITY        "Line Intensity"/*  */

/*** Clusters Table ***/
#define CR2RES_COL_XS               "xs"            /* pixel position */
#define CR2RES_COL_YS               "ys"            /* pixel position */
#define CR2RES_COL_CLUSTERS         "clusters"      /* cluster label */

/* SLIT_FUNC Table */
#define CR2RES_COL_SLIT_FUNC_SUFFIX "SLIT_FUNC"     /* Intensity */

/* POL_SPEC Table */
#define CR2RES_COL_POL_STOKES_SUFFIX        "STOKES"        /* TODO */
#define CR2RES_COL_POL_STOKES_ERROR_SUFFIX  "STOKES_ERR"    /* Error */
#define CR2RES_COL_POL_NULL_SUFFIX          "NULL"          /* TODO */
#define CR2RES_COL_POL_NULL_ERROR_SUFFIX    "NULL_ERR"      /* Error */
#define CR2RES_COL_POL_INTENS_SUFFIX        "INTENS"        /* TODO */
#define CR2RES_COL_POL_INTENS_ERROR_SUFFIX  "INTENS_ERR"    /* Error */

/* EXTRACT_1D Table */
#define CR2RES_COL_SPEC_SUFFIX      "SPEC"          /* Intensity */
#define CR2RES_COL_ERROR_SUFFIX     "ERR"           /* Error */
#define CR2RES_COL_WL_SUFFIX        "WL"            /* Wavelength */

/* SLIT_CURV Table */
#define CR2RES_COL_SLIT_CURV_SUFFIX "SLIT_CURV"     /* Polynomial */

/* SPLICED_1D Table */
#define CR2RES_COL_SPLICED_1D_SPEC  "SPLICED_1D_SPEC"  /* Intensity */
#define CR2RES_COL_SPLICED_1D_ERROR "SPLICED_1D_ERR"   /* Error */
#define CR2RES_COL_SPLICED_1D_WL    "SPLICED_1D_WL"    /* Wavelength */

/* EXTRACT_2D Table */
#define CR2RES_COL_POSITIONX_SUFFIX     "POSITIONX"     /* pixel position */
#define CR2RES_COL_POSITIONY_SUFFIX     "POSITIONY"     /* pixel position */
#define CR2RES_COL_SLIT_FRACTION_SUFFIX "SLIT_FRACTION" /* slit fraction */

/******************************************/

/*************************************/
/* Define here the PRO.TYPE keywords */

/* For  CR2RES_EMISSION_LINES_PROCATG */
/* Table with columns CR2RES_COL_EMISSION / CR2RES_COL_WAVELENGTH */
#define CR2RES_PROTYPE_CATALOG              "CATALOG"

/* For CR2RES_CAL_DETLIN_COEFFS_PROCATG */
#define CR2RES_DETLIN_COEFFS_PROTYPE        "DETLIN_COEFFS"

/* For  CR2RES_CAL_FLAT_BPM_PROCATG */
/*      CR2RES_CAL_DETLIN_BPM_PROCATG */
/*      CR2RES_CAL_DARK_BPM_PROCATG */
/*      CR2RES_UTIL_BPM_SPLIT_PROCATG */
/*      CR2RES_UTIL_NORM_BPM_PROCATG */
/* BPM Image with values CR2RES_BPM_DARK / CR2RES_BPM_FLAT... */
#define CR2RES_BPM_PROTYPE                  "BPM"

/* For CR2RES_CAL_DARK_MASTER_PROCATG */
#define CR2RES_MASTER_DARK_PROTYPE          "MASTER_DARK"

/* For  CR2RES_CAL_FLAT_EXTRACT_1D_PROCATG */
/*      CR2RES_UTIL_EXTRACT_1D_PROCATG */
/*      CR2RES_OBS_NODDING_EXTRACTA_PROCATG */
/*      CR2RES_OBS_NODDING_EXTRACTB_PROCATG */
/* Table with columns cr2res_dfs_SPEC_colname() */
/*                    cr2res_dfs_WAVELENGTH_colname() */
/*                    cr2res_dfs_SPEC_ERR_colname() */
#define CR2RES_EXTRACT_1D_PROTYPE           "EXTRACT_1D"

/* For CR2RES_OBS_POL_SPECA_PROCATG */
/*     CR2RES_OBS_POL_SPECB_PROCATG */
/* Table with columns cr2res_dfs_POL_WAVELENGTH_colname() */
/*                    cr2res_dfs_POL_STOKES_colname() */
/*                    cr2res_dfs_POL_STOKES_ERROR_colname() */
/*                    cr2res_dfs_POL_NULL_colname() */
/*                    cr2res_dfs_POL_NULL_ERROR_colname() */
/*                    cr2res_dfs_POL_INTENS_colname() */
/*                    cr2res_dfs_POL_INTENS_ERROR_colname() */
#define CR2RES_POL_SPEC_PROTYPE             "POL_SPEC"

/* For  CR2RES_OBS_2D_EXTRACT_PROCATG */
/* Table with columns cr2res_dfs_SPEC_colname() */
/*                    cr2res_dfs_SPEC_ERR_colname() */
/*                    cr2res_dfs_WAVELENGTH_colname() */
/*                    cr2res_dfs_POSITIONX_colname() */
/*                    cr2res_dfs_POSITIONY_colname() */
/*                    cr2res_dfs_SLIT_FRACTION_colname() */
#define CR2RES_EXTRACT_2D_PROTYPE           "EXTRACT_2D"

/* For  CR2RES_OBS_NODDING_COMBINEDA_PROCATG */
/*      CR2RES_OBS_NODDING_COMBINEDB_PROCATG */
#define CR2RES_COMBINED_PROTYPE             "COMBINED"

/* For CR2RES_UTIL_SPLICE_SPLICED_1D_PROCATG */
#define CR2RES_SPLICED_1D_PROTYPE           "SPLICED_1D"

/* For  CR2RES_CAL_FLAT_SLIT_MODEL_PROCATG */
/*      CR2RES_UTIL_SLIT_MODEL_PROCATG */
/*      CR2RES_OBS_NODDING_SLITMODELA_PROCATG */
/*      CR2RES_OBS_NODDING_SLITMODELB_PROCATG */
#define CR2RES_SLIT_MODEL_PROTYPE           "SLIT_MODEL"

/* For  CR2RES_CAL_FLAT_SLIT_FUNC_PROCATG */
/*      CR2RES_UTIL_SLIT_FUNC_PROCATG */
/*      CR2RES_OBS_NODDING_SLITFUNCA_PROCATG */
/*      CR2RES_OBS_NODDING_SLITFUNCB_PROCATG */
/* Table with columns from cr2res_dfs_SLIT_FUNC_colname() */
#define CR2RES_SLIT_FUNC_PROTYPE            "SLIT_FUNC"

/* For  CR2RES_CAL_FLAT_TW_PROCATG */
/*      CR2RES_CAL_FLAT_TW_MERGED_PROCATG */
/*      CR2RES_UTIL_TRACE_TW_PROCATG */
/*      CR2RES_UTIL_WAVE_TW_PROCATG */
/*      CR2RES_CAL_WAVE_TW_PROCATG */
/*      CR2RES_UTIL_SLIT_CURV_TW_PROCATG */
/* Table with Traces polynomials, orders, trace_nb wavelengths */
/*                  1 trace per Row */
#define CR2RES_TW_PROTYPE                   "TW"

/* For CR2RES_CAL_FLAT_MASTER_PROCATG */
/*     CR2RES_UTIL_MASTER_FLAT_PROCATG */
/* Master Flat image with values around 1 */
#define CR2RES_MASTER_FLAT_PROTYPE          "MASTER_FLAT"

/* For CR2RES_UTIL_CALIB_PROCATG */
#define CR2RES_CALIBRATED_PROTYPE           "CALIBRATED"

/* For  CR2RES_UTIL_SLIT_CURV_PROCATG  */
/* Table with columns from cr2res_dfs_SLIT_CURV_colname() */
#define CR2RES_SLIT_CURV_PROTYPE            "SLIT_CURV"

/* For  CR2RES_CALIB_COLLAPSED_PROCATG */
#define CR2RES_CALIB_COLLAPSED_PROTYPE      "CALIB_COLLAPSED"

/* For CR2RES_UTIL_WAVE_MAP_PROCATG */
/*     CR2RES_CAL_WAVE_MAP_PROCATG */
/*     CR2RES_UTIL_TRACE_MAP_WL_PROCATG */
#define CR2RES_WAVE_MAP_PROTYPE             "WAVE_MAP"
  
/* For CR2RES_UTIL_SLIT_CURV_MAP_PROCATG */
/*     CR2RES_UTIL_TRACE_MAP_SLIT_CURVE_PROCATG */
#define CR2RES_SLIT_CURV_MAP_PROTYPE        "SLIT_CURV_MAP"

/* For CR2RES_UTIL_TRACE_MAP_TRACE_PROCATG */
#define CR2RES_TRACE_MAP_PROTYPE            "TRACE_MAP"

/* For CR2RES_UTIL_WAVE_LINES_DIAGNOSTICS_PROCATG */
/*     CR2RES_CAL_WAVE_LINES_DIAGNOSTICS_PROCATG */
#define CR2RES_LINES_DIAGNOSTICS_PROTYPE    "LINES_DIAGNOSTICS"

/* For CR2RES_UTIL_WAVE_XCORR_PROCATG */
/* Col: IRPLIB_WLXCORR_COL_WAVELENGTH   */
/*      IRPLIB_WLXCORR_COL_CAT_INIT     */
/*      IRPLIB_WLXCORR_COL_CAT_FINAL    */
/*      IRPLIB_WLXCORR_COL_OBS          */
#define CR2RES_PROTYPE_XCORR        		"XCORR"

/*************************************/
/* Define here the PRO.CATG keywords */
/*************************************/
/* Produced by cr2res_cal_detlin */
#define CR2RES_CAL_DETLIN_COEFFS_PROCATG    "CAL_DETLIN_COEFFS"
#define CR2RES_CAL_DETLIN_BPM_PROCATG       "CAL_DETLIN_BPM"

/* Produced by cr2res_cal_dark */
#define CR2RES_CAL_DARK_MASTER_PROCATG      "CAL_DARK_MASTER"
#define CR2RES_CAL_DARK_BPM_PROCATG         "CAL_DARK_BPM"

/* Produced by cr2res_cal_flat */
#define CR2RES_CAL_FLAT_BPM_PROCATG         "CAL_FLAT_BPM"
#define CR2RES_CAL_FLAT_EXTRACT_1D_PROCATG  "CAL_FLAT_EXTRACT_1D"
#define CR2RES_CAL_FLAT_SLIT_MODEL_PROCATG  "CAL_FLAT_SLIT_MODEL"
#define CR2RES_CAL_FLAT_SLIT_FUNC_PROCATG   "CAL_FLAT_SLIT_FUNC"
#define CR2RES_CAL_FLAT_TW_PROCATG          "CAL_FLAT_TW"
#define CR2RES_CAL_FLAT_TW_MERGED_PROCATG   "CAL_FLAT_TW_MERGED"
#define CR2RES_CAL_FLAT_MASTER_PROCATG      "CAL_FLAT_MASTER"

/* Produced by cr2res_cal_wave */
#define CR2RES_CAL_WAVE_TW_PROCATG          "CAL_WAVE_TW"
#define CR2RES_CAL_WAVE_MAP_PROCATG         "CAL_WAVE_MAP"
#define CR2RES_CAL_WAVE_LINES_DIAGNOSTICS_PROCATG "CAL_WAVE_LINES_DIAGNOSTICS"

/* Produced by cr2res_obs_nodding */
#define CR2RES_OBS_NODDING_EXTRACTA_PROCATG     "OBS_NODDING_EXTRACTA"
#define CR2RES_OBS_NODDING_COMBINEDA_PROCATG    "OBS_NODDING_COMBINEDA"
#define CR2RES_OBS_NODDING_SLITFUNCA_PROCATG    "OBS_NODDING_SLITFUNCA"
#define CR2RES_OBS_NODDING_SLITMODELA_PROCATG   "OBS_NODDING_SLITMODELA"
#define CR2RES_OBS_NODDING_EXTRACTB_PROCATG     "OBS_NODDING_EXTRACTB"
#define CR2RES_OBS_NODDING_COMBINEDB_PROCATG    "OBS_NODDING_COMBINEDB"
#define CR2RES_OBS_NODDING_SLITFUNCB_PROCATG    "OBS_NODDING_SLITFUNCB"
#define CR2RES_OBS_NODDING_SLITMODELB_PROCATG   "OBS_NODDING_SLITMODELB"

/* Produced by cr2res_obs_2d */
#define CR2RES_OBS_2D_EXTRACT_PROCATG       "OBS_2D_EXTRACT"

/* Produced by cr2res_obs_pol */
#define CR2RES_OBS_POL_SPECA_PROCATG        "OBS_POL_SPECA"
#define CR2RES_OBS_POL_SPECB_PROCATG        "OBS_POL_SPECB"

/* Produced by cr2res_util_genlines */
#define CR2RES_EMISSION_LINES_PROCATG       "EMISSION_LINES"

/* Produced by cr2res_util_calib */
#define CR2RES_UTIL_CALIB_PROCATG           "UTIL_CALIB"

/* Produced by cr2res_util_bpm_split */
#define CR2RES_UTIL_BPM_SPLIT_PROCATG       "UTIL_BPM_SPLIT"

/* Produced by cr2res_util_trace */
#define CR2RES_UTIL_TRACE_TW_PROCATG        "UTIL_TRACE_TW"

/* Produced by cr2res_util_extract */
#define CR2RES_UTIL_SLIT_FUNC_PROCATG       "UTIL_SLIT_FUNC"
#define CR2RES_UTIL_SLIT_MODEL_PROCATG      "UTIL_SLIT_MODEL"
#define CR2RES_UTIL_EXTRACT_1D_PROCATG      "UTIL_EXTRACT_1D"

/* Produced by cr2res_util_normflat */
#define CR2RES_UTIL_MASTER_FLAT_PROCATG     "UTIL_MASTER_FLAT"
#define CR2RES_UTIL_NORM_BPM_PROCATG        "UTIL_NORM_BPM"

/* Produced by cr2res_util_slit_curv */
#define CR2RES_UTIL_SLIT_CURV_MAP_PROCATG   "UTIL_SLIT_CURV_MAP"
#define CR2RES_UTIL_SLIT_CURV_PROCATG       "UTIL_SLIT_CURV"
#define CR2RES_UTIL_SLIT_CURV_TW_PROCATG    "UTIL_SLIT_CURV_TW"

/* Produced by cr2res_util_wave */
#define CR2RES_UTIL_WAVE_TW_PROCATG         "UTIL_WAVE_TW"
#define CR2RES_UTIL_WAVE_MAP_PROCATG        "UTIL_WAVE_MAP"
#define CR2RES_UTIL_WAVE_XCORR_PROCATG      "UTIL_WAVE_XCORR"
#define CR2RES_UTIL_WAVE_LINES_DIAGNOSTICS_PROCATG "UTIL_WAVE_LINES_DIAGNOSTICS"

/* Produced by cr2res_util_trace_maps */
#define CR2RES_UTIL_TRACE_MAP_SLIT_CURVE_PROCATG "UTIL_TRACE_MAP_SLIT_CURVE"
#define CR2RES_UTIL_TRACE_MAP_WL_PROCATG    "UTIL_TRACE_MAP_WL"
#define CR2RES_UTIL_TRACE_MAP_TRACE_PROCATG "UTIL_TRACE_MAP_TRACE"

/* Produced by cr2res_util_splice */
#define CR2RES_UTIL_SPLICE_SPLICED_1D_PROCATG   "UTIL_SPLICE_SPLICED_1D"

/* Define here the DO.CATG keywords */
#define CR2RES_COMMAND_LINE                 "COMMAND_LINE"
#define CR2RES_EMISSION_LINES_TXT_RAW       "EMISSION_LINES_TXT"
#define CR2RES_DETLIN_RAW                   "DETLIN"
#define CR2RES_DARK_RAW                     "DARK"
#define CR2RES_FLAT_RAW                     "FLAT"
#define CR2RES_WAVE_RAW                     "WAVE"
#define CR2RES_OBS_NODDING_RAW              "OBS_NODDING"
#define CR2RES_OBS_2D_RAW                   "OBS_2D"
#define CR2RES_OBS_POL_RAW                  "OBS_POL"

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

cpl_error_code cr2res_dfs_set_groups(cpl_frameset *);
char * cr2res_dfs_SPEC_colname(int, int) ;
char * cr2res_dfs_WAVELENGTH_colname(int, int);
char * cr2res_dfs_SPEC_ERR_colname(int, int);
char * cr2res_dfs_SLIT_FUNC_colname(int, int) ;
char * cr2res_dfs_SLIT_CURV_colname(int, int) ;
char * cr2res_dfs_POSITIONX_colname(int, int) ;
char * cr2res_dfs_POSITIONY_colname(int, int) ;
char * cr2res_dfs_SLIT_FRACTION_colname(int, int) ;
char * cr2res_dfs_POL_WAVELENGTH_colname(int);
char * cr2res_dfs_POL_STOKES_colname(int) ;
char * cr2res_dfs_POL_STOKES_ERROR_colname(int) ;
char * cr2res_dfs_POL_NULL_colname(int) ;
char * cr2res_dfs_POL_NULL_ERROR_colname(int) ;
char * cr2res_dfs_POL_INTENS_colname(int) ;
char * cr2res_dfs_POL_INTENS_ERROR_colname(int) ;
char * cr2res_dfs_SPEC_colname_parse(
        const char  * colname,
        int         * order,
        int         * trace) ;

cpl_table * cr2res_dfs_create_lines_diagnostics_table(int nrows) ;

int cr2res_dfs_check_traces_table(const cpl_table * traces) ;

#endif
