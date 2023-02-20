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

/*** PHOTO_FLUX Table ***/
#define CR2RES_COL_STDNAME          "Std_Star_Name"     /* No Unit */
#define CR2RES_COL_RA               "Right_Ascension"   /* In degrees */
#define CR2RES_COL_DEC              "Declination"       /* In degrees */
#define CR2RES_COL_PHOTOFLUX        "Photospheric_Flux" /* In Jy */

/*** CATALOG Table ***/
#define CR2RES_COL_EMISSION         "Emission"      /* No Unit */
/*** CATALOG, CR2RES_OBS_NODDING_IDP_PROCATG and TRACE_WAVE Table ***/
#define CR2RES_COL_WAVELENGTH       "Wavelength"    /* In Nanometers */
/*** CATALOG and TRACE_WAVE Table ***/
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

/*** CR2RES_OBS_NODDING_IDP_PROCATG  ***/
#define CR2RES_IDP_COL_FLUX         "FLUX"       /* Intensity */
#define CR2RES_IDP_COL_ERR          "ERR"        /* Error */
#define CR2RES_IDP_COL_WAVE         "WAVE"       /* wavelength */
#define CR2RES_IDP_COL_QUAL         "QUAL"       /* Quality */
#define CR2RES_IDP_COL_ORDER        "ORDER"      /* Order Number */
#define CR2RES_IDP_COL_TRACE        "TRACE"      /* Order Number */
#define CR2RES_IDP_COL_DETEC        "DETEC"      /* Detector */
#define CR2RES_IDP_COL_XPOS         "XPOS"       /* Original X position */
#define CR2RES_IDP_COL_YPOS         "YPOS"       /* Original Y position */
#define CR2RES_IDP_COL_SLITFRAC     "SLITFRAC"   /* Original Y position */

/*** LINES_DIAGNOSTICS Table ***/
#define CR2RES_COL_MEASURED_LAMBDA  "Measured_WL"   /* In Nanometers */
#define CR2RES_COL_CATALOG_LAMBDA   "Catalog_WL"    /* In Nanometers */
#define CR2RES_COL_DELTA_LAMBDA     "Delta_WL"      /* In Nanometers */
#define CR2RES_COL_MEASURED_PIXEL   "Measured_Pix"  /* In Pixels */
#define CR2RES_COL_LINE_WIDTH       "Line_Width"    /* In Pixels */
#define CR2RES_COL_FIT_QUALITY      "Fit_Quality"   /*  */
#define CR2RES_COL_INTENSITY        "Line_Intensity"/*  */
#define CR2RES_COL_FPET_M           "FPET_Order"

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

/* THROUGHPUT Table */
#define CR2RES_COL_CONVERSION_SUFFIX    "CONVERSION"    /* Intensity */
#define CR2RES_COL_SENSITIVITY_SUFFIX   "SENSITIVITY"   /* Intensity */
#define CR2RES_COL_THROUGHPUT_SUFFIX    "THROUGHPUT"    /* Intensity */

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
/* Define here the DRS.TYPE keywords */

/* For  CR2RES_EMISSION_LINES_PROCATG */
/* Table with columns CR2RES_COL_EMISSION / CR2RES_COL_WAVELENGTH */
#define CR2RES_DRSTYPE_CATALOG              "CATALOG"

/* For  CR2RES_PHOTO_FLUX_PROCATG */
/* Col: CR2RES_COL_STDNAME      */
/*      CR2RES_COL_RA           */
/*      CR2RES_COL_DEC          */
/*      CR2RES_COL_PHOTOFLUX    */
#define CR2RES_PHOTO_FLUX_DRSTYPE           "PHOTO_FLUX"

/* For CR2RES_CAL_DETLIN_COEFFS_PROCATG */
#define CR2RES_DETLIN_COEFFS_DRSTYPE        "DETLIN_COEFFS"

/* For  CR2RES_CAL_FLAT_BPM_PROCATG */
/*      CR2RES_CAL_DETLIN_BPM_PROCATG */
/*      CR2RES_CAL_DARK_BPM_PROCATG */
/*      CR2RES_UTIL_BPM_MERGE_PROCATG */
/*      CR2RES_UTIL_BPM_SPLIT_PROCATG */
/*      CR2RES_UTIL_NORM_BPM_PROCATG */
/* BPM Image with values CR2RES_BPM_DARK / CR2RES_BPM_FLAT... */
#define CR2RES_BPM_DRSTYPE                  "BPM"

/* For CR2RES_CAL_DARK_MASTER_PROCATG */
#define CR2RES_MASTER_DARK_DRSTYPE          "MASTER_DARK"

/* For  CR2RES_CAL_FLAT_EXTRACT_1D_PROCATG */
/*      CR2RES_UTIL_EXTRACT_1D_PROCATG */
/*      CR2RES_UTIL_WAVE_EXTRACT_1D_PROCATG */
/*      CR2RES_CAL_WAVE_EXTRACT_1D_PROCATG */
/*      CR2RES_OBS_NODDING_EXTRACTA_PROCATG */
/*      CR2RES_OBS_NODDING_EXTRACTB_PROCATG */
/*      CR2RES_OBS_NODDING_EXTRACTC_PROCATG */
/*         CR2RES_OBS_STARING_EXTRACT_PROCATG */
/*      CR2RES_OBS_POL_EXTRACTA_PROCATG */
/*      CR2RES_OBS_POL_EXTRACTB_PROCATG */
/* Table with columns cr2res_dfs_SPEC_colname() */
/*                    cr2res_dfs_WAVELENGTH_colname() */
/*                    cr2res_dfs_SPEC_ERR_colname() */
#define CR2RES_EXTRACT_1D_DRSTYPE           "EXTRACT_1D"

/* For  CR2RES_OBS_NODDING_IDP_PROCATG */
/* Table with columns CR2RES_COL_WAVELENGTH */
/*                    CR2RES_COL_SPECTRUM */
/*                    CR2RES_COL_ERROR */
/*                    CR2RES_COL_QUALITY */
/*                    CR2RES_COL_ORDER */
/*                    CR2RES_COL_DETECTOR */
#define CR2RES_EXTRACT_1D_IDP_DRSTYPE       "EXTRACT_1D_IDP"

/* For  CR2RES_OBS_NODDING_THROUGHPUT_PROCATG */
/* Table with columns cr2res_dfs_WAVELENGTH_colname() */
/*                    cr2res_dfs_CONVERSION_colname() */
/*                    cr2res_dfs_THROUGHPUT_colname() */
/*                    cr2res_dfs_SENSITIVITY_colname() */
#define CR2RES_THROUGHPUT_DRSTYPE           "THROUGHPUT"

/* For CR2RES_OBS_POL_SPECA_PROCATG */
/*     CR2RES_OBS_POL_SPECB_PROCATG */
/* Table with columns cr2res_dfs_POL_WAVELENGTH_colname() */
/*                    cr2res_dfs_POL_STOKES_colname() */
/*                    cr2res_dfs_POL_STOKES_ERROR_colname() */
/*                    cr2res_dfs_POL_NULL_colname() */
/*                    cr2res_dfs_POL_NULL_ERROR_colname() */
/*                    cr2res_dfs_POL_INTENS_colname() */
/*                    cr2res_dfs_POL_INTENS_ERROR_colname() */
#define CR2RES_POL_SPEC_DRSTYPE             "POL_SPEC"

/* For  CR2RES_OBS_2D_EXTRACT_PROCATG */
/* Table with columns cr2res_dfs_SPEC_colname() */
/*                    cr2res_dfs_SPEC_ERR_colname() */
/*                    cr2res_dfs_WAVELENGTH_colname() */
/*                    cr2res_dfs_POSITIONX_colname() */
/*                    cr2res_dfs_POSITIONY_colname() */
/*                    cr2res_dfs_SLIT_FRACTION_colname() */
#define CR2RES_EXTRACT_2D_DRSTYPE           "EXTRACT_2D"

/* For  CR2RES_OBS_NODDING_COMBINEDA_PROCATG */
/*      CR2RES_OBS_NODDING_COMBINEDB_PROCATG */
#define CR2RES_COMBINED_DRSTYPE             "COMBINED"

/* For CR2RES_UTIL_SPLICE_SPLICED_1D_PROCATG */
#define CR2RES_SPLICED_1D_DRSTYPE           "SPLICED_1D"

/* For  CR2RES_CAL_FLAT_SLIT_MODEL_PROCATG */
/*      CR2RES_UTIL_SLIT_MODEL_PROCATG */
/*      CR2RES_OBS_NODDING_SLITMODELA_PROCATG */
/*      CR2RES_OBS_NODDING_SLITMODELB_PROCATG */
/*         CR2RES_OBS_STARING_SLITMODEL_PROCATG */
#define CR2RES_SLIT_MODEL_DRSTYPE           "SLIT_MODEL"

/* For  CR2RES_CAL_FLAT_SLIT_FUNC_PROCATG */
/*      CR2RES_UTIL_SLIT_FUNC_PROCATG */
/*      CR2RES_OBS_NODDING_SLITFUNCA_PROCATG */
/*      CR2RES_OBS_NODDING_SLITFUNCB_PROCATG */
/*         CR2RES_OBS_STARING_SLITFUNC_PROCATG */
/* Table with columns from cr2res_dfs_SLIT_FUNC_colname() */
#define CR2RES_SLIT_FUNC_DRSTYPE            "SLIT_FUNC"

/* For  CR2RES_CAL_FLAT_TW_PROCATG */
/*      CR2RES_CAL_FLAT_TW_MERGED_PROCATG */
/*      CR2RES_UTIL_TRACE_TW_PROCATG */
/*      CR2RES_UTIL_WAVE_TW_PROCATG */
/*      CR2RES_CAL_WAVE_TW_PROCATG */
/*      CR2RES_UTIL_SLIT_CURV_TW_PROCATG */
/*      CR2RES_OBS_NODDING_TWA_PROCATG */
/*      CR2RES_OBS_NODDING_TWB_PROCATG */
/*      CR2RES_OBS_POL_TWA_PROCATG */
/*      CR2RES_OBS_POL_TWB_PROCATG */
/* Table with Traces polynomials, orders, trace_nb wavelengths */
/*                  1 trace per Row */
#define CR2RES_TW_DRSTYPE                   "TW"

/* For CR2RES_CAL_FLAT_MASTER_PROCATG */
/*     CR2RES_UTIL_MASTER_FLAT_PROCATG */
/* Master Flat image with values around 1 */
#define CR2RES_MASTER_FLAT_DRSTYPE          "MASTER_FLAT"

/* For CR2RES_UTIL_CALIB_PROCATG */
/*     CR2RES_OBS_POL_CALIB_A_PROCATG */
/*     CR2RES_OBS_POL_CALIB_B_PROCATG */
/*     CR2RES_OBS_2D_CALIBRATED_PROCATG */
#define CR2RES_CALIBRATED_DRSTYPE           "CALIBRATED"

/* For  CR2RES_UTIL_SLIT_CURV_PROCATG  */
/* Table with columns from cr2res_dfs_SLIT_CURV_colname() */
#define CR2RES_SLIT_CURV_DRSTYPE            "SLIT_CURV"

/* For  CR2RES_CALIB_COLLAPSED_PROCATG */
#define CR2RES_CALIB_COLLAPSED_DRSTYPE      "CALIB_COLLAPSED"

/* For CR2RES_UTIL_WAVE_MAP_PROCATG */
/*     CR2RES_CAL_WAVE_MAP_PROCATG */
/*     CR2RES_UTIL_TRACE_MAP_WL_PROCATG */
#define CR2RES_WAVE_MAP_DRSTYPE             "WAVE_MAP"
  
/* For CR2RES_UTIL_SLIT_CURV_MAP_PROCATG */
/*     CR2RES_UTIL_TRACE_MAP_SLIT_CURVE_PROCATG */
#define CR2RES_SLIT_CURV_MAP_DRSTYPE        "SLIT_CURV_MAP"

/* For CR2RES_UTIL_TRACE_MAP_TRACE_PROCATG */
#define CR2RES_TRACE_MAP_DRSTYPE            "TRACE_MAP"

/* For CR2RES_UTIL_WAVE_LINES_DIAGNOSTICS_PROCATG */
/*     CR2RES_CAL_WAVE_LINES_DIAGNOSTICS_PROCATG */
#define CR2RES_LINES_DIAGNOSTICS_DRSTYPE    "LINES_DIAGNOSTICS"

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
#define CR2RES_CAL_WAVE_EXTRACT_1D_PROCATG  "CAL_WAVE_EXTRACT_1D"

/* Produced by cr2res_obs_nodding */
#define CR2RES_OBS_NODDING_EXTRACTA_PROCATG     "OBS_NODDING_EXTRACTA"
#define CR2RES_OBS_NODDING_COMBINEDA_PROCATG    "OBS_NODDING_COMBINEDA"
#define CR2RES_OBS_NODDING_SLITFUNCA_PROCATG    "OBS_NODDING_SLITFUNCA"
#define CR2RES_OBS_NODDING_SLITMODELA_PROCATG   "OBS_NODDING_SLITMODELA"
#define CR2RES_OBS_NODDING_EXTRACTB_PROCATG     "OBS_NODDING_EXTRACTB"
#define CR2RES_OBS_NODDING_COMBINEDB_PROCATG    "OBS_NODDING_COMBINEDB"
#define CR2RES_OBS_NODDING_SLITFUNCB_PROCATG    "OBS_NODDING_SLITFUNCB"
#define CR2RES_OBS_NODDING_SLITMODELB_PROCATG   "OBS_NODDING_SLITMODELB"
#define CR2RES_OBS_NODDING_EXTRACTC_PROCATG     "OBS_NODDING_EXTRACT_COMB"
#define CR2RES_OBS_NODDING_THROUGHPUT_PROCATG   "OBS_NODDING_THROUGHPUT"
#define CR2RES_OBS_NODDING_TWA_PROCATG          "OBS_NODDING_TWA"
#define CR2RES_OBS_NODDING_TWB_PROCATG          "OBS_NODDING_TWB"
#define CR2RES_OBS_NODDING_IDP_PROCATG          "OBS_NODDING_IDP"

/* Produced by cr2res_obs_staring */
#define CR2RES_OBS_STARING_EXTRACT_PROCATG         "OBS_STARING_EXTRACT"
#define CR2RES_OBS_STARING_SLITFUNC_PROCATG        "OBS_STARING_SLITFUNC"
#define CR2RES_OBS_STARING_SLITMODEL_PROCATG       "OBS_STARING_SLITMODEL"
#define CR2RES_OBS_STARING_IDP_PROCATG             "OBS_STARING_IDP"


/* Produced by cr2res_obs_2d */
#define CR2RES_OBS_2D_EXTRACT_PROCATG       "OBS_2D_EXTRACT"
#define CR2RES_OBS_2D_CALIBRATED_PROCATG    "OBS_2D_CALIBRATED"
#define CR2RES_OBS_2D_IDP_PROCATG           "OBS_2D_IDP"

/* Produced by cr2res_obs_pol */
#define CR2RES_OBS_POL_EXTRACTA_PROCATG     "OBS_POL_EXTRACTA"
#define CR2RES_OBS_POL_EXTRACTB_PROCATG     "OBS_POL_EXTRACTB"
#define CR2RES_OBS_POL_TWA_PROCATG          "OBS_POL_TWA"
#define CR2RES_OBS_POL_TWB_PROCATG          "OBS_POL_TWB"
#define CR2RES_OBS_POL_SPECA_PROCATG        "OBS_POL_SPECA"
#define CR2RES_OBS_POL_SPECB_PROCATG        "OBS_POL_SPECB"
#define CR2RES_OBS_POL_CALIB_A_PROCATG      "OBS_POL_CALIBA"
#define CR2RES_OBS_POL_CALIB_B_PROCATG      "OBS_POL_CALIBB"
#define CR2RES_OBS_POL_IDP_PROCATG          "OBS_POL_IDP"

/* Produced by cr2res_util_genlines */
#define CR2RES_EMISSION_LINES_PROCATG       "EMISSION_LINES"

/* Produced by cr2res_util_genstd */
#define CR2RES_PHOTO_FLUX_PROCATG           "PHOTO_FLUX"

/* Produced by cr2res_util_calib */
#define CR2RES_UTIL_CALIB_PROCATG           "UTIL_CALIB"

/* Produced by cr2res_util_bpm_merge */
#define CR2RES_UTIL_BPM_MERGE_PROCATG       "UTIL_BPM_MERGE"

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
#define CR2RES_UTIL_WAVE_LINES_DIAGNOSTICS_PROCATG "UTIL_WAVE_LINES_DIAGNOSTICS"
#define CR2RES_UTIL_WAVE_EXTRACT_1D_PROCATG "UTIL_WAVE_EXTRACT_1D"

/* Produced by cr2res_util_trace_maps */
#define CR2RES_UTIL_TRACE_MAP_SLIT_CURVE_PROCATG    "UTIL_TRACE_MAP_SLIT_CURVE"
#define CR2RES_UTIL_TRACE_MAP_WL_PROCATG            "UTIL_TRACE_MAP_WL"
#define CR2RES_UTIL_TRACE_MAP_TRACE_PROCATG         "UTIL_TRACE_MAP_TRACE"

/* Produced by cr2res_util_splice */
#define CR2RES_UTIL_SPLICE_SPLICED_1D_PROCATG   "UTIL_SPLICE_SPLICED_1D"

/* Define here the DO.CATG keywords */
#define CR2RES_EMISSION_LINES_TXT_RAW       "EMISSION_LINES_TXT"
#define CR2RES_LINES_SELECTION_TXT_RAW      "LINES_SELECTION_TXT"
#define CR2RES_PHOTO_FLUX_TXT_RAW           "PHOTO_FLUX_TXT"
#define CR2RES_DETLIN_DARK_RAW              "DETLIN_DARK"
#define CR2RES_DETLIN_LAMP_RAW              "DETLIN_LAMP"
#define CR2RES_DARK_RAW                     "DARK"
#define CR2RES_FLAT_RAW                     "FLAT"
#define CR2RES_WAVE_FPET_RAW                "WAVE_FPET"
#define CR2RES_WAVE_UNE_RAW                 "WAVE_UNE"
#define CR2RES_METROLOGY_RAW                "METROLOGY"
#define CR2RES_CAL_NODDING_OTHER_RAW        "CAL_NODDING_OTHER"
#define CR2RES_CAL_NODDING_JITTER_RAW       "CAL_NODDING_JITTER"
#define CR2RES_OBS_NODDING_OTHER_RAW        "OBS_NODDING_OTHER"
#define CR2RES_OBS_NODDING_JITTER_RAW       "OBS_NODDING_JITTER"
#define CR2RES_OBS_ASTROMETRY_OTHER_RAW     "OBS_ASTROMETRY_OTHER"
#define CR2RES_OBS_ASTROMETRY_JITTER_RAW    "OBS_ASTROMETRY_JITTER"
#define CR2RES_OBS_STARING_OTHER_RAW        "OBS_STARING_OTHER"
#define CR2RES_OBS_STARING_JITTER_RAW       "OBS_STARING_JITTER"
#define CR2RES_OBS_STARING_WAVE_SKY_RAW     "OBS_WAVE_SKY"
#define CR2RES_OBS_POLARIMETRY_OTHER_RAW    "OBS_POLARIMETRY_OTHER"
#define CR2RES_OBS_2D_OBJECT_RAW            "OBS_2D_OBJECT"
#define CR2RES_OBS_2D_SKY_RAW               "OBS_2D_SKY"

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/

cpl_error_code cr2res_dfs_set_groups(cpl_frameset *);
char * cr2res_dfs_SPEC_colname(int, int) ;
char * cr2res_dfs_THROUGHPUT_colname(int, int) ;
char * cr2res_dfs_SENSITIVITY_colname(int, int) ;
char * cr2res_dfs_CONVERSION_colname(int, int) ;
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
        const char  *   colname,
        int         *   order_idx,
        int         *   trace) ;

cpl_table * cr2res_dfs_create_lines_diagnostics_table(int nrows) ;

int cr2res_dfs_check_traces_table(const cpl_table * traces) ;

#endif
