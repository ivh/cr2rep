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

#ifndef CR2RES_PFITS_H
#define CR2RES_PFITS_H

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>
#include "cr2res_nodding.h"
#include "cr2res_utils.h"

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

/* RAW files header Keywords Names */
#define CR2RES_HEADER_INSTRUMENT        "INSTRUME"
#define CR2RES_HEADER_MJD_OBS           "MJD-OBS"
#define CR2RES_HEADER_EXTNAME           "EXTNAME"
#define CR2RES_HEADER_NAXIS1            "NAXIS1"
#define CR2RES_HEADER_NAXIS2            "NAXIS2"
#define CR2RES_HEADER_EXPNO             "ESO TPL EXPNO"
#define CR2RES_HEADER_NEXP              "ESO TPL NEXP"
#define CR2RES_HEADER_DECKER_POS        "ESO INS OPTI8 NO"
#define CR2RES_HEADER_NODPOS            "ESO SEQ NODPOS"
#define CR2RES_HEADER_DROT_POSANG       "ESO INS1 DROT POSANG"
#define CR2RES_HEADER_NODTHROW          "ESO SEQ NODTHROW"
#define CR2RES_HEADER_ARCFILE           "ARCFILE"
#define CR2RES_HEADER_MJDOBS            "MJD-OBS"
#define CR2RES_HEADER_WLEN_ID           "ESO INS WLEN ID"
#define CR2RES_HEADER_WLEN_BEGIN        "ESO INS WLEN BEGIN%d"
#define CR2RES_HEADER_WLEN_END          "ESO INS WLEN END%d"
#define CR2RES_HEADER_WLEN_CENY         "ESO INS WLEN CENY%d"
#define CR2RES_HEADER_WLEN_CWLEN        "ESO INS WLEN CWLEN"
#define CR2RES_HEADER_GRAT1_ZPORD       "ESO INS GRAT1 ZP_ORD"
#define CR2RES_HEADER_GRAT1_ORDER       "ESO INS GRAT1 ORDER"
#define CR2RES_HEADER_NDIT              "ESO DET NDIT"
#define CR2RES_HEADER_DIT               "ESO DET SEQ1 DIT"
#define CR2RES_HEADER_PROG_ID           "ESO OBS PROG ID"
#define CR2RES_HEADER_OBS_ID            "ESO OBS ID"
#define CR2RES_HEADER_DRS_TYPE          "ESO DRS TYPE"
#define CR2RES_HEADER_DRS_TMID          "ESO DRS TMID"

/* QC Parameter Names */
#define CR2RES_HEADER_QC_DARK_RON1          "ESO QC DARK RON1"
#define CR2RES_HEADER_QC_DARK_RON2          "ESO QC DARK RON2"
#define CR2RES_HEADER_QC_DARK_MEAN          "ESO QC DARK MEAN"
#define CR2RES_HEADER_QC_DARK_MEDIAN        "ESO QC DARK MEDIAN"
#define CR2RES_HEADER_QC_DARK_STDEV         "ESO QC DARK STDDEV"
#define CR2RES_HEADER_QC_DARK_NBAD          "ESO QC DARK NBAD"

#define CR2RES_HEADER_QC_DARK_RON1_AVG      "ESO QC DARK RON1 AVG"
#define CR2RES_HEADER_QC_DARK_RON1_RMS      "ESO QC DARK RON1 RMS"
#define CR2RES_HEADER_QC_DARK_RON2_AVG      "ESO QC DARK RON2 AVG"
#define CR2RES_HEADER_QC_DARK_RON2_RMS      "ESO QC DARK RON2 RMS"
#define CR2RES_HEADER_QC_DARK_MEAN_AVG      "ESO QC DARK MEAN AVG"
#define CR2RES_HEADER_QC_DARK_MEAN_RMS      "ESO QC DARK MEAN RMS"
#define CR2RES_HEADER_QC_DARK_MEDIAN_AVG    "ESO QC DARK MEDIAN AVG"
#define CR2RES_HEADER_QC_DARK_MEDIAN_RMS    "ESO QC DARK MEDIAN RMS"
#define CR2RES_HEADER_QC_DARK_STDEV_AVG     "ESO QC DARK STDEV AVG"
#define CR2RES_HEADER_QC_DARK_STDEV_RMS     "ESO QC DARK STDEV RMS"
#define CR2RES_HEADER_QC_DARK_NBAD_AVG      "ESO QC DARK NBAD AVG"
#define CR2RES_HEADER_QC_DARK_NBAD_RMS      "ESO QC DARK NBAD RMS"

#define CR2RES_HEADER_QC_DETLIN_NBAD        "ESO QC DETLIN NBBAD"
#define CR2RES_HEADER_QC_DETLIN_NBFAILED    "ESO QC DETLIN NBFAILED"
#define CR2RES_HEADER_QC_DETLIN_NBSUCCESS   "ESO QC DETLIN NBSUCCESS"
#define CR2RES_HEADER_QC_DETLIN_MEDIAN      "ESO QC DETLIN MEDIAN"
#define CR2RES_HEADER_QC_DETLIN_MINLEVEL    "ESO QC DETLIN MINLEVEL"
#define CR2RES_HEADER_QC_DETLIN_MAXLEVEL    "ESO QC DETLIN MAXLEVEL"
#define CR2RES_HEADER_QC_DETLIN_MEDA        "ESO QC DETLIN MEDA"
#define CR2RES_HEADER_QC_DETLIN_MEDB        "ESO QC DETLIN MEDB"
#define CR2RES_HEADER_QC_DETLIN_MEDC        "ESO QC DETLIN MEDC"
#define CR2RES_HEADER_QC_DETLIN_MEDQ        "ESO QC DETLIN MEDQ"

#define CR2RES_HEADER_QC_DETLIN_NBAD_AVG        "ESO QC DETLIN NBBAD AVG"
#define CR2RES_HEADER_QC_DETLIN_NBAD_RMS        "ESO QC DETLIN NBBAD RMS"
#define CR2RES_HEADER_QC_DETLIN_NBFAILED_AVG    "ESO QC DETLIN NBFAILED AVG"
#define CR2RES_HEADER_QC_DETLIN_NBFAILED_RMS    "ESO QC DETLIN NBFAILED RMS"
#define CR2RES_HEADER_QC_DETLIN_NBSUCCESS_AVG   "ESO QC DETLIN NBSUCCESS AVG"
#define CR2RES_HEADER_QC_DETLIN_NBSUCCESS_RMS   "ESO QC DETLIN NBSUCCESS RMS"
#define CR2RES_HEADER_QC_DETLIN_MEDIAN_AVG      "ESO QC DETLIN MEDIAN AVG"
#define CR2RES_HEADER_QC_DETLIN_MEDIAN_RMS      "ESO QC DETLIN MEDIAN RMS"
#define CR2RES_HEADER_QC_DETLIN_MINLEVEL_AVG    "ESO QC DETLIN MINLEVEL AVG"
#define CR2RES_HEADER_QC_DETLIN_MINLEVEL_RMS    "ESO QC DETLIN MINLEVEL RMS"
#define CR2RES_HEADER_QC_DETLIN_MAXLEVEL_AVG    "ESO QC DETLIN MAXLEVEL AVG"
#define CR2RES_HEADER_QC_DETLIN_MAXLEVEL_RMS    "ESO QC DETLIN MAXLEVEL RMS"
#define CR2RES_HEADER_QC_DETLIN_MEDA_AVG    "ESO QC DETLIN MEDA AVG"
#define CR2RES_HEADER_QC_DETLIN_MEDA_RMS    "ESO QC DETLIN MEDA RMS"
#define CR2RES_HEADER_QC_DETLIN_MEDB_AVG    "ESO QC DETLIN MEDB AVG"
#define CR2RES_HEADER_QC_DETLIN_MEDB_RMS    "ESO QC DETLIN MEDB RMS"
#define CR2RES_HEADER_QC_DETLIN_MEDC_AVG    "ESO QC DETLIN MEDC AVG"
#define CR2RES_HEADER_QC_DETLIN_MEDC_RMS    "ESO QC DETLIN MEDC RMS"
#define CR2RES_HEADER_QC_DETLIN_MEDQ_AVG    "ESO QC DETLIN MEDQ AVG"
#define CR2RES_HEADER_QC_DETLIN_MEDQ_RMS    "ESO QC DETLIN MEDQ RMS"

#define CR2RES_HEADER_QC_FLAT_MEAN          "ESO QC FLAT MEAN"
#define CR2RES_HEADER_QC_FLAT_MEDIAN        "ESO QC FLAT MEDIAN"
#define CR2RES_HEADER_QC_FLAT_FLUX          "ESO QC FLAT EFFICIENCY"
#define CR2RES_HEADER_QC_FLAT_RMS           "ESO QC FLAT RMS"
#define CR2RES_HEADER_QC_FLAT_S2N           "ESO QC FLAT S2N"
#define CR2RES_HEADER_QC_FLAT_NBBAD         "ESO QC FLAT NBBAD"
#define CR2RES_HEADER_QC_FLAT_ORDERPOS      "ESO QC FLAT ORDERPOS"
#define CR2RES_HEADER_QC_FLAT_CENTERY       "ESO QC FLAT TRACE CENTERY"

#define CR2RES_HEADER_QC_FLAT_MEAN_AVG      "ESO QC FLAT MEAN AVG"
#define CR2RES_HEADER_QC_FLAT_MEAN_RMS      "ESO QC FLAT MEAN RMS"
#define CR2RES_HEADER_QC_FLAT_MEDIAN_AVG    "ESO QC FLAT MEDIAN AVG"
#define CR2RES_HEADER_QC_FLAT_MEDIAN_RMS    "ESO QC FLAT MEDIAN RMS"
#define CR2RES_HEADER_QC_FLAT_FLUX_AVG      "ESO QC FLAT EFFICIENCY AVG"
#define CR2RES_HEADER_QC_FLAT_FLUX_RMS      "ESO QC FLAT EFFICIENCY RMS"
#define CR2RES_HEADER_QC_FLAT_RMS_AVG       "ESO QC FLAT RMS AVG"
#define CR2RES_HEADER_QC_FLAT_RMS_RMS       "ESO QC FLAT RMS RMS"
#define CR2RES_HEADER_QC_FLAT_S2N_AVG       "ESO QC FLAT S2N AVG"
#define CR2RES_HEADER_QC_FLAT_S2N_RMS       "ESO QC FLAT S2N RMS"
#define CR2RES_HEADER_QC_FLAT_NBBAD_AVG     "ESO QC FLAT NBBAD AVG"
#define CR2RES_HEADER_QC_FLAT_NBBAD_RMS     "ESO QC FLAT NBBAD RMS"
#define CR2RES_HEADER_QC_FLAT_CENTERY_AVG   "ESO QC FLAT CENTERY AVG"
#define CR2RES_HEADER_QC_FLAT_CENTERY_RMS   "ESO QC FLAT CENTERY RMS"

#define CR2RES_HEADER_QC_WAVE_BESTXCORR     "ESO QC WAVE BESTXCORR"
#define CR2RES_HEADER_QC_WAVE_CENTWL        "ESO QC WAVE CENTWL"
#define CR2RES_HEADER_QC_WAVE_DISPWL        "ESO QC WAVE DISPWL"
#define CR2RES_HEADER_QC_WAVE_RESOL_FWHM    "ESO QC WAVE RESOLFWHM"
#define CR2RES_HEADER_QC_WAVE_POS           "ESO QC WAVE RESOLFWHM POS"
#define CR2RES_HEADER_QC_WAVE_RESOL         "ESO QC WAVE RESOL"
#define CR2RES_HEADER_QC_WAVE_LAMP_EFFIC    "ESO QC WAVE LAMP EFFIC"

#define CR2RES_HEADER_QC_OVEREXPOSED        "ESO QC OVEREXPOSED"
#define CR2RES_HEADER_QC_SIGNAL             "ESO QC SIGNAL"
#define CR2RES_HEADER_QC_STANDARD_FLUX      "ESO QC STANDARD FLUX"
#define CR2RES_HEADER_QC_SNR                "ESO QC SNR%d"
#define CR2RES_HEADER_QC_SLITFWHM_ORDER     "ESO QC SLITFWHM%d"
#define CR2RES_HEADER_QC_SLITFWHM_MED       "ESO QC SLITFWHM MED"
#define CR2RES_HEADER_QC_REAL_ORDER         "ESO QC REALORDER%d"
#define CR2RES_HEADER_QC_THROUGHPUT         "ESO QC THROUGHPUT"
#define CR2RES_HEADER_QC_TILT               "ESO QC TILT%d"
#define CR2RES_HEADER_QC_TILT_GLOBAL        "ESO QC TILT GLOBAL"
#define CR2RES_HEADER_QC_FPI_CONTRAST       "ESO QC FPI CONTRAST"
#define CR2RES_HEADER_QC_FPI_SEPARATION     "ESO QC FPI SEPARATION"
#define CR2RES_HEADER_QC_UNE_FLUX           "ESO QC UNE FLUX"

/*-----------------------------------------------------------------------------
                                   Functions prototypes
 -----------------------------------------------------------------------------*/

cr2res_nodding_pos cr2res_pfits_get_nodding_pos(const cpl_propertylist * plist);

const char * cr2res_pfits_get_procatg(const cpl_propertylist *) ;
const char * cr2res_pfits_get_drstype(const cpl_propertylist *) ;
const char * cr2res_pfits_get_protype(const cpl_propertylist *) ;
const char * cr2res_pfits_get_wlen_id(const cpl_propertylist *) ;
const char * cr2res_pfits_get_arcfile(const cpl_propertylist *) ;
const char * cr2res_pfits_get_progid(const cpl_propertylist *) ;


double cr2res_pfits_get_ra(const cpl_propertylist *) ;
double cr2res_pfits_get_drot_posang(const cpl_propertylist *) ;
double cr2res_pfits_get_dec(const cpl_propertylist *) ;
double cr2res_pfits_get_nodthrow(const cpl_propertylist *) ;
double cr2res_pfits_get_dit(const cpl_propertylist *) ;
double cr2res_pfits_get_cwlen(const cpl_propertylist * plist) ;
double cr2res_pfits_get_mjd_obs(const cpl_propertylist * plist) ;
double cr2res_pfits_get_wstrt(const cpl_propertylist * plist, int order_idx) ;
double cr2res_pfits_get_wend(const cpl_propertylist * plist, int order_idx) ;

int cr2res_pfits_get_naxis1(const cpl_propertylist * plist) ;
int cr2res_pfits_get_naxis2(const cpl_propertylist * plist) ;
int cr2res_pfits_get_expno(const cpl_propertylist * plist) ;
int cr2res_pfits_get_nexp(const cpl_propertylist * plist) ;
int cr2res_pfits_get_ndit(const cpl_propertylist * plist) ;
int cr2res_pfits_get_obs_id(const cpl_propertylist * plist) ;
int cr2res_pfits_get_order_zp(const cpl_propertylist * plist) ;
int cr2res_pfits_get_order(const cpl_propertylist * plist) ;
int cr2res_pfits_get_order_idx(const cpl_propertylist * plist,double yposition);

cr2res_decker cr2res_pfits_get_decker_position(const cpl_propertylist * plist) ;

#endif
