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
#define CR2RES_HEADER_NODTHROW          "ESO SEQ NODTHROW"
#define CR2RES_HEADER_ARCFILE           "ARCFILE"
#define CR2RES_HEADER_MJDOBS            "MJD-OBS"
#define CR2RES_HEADER_WLEN_ID           "ESO INS WLEN ID"
#define CR2RES_HEADER_WLEN_BEGIN        "ESO INS WLEN BEGIN%d"
#define CR2RES_HEADER_WLEN_END          "ESO INS WLEN END%d"
#define CR2RES_HEADER_WLEN_CENY         "ESO INS WLEN CENY%d"
#define CR2RES_HEADER_GRAT1_ZPORD       "ESO INS GRAT1 ZP_ORD"
#define CR2RES_HEADER_GRAT1_ORDER       "ESO INS GRAT1 ORDER"
#define CR2RES_HEADER_NDIT              "ESO DET NDIT"
#define CR2RES_HEADER_DIT               "ESO DET SEQ1 DIT"
#define CR2RES_HEADER_PROG_ID           "ESO OBS PROG ID"
#define CR2RES_HEADER_OBS_ID            "ESO OBS ID"


/* QC Parameter Names */
#define CR2RES_HEADER_QC_DARK_RON1          "ESO QC DARK RON1"
#define CR2RES_HEADER_QC_DARK_RON2          "ESO QC DARK RON2"
#define CR2RES_HEADER_QC_DARK_MEAN          "ESO QC DARK MEAN"
#define CR2RES_HEADER_QC_DARK_MEDIAN        "ESO QC DARK MEDIAN"
#define CR2RES_HEADER_QC_DARK_STDEV         "ESO QC DARK STDDEV"
#define CR2RES_HEADER_QC_DARK_NBAD          "ESO QC DARK NBAD"

#define CR2RES_HEADER_QC_DETLIN_NBAD        "ESO QC DETLIN NBBAD"
#define CR2RES_HEADER_QC_DETLIN_NBFAILED    "ESO QC DETLIN NBFAILED"
#define CR2RES_HEADER_QC_DETLIN_NBSUCCESS   "ESO QC DETLIN NBSUCCESS"
#define CR2RES_HEADER_QC_DETLIN_FITQUALITY  "ESO QC DETLIN FITQUALITY"
#define CR2RES_HEADER_QC_DETLIN_MEDIAN      "ESO QC DETLIN MEDIAN"
#define CR2RES_HEADER_QC_DETLIN_GAIN        "ESO QC DETLIN GAIN"
#define CR2RES_HEADER_QC_DETLIN_MINLEVEL    "ESO QC DETLIN MINLEVEL"
#define CR2RES_HEADER_QC_DETLIN_MAXLEVEL    "ESO QC DETLIN MAXLEVEL"

#define CR2RES_HEADER_QC_FLAT_ORDERPOS      "ESO QC FLAT ORDERPOS"
#define CR2RES_HEADER_QC_FLAT_MEAN          "ESO QC FLAT MEAN"
#define CR2RES_HEADER_QC_FLAT_MEDIAN        "ESO QC FLAT MEDIAN"
#define CR2RES_HEADER_QC_FLAT_FLUX          "ESO QC FLAT FLUX"
#define CR2RES_HEADER_QC_FLAT_RMS           "ESO QC FLAT RMS"
#define CR2RES_HEADER_QC_FLAT_S2N           "ESO QC FLAT S2N"
#define CR2RES_HEADER_QC_FLAT_OVEREXPOSED   "ESO QC FLAT OVEREXPOSED"
#define CR2RES_HEADER_QC_FLAT_TRACE_CENTERY "ESO QC FLAT TRACE CENTERY"
#define CR2RES_HEADER_QC_FLAT_NBBAD         "ESO QC FLAT NBBAD"

#define CR2RES_HEADER_QC_WAVE_BESTXCORR     "ESO QC WAVE BESTXCORR"
#define CR2RES_HEADER_QC_WAVE_CENTWL        "ESO QC WAVE CENTWL"
#define CR2RES_HEADER_QC_WAVE_RESOL_FWHM    "ESO QC WAVE RESOL FWHM"
#define CR2RES_HEADER_QC_WAVE_RESOL         "ESO QC WAVE RESOL"
#define CR2RES_HEADER_QC_WAVE_DISPWL        "ESO QC WAVE DISPWL"
#define CR2RES_HEADER_QC_WAVE_LAMP_EFFIC    "ESO QC WAVE LAMP EFFIC"
#define CR2RES_HEADER_QC_WAVE_OVEREXPOSED   "ESO QC WAVE OVEREXPOSED"

#define CR2RES_HEADER_QC_SIGNAL             "ESO QC SIGNAL"
#define CR2RES_HEADER_QC_SNR                "ESO QC SNR%d"
#define CR2RES_HEADER_QC_SLITFWHM_ORDER     "ESO QC SLITFWHM%d"
#define CR2RES_HEADER_QC_SLITFWHM_MED       "ESO QC SLITFWHM MED"
#define CR2RES_HEADER_QC_REAL_ORDER         "ESO QC REALORDER%d"

/*-----------------------------------------------------------------------------
                                   Functions prototypes
 -----------------------------------------------------------------------------*/

cr2res_nodding_pos cr2res_pfits_get_nodding_pos(const cpl_propertylist * plist);

const char * cr2res_pfits_get_procatg(const cpl_propertylist *) ;
const char * cr2res_pfits_get_protype(const cpl_propertylist *) ;
const char * cr2res_pfits_get_wlen_id(const cpl_propertylist *) ;
const char * cr2res_pfits_get_arcfile(const cpl_propertylist *) ;
const char * cr2res_pfits_get_progid(const cpl_propertylist *) ;


double cr2res_pfits_get_ra(const cpl_propertylist *) ;
double cr2res_pfits_get_dec(const cpl_propertylist *) ;
double cr2res_pfits_get_nodthrow(const cpl_propertylist *) ;
double cr2res_pfits_get_dit(const cpl_propertylist *) ;
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
