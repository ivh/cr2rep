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

#include <cpl.h>
#include <math.h>
#include "cr2res_pfits.h"
#include "cr2res_utils.h"
#include "cr2res_io.h"

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_pfits     FITS header protected access
 *
 */
/*----------------------------------------------------------------------------*/
/**@{*/

/*-----------------------------------------------------------------------------
                              Function codes
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @brief    find out the nodding position
  @param    plist       property list to read from
  @return   
 */
/*----------------------------------------------------------------------------*/
cr2res_nodding_pos cr2res_pfits_get_nodding_pos(const cpl_propertylist * plist)
{
    const char  *   sval ;
    sval = cpl_propertylist_get_string(plist, CR2RES_HEADER_NODPOS);

    if (sval==NULL) {
        cpl_error_reset() ;
        return CR2RES_NODDING_NONE ;
    }
    if (sval[0] == 'A') return CR2RES_NODDING_A ;
    if (sval[0] == 'B') return CR2RES_NODDING_B ;
    return CR2RES_NODDING_NONE ; 
}

/*----------------------------------------------------------------------------*/
/**
  @brief    find out the PRO.CATG
  @param    plist       property list to read from
  @return   pointer to statically allocated character string
 */
/*----------------------------------------------------------------------------*/
const char * cr2res_pfits_get_procatg(const cpl_propertylist * plist)
{
    return (const char *) cpl_propertylist_get_string(plist, CPL_DFS_PRO_CATG);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    find out the PRO.TYPE
  @param    plist       property list to read from
  @return   pointer to statically allocated character string
 */
/*----------------------------------------------------------------------------*/
const char * cr2res_pfits_get_protype(const cpl_propertylist * plist)
{
    return (const char *) cpl_propertylist_get_string(plist, CPL_DFS_PRO_TYPE);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    find out the arcfile
  @param    plist       property list to read from
  @return   pointer to statically allocated character string
 */
/*----------------------------------------------------------------------------*/
const char * cr2res_pfits_get_arcfile(const cpl_propertylist * plist)
{
    return (const char *) cpl_propertylist_get_string(plist, 
            CR2RES_HEADER_ARCFILE);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    find out the PROG ID
  @param    plist       property list to read from
  @return   pointer to statically allocated character string
 */
/*----------------------------------------------------------------------------*/
const char * cr2res_pfits_get_progid(const cpl_propertylist * plist)
{
    return (const char *) cpl_propertylist_get_string(plist,
            CR2RES_HEADER_PROG_ID);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    find out the Setting
  @param    plist       property list to read from
  @return   pointer to statically allocated character string
 */
/*----------------------------------------------------------------------------*/
const char * cr2res_pfits_get_wlen_id(const cpl_propertylist * plist)
{
    return (const char *) cpl_propertylist_get_string(plist, 
            CR2RES_HEADER_WLEN_ID);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    find out the DEC
  @param    plist       property list to read from
  @return   the requested value
 */
/*----------------------------------------------------------------------------*/
double cr2res_pfits_get_dec(const cpl_propertylist * plist)
{
    return cpl_propertylist_get_double(plist, "DEC") ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    find out the RA
  @param    plist       property list to read from
  @return   the requested value
 */
/*----------------------------------------------------------------------------*/
double cr2res_pfits_get_ra(const cpl_propertylist * plist)
{
    return cpl_propertylist_get_double(plist, "RA") ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    find out the DROT POSANG value
  @param    plist       property list to read from
  @return   the requested value
 */
/*----------------------------------------------------------------------------*/
double cr2res_pfits_get_drot_posang(const cpl_propertylist * plist)
{
    return cpl_propertylist_get_double(plist, CR2RES_HEADER_DROT_POSANG)  ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    find out the NODTHROW value
  @param    plist       property list to read from
  @return   the requested value
 */
/*----------------------------------------------------------------------------*/
double cr2res_pfits_get_nodthrow(const cpl_propertylist * plist)
{
    return cpl_propertylist_get_double(plist, CR2RES_HEADER_NODTHROW)  ;
}


/*----------------------------------------------------------------------------*/
/**
  @brief    find out the DIT value
  @param    plist       property list to read from
  @return   the requested value
 */
/*----------------------------------------------------------------------------*/
double cr2res_pfits_get_dit(const cpl_propertylist * plist)
{
    double val;
    val = cpl_propertylist_get_double(plist, CR2RES_HEADER_DIT)  ;
    if (cpl_error_get_code() == CPL_ERROR_TYPE_MISMATCH){
        cpl_error_reset();
        val = (double) cpl_propertylist_get_int(plist, CR2RES_HEADER_DIT);
    }
    return val;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    find out the MJD-OBS value 
  @param    plist       property list to read from
  @return   the requested value
 */
/*----------------------------------------------------------------------------*/
double cr2res_pfits_get_mjd_obs(const cpl_propertylist * plist)
{
    return cpl_propertylist_get_double(plist, CR2RES_HEADER_MJDOBS)  ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    find out the CWLEN value 
  @param    plist       property list to read from
  @return   the requested value
 */
/*----------------------------------------------------------------------------*/
double cr2res_pfits_get_cwlen(const cpl_propertylist * plist)
{
    return cpl_propertylist_get_double(plist, CR2RES_HEADER_WLEN_CWLEN)  ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    find out the Start wavelength for an order_idx (current detector)
  @param    plist       property list to read from
  @param    order_idx   Order INDEX (-49 -> 50)
  @return   the requested value
 */
/*----------------------------------------------------------------------------*/
double cr2res_pfits_get_wstrt(const cpl_propertylist * plist, int order_idx)
{
    char    *   key_name ;
    int         order_idxp ;
    double      val  ;

    /* Check entries */
    if (plist == NULL) return -1.0 ;

    /* Conversion order_idx -> order_idxp */
    if ((order_idxp = cr2res_io_convert_order_idx_to_idxp(order_idx)) < 0) 
        return -1.0 ;

    /* Create key name */
    key_name = cpl_sprintf(CR2RES_HEADER_WLEN_BEGIN, order_idxp) ;

    /* Get the value */
    val = cpl_propertylist_get_double(plist, key_name) ;

    cpl_free(key_name) ;
    return val ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    find out the End wavelength for an order_idx (current detector)
  @param    plist       property list to read from
  @param    order_idx   Order INDEX (-49 -> 50)
  @return   the requested value
 */
/*----------------------------------------------------------------------------*/
double cr2res_pfits_get_wend(const cpl_propertylist * plist, int order_idx)
{
    char    *   key_name ;
    int         order_idxp ;
    double      val  ;

    /* Check entries */
    if (plist == NULL) return -1.0 ;

    /* Conversion order_idx -> order_idxp */
    if ((order_idxp = cr2res_io_convert_order_idx_to_idxp(order_idx)) < 0) 
        return -1.0 ;

    /* Create key name */
    key_name = cpl_sprintf(CR2RES_HEADER_WLEN_END, order_idxp) ;

    /* Get the value */
    val = cpl_propertylist_get_double(plist, key_name) ;

    cpl_free(key_name) ;
    return val ;
}

/*----------------------------------------------------------------------------*/
  /*
  @brief    find out the NAXIS1 value
  @param    plist       property list to read from
  @return   the requested value
 */
/*----------------------------------------------------------------------------*/
int cr2res_pfits_get_naxis1(const cpl_propertylist * plist)
{
    return cpl_propertylist_get_int(plist, CR2RES_HEADER_NAXIS1)  ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    find out the NAXIS2 value
  @param    plist       property list to read from
  @return   the requested value
 */
/*----------------------------------------------------------------------------*/
int cr2res_pfits_get_naxis2(const cpl_propertylist * plist)
{
    return cpl_propertylist_get_int(plist, CR2RES_HEADER_NAXIS2)  ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    find out the NEXP value
  @param    plist       property list to read from
  @return   the requested value
 */
/*----------------------------------------------------------------------------*/
int cr2res_pfits_get_nexp(const cpl_propertylist * plist)
{
    return cpl_propertylist_get_int(plist, CR2RES_HEADER_NEXP)  ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    find out the EXPNO value
  @param    plist       property list to read from
  @return   the requested value
 */
/*----------------------------------------------------------------------------*/
int cr2res_pfits_get_expno(const cpl_propertylist * plist)
{
    return cpl_propertylist_get_int(plist, CR2RES_HEADER_EXPNO)  ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    find out the NDIT value
  @param    plist       property list to read from
  @return   the requested value
 */
/*----------------------------------------------------------------------------*/
int cr2res_pfits_get_ndit(const cpl_propertylist * plist)
{
    return cpl_propertylist_get_int(plist, CR2RES_HEADER_NDIT)  ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    find out the OBS ID value 
  @param    plist       property list to read from
  @return   the requested value
 */
/*----------------------------------------------------------------------------*/
int cr2res_pfits_get_obs_id(const cpl_propertylist * plist)
{
    return cpl_propertylist_get_int(plist, CR2RES_HEADER_OBS_ID)  ;
}

/*----------------------------------------------------------------------------*/
  /*
  @brief    find out the ORDER value
  @param    plist       property list to read from
  @return   the requested value
 */
/*----------------------------------------------------------------------------*/
int cr2res_pfits_get_order(const cpl_propertylist * plist)
{
    return cpl_propertylist_get_int(plist, CR2RES_HEADER_GRAT1_ORDER)  ;
}

/*----------------------------------------------------------------------------*/
  /*
  @brief    find out the ZP_ORD value
  @param    plist       property list to read from
  @return   the requested value
 */
/*----------------------------------------------------------------------------*/
int cr2res_pfits_get_order_zp(const cpl_propertylist * plist)
{
    return cpl_propertylist_get_int(plist, CR2RES_HEADER_GRAT1_ZPORD)  ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    find out the order_idx closest to the passed y position
  @param    plist       property list to read from
  @param    yposition   Y position
  @return   the order_idx (-49 -> 50)
 */
/*----------------------------------------------------------------------------*/
int cr2res_pfits_get_order_idx(
        const cpl_propertylist  *   plist,
        double                      yposition)
{
    char    *   key_name ;
    int         i, order_idxp;
    int         order_idx = -1000 ;
    int         min_order_idx = -49 ;
    int         max_order_idx =  50 ;
    double      ycen, curr_diff;
    double      best_diff = CR2RES_DETECTOR_SIZE;

    /* Check entries */
    if (plist == NULL) return -1 ;
    if (yposition < 1) return -1 ;
    if (cpl_error_get_code() != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Error already set - abort") ;
        return -1 ;
    }

    for (i=min_order_idx ; i <= max_order_idx ; i++) {
        order_idxp = cr2res_io_convert_order_idx_to_idxp(i);
        key_name = cpl_sprintf(CR2RES_HEADER_WLEN_CENY, order_idxp) ;
        ycen = cpl_propertylist_get_double(plist, key_name);
        cpl_free(key_name) ;
        if (cpl_error_get_code() != CPL_ERROR_NONE) {
            cpl_error_reset();
            continue;
        }
        curr_diff = fabs(yposition - ycen );
        if (curr_diff < best_diff){
               best_diff = curr_diff;
               order_idx = i;
        }
    }

    cpl_msg_debug(__func__, 
        "Order %d identified with %.1f pix difference from expectation",
        order_idx, best_diff);
    if (best_diff > 80.0)
        cpl_msg_warning(__func__,
                "Order %d identified with large difference of %.1f pix",
                order_idx, best_diff);
    return order_idx ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    find out the decker position
  @param    plist       property list to read from
  @return   0, 1 or 2 or -1 in error case
 */
/*----------------------------------------------------------------------------*/
cr2res_decker cr2res_pfits_get_decker_position(const cpl_propertylist * plist)
{
    int         decker_value ;
    decker_value = cpl_propertylist_get_int(plist, CR2RES_HEADER_DECKER_POS);
    if (cpl_error_get_code() != CPL_ERROR_NONE) {
        cpl_error_reset();
        return CR2RES_DECKER_INVALID ;
    }
    /* TODO: enter proper encoder values */
    if (decker_value == 3) return CR2RES_DECKER_NONE ;
    if (decker_value == 1) return CR2RES_DECKER_1_3 ;
    if (decker_value == 2) return CR2RES_DECKER_2_4 ;
    return CR2RES_DECKER_INVALID ;
}

/**@}*/
