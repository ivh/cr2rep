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

/* TODO */
/*----------------------------------------------------------------------------*/
/**
  @brief    find out the order number closest to the passed y position
  @param    plist       property list to read from
  @param    yposition   Y position
  @param    detector    Detector number 1-3
  @return   the requested value
 */
/*----------------------------------------------------------------------------*/
int cr2res_pfits_get_order(
            const cpl_propertylist * plist,
            double yposition)
{
    /* Check entries */
    if (plist == NULL) return -1 ;
    if (yposition < 1) return -1 ;

    char    *   key_name ;
    int         i, best_number;
    double      curr_diff;
    double      best_diff = 2048.0;

    for (i=0; i<9; i++) { // NOTE: There are never more than 9 orders!
        key_name = cpl_sprintf("ESO INS WLEN CENY%d",i);
        curr_diff = fabs(yposition -
                        cpl_propertylist_get_double(plist, key_name));
        if (curr_diff < best_diff){
               best_diff = curr_diff;
               best_number = i;
        }
    }
    if (best_diff > 50.0)
        cpl_msg_warning(__func__,
                "Order %d identified with large difference of %.1f pix",  
                best_number, best_diff);
    cpl_free(key_name) ;
    return best_number ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    find out the Min wavelength for an order
  @param    plist       property list to read from
  @return   the requested value
 */
/*----------------------------------------------------------------------------*/
double cr2res_pfits_get_wmin(const cpl_propertylist * plist, int order)
{
    char    *   key_name ;
    double      val  ;

    /* Check entries */
    if (plist == NULL) return -1.0 ;
    if (order < 0) return -1.0 ;

    /* Create key name */
    key_name = cpl_sprintf("ESO INS WLEN STRT%d", order) ;

    /* Get the value */
    val = cpl_propertylist_get_double(plist, key_name) ;

    cpl_free(key_name) ;
    return val ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    find out the Max wavelength for an order
  @param    plist       property list to read from
  @return   the requested value
 */
/*----------------------------------------------------------------------------*/
double cr2res_pfits_get_wmax(const cpl_propertylist * plist, int order)
{
    char    *   key_name ;
    double      val  ;

    /* Check entries */
    if (plist == NULL) return -1.0 ;
    if (order < 0) return -1.0 ;

    /* Create key name */
    key_name = cpl_sprintf("ESO INS WLEN END%d", order) ;

    /* Get the value */
    val = cpl_propertylist_get_double(plist, key_name) ;

    cpl_free(key_name) ;
    return val ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    find out the Y pos of an order
  @param    plist       property list to read from
  @return   the requested value
 */
/*----------------------------------------------------------------------------*/
double cr2res_pfits_get_ypos(const cpl_propertylist * plist, int order)
{
    char    *   key_name ;
    double      val  ;

    /* Check entries */
    if (plist == NULL) return -1.0 ;
    if (order < 0) return -1.0 ;

    /* Create key name */
    key_name = cpl_sprintf("ESO INS WLEN CENY%d", order) ;

    /* Get the value */
    val = cpl_propertylist_get_double(plist, key_name) ;

    cpl_free(key_name) ;
    return val ;
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
    return (const char *) cpl_propertylist_get_string(plist, "ARCFILE");
}



/**@}*/
