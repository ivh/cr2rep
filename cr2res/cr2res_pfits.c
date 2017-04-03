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

/*----------------------------------------------------------------------------*/
/**
  @brief    find out the Min wavelength for an order
  @param    plist       property list to read from
  @return   the requested value
 */
/*----------------------------------------------------------------------------*/
double kmos_pfits_get_wmin(const cpl_propertylist * plist, int order)
{
    char    *   key_name ;
    double      val  ;

    /* Check entries */
    if (plist == NULL) return -1.0 ;
    if (order < 0) return -1.0 ;

    /* Create key name */
    key_name = cpl_sprintf("ESO WMIN_%02d", order) ;

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
double kmos_pfits_get_wmax(const cpl_propertylist * plist, int order)
{
    char    *   key_name ;
    double      val  ;

    /* Check entries */
    if (plist == NULL) return -1.0 ;
    if (order < 0) return -1.0 ;

    /* Create key name */
    key_name = cpl_sprintf("ESO WMAX_%02d", order) ;

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
double kmos_pfits_get_ypos(const cpl_propertylist * plist, int order)
{
    char    *   key_name ;
    double      val  ;

    /* Check entries */
    if (plist == NULL) return -1.0 ;
    if (order < 0) return -1.0 ;

    /* Create key name */
    key_name = cpl_sprintf("ESO YPOS_%02d", order) ;

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
