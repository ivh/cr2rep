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

#include <math.h>
#include <string.h>

#include <cpl.h>
#include "cr2res_wave.h"
#include "cr2res_io.h"
#include "cr2res_pfits.h"

/*-----------------------------------------------------------------------------
                                   Defines
 -----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
                                Functions prototypes
 -----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_wave        Wavelength Calibration
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief   Top level function that takes spectrum, returns solution.
  @param    spectrum        Input spectrum: arc lamp, etalon etc.
  @param    initial_guess   Starting wavelength solution
  @param    catalog         Line catalog, wavelengths, strengths
  @param    template        Template spectrum (flux, wavelengths)
  @return  Wavelength solution, i.e. polynomial that translates pixel
            values to wavelength.
 */
/*----------------------------------------------------------------------------*/
cpl_polynomial * cr2res_wave(
        cpl_vector      *   spectrum,
        cpl_polynomial  *   initial_guess,
        cpl_table       *   catalog,
        cpl_bivector    *   template)
{
    return NULL ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief  Find solution by cross-correlating template spectrum
  @param    spectrum        Input spectrum
  @param    initial_guess   Starting wavelength solution
  @param    template        Template spectrum (flux, wavelengths)
  @return  Wavelength solution, i.e. polynomial that translates pixel
            values to wavelength.

    TODO: Summarize method
 */
/*----------------------------------------------------------------------------*/
cpl_polynomial * cr2res_wave_xcorr(
        cpl_vector      *   spectrum,
        cpl_polynomial  *   initial_guess,
        cpl_bivector    *   template)
{
    return NULL ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief   Find solution by finding lines and fitting
  @param    spectrum        Input spectrum
  @param    initial_guess   Starting wavelength solution
  @param    catalog         Line catalog, wavelengths, strengths
  @param    template        Template spectrum (flux, wavelengths)
  @return  Wavelength solution, i.e. polynomial that translates pixel
            values to wavelength.
 */
/*----------------------------------------------------------------------------*/
cpl_polynomial * cr2res_wave_line_fitting(
        cpl_vector      *   spectrum,
        cpl_polynomial  *   initial_guess,
        cpl_table       *   catalog)
{
    return NULL ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief   Find solution from etalon
  @param    spectrum        Input spectrum: etalon
  @param    initial_guess   Starting wavelength solution
  @return  Wavelength solution, i.e. polynomial that translates pixel
            values to wavelength.

    This function uses the intrinsic property of the etalon spectrum,
    that the lines are *supposed* to be equi-spaced, to refine the wavelength
    solution. The input solution needs to be good enough for the zero-point
    since the etalon spectrum carries no information about the absolute
    wavelength scale.

    The method involves these steps:
    * Identify lines (thresholding)
    * Determine line centers (ceter of gravity or gauss fit, needs testing)
    * Subtract x-coods from subsequent lines, i.e. measure d many times
        d is the distance between fringes in pixels, assumed to be constant
    * Determine mis-counts in fringes, by looking at outliers in d-distribution,
        re-count with e.g. half-distance between two fringes, if one is missing.
    * Fit the d-distibution to determine "true d"
    * Use d to calculate new x-coodinates, use input-solution to translate to
        wavelength
    * Fit measured x-coords to new lamdas to get new solution.
 */
/*----------------------------------------------------------------------------*/
cpl_polynomial * cr2res_wave_etalon(
        cpl_vector      *   spectrum,
        cpl_polynomial  *   initial_guess)
{
    return NULL ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief
  @param
  @return
 */
/*----------------------------------------------------------------------------*/
cpl_vector * cr2res_wave_line_detection(
        cpl_vector      *   spectrum)
{
    return NULL ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief
  @param
  @return
 */
/*----------------------------------------------------------------------------*/
cpl_vector * cr2res_wave_gen_spectrum(
        cpl_table       *   catalog,
        cpl_polynomial  *   initial_guess)
{
    return NULL ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Compute the wavelength polynomial from boundaries
  @param    wmin    First pixel wavelength
  @param    wmax    Last pixel wavelength
  @return   the polynomial or NULL in error case

  The returned polynomial must be deallocated with cpl_polynomial_delete()
 */
/*----------------------------------------------------------------------------*/
cpl_array * cr2res_wave_get_estimate(
        const char	*	filename,
        int				detector,
        int             order)
{
    double                  wmin, wmax ;
    cpl_array           *   wl ;
    cpl_propertylist    *   plist ;
    int                     wished_ext_nb ;

    /* Check Entries */
    if (filename == NULL) return NULL ;
    if (order < 0 || detector <= 0 || detector > CR2RES_NB_DETECTORS)
        return NULL ;

    /* Load the propertylist */
    wished_ext_nb = cr2res_io_get_ext_idx(filename, detector) ;
    plist = cpl_propertylist_load(filename, wished_ext_nb) ;

    /* Get the values for this order */
    wmin = kmos_pfits_get_wmin(plist, order) ;
    wmax = kmos_pfits_get_wmax(plist, order) ;
    cpl_propertylist_delete(plist) ;
    if (cpl_error_get_code() != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot Get the WMIN/WMAX from the header") ;
        return NULL ;
    }

    /* Create the array */
    wl = cpl_array_new(2, CPL_TYPE_DOUBLE) ;
    cpl_array_set(wl, 0, wmin) ;
    cpl_array_set(wl, 1, wmax) ;

    return wl ;
}

/**@}*/
