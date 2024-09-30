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

#include <string.h>
#include <math.h>
#include <cpl.h>

#include "cr2res_photom.h"
#include "cr2res_dfs.h"
#include "cr2res_pfits.h"
#include "cr2res_utils.h"
#include "cr2res_extract.h"

/*-----------------------------------------------------------------------------
                            Functions Prototypes
 -----------------------------------------------------------------------------*/

static double cr2res_photom_throughput_avg(
        const cpl_vector    *   vec,
        const cpl_vector    *   wl,
        const char          *   setting) ;

static cpl_vector * cr2res_photom_conversion(
        const cpl_bivector  *   star,
        const cpl_bivector  *   photflux,
        double                  exptime) ;

static cpl_vector * cr2res_photom_throughput(
        const cpl_vector    *   conversion,
        const cpl_vector    *   wl,
        double                  gain) ;

static cpl_vector * cr2res_photom_sensitivity(
        const cpl_vector    *   conversion,
        const cpl_vector    *   sigma,
        double                  exptime) ;

static int cr2res_photom_star_find(
        const cpl_vector    *   v_ra, 
        const cpl_vector    *   v_dec,
        double                  ra, 
        double                  dec, 
        double                  maxdist, 
        double              *   pdist) ;

static double cr2res_photom_great_circle_dist(
        double      ra1, 
        double      dec1,
        double      ra2, 
        double      dec2) ;

/*----------------------------------------------------------------------------*/
/**
 * @defgroup cr2res_photom      Photometry
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Compute the conversion/throughput/sensitivity high level
  @param    extr        the extracted spectrum table
  @param    std_star_tab The STD STARS table
  @param    setting     The setting
  @param    ra          RA
  @param    dec         DEC
  @param    gain        The gain
  @param    exptime     The exposure time
  @param    display         Flag to allow display
  @param    display_order   Order to display
  @param    display_trace   Trace to display 
  @param    throughput  [out]   the throughput table
  @param    ext_plist   [out]   the QCs
  @return   0 if ok, -1 otherwise
 */
/*----------------------------------------------------------------------------*/
int cr2res_photom_engine(
        const cpl_table     *   extr,
        const char          *   std_star_file,
        const char          *   setting,
        double                  ra,
        double                  dec,
        double                  gain,
        double                  exptime,
        int                     display,
        int                     display_order,
        int                     display_trace,
        cpl_table           **  throughput,
        cpl_propertylist    **  ext_plist)
{
    cpl_table       *   std_star_tab ;
    cpl_bivector    *   std_star_biv ;
    cpl_table       *   throughput_loc ;
    cpl_array       *   col_names ;
    char            *   cwname ;
    char            *   csname ;
    char            *   ccname ;
    char            *   ctname ;
    cpl_bivector    *   spec_biv ;
    cpl_bivector    *   spec_err_biv ;
    cpl_vector      *   conversion_vec ;
    cpl_vector      *   throughput_vec ;
    cpl_vector      *   sensitivity_vec ;
    double              throughput_avg_val, throughput_center ;
    int                 throughput_avg_nb ;
    cpl_propertylist *  plist ;
    const double    *   pdata ;
    cpl_size            j, ncols ;
    int                 order, trace_nb ;

    /* Test entries */
    if (extr==NULL || std_star_file==NULL || setting == NULL || 
            throughput==NULL || ext_plist == NULL) return -1 ;

    /* Initialise */
    throughput_avg_val = 0.0 ;
    throughput_avg_nb = 0 ;

    /* Load the std stars table */
    if ((std_star_tab = cpl_table_load(std_star_file, 1, 0)) == NULL) {
        cpl_msg_error(__func__, "Cannot load the std star table") ;
        return -1 ;
    }

    /* Search for the proper star */
    if ((std_star_biv = cr2res_photom_conv_get_star(std_star_tab, ra, dec)) == 
            NULL) {
        cpl_msg_error(__func__, "Cannot find the star flux") ;
        cpl_table_delete(std_star_tab) ;
        return -1 ;
    }
    cpl_table_delete(std_star_tab) ;

    /* Create the output table */
    throughput_loc = cpl_table_new(cpl_table_get_nrow(extr)) ;

    /* Loop on the extracted spectra */
	col_names = cpl_table_get_column_names(extr);
	ncols = cpl_table_get_ncol(extr) ;
	for (j=0 ; j<ncols ; j++) {
        const char *col_name;
        char *col_type;
        col_name = cpl_array_get_string(col_names, j);
		col_type = cr2res_dfs_SPEC_colname_parse(col_name, &order, &trace_nb) ;
		if (col_type != NULL && !strcmp(col_type, CR2RES_COL_SPEC_SUFFIX)) {
            /* Get this extracted spectrum */
			if (cr2res_extract_EXTRACT1D_get_spectrum(extr, order, trace_nb,
						&spec_biv, &spec_err_biv)) {
                cpl_msg_error(__func__, "Cannot get the extracted spectrum") ;
                cpl_free(col_type) ;
                continue ;
            }
            cpl_msg_info(__func__, 
                    "Conversion/Sensitivity/Throughput for order/trace: %d/%d",
                    order, trace_nb) ;

            /* Conversion */
            if ((conversion_vec = cr2res_photom_conversion(spec_biv, 
                            std_star_biv, exptime))==NULL) {
                cpl_msg_warning(__func__, 
                        "Cannot compute the conversion factor");
                cpl_error_reset() ;
                cpl_bivector_delete(spec_biv) ;
                cpl_bivector_delete(spec_err_biv) ;
                cpl_free(col_type) ;
                continue ;
            }

            /* Throughput */
            if ((throughput_vec = cr2res_photom_throughput(conversion_vec,
                            cpl_bivector_get_x_const(spec_biv), gain))==NULL) {
                cpl_msg_warning(__func__, "Cannot compute the throughput");
                cpl_error_reset() ;
                cpl_bivector_delete(spec_biv) ;
                cpl_bivector_delete(spec_err_biv) ;
                cpl_vector_delete(conversion_vec) ;
                if (col_type != NULL) cpl_free(col_type) ;
                continue ;
            }

            /* Sensitivity */
            if ((sensitivity_vec = cr2res_photom_sensitivity(conversion_vec, 
                            cpl_bivector_get_y_const(spec_err_biv), 
                            exptime))==NULL) {
                cpl_msg_warning(__func__, "Cannot compute the sensitivity");
                cpl_error_reset() ;
                cpl_bivector_delete(spec_biv) ;
                cpl_bivector_delete(spec_err_biv) ;
                cpl_vector_delete(conversion_vec) ;
                cpl_vector_delete(throughput_vec) ;
                if (col_type != NULL) cpl_free(col_type) ;
                continue ;
            }

            /* Add  to the table */
			/* Create/Fill WAVELENGTH column */
			cwname = cr2res_dfs_WAVELENGTH_colname(order, trace_nb) ;
        	cpl_table_new_column(throughput_loc, cwname, CPL_TYPE_DOUBLE);
            pdata = cpl_bivector_get_x_data_const(spec_biv) ;
			cpl_table_copy_data_double(throughput_loc, cwname, pdata) ;

            /* Create/Fill CONVERSION column */
			ccname = cr2res_dfs_CONVERSION_colname(order, trace_nb) ;
        	cpl_table_new_column(throughput_loc, ccname, CPL_TYPE_DOUBLE);
            pdata = cpl_vector_get_data_const(conversion_vec) ;
			cpl_table_copy_data_double(throughput_loc, ccname, pdata) ;

            /* Create/Fill SENSITIVITY column */
			csname = cr2res_dfs_SENSITIVITY_colname(order, trace_nb) ;
        	cpl_table_new_column(throughput_loc, csname, CPL_TYPE_DOUBLE);
            pdata = cpl_vector_get_data_const(sensitivity_vec) ;
			cpl_table_copy_data_double(throughput_loc, csname, pdata) ;

            /* Create/Fill THROUGHPUT column */
			ctname = cr2res_dfs_THROUGHPUT_colname(order, trace_nb) ;
        	cpl_table_new_column(throughput_loc, ctname, CPL_TYPE_DOUBLE);
            pdata = cpl_vector_get_data_const(throughput_vec) ;
			cpl_table_copy_data_double(throughput_loc, ctname, pdata) ;

            /* Compute QCs */
            throughput_center = cr2res_photom_throughput_avg(throughput_vec,
                        cpl_bivector_get_x_const(spec_biv), setting) ;
            if (fabs(throughput_center+1.0) > 1e-3) {
                throughput_avg_nb ++ ;
                throughput_avg_val += throughput_center ;
            }

            cpl_vector_delete(conversion_vec) ;
            cpl_vector_delete(sensitivity_vec) ;
            cpl_vector_delete(throughput_vec) ;
            cpl_bivector_delete(spec_biv) ;
            cpl_bivector_delete(spec_err_biv) ;
            
            /* Plot on request */
            if (display && display_trace==trace_nb && display_order==order) {
                cpl_plot_column(
"set grid;set xlabel 'Wavelength (nm)';set ylabel 'Conversion (ADU/sec/Jy)';",
                    "t 'Conversion factor' w lines", "",
                    throughput_loc, cwname, ccname) ;
                cpl_plot_column(
"set grid;set xlabel 'Wavelength (nm)';set ylabel 'Sensitivity (Jy/10sig/hour)';",
                    "t 'Sensitivity' w lines", "", 
                    throughput_loc, cwname, csname) ;
                cpl_plot_column(
"set grid;set xlabel 'Wavelength (nm)';set ylabel 'Throughput (e-/photons)';",
                    "t 'Throughput' w lines", "", 
                    throughput_loc, cwname, ctname) ;
            }
			cpl_free(cwname) ;
			cpl_free(ccname) ;
			cpl_free(csname) ;
			cpl_free(ctname) ;
        }
		if (col_type != NULL) cpl_free(col_type) ;
	}
	cpl_array_delete(col_names) ;
    cpl_bivector_delete(std_star_biv) ;

    plist = cpl_propertylist_new() ;
    if (throughput_avg_nb > 0) {
        cpl_propertylist_append_double(plist, CR2RES_HEADER_QC_THROUGHPUT, 
                throughput_avg_val / throughput_avg_nb) ;
    } 
    *throughput = throughput_loc ;
    *ext_plist = plist ;
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Get the wished std star photospheric flux
  @param    tab     the std stars table
  @param    ra      RA
  @param    dec     DEC
  @return   the bivector with the flux or NULL in error case
 */
/*----------------------------------------------------------------------------*/
cpl_bivector * cr2res_photom_conv_get_star(
        const cpl_table     *   tab,
        double                  ra,
        double                  dec)
{
    cpl_vector          *   ras ;
    cpl_vector          *   decs ;
    double                  maxdist, dist ;
    int                     ind ; 
    cpl_array           *   wave ;
    cpl_array           *   flux ;
    cpl_vector          *   wave_vec ;
    cpl_vector          *   flux_vec ;
    cpl_bivector        *   flux_biv ;
    
    /* Test entries */
    if (tab == NULL) return NULL ;

    /* Initialise CR2RES_PHOTOM_STARSDIST_ARCSEC  arcsec in degrees*/
    maxdist = CR2RES_PHOTOM_STARSDIST_ARCSEC / 3600.0 ;
    
    /* Get the RAs and DECs */
    ras = cpl_vector_wrap(cpl_table_get_nrow(tab), 
            cpl_table_get_data_double((cpl_table*)tab, CR2RES_COL_RA)) ;
    decs = cpl_vector_wrap(cpl_table_get_nrow(tab), 
            cpl_table_get_data_double((cpl_table*)tab, CR2RES_COL_DEC)) ;

    /* Find the star */
    cpl_msg_info(__func__, "Search a star around RA=%g, DEC=%g in the catalog",
            ra, dec) ;
    if ((ind = cr2res_photom_star_find(ras, decs, ra, dec, maxdist,
                    NULL)) < 0) {
        cpl_msg_error(__func__, "No star found") ;
        cpl_vector_unwrap(ras) ;
        cpl_vector_unwrap(decs) ;
        return NULL ;
    }
    dist = 3600 * cr2res_photom_great_circle_dist(cpl_vector_get(ras,
                ind), cpl_vector_get(decs, ind), ra, dec) ;
    cpl_msg_info(__func__, "Star found: %s at RA=%g, DEC=%g at %g arcsec", 
            cpl_table_get_string(tab, CR2RES_COL_STDNAME, ind),
            cpl_vector_get(ras, ind), cpl_vector_get(decs, ind),
            dist) ;
    cpl_vector_unwrap(ras) ;
    cpl_vector_unwrap(decs) ;

    /* Load the star */
    wave=cpl_array_duplicate(cpl_table_get_array(tab,CR2RES_COL_PHOTOFLUX,0));
    flux=cpl_array_duplicate(cpl_table_get_array(tab,CR2RES_COL_PHOTOFLUX,ind));
    wave_vec = cpl_vector_wrap(cpl_array_get_size(wave), 
            cpl_array_get_data_double(wave)) ;
    flux_vec = cpl_vector_wrap(cpl_array_get_size(flux), 
            cpl_array_get_data_double(flux)) ;
    cpl_array_unwrap(wave) ;
    cpl_array_unwrap(flux) ;
    flux_biv=cpl_bivector_wrap_vectors(wave_vec, flux_vec) ;
    return flux_biv ;
}

/**@}*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Compute the Central avg value of a vector
  @param    vec     Input vector
  @return   the avg or -1.0 in error case
    200 central values 
 */
/*----------------------------------------------------------------------------*/
static double cr2res_photom_throughput_avg(
        const cpl_vector    *   vec,
        const cpl_vector    *   wl,
        const char          *   setting)
{
    const double    *       pvec ;
    const double    *       pwl ;
    double                  sum, wmin, wmax ; 
    cpl_size                nelem, i ;
    int                     nsum ;

    /* Check entries */
    if (vec == NULL || wl == NULL) return -1.0 ;
    nelem = cpl_vector_get_size(vec) ;
    if (cpl_vector_get_size(wl) != nelem) return -1.0 ;

    /* Initialize */
    wmin = wmax = -1.0 ;

    /* Check the setting */
    if (!strcmp(setting, "Y1028")) {
        wmin = 1045.70 ;
        wmax = 1047.00 ;
    } else if (!strcmp(setting, "J1228")) {
        wmin = 1245.50 ;
        wmax = 1246.40 ;
    } else if (!strcmp(setting, "H1567")) {
        wmin = 1695.20 ;
        wmax = 1696.20 ;
    } else if (!strcmp(setting, "K2148")) {
        wmin = 2312.60 ;
        wmax = 2313.30 ;
    } else if (!strcmp(setting, "L3377")) {
        wmin = 3828.70 ;
        wmax = 3829.30 ;
    } else if (!strcmp(setting, "M4266")) {
        wmin = 3953.00 ;
        wmax = 3953.90 ;
    } 
    if (wmin < 0.0 || wmax < 0.0) return -1.0 ;
    
    cpl_msg_info(__func__, "Compute the QC.THROUGHPUT between %g and %g nm",
            wmin, wmax) ;

    pvec = cpl_vector_get_data_const(vec) ;
    pwl = cpl_vector_get_data_const(wl) ;
    sum = 0.0 ;
    nsum = 0 ;

    for (i=0 ; i<nelem ; i++) {
        if (pwl[i]>wmin && pwl[i]<wmax) {
            sum += pvec[i] ;
            nsum++ ;
        }
    }
    if (nsum > 0) return sum/nsum ;
    else return -1.0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Compute the conversion signal
  @param    star        the star spectrum (ADU)
  @param    photflux    the star photospheric flux (from a catalog)
  @param    exptime     For the conversion to ADU/sec
  @return   the conversion signal or NULL in error case

  The returned conversion signal is a vector with the size of the star
  bivector. It corresponds, for each wavelength of the star bivector, of
  the star value divided by the photflux value (interpolated).

  conversion = star_spec / model

  star_spec in ADU/sec
  model in Jy
  conversion in ADU/sec/Jy
 */
/*----------------------------------------------------------------------------*/
static cpl_vector * cr2res_photom_conversion(
        const cpl_bivector  *   star,
        const cpl_bivector  *   photflux,
        double                  exptime)
{
    cpl_vector      *   conversion ;
    cpl_bivector    *   model ;
    double          *   pmodel ;
    int                 i ;

    /* Check entries */
    if (star == NULL) return NULL ;
    if (photflux == NULL) return NULL ;
    if (fabs(exptime) < 1e-3) return NULL ;
    
    /* Resample the model */
    model = cpl_bivector_duplicate(star) ;
    if (cpl_bivector_interpolate_linear(model, photflux) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot interpolate the signal") ;
        cpl_bivector_delete(model) ;
        return NULL ;
    }

    /* Set to 1 the values that are 0 */
    pmodel = cpl_vector_get_data(cpl_bivector_get_y(model)) ;
    for (i=0 ; i<cpl_bivector_get_size(model) ; i++) {
        if (fabs(pmodel[i]) == 0.0) pmodel[i] = 1.0 ;
    }

    /* Create the conversion */
    conversion = cpl_vector_duplicate(cpl_bivector_get_y_const(star)) ;
    cpl_vector_divide_scalar(conversion, exptime) ;

    if (cpl_vector_divide(conversion, 
                cpl_bivector_get_y(model)) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot divide the star by the model") ;
        cpl_bivector_delete(model) ;
        return NULL ;
    }
    cpl_bivector_delete(model) ;
    return conversion ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Compute the throughput signal
  @param    conversion  the conversion in ADU/sec/Jy
  @param    wl          the wavelength in nanometers
  @param    gain        the detector gain
  @return   the throughput signal or NULL in error case

  The returned throughput signal is a vector with the size of the input
  vectors.

  throughput = (conversion*gain) / ((1e-32*S*Delta_lam*1e6)/(wl*h))

  Delta_lam is the wl/pix in m./pix
  wl is the wavelength in m.
  h = 6.626e-34
  S = PI * 4^2 m^2 
  gain = *** / 0.9  phot/sec/(ADU/sec) 
 */
/*----------------------------------------------------------------------------*/
static cpl_vector * cr2res_photom_throughput(
        const cpl_vector    *   conversion,
        const cpl_vector    *   wl,
        double                  gain)
{
    double          h, S ;
    cpl_vector  *   throughput ;
    cpl_vector  *   wave_m ;
    cpl_vector  *   dwave_m ;
    double      *   pdwave_m ;
    int             i ;

    /* Check entries */
    if (conversion == NULL) return NULL ;
    if (wl == NULL) return NULL ;
    
    /* Initialise  */
    h = 6.626e-34 ;
    S = 16 * 3.141592653589793238462 ;
 
    /* Wavelength in meters */
    wave_m = cpl_vector_duplicate(wl) ;
    cpl_vector_multiply_scalar(wave_m, 1e-9) ;

    /* Wavelength per pixels */
    dwave_m = cpl_vector_duplicate(wave_m) ;
    pdwave_m = cpl_vector_get_data(dwave_m) ;
    for (i=cpl_vector_get_size(dwave_m)-1 ; i>0 ; i--) 
        pdwave_m[i] -= pdwave_m[i-1] ;
    pdwave_m[0] = pdwave_m[1] ;

    /* Compute the throughput */
    throughput = cpl_vector_duplicate(conversion) ;
    cpl_vector_multiply_scalar(throughput, gain) ;
    cpl_vector_multiply_scalar(throughput, h) ;
    cpl_vector_multiply(throughput, wave_m) ;
    cpl_vector_divide_scalar(throughput, 1e-32) ;
    cpl_vector_divide_scalar(throughput, S) ;
    cpl_vector_divide(throughput, dwave_m) ;
    cpl_vector_divide_scalar(throughput, 1e6) ;

    cpl_vector_delete(dwave_m) ;
    cpl_vector_delete(wave_m) ;
    return throughput ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Compute the sensitivity
  @param    conversion  the conversion signal
  @param    sigma       the sigma (error) on the extraction
  @param    exptime     the exposure time
  @return   the sensitivity or NULL in error case

  The input vectors must have the same size.

  sensitivity = (10 * sigma * sqrt(exptime/3600)) / conversion
  
  sigma in ADU/sec
  conversion in ADU/sec/Jy
  exptime in sec
  sensitivity in Jy/10sigma/hour
 */
/*----------------------------------------------------------------------------*/
static cpl_vector * cr2res_photom_sensitivity(
        const cpl_vector    *   conversion,
        const cpl_vector    *   sigma,
        double                  exptime)
{
    cpl_vector      *   sensitivity ;
    double              factor ;
    const double    *   pconversion ;
    double          *   psensitivity ;
    int                 i ;

    /* Check entries */
    if (conversion == NULL) return NULL ;
    if (sigma == NULL) return NULL ;
    if (exptime < 0) return NULL ;

    /* Get access  to data */
    pconversion = cpl_vector_get_data_const(conversion) ;
    sensitivity = cpl_vector_duplicate(sigma) ;
    psensitivity = cpl_vector_get_data(sensitivity) ;

    /* Compute sensitivity */
    for (i=0 ; i<cpl_vector_get_size(conversion) ; i++) {
        if (fabs(pconversion[i])<1e-5) psensitivity[i] = 0.0 ;
        else psensitivity[i] /= pconversion[i] ;
    }
    for (i=0 ; i<11 ; i++) psensitivity[i] = 0 ;

    /* Compute the factor */
    factor = sqrt(exptime/3600) * 10 ;

    /* Multiply by the factor */
    if (cpl_vector_multiply_scalar(sensitivity, factor) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot multiply by the factor") ;
        cpl_vector_delete(sensitivity) ;
        return NULL ;
    }
    return sensitivity ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Find the index of the star closest to (ra, dec)
  @param    v_ra     Vector of Right Ascensions (ra, dec) [degree]
  @param    v_dec    Vector of Declinations (ra, dec)  [degree]
  @param    ra       Right Ascension of position to search for
  @param    dec      Declination of ditto
  @param    maxdist  Maximum acceptable (non-negative) distance [degree]
  @param    pdist    Actual distance (if not NULL)
  @return   The index (starting from zero), or negative on error.

  The two vectors must be of identical length.
  pdist may be NULL.
  It is an error if no star is found.
 */
/*----------------------------------------------------------------------------*/
static int cr2res_photom_star_find(
        const cpl_vector    *   v_ra, 
        const cpl_vector    *   v_dec,
        double                  ra, 
        double                  dec, 
        double                  maxdist, 
        double              *   pdist)
{
    int                 nra, ndec, minind ;
    double              dmin;
    const double    *   pv_ra ;
    const double    *   pv_dec ;
    int                 i ;

    /* Test inputs */
    if (v_ra == NULL || v_dec == NULL) return -1 ;

    /* Initialise */
    dmin = 0 ;
    minind = -1 ;
    nra = cpl_vector_get_size(v_ra) ;
    ndec = cpl_vector_get_size(v_dec) ;
    pv_ra = cpl_vector_get_data_const(v_ra) ;
    pv_dec = cpl_vector_get_data_const(v_dec) ;

    /* Test inputs */
    if (nra != ndec) return -1 ;
    if (maxdist < 0) return -1 ;
    
    /* Find the index of the star closest to the given coordinate */
    for (i=0 ; i < nra ; i++) {
        double gdist;
        /* Get the distance */
        gdist = cr2res_photom_great_circle_dist(pv_ra[i], pv_dec[i], ra, dec);
        if (i == 0 || gdist < dmin) {
            minind = i;
            dmin = gdist;
        }
    }
    if (pdist != NULL) *pdist = dmin;

    /* Check that it is close enough */
    if (dmin > maxdist) return -1 ;
    return minind ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief   Compute the great-circle distance between two points on a sphere
  @param   ra1    Right ascension of first point [degrees]
  @param   dec1   Declination of first point [degrees]
  @param   ra2    Right ascension of second point [degrees]
  @param   dec2   Declination of second point [degrees]
  @return  Non-negative distance [degrees].
  @see     http://en.wikipedia.org/wiki/Great-circle_distance (on 2005-10-23)
 */
/*----------------------------------------------------------------------------*/
static double cr2res_photom_great_circle_dist(
        double      ra1, 
        double      dec1,
        double      ra2, 
        double      dec2)
{
  /* Convert all input from degrees to radian - and back for the result */
  const double dra  = sin( atan(1.0)/45.0 * (ra2  - ra1 )/2.0 );
  const double ddec = sin( atan(1.0)/45.0 * (dec2 - dec1)/2.0 );

  dec1 *= atan(1.0)/45.0;
  dec2 *= atan(1.0)/45.0;

  return 2.0 * asin(sqrt( ddec*ddec + cos(dec1)*cos(dec2)*dra*dra))
      * 45.0/atan(1.0);
}
