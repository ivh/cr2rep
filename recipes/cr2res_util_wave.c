
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

#include <locale.h>
#include <string.h>

#include <cpl.h>
#include <math.h>
#include "hdrl.h"

#include "cr2res_utils.h"
#include "cr2res_pfits.h"
#include "cr2res_dfs.h"
#include "cr2res_io.h"
#include "cr2res_extract.h"
#include "cr2res_trace.h"
#include "cr2res_wave.h"

/*-----------------------------------------------------------------------------
                                Define
 -----------------------------------------------------------------------------*/

#define RECIPE_STRING "cr2res_util_wave"

/*-----------------------------------------------------------------------------
                             Plugin registration
 -----------------------------------------------------------------------------*/

int cpl_plugin_get_info(cpl_pluginlist * list);

/*-----------------------------------------------------------------------------
                            Private function prototypes
 -----------------------------------------------------------------------------*/

static int cr2res_util_wave_create(cpl_plugin *);
static int cr2res_util_wave_exec(cpl_plugin *);
static int cr2res_util_wave_destroy(cpl_plugin *);
static int cr2res_util_wave(cpl_frameset *, const cpl_parameterlist *);

/*-----------------------------------------------------------------------------
                            Static variables
 -----------------------------------------------------------------------------*/

static char cr2res_util_wave_description[] =
"The utility expects 3 files as input:\n"
"   * extracted.fits " CR2RES_COMMAND_LINE "\n"
"   * trace_wave.fits " CR2RES_COMMAND_LINE "\n"
"   * static_calib.fits (optional) " CR2RES_COMMAND_LINE "\n"
"Four different methods are offered by the utility. They are controlled\n"
"by the --method parameter:\n"
"   XCORR:  Cross Correlation with a emission lines catalog (default)\n"
"   LINE1D: Line identification and fitting for each 1D spectra\n"
"   LINE2D: Line identification and fitting for all 1D spectra at once\n"
"   ETALON: Does not require any static calibration file.\n"
"The recipe produces the following products:\n"
"   * TRACE_WAVE\n"
"   * WAVE_MAP\n"
"\n";

/*-----------------------------------------------------------------------------
                                Function code
 -----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/**
  @defgroup cr2res_util_wave    Wavelength Calibration
 */
/*----------------------------------------------------------------------------*/

/**@{*/

/*----------------------------------------------------------------------------*/
/**
  @brief    Build the list of available plugins, for this module.
  @param    list    the plugin list
  @return   0 if everything is ok, 1 otherwise
  @note     Only this function is exported

  Create the recipe instance and make it available to the application using the
  interface.
 */
/*----------------------------------------------------------------------------*/
int cpl_plugin_get_info(cpl_pluginlist * list)
{
    cpl_recipe  *   recipe = cpl_calloc(1, sizeof *recipe );
    cpl_plugin  *   plugin = &recipe->interface;

    if (cpl_plugin_init(plugin,
                    CPL_PLUGIN_API,
                    CR2RES_BINARY_VERSION,
                    CPL_PLUGIN_TYPE_RECIPE,
                    "cr2res_util_wave",
                    "Wavelength Calibration",
                    cr2res_util_wave_description,
                    "Thomas Marquart, Yves Jung",
                    PACKAGE_BUGREPORT,
                    cr2res_get_license(),
                    cr2res_util_wave_create,
                    cr2res_util_wave_exec,
                    cr2res_util_wave_destroy)) {
        cpl_msg_error(cpl_func, "Plugin initialization failed");
        (void)cpl_error_set_where(cpl_func);
        return 1;
    }

    if (cpl_pluginlist_append(list, plugin)) {
        cpl_msg_error(cpl_func, "Error adding plugin to list");
        (void)cpl_error_set_where(cpl_func);
        return 1;
    }

    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Setup the recipe options
  @param    plugin  the plugin
  @return   0 if everything is ok

  Defining the command-line/configuration parameters for the recipe.
 */
/*----------------------------------------------------------------------------*/
static int cr2res_util_wave_create(cpl_plugin * plugin)
{
    cpl_recipe    * recipe;
    cpl_parameter * p;

    /* Check that the plugin is part of a valid recipe */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else
        return -1;

    /* Create the parameters list in the cpl_recipe object */
    recipe->parameters = cpl_parameterlist_new();

    /* Fill the parameters list */
    p = cpl_parameter_new_value("cr2res.cr2res_util_wave.detector",
            CPL_TYPE_INT, "Only reduce the specified detector",
            "cr2res.cr2res_util_wave", 0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "detector");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_wave.order",
            CPL_TYPE_INT, "Only reduce the specified order",
            "cr2res.cr2res_util_wave", -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "order");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_wave.trace_nb",
            CPL_TYPE_INT, "Only reduce the specified trace number",
            "cr2res.cr2res_util_wave", -1);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "trace_nb");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_wave.method",
            CPL_TYPE_STRING, "Data Type (XCORR / LINE1D / LINE2D / ETALON)",
            "cr2res.cr2res_util_wave", "XCORR");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "method");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_wave.wl_shift",
            CPL_TYPE_DOUBLE, "Wavelength shift (nm) to apply to the guess",
            "cr2res.cr2res_util_wave", 0.0);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "wl_shift");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_wave.wl_est",
            CPL_TYPE_STRING, "Estimated wavelength start and end",
            "cr2res.cr2res_util_wave", "-1.0, -1.0");
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "wl_est");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_wave.degree",
            CPL_TYPE_INT, "Wavelength Polynomial degree",
            "cr2res.cr2res_util_wave", 3);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "degree");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);

    p = cpl_parameter_new_value("cr2res.cr2res_util_wave.display",
            CPL_TYPE_BOOL, "Flag for display",
            "cr2res.cr2res_util_wave", FALSE);
    cpl_parameter_set_alias(p, CPL_PARAMETER_MODE_CLI, "display");
    cpl_parameter_disable(p, CPL_PARAMETER_MODE_ENV);
    cpl_parameterlist_append(recipe->parameters, p);
    return 0;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Execute the plugin instance given by the interface
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_util_wave_exec(cpl_plugin * plugin)
{
    cpl_recipe  *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1;

    return cr2res_util_wave(recipe->frames, recipe->parameters);
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Destroy what has been created by the 'create' function
  @param    plugin  the plugin
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_util_wave_destroy(cpl_plugin * plugin)
{
    cpl_recipe *recipe;

    /* Get the recipe out of the plugin */
    if (cpl_plugin_get_type(plugin) == CPL_PLUGIN_TYPE_RECIPE)
        recipe = (cpl_recipe *)plugin;
    else return -1 ;

    cpl_parameterlist_delete(recipe->parameters);
    return 0 ;
}

/*----------------------------------------------------------------------------*/
/**
  @brief    Interpret the command line options and execute the data processing
  @param    frameset   the frames list
  @param    parlist    the parameters list
  @return   0 if everything is ok
 */
/*----------------------------------------------------------------------------*/
static int cr2res_util_wave(
        cpl_frameset            *   frameset,
        const cpl_parameterlist *   parlist)
{
    const cpl_parameter *   param;
    int                     reduce_det, reduce_order, reduce_trace,
                            degree, display ;
    cpl_frame           *   fr ;
    const char          *   sval ;
    cr2res_wavecal_type     wavecal_type ;
    const char          *   trace_wave_file ;
    const char          *   extracted_file ;
    const char          *   static_calib_file ;
    char                *   out_file;
    cpl_array           *   wl_array ;
    cpl_array           *   wl_err_array ;
    cpl_table           *   out_trace_wave[CR2RES_NB_DETECTORS] ;
    hdrl_image          *   out_wave_map[CR2RES_NB_DETECTORS] ;
    cpl_propertylist    *   ext_plist[CR2RES_NB_DETECTORS] ;
    cpl_table           *   trace_wave_table ;
    cpl_table           *   extracted_table ;
    cpl_bivector        *   extracted_bivec ;
    cpl_bivector        *   extracted_err_bivec ;

    // For LINE2D
    cpl_bivector        **  spectra_2D;
    cpl_bivector        **  spectra_err_2D;
    cpl_polynomial      **  wavesol_init_2D;
    const cpl_array           **  wavesol_error_2D;
    int                  *  orders;
    int                     nb_orders;
    cpl_size                degree_2d[2];
    cpl_size                zero = 0;
    cpl_polynomial      *   tmp_sol, * collapsed_sol;

    cpl_polynomial      *   wave_guess ;
    const cpl_array     *   wave_error_guess ;
    cpl_polynomial      *   wave_sol ;
    double                  wstart, wend, wl_shift ;
    int                     det_nr, nb_traces, trace_id, order, i ;

    /* Needed for sscanf() */
    setlocale(LC_NUMERIC, "C");

    /* Initialise */
    wavecal_type = CR2RES_UNSPECIFIED ;
    wstart = wend = -1.0 ;
    wl_shift = 0.0 ;
    wl_err_array = NULL ;

    /* RETRIEVE INPUT PARAMETERS */
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_wave.degree");
    degree = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_wave.detector");
    reduce_det = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_wave.order");
    reduce_order = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_wave.trace_nb");
    reduce_trace = cpl_parameter_get_int(param);
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_wave.method");
    sval = cpl_parameter_get_string(param) ;
    if (!strcmp(sval, "XCORR"))         wavecal_type = CR2RES_XCORR ;
    else if (!strcmp(sval, "LINE1D"))   wavecal_type = CR2RES_LINE1D ;
    /* TODO : add support for 2D */
    else if (!strcmp(sval, "LINE2D"))   wavecal_type = CR2RES_LINE2D ;
    else if (!strcmp(sval, "ETALON"))   wavecal_type = CR2RES_ETALON ;
    else {
        cpl_msg_error(__func__, "Invalid Data Type specified");
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1;
    }
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_wave.wl_est");
    sval = cpl_parameter_get_string(param) ;
     if (sscanf(sval, "%lg,%lg", &wstart, &wend) != 2) {
        return -1 ;
    }
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_wave.wl_shift");
    wl_shift = cpl_parameter_get_double(param) ;
    param = cpl_parameterlist_find_const(parlist,
            "cr2res.cr2res_util_wave.display");
    display = cpl_parameter_get_bool(param) ;

    /* Check Parameters */
    if (degree < 0) {
        cpl_msg_error(__func__, "The degree needs to be >= 0");
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Identify the RAW and CALIB frames in the input frameset */
    if (cr2res_dfs_set_groups(frameset) != CPL_ERROR_NONE) {
        cpl_msg_error(__func__, "Cannot identify RAW and CALIB frames") ;
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Get Inputs */
    fr = cpl_frameset_get_position(frameset, 0);
    extracted_file = cpl_frame_get_filename(fr) ;
    fr = cpl_frameset_get_position(frameset, 1);
    trace_wave_file = cpl_frame_get_filename(fr) ;
    if (cpl_frameset_get_size(frameset) > 2) {
        fr = cpl_frameset_get_position(frameset, 2);
        static_calib_file = cpl_frame_get_filename(fr) ;
    } else {
        static_calib_file = NULL ;
    }
    if (trace_wave_file==NULL || extracted_file==NULL) {
        cpl_msg_error(__func__, "The utility needs at least 2 files as input");
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Get Data Type */
    if (wavecal_type == CR2RES_UNSPECIFIED) {
        /* Get the wavecal_type from the input extracted spectrum possible */
        cpl_msg_error(__func__, "Please use the --data_type option") ;
        return -1 ;
    }
    if ((wavecal_type == CR2RES_XCORR || wavecal_type == CR2RES_LINE1D || 
                wavecal_type == CR2RES_LINE2D) && static_calib_file == NULL) {
        cpl_msg_error(__func__,
                "The catalog file is needed for XCORR/LINE1D/LINE2D");
        cpl_error_set(__func__, CPL_ERROR_ILLEGAL_INPUT) ;
        return -1 ;
    }

    /* Loop over the detectors */
    for (det_nr=1 ; det_nr<=CR2RES_NB_DETECTORS ; det_nr++) {

        /* Initialise */
        out_trace_wave[det_nr-1] = NULL ;
        out_wave_map[det_nr-1] = NULL ;
        ext_plist[det_nr-1] = NULL ;

        /* Store the extenѕion header for product saving */
        ext_plist[det_nr-1] = cpl_propertylist_load(extracted_file,
                cr2res_io_get_ext_idx(extracted_file, det_nr, 1)) ;

        /* Compute only one detector */
        if (reduce_det != 0 && det_nr != reduce_det) continue ;

        cpl_msg_info(__func__, "Process detector number %d", det_nr) ;
        cpl_msg_indent_more() ;

        /* Load the TRACE_WAVE table of this detector */
        cpl_msg_info(__func__, "Load the TRACE_WAVE table") ;
        if ((trace_wave_table = cr2res_io_load_TRACE_WAVE(trace_wave_file,
                        det_nr)) == NULL) {
            cpl_msg_error(__func__,"Failed to load table - skip detector");
            cpl_error_reset() ;
            cpl_msg_indent_less() ;
            continue ;
        }
        nb_traces = cpl_table_get_nrow(trace_wave_table) ;

        /* Load the EXTRACT1D table of this detector */
        cpl_msg_info(__func__, "Load the EXTRACT1D table") ;
        if ((extracted_table = cr2res_io_load_EXTRACT_1D(extracted_file,
                        det_nr)) == NULL) {
            cpl_msg_error(__func__,"Failed to load table - skip detector");
            cpl_error_reset() ;
            cpl_msg_indent_less() ;
            continue ;
        }

        /* Create Output table for this detector */
        out_trace_wave[det_nr-1] = cpl_table_duplicate(trace_wave_table) ;

        /* Clear the Wavelength column */
        cpl_table_erase_column(out_trace_wave[det_nr-1], CR2RES_COL_WAVELENGTH);
        cpl_table_new_column_array(out_trace_wave[det_nr-1],
                CR2RES_COL_WAVELENGTH, CPL_TYPE_DOUBLE, degree+1) ;

        /* Clear the Wavelength Error column */
        cpl_table_erase_column(out_trace_wave[det_nr-1], 
                CR2RES_COL_WAVELENGTH_ERROR) ;
        cpl_table_new_column_array(out_trace_wave[det_nr-1],
                CR2RES_COL_WAVELENGTH_ERROR, CPL_TYPE_DOUBLE, 2) ;


        if (wavecal_type == CR2RES_LINE2D){
            spectra_2D = cpl_malloc(nb_traces * sizeof(cpl_bivector *));
            spectra_err_2D = cpl_malloc(nb_traces * sizeof(cpl_bivector *));
            wavesol_init_2D = cpl_malloc(nb_traces * sizeof(cpl_polynomial *));
            wavesol_error_2D = cpl_malloc(nb_traces * sizeof(cpl_array *));
            orders = cr2res_trace_get_order_numbers(trace_wave_table, &nb_orders);
        }

        /* Loop over the traces spectra */
        for (i=0 ; i<nb_traces ; i++) {
            /* Initialise */

            /* Get Order and trace id */
            order = cpl_table_get(trace_wave_table, CR2RES_COL_ORDER, i, NULL) ;
            trace_id = cpl_table_get(trace_wave_table, CR2RES_COL_TRACENB, i,
                    NULL) ;

            /* Check if this order needs to be skipped */
            if (reduce_order > -1 && order != reduce_order) {
                continue ;
            }

            /* Check if this trace needs to be skipped */
            if (reduce_trace > -1 && trace_id != reduce_trace) {
                continue ;
            }

            cpl_msg_info(__func__, "Process Order %d/Trace %d",order,trace_id) ;
            cpl_msg_indent_more() ;

            /* Get the extracted spectrum */
            if (cr2res_extract_EXTRACT1D_get_spectrum(
                    extracted_table, order, trace_id,
                    &extracted_bivec, &extracted_err_bivec)) {
                cpl_msg_error(__func__, "Cannot get the extracted spectrum") ;
                cpl_msg_indent_less() ;
                continue ;
            }

            /* Get the wavelength guess */
            if (wstart > 0.0 && wend > 0.0) {
                wave_guess = cr2res_wlestimate_compute(wstart, wend) ;
                wave_error_guess = NULL;
            } else {
                if ((wave_guess=cr2res_get_trace_wave_poly(trace_wave_table,
                            CR2RES_COL_WAVELENGTH, order, trace_id)) == NULL) {
                    cpl_msg_error(__func__, "Cannot get the WL guess") ;
                    cpl_msg_indent_less();
                    cpl_bivector_delete(extracted_bivec);
                    cpl_bivector_delete(extracted_err_bivec);
                    continue ;
                }
                if ((wave_error_guess = cpl_table_get_array(trace_wave_table, 
                        CR2RES_COL_WAVELENGTH_ERROR, 
                        cr2res_get_trace_table_index(trace_wave_table, 
                            order, trace_id))) == NULL) {
                    cpl_msg_error(__func__, "Cannot get the WL ERROR guess") ;
                    cpl_msg_indent_less() ;
                    cpl_bivector_delete(extracted_bivec) ;
                    cpl_bivector_delete(extracted_err_bivec) ;
                    cpl_polynomial_delete(wave_guess) ;
                    continue ;
                }
            }

            /* Apply the shift */
            if (fabs(wl_shift) > 1e-3) {
                const cpl_size power = 0;
                cpl_polynomial_set_coeff(wave_guess, &power,
                        cpl_polynomial_get_coeff(wave_guess, &power)+wl_shift) ;
            }

            if (wavecal_type == CR2RES_LINE2D){
                spectra_2D[i] = extracted_bivec;
                spectra_err_2D[i] = extracted_err_bivec;
                wavesol_init_2D[i] = wave_guess;
                wavesol_error_2D[i] = wave_error_guess;
            }
            else{
                /* Call the Wavelength Calibration */
                if ((wave_sol = cr2res_wave_1d(extracted_bivec,
                                extracted_err_bivec, wave_guess, 
                                wave_error_guess, wavecal_type, static_calib_file, 
                                degree, display, &wl_err_array)) == NULL) {
                    cpl_msg_error(__func__, "Cannot calibrate in Wavelength") ;
                    cpl_polynomial_delete(wave_guess) ;
                    cpl_bivector_delete(extracted_bivec) ;
                    cpl_bivector_delete(extracted_err_bivec) ;
                    cpl_error_reset() ;
                    cpl_msg_indent_less() ;
                    continue ;
                }
                cpl_bivector_delete(extracted_bivec) ;
                cpl_bivector_delete(extracted_err_bivec) ;
                cpl_polynomial_delete(wave_guess) ;

                /* Store the Solution in the table */
                wl_array = cr2res_convert_poly_to_array(wave_sol, degree+1) ;
                cpl_polynomial_delete(wave_sol);
                if (wl_array != NULL) {
                    cpl_table_set_array(out_trace_wave[det_nr-1],
                            CR2RES_COL_WAVELENGTH, i, wl_array);
                    cpl_array_delete(wl_array) ;
                }
                if (wl_err_array != NULL) {
                    cpl_table_set_array(out_trace_wave[det_nr-1],
                            CR2RES_COL_WAVELENGTH_ERROR, i, wl_err_array);
                    cpl_array_delete(wl_err_array) ;
                }
            }
            cpl_msg_indent_less() ;
        }

        if (wavecal_type == CR2RES_LINE2D){
            degree_2d[0] = degree;
            degree_2d[1] = degree;
            if ((wave_sol = cr2res_wave_2d(spectra_2D, spectra_err_2D, wavesol_init_2D, wavesol_error_2D, wavecal_type,
                    static_calib_file, orders, nb_orders, degree_2d, display, &wl_err_array)) == NULL){
                        cpl_msg_error(__func__, "Cannot calibrate in Wavelength") ;
                        for(cpl_size i = 0; i < nb_traces; i++)
                        {
                            cpl_bivector_delete(spectra_2D[i]);
                            cpl_bivector_delete(spectra_err_2D[i]);
                            cpl_polynomial_delete(wavesol_init_2D[i]);
                        }
                        cpl_free(spectra_2D);
                        cpl_free(spectra_err_2D);
                        cpl_free(wavesol_init_2D);
                        cpl_free(wavesol_error_2D);
                        cpl_free(orders);
                        cpl_table_delete(trace_wave_table) ;
                        cpl_table_delete(extracted_table) ;
                        continue;
                    }

            /* Store the Solution in the table */
            tmp_sol = cpl_polynomial_new(1);

            for(cpl_size i = 0; i < nb_traces; i++)
            {
                cpl_polynomial_set_coeff(tmp_sol, &zero, i);
                collapsed_sol = cpl_polynomial_extract(wave_sol, 1, tmp_sol);
                wl_array = cr2res_convert_poly_to_array(collapsed_sol, degree+1) ;
                cpl_polynomial_delete(collapsed_sol);
                if (wl_array != NULL) {
                    cpl_table_set_array(out_trace_wave[det_nr-1],
                            CR2RES_COL_WAVELENGTH, i, wl_array);
                    cpl_array_delete(wl_array) ;
                }
                if (wl_err_array != NULL) {
                    cpl_table_set_array(out_trace_wave[det_nr-1],
                            CR2RES_COL_WAVELENGTH_ERROR, i, wl_err_array);
                    cpl_array_delete(wl_err_array) ;
                }
            }
            // Deallocate 2D structures
            cpl_polynomial_delete(tmp_sol);
            cpl_polynomial_delete(wave_sol);
            
            for(cpl_size i = 0; i < nb_traces; i++)
            {
                cpl_bivector_delete(spectra_2D[i]);
                cpl_bivector_delete(spectra_err_2D[i]);
                cpl_polynomial_delete(wavesol_init_2D[i]);
            }
            cpl_free(spectra_2D);
            cpl_free(spectra_err_2D);
            cpl_free(wavesol_init_2D);
            cpl_free(wavesol_error_2D);
            cpl_free(orders);
        }

        cpl_table_delete(trace_wave_table) ;
        cpl_table_delete(extracted_table) ;

        /* Generate the Wave Map */
        out_wave_map[det_nr-1] =
            cr2res_wave_gen_wave_map(out_trace_wave[det_nr-1]) ;

        /* Deallocate */
        cpl_msg_indent_less() ;
    }

    /* Save the new trace_wave table */
    out_file = cpl_sprintf("%s_wave.fits",
            cr2res_get_base_name(cr2res_get_root_name(extracted_file)));
    cr2res_io_save_TRACE_WAVE(out_file, frameset, frameset, parlist, 
            out_trace_wave, NULL, ext_plist, 
            CR2RES_UTIL_WAVE_TRACE_WAVE_PROCATG, RECIPE_STRING) ;
    cpl_free(out_file);

    /* Save the Wave Map */
    out_file = cpl_sprintf("%s_wave_map.fits",
            cr2res_get_base_name(cr2res_get_root_name(extracted_file)));
    cr2res_io_save_WAVE_MAP(out_file, frameset, frameset, parlist, out_wave_map,
            NULL, ext_plist, CR2RES_UTIL_WAVE_MAP_PROCATG, RECIPE_STRING) ;
    cpl_free(out_file);

    /* Free and return */
    for (i=0 ; i<CR2RES_NB_DETECTORS ; i++) {
        if (ext_plist[i] != NULL)
            cpl_propertylist_delete(ext_plist[i]) ;
        if (out_trace_wave[i] != NULL)
            cpl_table_delete(out_trace_wave[i]) ;
        if (out_wave_map[i] != NULL) {
            hdrl_image_delete(out_wave_map[i]) ;
        }
    }
    return (int)cpl_error_get_code();
}

/**@}*/
