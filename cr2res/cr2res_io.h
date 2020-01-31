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

#ifndef CR2RES_IO_H
#define CR2RES_IO_H

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>
#include "hdrl.h"

#include "cr2res_utils.h"

/*-----------------------------------------------------------------------------
                                   Functions prototypes
 -----------------------------------------------------------------------------*/

cpl_frame * cr2res_io_find_SLIT_MODEL(
        const cpl_frameset  *   in,
        const char          *   setting_id,
        cr2res_decker           cr2res_decker) ;
const cpl_frame * cr2res_io_find_TRACE_WAVE(const cpl_frameset * in) ;
cpl_frameset * cr2res_io_find_TRACE_WAVE_all(const cpl_frameset * in) ;
cpl_frameset * cr2res_io_find_EXTRACT_1D_all(const cpl_frameset * in) ;
const cpl_frame * cr2res_io_find_BPM(const cpl_frameset * in) ;
cpl_frameset * cr2res_io_find_BPM_all(const cpl_frameset * in) ;

cpl_vector * cr2res_io_read_dits(const cpl_frameset * in) ;

cr2res_decker * cr2res_io_read_decker_positions(const cpl_frameset * in) ;

cpl_frameset * cr2res_io_extract_decker_frameset(
        const cpl_frameset  *   in,
        const char          *   tag,
        cr2res_decker           decker) ;

int cr2res_io_convert_order_to_idx(int order) ;
int cr2res_io_convert_idx_to_order(int order_idx) ;

char * cr2res_io_create_extname(
        int             detector,
        int             data) ;

int cr2res_io_get_ext_idx(
        const char  *   filename,
        int             detector,
        int             data) ;

hdrl_image * cr2res_io_load_image(
        const char  *   in,
        int             detector) ;

hdrl_imagelist * cr2res_io_load_image_list(
        const char  *   in,
        int             detector) ;

hdrl_imagelist * cr2res_io_load_image_list_from_set(
        const cpl_frameset  *   in,
        int                     detector) ;

cpl_table * cr2res_load_table(
        const char  *   in,
        int             det_nr,
        int             pmin,
        int             pmax) ;

cpl_bivector * cr2res_io_load_EMISSION_LINES(
        const char  *   filename) ;

cpl_image * cr2res_io_load_BPM(
        const char  *   filename,
        int             detector,
        int             data) ;

hdrl_image * cr2res_io_load_MASTER_DARK(
        const char  *   filename,
        int             detector) ;

hdrl_imagelist * cr2res_io_load_DETLIN_COEFFS(
        const char  *   filename,
        int             detector) ;

hdrl_image * cr2res_io_load_MASTER_FLAT(
        const char  *   filename,
        int             detector) ;

cpl_table * cr2res_io_load_TRACE_WAVE(
        const char  *   filename,
        int             detector);

hdrl_image * cr2res_io_load_SLIT_MODEL(
        const char  *   filename,
        int             detector) ;

hdrl_image * cr2res_io_load_TRACE_MAP(
        const char  *   filename,
        int             detector) ;

hdrl_image * cr2res_io_load_WAVE_MAP(
        const char  *   filename,
        int             detector) ;

hdrl_image * cr2res_io_load_SLIT_CURV_MAP(
        const char  *   filename,
        int             detector) ;

cpl_table * cr2res_io_load_SLIT_CURV(
        const char  *   filename,
        int             detector);

cpl_table * cr2res_io_load_EXTRACT_1D(
        const char  *   filename,
        int             detector);

cpl_table * cr2res_io_load_SPLICED_1D(
        const char  *   filename,
        int             detector);

cpl_table * cr2res_io_load_EXTRACT_2D(
        const char  *   filename,
        int             detector);

int cr2res_io_save_PHOTO_FLUX(
        const char              *   filename,
        cpl_table               *   out_table,
        const cpl_parameterlist *   parlist,
        cpl_frameset            *   set,
        const char              *   recipe) ;

int cr2res_io_save_EMISSION_LINES(
        const char              *   filename,
        cpl_table               *   out_table,
        const cpl_parameterlist *   parlist,
        cpl_frameset            *   set,
        const char              *   recipe,
        const char              *   setting_string) ;

int cr2res_io_save_MASTER_DARK(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  master_darks,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe) ;

int cr2res_io_save_DETLIN_COEFFS(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        hdrl_imagelist          **  coeffs,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe) ;

int cr2res_io_save_BPM(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        cpl_image               **  bpms,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe) ;

int cr2res_io_save_MASTER_FLAT(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  master_flats,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe) ;

int cr2res_io_save_CALIBRATED(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  calib_collapsed,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe) ;

int cr2res_io_save_TRACE_WAVE(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe) ;

int cr2res_io_save_LINES_DIAGNOSTICS(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe) ;

int cr2res_io_save_EXTRACT_1D(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe) ;

int cr2res_io_save_THROUGHPUT(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe) ;

int cr2res_io_save_SLIT_FUNC(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  slit_func,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe) ;

int cr2res_io_save_COMBINED(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  data,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe) ;

int cr2res_io_save_SLIT_MODEL(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  data,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe);

int cr2res_io_save_TRACE_MAP(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  data,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe) ;

int cr2res_io_save_WAVE_MAP(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  data,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe) ;

int cr2res_io_save_SLIT_CURV_MAP(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  data,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe) ;

int cr2res_io_save_SLIT_CURV(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe) ;

int cr2res_io_save_SPLICED_1D(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        cpl_table               *   spliced_1d,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        *   ext_plist,
        const char              *   procatg,
        const char              *   recipe) ;

int cr2res_io_save_EXTRACT_2D(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe) ;

int cr2res_io_save_POL_SPEC(
        const char              *   filename,
        cpl_frameset            *   allframes,
        cpl_frameset            *   inframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        cpl_propertylist        **  ext_plist,
        const char              *   procatg,
        const char              *   recipe) ;

#endif

