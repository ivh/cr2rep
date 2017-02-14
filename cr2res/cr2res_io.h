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

cpl_image * cr2res_io_load_MASTER_DARK(
        const char  *   filename,
        int             detector,
        int             data);

cpl_image * cr2res_io_load_MASTER_BPM(
        const char  *   filename,
        int             detector);

cpl_imagelist * cr2res_io_load_DETLIN_COEFFS(
        const char  *   filename,
        int             detector);

cpl_image * cr2res_io_load_MASTER_FLAT(
        const char  *   filename,
        int             detector,
        int             data);

cpl_table * cr2res_io_load_TRACE_OPEN(
        const char  *   filename,
        int             detector);

cpl_table * cr2res_io_load_TRACE_DECKER(
        const char  *   filename,
        int             detector,
        cr2res_decker * decker_type);

cpl_table * cr2res_io_load_BLAZE(
        const char  *   filename,
        int             detector);

cpl_image * cr2res_io_load_BLAZE_IMAGE(
        const char  *   filename,
        int             detector);

cpl_image * cr2res_io_load_SLIT_MODEL(
        const char  *   filename,
        int             detector,
        int             data);

cpl_table * cr2res_io_load_SLIT_FUNC(
        const char  *   filename,
        int             detector) ;

cpl_image * cr2res_io_load_WAVE_MAP(
        const char  *   filename,
        int             detector);

cpl_table * cr2res_io_load_WAVE_SUB_ORDER(
        const char  *   filename,
        int             detector);

cpl_image * cr2res_io_load_SLITPOS_MAP(
        const char  *   filename,
        int             detector);

cpl_image * cr2res_io_load_TILT_MAP(
        const char  *   filename,
        int             detector);

cpl_table * cr2res_io_load_TILT_POLY(
        const char  *   filename,
        int             detector);

cpl_table * cr2res_io_load_EXTRACT_1D(
        const char  *   filename,
        int             detector);

cpl_table * cr2res_io_load_SPLICED_1D(
        const char  *   filename,
        int             detector);

cpl_table * cr2res_io_load_EXTRACT_2D(
        const char  *   filename);

cpl_table * cr2res_io_load_EXTRACT_POL(
        const char  *   filename,
        int             detector);

int cr2res_io_save_MASTER_DARK(
        cpl_frameset            *   allframes,
        const char              *   filename,
        cpl_frameset            *   used_frames,
        const cpl_parameterlist *   parlist,
        hdrl_image              **  master_darks,
        const cpl_propertylist  *   qc_list,
        const char              *   procatg,
        const char              *   recipe) ;

int cr2res_io_save_MASTER_BPM(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_imagelist           *   data,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id);

int cr2res_io_save_DETLIN_COEFFS(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_imagelist           **  data,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id);

int cr2res_io_save_MASTER_FLAT(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_imagelist           *   data,
        cpl_imagelist           *   errors,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id);

int cr2res_io_save_TRACE_OPEN(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id);

int cr2res_io_save_TRACE_DECKER(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id);

int cr2res_io_save_BLAZE(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id);

int cr2res_io_save_BLAZE_IMAGE(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_imagelist           *   data,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id);

int cr2res_io_save_SLIT_MODEL(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_imagelist           *   data,
        cpl_imagelist           *   errors,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id);

int cr2res_io_save_SLIT_FUNC(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id);

int cr2res_io_save_WAVE_MAP(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_imagelist           *   data,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id);

int cr2res_io_save_WAVE_SUB_ORDER(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id);

int cr2res_io_save_SLITPOS_MAP(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_imagelist           *   data,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id);

int cr2res_io_save_TILT_MAP(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_imagelist           *   data,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id);

int cr2res_io_save_TILT_POLY(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id);

int cr2res_io_save_EXTRACT_1D(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id);

int cr2res_io_save_SPLICED_1D(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id);

int cr2res_io_save_EXTRACT_2D(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               *   table,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id);

int cr2res_io_save_EXTRACT_POL(
        const char              *   filename,
        cpl_frameset            *   allframes,
        const cpl_parameterlist *   parlist,
        cpl_table               **  tables,
        const cpl_propertylist  *   qc_list,
        const char              *   recipe,
        const char              *   pipe_id);

#endif

