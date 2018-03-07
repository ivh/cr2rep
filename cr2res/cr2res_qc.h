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

#ifndef CR2RES_QC_H
#define CR2RES_QC_H

/*-----------------------------------------------------------------------------
                                   Includes
 -----------------------------------------------------------------------------*/

#include <cpl.h>

/*-----------------------------------------------------------------------------
                           Functions prototypes
 -----------------------------------------------------------------------------*/

int cr2res_qc_count_badpix(cpl_image * bpm, int type);

double cr2res_qc_read_out_noise(cpl_image * im1, cpl_image * im2);
double cr2res_qc_dark_mean(cpl_image * dark);
double cr2res_qc_dark_median(cpl_image * dark);
double cr2res_qc_dark_stddev(cpl_image * );

int cr2res_qc_trace_count_orders(cpl_table * tracewave);
int cr2res_qc_trace_count_traces(cpl_table * tracewave);
int cr2res_qc_trace_get_ypos(cpl_table * tracewave, int order, int trace);

double cr2res_qc_wave_zeropoint(cpl_table * tracewave, int order, int trace);

#endif

