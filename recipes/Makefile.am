## Process this file with automake to produce Makefile.in

##   This file is part of the CR2RES Pipeline
##   Copyright (C) 2002,2003 European Southern Observatory

##   This library is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.

##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.

##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02111-1307  USA

AUTOMAKE_OPTIONS = 1.8 foreign

DISTCLEANFILES = *~


if MAINTAINER_MODE

MAINTAINERCLEANFILES = Makefile.in

endif


AM_CPPFLAGS = $(all_includes)

noinst_HEADERS =

plugin_LTLIBRARIES = cr2res_cal_detlin.la \
                     cr2res_cal_dark.la \
                     cr2res_cal_flat.la \
                     cr2res_cal_wave.la \
                     cr2res_obs_nodding.la \
                     cr2res_obs_staring.la \
                     cr2res_obs_2d.la \
                     cr2res_obs_pol.la \
                     cr2res_util_bpm_merge.la \
                     cr2res_util_bpm_split.la \
                     cr2res_util_calib.la \
                     cr2res_util_extract.la \
                     cr2res_util_genlines.la \
                     cr2res_util_genstd.la \
                     cr2res_util_plot.la \
                     cr2res_util_normflat.la \
                     cr2res_util_slit_curv.la \
                     cr2res_util_splice.la \
                     cr2res_util_trace.la \
                     cr2res_util_trace_map.la \
                     cr2res_util_wave.la 

cr2res_cal_detlin_la_SOURCES = cr2res_cal_detlin.c
cr2res_cal_detlin_la_LIBADD = $(HDRL_LIBS) $(LIBCR2RES)
cr2res_cal_detlin_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version -no-undefined
cr2res_cal_detlin_la_DEPENDENCIES = $(LIBCR2RES)

cr2res_cal_dark_la_SOURCES = cr2res_cal_dark.c
cr2res_cal_dark_la_LIBADD = $(HDRL_LIBS) $(LIBCR2RES)
cr2res_cal_dark_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version -no-undefined
cr2res_cal_dark_la_DEPENDENCIES = $(LIBCR2RES)

cr2res_cal_flat_la_SOURCES = cr2res_cal_flat.c
cr2res_cal_flat_la_LIBADD = $(HDRL_LIBS) $(LIBCR2RES)
cr2res_cal_flat_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version -no-undefined
cr2res_cal_flat_la_DEPENDENCIES = $(LIBCR2RES)

cr2res_cal_wave_la_SOURCES = cr2res_cal_wave.c
cr2res_cal_wave_la_LIBADD = $(HDRL_LIBS) $(LIBCR2RES)
cr2res_cal_wave_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version -no-undefined
cr2res_cal_wave_la_DEPENDENCIES = $(LIBCR2RES)

cr2res_obs_nodding_la_SOURCES = cr2res_obs_nodding.c
cr2res_obs_nodding_la_LIBADD = $(HDRL_LIBS) $(LIBCR2RES)
cr2res_obs_nodding_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version -no-undefined
cr2res_obs_nodding_la_DEPENDENCIES = $(LIBCR2RES)

cr2res_obs_staring_la_SOURCES = cr2res_obs_staring.c
cr2res_obs_staring_la_LIBADD = $(HDRL_LIBS) $(LIBCR2RES)
cr2res_obs_staring_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version -no-undefined
cr2res_obs_staring_la_DEPENDENCIES = $(LIBCR2RES)

cr2res_obs_2d_la_SOURCES = cr2res_obs_2d.c
cr2res_obs_2d_la_LIBADD = $(HDRL_LIBS) $(LIBCR2RES)
cr2res_obs_2d_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version -no-undefined
cr2res_obs_2d_la_DEPENDENCIES = $(LIBCR2RES)

cr2res_obs_pol_la_SOURCES = cr2res_obs_pol.c
cr2res_obs_pol_la_LIBADD = $(HDRL_LIBS) $(LIBCR2RES)
cr2res_obs_pol_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version -no-undefined
cr2res_obs_pol_la_DEPENDENCIES = $(LIBCR2RES)

cr2res_util_bpm_merge_la_SOURCES = cr2res_util_bpm_merge.c
cr2res_util_bpm_merge_la_LIBADD = $(HDRL_LIBS) $(LIBCR2RES)
cr2res_util_bpm_merge_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version -no-undefined
cr2res_util_bpm_merge_la_DEPENDENCIES = $(LIBCR2RES)

cr2res_util_bpm_split_la_SOURCES = cr2res_util_bpm_split.c
cr2res_util_bpm_split_la_LIBADD = $(HDRL_LIBS) $(LIBCR2RES)
cr2res_util_bpm_split_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version -no-undefined
cr2res_util_bpm_split_la_DEPENDENCIES = $(LIBCR2RES)

cr2res_util_calib_la_SOURCES = cr2res_util_calib.c
cr2res_util_calib_la_LIBADD = $(HDRL_LIBS) $(LIBCR2RES)
cr2res_util_calib_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version -no-undefined
cr2res_util_calib_la_DEPENDENCIES = $(LIBCR2RES)

cr2res_util_extract_la_SOURCES = cr2res_util_extract.c
cr2res_util_extract_la_LIBADD = $(HDRL_LIBS) $(LIBCR2RES)
cr2res_util_extract_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version -no-undefined
cr2res_util_extract_la_DEPENDENCIES = $(LIBCR2RES)

cr2res_util_genlines_la_SOURCES = cr2res_util_genlines.c
cr2res_util_genlines_la_LIBADD = $(HDRL_LIBS) $(LIBCR2RES)
cr2res_util_genlines_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version -no-undefined
cr2res_util_genlines_la_DEPENDENCIES = $(LIBCR2RES)

cr2res_util_genstd_la_SOURCES = cr2res_util_genstd.c
cr2res_util_genstd_la_LIBADD = $(HDRL_LIBS) $(LIBCR2RES)
cr2res_util_genstd_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version -no-undefined
cr2res_util_genstd_la_DEPENDENCIES = $(LIBCR2RES)

cr2res_util_plot_la_SOURCES = cr2res_util_plot.c
cr2res_util_plot_la_LIBADD = $(HDRL_LIBS) $(LIBCR2RES)
cr2res_util_plot_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version -no-undefined
cr2res_util_plot_la_DEPENDENCIES = $(LIBCR2RES)

cr2res_util_normflat_la_SOURCES = cr2res_util_normflat.c
cr2res_util_normflat_la_LIBADD = $(HDRL_LIBS) $(LIBCR2RES)
cr2res_util_normflat_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version -no-undefined
cr2res_util_normflat_la_DEPENDENCIES = $(LIBCR2RES)

cr2res_util_slit_curv_la_SOURCES = cr2res_util_slit_curv.c
cr2res_util_slit_curv_la_LIBADD = $(HDRL_LIBS) $(LIBCR2RES)
cr2res_util_slit_curv_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version -no-undefined
cr2res_util_slit_curv_la_DEPENDENCIES = $(LIBCR2RES)

cr2res_util_splice_la_SOURCES = cr2res_util_splice.c
cr2res_util_splice_la_LIBADD = $(HDRL_LIBS) $(LIBCR2RES)
cr2res_util_splice_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version -no-undefined
cr2res_util_splice_la_DEPENDENCIES = $(LIBCR2RES)

cr2res_util_trace_la_SOURCES = cr2res_util_trace.c
cr2res_util_trace_la_LIBADD = $(HDRL_LIBS) $(LIBCR2RES)
cr2res_util_trace_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version -no-undefined
cr2res_util_trace_la_DEPENDENCIES = $(LIBCR2RES)

cr2res_util_trace_map_la_SOURCES = cr2res_util_trace_map.c
cr2res_util_trace_map_la_LIBADD = $(HDRL_LIBS) $(LIBCR2RES)
cr2res_util_trace_map_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version -no-undefined
cr2res_util_trace_map_la_DEPENDENCIES = $(LIBCR2RES)

cr2res_util_wave_la_SOURCES = cr2res_util_wave.c
cr2res_util_wave_la_LIBADD = $(HDRL_LIBS) $(LIBCR2RES)
cr2res_util_wave_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version -no-undefined
cr2res_util_wave_la_DEPENDENCIES = $(LIBCR2RES)

