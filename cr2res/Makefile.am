## Process this file with automake to produce Makefile.in

##   This file is part of the CR2RE Pipeline
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

SUBDIRS = . tests

DISTCLEANFILES = *~

if MAINTAINER_MODE

MAINTAINERCLEANFILES = Makefile.in

endif

AM_CPPFLAGS = -DCX_LOG_DOMAIN=\"Cr2resLib\" $(all_includes)

noinst_HEADERS =    cr2res_bpm.h \
                    cr2res_calib.h \
                    cr2res_detlin.h \
                    cr2res_dfs.h \
                    cr2res_extract.h \
                    cr2res_etalon.h \
                    cr2res_flat.h \
                    cr2res_idp.h \
                    cr2res_io.h \
                    cr2res_io_eop_table.h \
                    cr2res_nodding.h \
                    cr2res_pfits.h \
                    cr2res_photom.h \
                    cr2res_qc.h \
                    cr2res_qc_lines.h \
                    cr2res_slit_curv.h \
                    cr2res_splice.h \
                    cr2res_trace.h \
                    cr2res_utils.h \
                    cr2res_wave.h \
                    cr2res_pol.h

pkginclude_HEADERS =

privatelib_LTLIBRARIES = libcr2res.la

libcr2res_la_SOURCES =  cr2res_bpm.c \
                        cr2res_calib.c \
                        cr2res_detlin.c \
                        cr2res_dfs.c \
                        cr2res_extract.c \
                        cr2res_etalon.c \
                        cr2res_flat.c \
                        cr2res_idp.c \
                        cr2res_io.c \
                        cr2res_nodding.c \
                        cr2res_pfits.c \
                        cr2res_photom.c \
                        cr2res_qc.c \
                        cr2res_slit_curv.c \
                        cr2res_splice.c \
                        cr2res_trace.c \
                        cr2res_utils.c \
                        cr2res_wave.c \
                        cr2res_pol.c

libcr2res_la_LDFLAGS = $(CPL_LDFLAGS) $(HDRL_LDFLAGS) -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)
libcr2res_la_LIBADD = $(LIBIRPLIB) $(HDRL_LIBS) $(LIBCPLDFS) $(LIBCPLUI) $(LIBCPLDRS) $(LIBCPLCORE)
libcr2res_la_DEPENDENCIES = $(LIBIRPLIB)

