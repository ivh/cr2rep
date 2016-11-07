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

DISTCLEANFILES = *~


if MAINTAINER_MODE

MAINTAINERCLEANFILES = Makefile.in

endif


AM_CPPFLAGS = $(all_includes)

noinst_HEADERS = 

plugin_LTLIBRARIES = cr2res_bias.la \
                     cr2res_cal_dark.la \
                     cr2res_etalon.la \
                     cr2res_test.la \
                     cr2res_trace.la 


cr2res_cal_dark_la_SOURCES = cr2res_cal_dark.c 
cr2res_cal_dark_la_LIBADD = $(HDRL_LIBS) $(LIBCR2RE)
cr2res_cal_dark_la_LDFLAGS = $(HDRL_LDFLAGS) -module -avoid-version -no-undefined
cr2res_cal_dark_la_DEPENDENCIES = $(LIBCR2RE)

