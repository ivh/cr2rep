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


AM_CPPFLAGS = -DCX_LOG_DOMAIN=\"Cr2reLib\" $(all_includes)

noinst_HEADERS =    cr2re_utils.h \
                    cr2re_cluster.h \
                    cr2re_etalon.h \
                    cr2re_pfits.h \
                    cr2re_dfs.h

pkginclude_HEADERS =

privatelib_LTLIBRARIES = libcr2re.la

libcr2re_la_SOURCES =    cr2re_utils.c \
                                cr2re_pfits.c \
                                cr2re_cluster.c \
                                cr2re_etalon.c \
                                cr2re_dfs.c

libcr2re_la_LDFLAGS = $(CPL_LDFLAGS) $(XXCLIPM_LDFLAGS) -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)
libcr2re_la_LIBADD = $(LIBCPLDFS) $(LIBCPLUI) $(LIBCPLDRS) $(LIBCPLCORE) $(LIBXXCLIPM)
libcr2re_la_DEPENDENCIES =
