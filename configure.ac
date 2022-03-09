# Process this file with autoconf to produce a configure script.

AC_INIT([CR2RES Instrument Pipeline], [1.1.2], [https://support.eso.org/], [cr2re])
AC_PREREQ([2.59])

AC_CONFIG_SRCDIR([Makefile.am])
AC_CONFIG_AUX_DIR([admin])

AC_CANONICAL_HOST
AC_CANONICAL_TARGET

AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config.h])

AM_MAINTAINER_MODE

CR2RES_SET_PREFIX(${PIPE_HOME:-/usr/local})

# Immediately before every release do:
#-------------------------------------
#   if (the interface is totally unchanged from previous release)
#       REVISION++;
#   else {
#       /* interfaces have been added, removed or changed */
#       REVISION = 0;
#       CURRENT++;
#       if (any interfaces have been _added_ since last release)
#           AGE++;
#       if (any interfaces have been _removed_ or incompatibly changed)
#           AGE = 0;
#   }

CR2RES_SET_VERSION_INFO([$VERSION])

# Checks for programs.
AC_PROG_CC
AC_PROG_CPP
AM_PROG_CC_C_O
AC_PROG_LN_S

ESO_ENABLE_DEBUG(no)
ESO_ENABLE_STRICT(no)
# probably unecessary, was required for some ancient HP/UX machines
ESO_PROG_CC_FLAG([fno-builtin], [CFLAGS=" -fno-builtin $CFLAGS"])
ESO_PROG_CC_FLAG([std=c99], [CFLAGS="$CFLAGS -std=c99"])

ESO_CHECK_DOCTOOLS

AC_ENABLE_STATIC(no)
AC_ENABLE_SHARED(yes)

AC_PROG_LIBTOOL
AC_SUBST(LIBTOOL_DEPS)

# Checks for libraries.
AC_CHECK_LIB(m, pow, [LIBS="$LIBS -lm"])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE

# Checks for library functions.
AC_CHECK_FUNCS([floor pow sqrt isinf isnan])
ESO_FUNC_STRDUP

# Check for CPL presence and usability
CPL_CHECK_LIBS

# setup hdrl folder
HDRL_CHECK([hdrl])

# Check for <xx>clipm
# - can be called multiple times
# - replace the argument "xxclipm" by your library name,
#   e.g. with KMOS: XXCLIPM_CHECK_LIB(kmclipm)
# Comment in if required!
# XXCLIPM_CHECK_LIB(xxclipm)

# call after CPL_CHECK_LIBS and XXCLIPM_CHECK_LIB
CR2RES_SET_PATHS 
CR2RES_CREATE_SYMBOLS

AM_WITH_DMALLOC

AC_CONFIG_SUBDIRS([hdrl])

AC_CONFIG_FILES(Makefile
                doxygen/Doxyfile
                cr2res/Makefile
                cr2res/tests/Makefile
                tools/Makefile
                reflex/Makefile
                reflex/cr2res.xml
                recipes/Makefile
                regtests/Makefile
                regtests/tests/Makefile
                irplib/Makefile
                irplib/tests/Makefile)
AC_OUTPUT
