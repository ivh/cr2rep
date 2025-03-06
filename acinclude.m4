# CR2RES_SET_PREFIX(PREFIX)
#---------------------------
AC_DEFUN([CR2RES_SET_PREFIX],
[
    unset CDPATH
    # make $PIPE_HOME the default for the installation
    AC_PREFIX_DEFAULT($1)

    if test "x$prefix" = "xNONE"; then
        prefix=$ac_default_prefix
        ac_configure_args="$ac_configure_args --prefix $prefix"
    fi

    if test "x$exec_prefix" = "xNONE"; then
        exec_prefix=$prefix
    fi

])


# CR2RES_SET_VERSION_INFO(VERSION, [CURRENT], [REVISION], [AGE])
#----------------------------------------------------------------
# Setup various version information, especially the libtool versioning
AC_DEFUN([CR2RES_SET_VERSION_INFO],
[
    cr2res_version=`echo "$1" | sed -e 's/[[a-z,A-Z]].*$//'`

    cr2res_major_version=`echo "$cr2res_version" | \
        sed 's/\([[0-9]]*\).\(.*\)/\1/'`
    cr2res_minor_version=`echo "$cr2res_version" | \
        sed 's/\([[0-9]]*\).\([[0-9]]*\)\(.*\)/\2/'`
    cr2res_micro_version=`echo "$cr2res_version" | \
        sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\3/'`

    if test -z "$cr2res_major_version"; then cr2res_major_version=0
    fi

    if test -z "$cr2res_minor_version"; then cr2res_minor_version=0
    fi

    if test -z "$cr2res_micro_version"; then cr2res_micro_version=0
    fi

    CR2RES_VERSION="$cr2res_version"
    CR2RES_MAJOR_VERSION=$cr2res_major_version
    CR2RES_MINOR_VERSION=$cr2res_minor_version
    CR2RES_MICRO_VERSION=$cr2res_micro_version

    if test -z "$4"; then CR2RES_INTERFACE_AGE=0
    else CR2RES_INTERFACE_AGE="$4"
    fi

    CR2RES_BINARY_AGE=`expr 100 '*' $CR2RES_MINOR_VERSION + $CR2RES_MICRO_VERSION`
    CR2RES_BINARY_VERSION=`expr 10000 '*' $CR2RES_MAJOR_VERSION + \
                          $CR2RES_BINARY_AGE`

    AC_SUBST(CR2RES_VERSION)
    AC_SUBST(CR2RES_MAJOR_VERSION)
    AC_SUBST(CR2RES_MINOR_VERSION)
    AC_SUBST(CR2RES_MICRO_VERSION)
    AC_SUBST(CR2RES_INTERFACE_AGE)
    AC_SUBST(CR2RES_BINARY_VERSION)
    AC_SUBST(CR2RES_BINARY_AGE)

    AC_DEFINE_UNQUOTED(CR2RES_MAJOR_VERSION, $CR2RES_MAJOR_VERSION,
                       [CR2RES major version number])
    AC_DEFINE_UNQUOTED(CR2RES_MINOR_VERSION, $CR2RES_MINOR_VERSION,
                       [CR2RES minor version number])
    AC_DEFINE_UNQUOTED(CR2RES_MICRO_VERSION, $CR2RES_MICRO_VERSION,
                       [CR2RES micro version number])
    AC_DEFINE_UNQUOTED(CR2RES_INTERFACE_AGE, $CR2RES_INTERFACE_AGE,
                       [CR2RES interface age])
    AC_DEFINE_UNQUOTED(CR2RES_BINARY_VERSION, $CR2RES_BINARY_VERSION,
                       [CR2RES binary version number])
    AC_DEFINE_UNQUOTED(CR2RES_BINARY_AGE, $CR2RES_BINARY_AGE,
                       [CR2RES binary age])

    ESO_SET_LIBRARY_VERSION([$2], [$3], [$4])
])


# CR2RES_SET_PATHS
#------------------
# Define auxiliary directories of the installed directory tree.
AC_DEFUN([CR2RES_SET_PATHS],
[

    if test -z "$plugindir"; then
        plugindir='${libdir}/esopipes-plugins/${PACKAGE}-${VERSION}'
    fi

    if test -z "$privatelibdir"; then
        privatelibdir='${libdir}/${PACKAGE}-${VERSION}'
    fi

    if test -z "$apidocdir"; then
        apidocdir='${datadir}/doc/esopipes/${PACKAGE}-${VERSION}/html'
    fi

    if test -z "$pipedocsdir"; then
        pipedocsdir='${datadir}/doc/esopipes/${PACKAGE}-${VERSION}'
    fi

    if test -z "$configdir"; then
        configdir='${datadir}/esopipes/${PACKAGE}-${VERSION}/config'
    fi

    if test -z "$wkfextradir"; then
        wkfextradir='${datadir}/esopipes/${PACKAGE}-${VERSION}/reflex'
    fi

    if test -z "$wkfcopydir"; then
        wkfcopydir='${datadir}/reflex/workflows/${PACKAGE}-${VERSION}'
    fi

    if test -z "$workflowdir"; then
        workflowdir='${datadir}/esopipes/workflows/${PACKAGE}-${VERSION}/crires'
    fi

    if test -z "$reportsdir"; then
       reportsdir='${datadir}/esopipes/reports/${PACKAGE}-${VERSION}'
    fi

    AC_SUBST(plugindir)
    AC_SUBST(privatelibdir)
    AC_SUBST(apidocdir)
    AC_SUBST(pipedocsdir)
    AC_SUBST(configdir)
    AC_SUBST(wkfextradir)
    AC_SUBST(wkfcopydir)
    AC_SUBST(workflowdir)
    AC_SUBST(reportsdir)



    # Define a preprocesor symbol for the plugin search paths

    AC_DEFINE_UNQUOTED(CR2RES_PLUGIN_DIR, "${PACKAGE}/plugins",
                       [Plugin directory tree prefix])

    eval plugin_dir="$plugindir"
    plugin_path=`eval echo $plugin_dir | \
                sed -e "s/\/${PACKAGE}-${VERSION}.*$//"`

    AC_DEFINE_UNQUOTED(CR2RES_PLUGIN_PATH, "$plugin_path",
                       [Absolute path to the plugin directory tree])

])


# CR2RES_CREATE_SYMBOLS
#-----------------------
# Define include and library related makefile symbols
AC_DEFUN([CR2RES_CREATE_SYMBOLS],
[

    # Symbols for package include file and library search paths

    CR2RES_INCLUDES='-I$(top_srcdir)/cr2res -I$(top_srcdir)/hdrl -I$(top_srcdir)/irplib'
    CR2RES_LDFLAGS='-L$(top_builddir)/cr2res'

    # Library aliases

    LIBCR2RES='$(top_builddir)/cr2res/libcr2res.la'
    LIBIRPLIB='$(top_builddir)/irplib/libirplib.la'
    LIBHDRL='$(top_builddir)/hdrl/libhdrl.la'

    # Substitute the defined symbols

    AC_SUBST(CR2RES_INCLUDES)
    AC_SUBST(CR2RES_LDFLAGS)

    AC_SUBST(LIBCR2RES)
    AC_SUBST(LIBIRPLIB)
    AC_SUBST(LIBHDRL)

    # Check for CPL and user defined libraries
    AC_REQUIRE([CPL_CHECK_LIBS])
    AC_REQUIRE([ESO_CHECK_EXTRA_LIBS])

    all_includes='$(CR2RES_INCLUDES) $(CPL_INCLUDES) $(EXTRA_INCLUDES)'
    all_ldflags='$(CR2RES_LDFLAGS) $(CPL_LDFLAGS) $(EXTRA_LDFLAGS)'

    AC_SUBST(all_includes)
    AC_SUBST(all_ldflags)
])
