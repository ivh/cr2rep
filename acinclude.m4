# CR2RE_SET_PREFIX(PREFIX)
#---------------------------
AC_DEFUN([CR2RE_SET_PREFIX],
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


# CR2RE_SET_VERSION_INFO(VERSION, [CURRENT], [REVISION], [AGE])
#----------------------------------------------------------------
# Setup various version information, especially the libtool versioning
AC_DEFUN([CR2RE_SET_VERSION_INFO],
[
    cr2re_version=`echo "$1" | sed -e 's/[[a-z,A-Z]].*$//'`

    cr2re_major_version=`echo "$cr2re_version" | \
        sed 's/\([[0-9]]*\).\(.*\)/\1/'`
    cr2re_minor_version=`echo "$cr2re_version" | \
        sed 's/\([[0-9]]*\).\([[0-9]]*\)\(.*\)/\2/'`
    cr2re_micro_version=`echo "$cr2re_version" | \
        sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\3/'`

    if test -z "$cr2re_major_version"; then cr2re_major_version=0
    fi

    if test -z "$cr2re_minor_version"; then cr2re_minor_version=0
    fi

    if test -z "$cr2re_micro_version"; then cr2re_micro_version=0
    fi

    CR2RE_VERSION="$cr2re_version"
    CR2RE_MAJOR_VERSION=$cr2re_major_version
    CR2RE_MINOR_VERSION=$cr2re_minor_version
    CR2RE_MICRO_VERSION=$cr2re_micro_version

    if test -z "$4"; then CR2RE_INTERFACE_AGE=0
    else CR2RE_INTERFACE_AGE="$4"
    fi

    CR2RE_BINARY_AGE=`expr 100 '*' $CR2RE_MINOR_VERSION + $CR2RE_MICRO_VERSION`
    CR2RE_BINARY_VERSION=`expr 10000 '*' $CR2RE_MAJOR_VERSION + \
                          $CR2RE_BINARY_AGE`

    AC_SUBST(CR2RE_VERSION)
    AC_SUBST(CR2RE_MAJOR_VERSION)
    AC_SUBST(CR2RE_MINOR_VERSION)
    AC_SUBST(CR2RE_MICRO_VERSION)
    AC_SUBST(CR2RE_INTERFACE_AGE)
    AC_SUBST(CR2RE_BINARY_VERSION)
    AC_SUBST(CR2RE_BINARY_AGE)

    AC_DEFINE_UNQUOTED(CR2RE_MAJOR_VERSION, $CR2RE_MAJOR_VERSION,
                       [CR2RE major version number])
    AC_DEFINE_UNQUOTED(CR2RE_MINOR_VERSION, $CR2RE_MINOR_VERSION,
                       [CR2RE minor version number])
    AC_DEFINE_UNQUOTED(CR2RE_MICRO_VERSION, $CR2RE_MICRO_VERSION,
                       [CR2RE micro version number])
    AC_DEFINE_UNQUOTED(CR2RE_INTERFACE_AGE, $CR2RE_INTERFACE_AGE,
                       [CR2RE interface age])
    AC_DEFINE_UNQUOTED(CR2RE_BINARY_VERSION, $CR2RE_BINARY_VERSION,
                       [CR2RE binary version number])
    AC_DEFINE_UNQUOTED(CR2RE_BINARY_AGE, $CR2RE_BINARY_AGE,
                       [CR2RE binary age])

    ESO_SET_LIBRARY_VERSION([$2], [$3], [$4])
])


# CR2RE_SET_PATHS
#------------------
# Define auxiliary directories of the installed directory tree.
AC_DEFUN([CR2RE_SET_PATHS],
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
       configdir='${datadir}/${PACKAGE}/config'
    fi

    if test -z "$wkfextradir"; then
        wkfextradir='${datadir}/esopipes/${PACKAGE}-${VERSION}/reflex'
    fi

    if test -z "$wkfcopydir"; then
        wkfcopydir='${datadir}/reflex/workflows/${PACKAGE}-${VERSION}'
    fi

    AC_SUBST(plugindir)
    AC_SUBST(privatelibdir)
    AC_SUBST(apidocdir)
    AC_SUBST(pipedocsdir)
    AC_SUBST(configdir)
    AC_SUBST(wkfextradir)
    AC_SUBST(wkfcopydir)


    # Define a preprocesor symbol for the plugin search paths

    AC_DEFINE_UNQUOTED(CR2RE_PLUGIN_DIR, "${PACKAGE}/plugins",
                       [Plugin directory tree prefix])

    eval plugin_dir="$plugindir"
    plugin_path=`eval echo $plugin_dir | \
                sed -e "s/\/${PACKAGE}-${VERSION}.*$//"`

    AC_DEFINE_UNQUOTED(CR2RE_PLUGIN_PATH, "$plugin_path",
                       [Absolute path to the plugin directory tree])

])


# CR2RE_CREATE_SYMBOLS
#-----------------------
# Define include and library related makefile symbols
AC_DEFUN([CR2RE_CREATE_SYMBOLS],
[

    # Symbols for package include file and library search paths

    CR2RE_INCLUDES='-I$(top_srcdir)/cr2re'
    CR2RE_LDFLAGS='-L$(top_builddir)/cr2re'

    # Library aliases

    LIBCR2RE='$(top_builddir)/cr2re/libcr2re.la'
    LIBIRPLIB='$(top_builddir)/irplib/libirplib.la'
    LIBHDRL='$(top_builddir)/hdrl/libhdrl.la'

    # Substitute the defined symbols

    AC_SUBST(CR2RE_INCLUDES)
    AC_SUBST(CR2RE_LDFLAGS)

    AC_SUBST(LIBCR2RE)
    AC_SUBST(LIBHDRL)

    # Check for CPL and user defined libraries
    AC_REQUIRE([CPL_CHECK_LIBS])
    AC_REQUIRE([ESO_CHECK_EXTRA_LIBS])

    all_includes='$(CR2RE_INCLUDES) $(CPL_INCLUDES) $(EXTRA_INCLUDES) $(XXCLIPM_INCLUDES)'
    all_ldflags='$(CR2RE_LDFLAGS) $(CPL_LDFLAGS) $(EXTRA_LDFLAGS) $(XXCLIPM_LDFLAGS)'

    AC_SUBST(all_includes)
    AC_SUBST(all_ldflags)
])
