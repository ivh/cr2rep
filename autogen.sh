#!/bin/sh

bt_usage ()
{
    echo "Usage: $1 $2"
}

bt_libtoolize ()
{

  uses_libltdl=0
  libltdl_directory="$1"

  egrep '^ *AC_LIBLTDL_(CONVENIENCE|INSTALLABLE)' configure.ac >/dev/null 2>&1

  if test "$?" = 0; then
    uses_libltdl=1
  else
    egrep '^ *LT_CONFIG_LTDL_DIR' configure.ac >/dev/null 2>&1

    if test "$?" != 0 && test -d "$libltdl_directory"; then
        uses_libltdl=1
    fi
  fi

  if test $uses_libltdl = 1; then
    libtoolize --force --copy --ltdl >/dev/null 2>&1
  fi
  return 0

}

bt_script=`basename $0`
bt_options="hvW"

set -- `getopt $bt_options $*`

for opt in $*; do
    case $opt in
    -h) bt_usage $bt_script "[-$bt_options]"
        exit 0
        ;;

    -v) bt_args="-v"
        shift
        ;;

    -W) bt_warnings=yes
        bt_args="-v"
        shift
        ;;

    --) shift
        break
        ;;
    esac
done


bt_cmd="`which autoreconf 2>/dev/null | egrep '^/' | head -n 1`"
bt_args="$bt_args -if"

if test ! -x $bt_cmd; then
    echo "Error: 'autoreconf' not found! Check your autoconf installation!"
    exit 1
fi

echo "Bootstrapping build tree in \`$PWD'"
echo "Please wait ..."

bootstrap="$bt_cmd $bt_args"

if test -f configure; then
    rm -f configure
fi


# The following is only required for libtool versions which do not support
# the LT_CONFIG_LTDL_DIR macro, i.e. all versions prior to libtool 2.0
# in this case we have to update the ltdl files by explicitly calling
# libtoolize with the appropriate options. Libtool 2.x and later is properly
# supported by autoreconf and the workaround below is not needed.

bt_libtool_version="`libtool --version 2>/dev/null | head -n 1 |
                        sed -e 's/^[^0-9]*//'`"
bt_libtool_major="`echo $bt_libtool_version | cut -d . -f 1`"
bt_libtool_minor="`echo $bt_libtool_version | cut -d . -f 2`"

test -z "$bt_libtool_major" && bt_libtool_major=0
test -z "$bt_libtool_minor" && bt_libtool_minor=0

bt_libtoolize "libltdl"


if test x"$bt_warnings" = xyes; then
    eval $bootstrap 2>&1 | tee autogen.log | egrep '^autoreconf:'
    update="`grep 'autoupdate' autogen.log | wc -l`"

    if test x$update != x0; then
        echo -n "Warning: obsolete constructs found. Running 'autoupdate' "
        echo "is recommended! Check logfile './autogen.log' for details."
    fi

else
    eval $bootstrap >/dev/null 2>&1
fi

if test -f configure; then
    echo ""
    echo "Don't forget to run ./configure"
    echo "If you haven't done so in a while, run ./configure --help"
    echo ""
fi

exit 0
