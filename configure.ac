#  configure.ac   for VFGEN (https://github.com/WarrenWeckesser/vfgen)
#  Copyright (c) 2008-2022 Warren Weckesser
#
#  Process this file with autoconf to produce a configure script.
#

AC_PREREQ([2.71])
AC_INIT([vfgen],[m4_esyscmd_s(cat version.txt)],[warren DOT weckesser AT gmail])

AC_MSG_NOTICE([[
------------------------------------
VFGEN Configure Script
------------------------------------
Basic environment check...
------------------------------------]])

AM_INIT_AUTOMAKE([foreign 1.8])

AC_MSG_NOTICE([[
------------------------------------
Checking your C++ compiler...
------------------------------------]])

AC_PROG_CXX
AC_LANG(C++)

AC_MSG_NOTICE([[
------------------------------------
Checking for MXML and GINAC libs...
------------------------------------]])

#
# Check for pkg-config first
#
AC_CHECK_PROG(HAVE_PKGCONFIG, pkg-config, yes, no)
dnl Give error and exit if we don't have pkgconfig
if test "x$HAVE_PKGCONFIG" = "xno"; then
    AC_MSG_ERROR(dnl
[The pkg-config script could not be found. Make sure it is
in your path, or set the PKG_CONFIG environment variable
to the full path to pkg-config.
Or see http://www.freedesktop.org/software/pkgconfig to get pkg-config.])
fi

#
# Check for ginac and mxml (Mini-XML)
#

PKG_CHECK_MODULES([MXML], [mxml], HAVE_MXML=yes,HAVE_MXML=no)
if test "x$HAVE_MXML" = "xno"; then
    AC_MSG_ERROR(dnl
[The Mini-XML (http://www.minixml.org) library was not found.])
fi

PKG_CHECK_MODULES([GINAC], [ginac], HAVE_GINAC=yes,HAVE_GINAC=no)
if test "x$HAVE_GINAC" = "xno"; then
    AC_MSG_ERROR(dnl
[The GiNaC (http://www.ginac.de) library was not found.])
fi

# If the MXML and GINAC CFLAGS options are the same, we don't need both...
if test "$GINAC_CFLAGS" = "$MXML_CFLAGS"; then
    MXML_CFLAGS=""
fi

AC_MSG_NOTICE([[
------------------------------------
Creating the Makefile...
------------------------------------]])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT

AC_MSG_NOTICE([[
------------------------------------
Now run
$ make
$ make install
to compile and install VFGEN
------------------------------------]])
