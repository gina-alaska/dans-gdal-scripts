# Copyright 2006 University of Alaska,
# Geographic Information Network of Alaska
# All Rights Reserved
# UAF Confidential Information

AC_PREREQ(2.59)
AC_INIT([dans-gdal-scripts], [0.06], [dstahlke@gi.alaska.edu])
AC_CONFIG_SRCDIR([src/common.h])
AC_CONFIG_HEADER([config.h])

AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AM_PROG_CC_C_O

# beginning of GDAL stuff #####################

AC_ARG_WITH(gdal,
[  --with-gdal=ARG  Include gdal support (ARG=path to gdal-config)],,)

if test "z$with_gdal" = "zyes" ; then
	AC_PATH_PROG([with_gdal], [gdal-config])
fi

if test "z$with_gdal" = "z"; then
	AC_PATH_PROG([with_gdal], [gdal-config])
fi

if test "z$with_gdal" = "zno" -o "z$with_gdal" = "z"; then
	AC_MSG_ERROR(gdal-config not found)
fi

GDAL_CONFIG=$with_gdal
AC_SUBST(GDAL_CFLAGS)
GDAL_CFLAGS=`${GDAL_CONFIG} --cflags`
AC_SUBST(GDAL_LIBS)
GDAL_LIBS=`${GDAL_CONFIG} --libs`

# end of GDAL stuff #####################

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([inttypes.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([memset sqrt strtol])

AC_CHECK_LIB([m], [pow], , AC_MSG_ERROR([math library is required]))
AC_CHECK_LIB([proj], [pj_init], , AC_MSG_ERROR([proj4 library is required]))

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
