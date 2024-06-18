#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([ivar], [1.4.3], [gkarthik@scripps.edu])
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_AWK
AC_PROG_CXX
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_CPP
AC_PROG_MKDIR_P

LIBS="-lhts -lp -lpthread"

# introduce the optional configure parameter for a non-standard install prefix of XXX
AC_ARG_WITH([hts],
    [AS_HELP_STRING([--with-hts=prefix],
        [use this to specify a non-standard install prefix of the htslib])],
    [HTSPATHSET=1],
    [HTSPATHSET=0])

# if optional parameter used, extend path flags for compliler and linker
if test $HTSPATHSET = 1 ; then
    # extend the compiler and linker flags according to the path set
    AM_CXXFLAGS="$AM_CXXFLAGS -I$with_hts/include"
    AM_LDFLAGS="$AM_LDFLAGS -L$with_hts/lib"
    AC_SUBST([LIBS])
    AC_SUBST([AM_CXXFLAGS])
    AC_SUBST([AM_LDFLAGS])
fi

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lhts':
AC_CHECK_LIB([hts], [main])
# FIXME: Replace `main' with a function in `-lpthread':
AC_CHECK_LIB([pthread], [main])
# FIXME: Replace `main' with a function in `-lz':
AC_CHECK_LIB([z], [main])

# Checks for header files.
AC_CHECK_HEADERS([unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT8_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memmove])

AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])
AC_OUTPUT
