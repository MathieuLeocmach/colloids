AC_INIT([colloids-tools], [1.0], [mathieu@iis.u-tokyo.ac.jp])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
AC_PROG_CXX
AC_LANG([C++])
AC_CONFIG_SRCDIR([rdf.cpp])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile])
LT_INIT

CXXFLAGS="-O3 -Wall"

AC_OPENMP
if test "$OPENMP_CXXFLAGS" != ""; then
	CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
	CPPFLAGS="$CPPFLAGS -Dcimg_use_openmp"
fi

AC_ARG_WITH(with-voro-src, [  --with-voro-src=DIR   Voro++ source files are in DIR])
AM_CONDITIONAL([WANT_VORO], test $with_voro_src)
if test $with_voro_src; then 	
	AC_SUBST(VORO_SRC, $with_voro_src)
fi
AC_OUTPUT
