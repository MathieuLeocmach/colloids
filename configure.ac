AC_INIT([colloids], [1.0], [mathieu@iis.u-tokyo.ac.jp])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
AC_PROG_CXX
AC_LANG([C++])
AC_CONFIG_SRCDIR([lib/particles.hpp])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SUBDIRS([mains graphic])
AC_CONFIG_FILES([Makefile lib/Makefile])
LT_INIT

CXXFLAGS="-O3 -Wall"
AX_BOOST_BASE([1.35], , AC_MSG_ERROR('boost >1.35 is needed to calculate spherical harmonics'))
AX_BOOST_PROGRAM_OPTIONS

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
