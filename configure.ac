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

dnl get the flags
CXXFLAGS="${CXXFLAGS=}"
dnl this macro is used to get the arguments supplied
dnl to the configure script (./configure --enable-debug)
dnl Check if we have enable debug support.
AC_MSG_CHECKING(whether to enable debugging)
debug_default="no"
AC_ARG_ENABLE(debug, [  --enable-debug=[no/yes] turn on debugging
                       [default=$debug_default]],, enable_debug=$debug_default)
dnl Yes, shell scripts can be used
if test "x$enable_debug" = "xyes"; then
         CXXFLAGS="$CXXFLAGS -g -DDEBUG -mtune=native"
AC_MSG_RESULT(yes)
else
CXXFLAGS="$CXXFLAGS -O3 -Wall -march=native -DNDEBUG"
AC_MSG_RESULT(no)
fi

AX_BOOST_BASE([1.35], , AC_MSG_ERROR('boost >1.35 is needed to calculate spherical harmonics'))
AX_BOOST_PROGRAM_OPTIONS

AX_COUNT_CPUS
CPPFLAGS="$CPPFLAGS -DTRACKER_N_THREADS=$CPU_COUNT"

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
