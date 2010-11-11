AC_INIT([colloids], [1.0], [mathieu@iis.u-tokyo.ac.jp])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([subdir-objects foreign -Wall -Werror])
AC_PROG_CXX
AC_LANG([C++])
AC_CONFIG_SRCDIR([lib/particles.hpp])
AC_CONFIG_HEADERS([config.h])
#AC_CONFIG_SUBDIRS([mains graphic])
AC_CONFIG_FILES([Makefile])
LT_INIT
AM_MAINTAINER_MODE([disable])
#AX_MAINTAINER_MODE_AUTO_SILENT

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
         CXXFLAGS="$CXXFLAGS -g -DDEBUG"
AC_MSG_RESULT(yes)
else
CXXFLAGS="$CXXFLAGS -O3 -Wall -DNDEBUG"
AC_MSG_RESULT(no)
fi

#Try to recognize the architecture
AX_GCC_ARCHFLAG([yes], [CXXFLAGS="$CXXFLAGS $ax_cv_gcc_archflag"])

AX_BOOST_BASE([1.35], , AC_MSG_ERROR('boost >1.35 is needed to calculate spherical harmonics'))
AX_BOOST_PROGRAM_OPTIONS

AC_CHECK_HEADER([CImg.h], , AC_MSG_ERROR('CImg >1.32 is needed'))

AX_COUNT_CPUS
CPPFLAGS="$CPPFLAGS -DTRACKER_N_THREADS=$CPU_COUNT"

AC_OPENMP
if test "$OPENMP_CXXFLAGS" != ""; then
	CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
	CPPFLAGS="$CPPFLAGS -Dcimg_use_openmp"
fi

# Checks for pthread.
AC_CHECK_HEADERS([pthread.h], [
    AC_CHECK_LIB(pthread, pthread_create, [
	LIBS="$LIBS -lpthread"
	AC_DEFINE(HAVE_LIBNCURSES, 1, [Define to 1 if you have the 'pthread' library (-lpthread).])
    ])
])

if test "$PTHREAD_LIBS" != ""; then
	LIBS="$PTHREAD_LIBS $LIBS"
	CXXFLAGS="$CXXFLAGS $PTHREAD_CXXFLAGS"
	CC="$PTHREAD_CC"
fi

#FFTW
AC_CHECK_HEADER([fftw3.h], , AC_MSG_ERROR('FFTW >3 is needed'))
AC_CHECK_LIB([fftw3f], [fftwf_free])
AC_CHECK_LIB([fftw3f_threads], [fftwf_init_threads])

#Display Unix
AC_CHECK_LIB([X11], [XOpenDisplay])
#Display Windows
AC_CHECK_LIB([gdi32],[GetPixel])
#Math library
AC_CHECK_LIB([m], [sqrt])

AC_ARG_WITH(voro-src, [  --with-voro-src=DIR   Voro++ source files are in DIR])
if test $with_voro_src; then 	
	AC_CHECK_FILE(
		[$with_voro_src/voro++.cc],
		[
			voro=true
			AC_SUBST(VORO_SRC, $with_voro_src)
			AC_SUBST(binvoro, "cgVoro periodic_cgVoro")
		],
		[AC_MSG_ERROR("No voro++.cc at the indicated path")]
	)
fi
AM_CONDITIONAL([WANT_VORO],[test x$voro = xtrue])

AC_OUTPUT
