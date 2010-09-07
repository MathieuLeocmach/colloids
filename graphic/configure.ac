AC_INIT([colloids-graphic], [1.0], [mathieu@iis.u-tokyo.ac.jp])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
AC_PROG_CXX
AC_LANG([C++])
AC_CONFIG_SRCDIR([lifFile.hpp])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SUBDIRS([mains])
AC_CONFIG_FILES([Makefile tinyxml/Makefile])
LT_INIT

CXXFLAGS="-O3 -Wall"
AX_BOOST_BASE([1.35], , AC_MSG_ERROR('boost >1.35 is needed to calculate spherical harmonics'))
AX_BOOST_PROGRAM_OPTIONS
AC_CHECK_HEADER([CImg.h], , AC_MSG_ERROR('CImg >1.32 is needed'))

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

AC_CHECK_LIB([fftw3f], [fftwf_free])
AC_CHECK_LIB([fftw3f_threads], [fftwf_init_threads])
AC_CHECK_LIB([X11], [XOpenDisplay])
AC_CHECK_LIB([gdi32],[GetPixel])
AC_CHECK_LIB([m], [sqrt])

AC_OUTPUT
