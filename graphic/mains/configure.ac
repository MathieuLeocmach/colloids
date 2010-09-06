AC_INIT([colloids-graphic-tools], [1.0], [mathieu@iis.u-tokyo.ac.jp])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
AC_PROG_CXX
AC_LANG([C++])
AC_CONFIG_SRCDIR([tracker.cpp])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile])
LT_INIT

CXXFLAGS="-O3 -Wall"
AX_BOOST_BASE([1.35], , AC_MSG_ERROR('boost >1.35 is needed to calculate spherical harmonics'))
AX_BOOST_PROGRAM_OPTIONS
AC_CHECK_HEADER([CImg.h], , AC_MSG_ERROR('CImg >1.32 is needed'))
AC_CHECK_HEADER([fftw3.h], , AC_MSG_ERROR('FFTW >3 is needed'))

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

AC_CHECK_LIB([fftw3f], [fftwf_free])
AC_CHECK_LIB([X11], [XOpenDisplay])
AC_CHECK_LIB([gdi32],[GetPixel])

AC_OUTPUT
