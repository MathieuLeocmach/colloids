ACLOCAL_AMFLAGS = -I m4
SUBDIRS = tinyxml . mains
lib_LTLIBRARIES = libcolloids-graphic.la
LDADD = ../lib/libcolloids.la tinyxml/libtinyxml.a
AM_CPPFLAGS = -I$(srcdir)/../lib -DTIXML_USE_STL
libcolloids_graphic_la_SOURCES = lifFile.cpp lifFile.hpp lifTracker.cpp lifTracker.hpp radiiTracker.cpp radiiTracker.hpp serieTracker.cpp serieTracker.hpp tracker.cpp tracker.hpp
include_HEADERS = lifFile.hpp lifTracker.hpp radiiTracker.hpp serieTracker.hpp tracker.hpp
