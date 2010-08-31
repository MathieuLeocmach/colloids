colloidsconfdir = $(HOME)/.colloids

ACLOCAL_AMFLAGS = -I m4
LDADD = ../../lib/libcolloids.la ../libcolloids-graphic.la ../tinyxml/libtinyxml.a
AM_CPPFLAGS = -I$(srcdir)/../../lib -I$(srcdir)/.. -DTIXML_USE_STL -DINSTAL_PATH="\"$(colloidsconfdir)/\""

bin_PROGRAMS = aquireWisdom tracker
colloidsconf_DATA = tracker.ini

aquireWisdom_SOURCES = aquireWisdom.cpp
tracker_SOURCES = tracker.cpp

tracker_LDFLAGS = $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB)
