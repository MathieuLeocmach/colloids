ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES = libcolloids.la libcolloids-graphic.la

include_HEADERS = lib/boo_data.hpp lib/fields.hpp lib/particles.hpp lib/dynamicClusters.hpp lib/files_series.hpp lib/periodic.hpp lib/dynamicParticles.hpp lib/index.hpp lib/traj.hpp graphic/lifFile.hpp graphic/lifTracker.hpp graphic/radiiTracker.hpp graphic/serieTracker.hpp graphic/tracker.hpp

libcolloids_la_SOURCES = lib/boo_data.cpp lib/fields.cpp lib/particles.cpp lib/dynamicClusters.cpp lib/files_series.cpp lib/periodic.cpp lib/dynamicParticles.cpp lib/index.cpp lib/traj.cpp lib/boo_data.hpp lib/fields.hpp lib/particles.hpp lib/dynamicClusters.hpp lib/files_series.hpp lib/periodic.hpp lib/dynamicParticles.hpp lib/index.hpp lib/traj.hpp lib/RStarTree/RStarBoundingBox.h lib/RStarTree/RStarTree.h lib/RStarTree/RStarVisitor.h

LDADD = libcolloids.la
AM_CPPFLAGS = -I$(srcdir)/lib -I$(srcdir)/graphic -DTIXML_USE_STL 
bin_PROGRAMS = bonds boo boo_flip cage2vtk cutter dat2vtk dhlngb drift dynamics g6 ISF linkboo linker lostngb MSD percolation rdf sp5c totalRdf traj2vtk periodic_rdf periodic_boo periodic_g6 $(binvoro) aquireWisdom tracker

EXTRA_PROGRAMS = cgVoro periodic_cgVoro
periodic_boo_CPPFLAGS = $(AM_CPPFLAGS) -Duse_periodic
periodic_rdf_CPPFLAGS = $(AM_CPPFLAGS) -Duse_periodic
periodic_g6_CPPFLAGS = $(AM_CPPFLAGS) -Duse_periodic
cgVoro_CPPFLAGS = $(AM_CPPFLAGS) -I$(VORO_SRC)
periodic_cgVoro_CPPFLAGS = $(cgVoro_CPPFLAGS) -Duse_periodic

bonds_SOURCES = mains/bonds.cpp
boo_SOURCES = mains/boo.cpp
periodic_boo_SOURCES = mains/boo.cpp
boo_flip_SOURCES = mains/boo_flip.cpp
cage2vtk_SOURCES = mains/cage2vtk.cpp
cgVoro_SOURCES = mains/cgVoro.cpp
periodic_cgVoro_SOURCES = mains/cgVoro.cpp
cutter_SOURCES = mains/cutter.cpp
dat2vtk_SOURCES = mains/dat2vtk.cpp
dhlngb_SOURCES = mains/dhlngb.cpp
drift_SOURCES = mains/drift.cpp
dynamics_SOURCES = mains/dynamics.cpp
g6_SOURCES = mains/g6.cpp
periodic_g6_SOURCES = mains/g6.cpp
ISF_SOURCES = mains/ISF.cpp
linkboo_SOURCES = mains/linkboo.cpp
linker_SOURCES = mains/linker.cpp
lostngb_SOURCES = mains/lostngb.cpp
MSD_SOURCES = mains/MSD.cpp
percolation_SOURCES = mains/percolation.cpp
rdf_SOURCES = mains/rdf.cpp
periodic_rdf_SOURCES = mains/rdf.cpp
sp5c_SOURCES = mains/sp5c.cpp
totalRdf_SOURCES = mains/totalRdf.cpp
traj2vtk_SOURCES = mains/traj2vtk.cpp

cutter_LDFLAGS = $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB)


libcolloids_graphic_la_SOURCES = graphic/lifFile.cpp graphic/lifFile.hpp graphic/lifTracker.cpp graphic/lifTracker.hpp graphic/radiiTracker.cpp graphic/radiiTracker.hpp graphic/serieTracker.cpp graphic/serieTracker.hpp graphic/tracker.cpp graphic/tracker.hpp graphic/tinyxml/tinystr.h graphic/tinyxml/tinyxmlerror.cpp graphic/tinyxml/tinyxmlparser.cpp graphic/tinyxml/tinystr.cpp graphic/tinyxml/tinyxml.cpp graphic/tinyxml/tinyxml.h

LDADD += libcolloids-graphic.la

aquireWisdom_SOURCES = graphic/mains/aquireWisdom.cpp
tracker_SOURCES = graphic/mains/tracker.cpp

tracker_LDFLAGS = $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIB)
