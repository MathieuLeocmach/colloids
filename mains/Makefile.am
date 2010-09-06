ACLOCAL_AMFLAGS = -I m4
LDADD = ../lib/libcolloids.la
AM_CPPFLAGS = -I$(srcdir)/../lib 

bin_PROGRAMS = bonds boo boo_flip cage2vtk cutter dat2vtk drift dynamics g6 ISF linkboo linker lostngb MSD percolation rdf sp5c totalRdf traj2vtk periodic_rdf periodic_boo $(binvoro)

EXTRA_PROGRAMS = cgVoro periodic_cgVoro

periodic_boo_CPPFLAGS = $(AM_CPPFLAGS) -Duse_periodic
periodic_rdf_CPPFLAGS = $(AM_CPPFLAGS) -Duse_periodic
cgVoro_CPPFLAGS = $(AM_CPPFLAGS) -I$(VORO_SRC)
periodic_cgVoro_CPPFLAGS = $(cgVoro_CPPFLAGS) -Duse_periodic

bonds_SOURCES = bonds.cpp
boo_SOURCES = boo.cpp
periodic_boo_SOURCES = boo.cpp
boo_flip_SOURCES = boo_flip.cpp
cage2vtk_SOURCES = cage2vtk.cpp
cgVoro_SOURCES = cgVoro.cpp
periodic_cgVoro_SOURCES = cgVoro.cpp
cutter_SOURCES = cutter.cpp
dat2vtk_SOURCES = dat2vtk.cpp
drift_SOURCES = drift.cpp
dynamics_SOURCES = dynamics.cpp
g6_SOURCES = g6.cpp
ISF_SOURCES = ISF.cpp
linkboo_SOURCES = linkboo.cpp
linker_SOURCES = linker.cpp
lostngb_SOURCES = lostngb.cpp
MSD_SOURCES = MSD.cpp
percolation_SOURCES = percolation.cpp
rdf_SOURCES = rdf.cpp
periodic_rdf_SOURCES = rdf.cpp
sp5c_SOURCES = sp5c.cpp
totalRdf_SOURCES = totalRdf.cpp
traj2vtk_SOURCES = traj2vtk.cpp
