#### Basic compiler flags ######################################################

CXX=g++
DEBUGFLAGS=-g -Wall
GPPFLAGS=-O2 -fopenmp
CXXFLAGS=$(DEBUGFLAGS)
#CXXFLAGS=$(GPPFLAGS)

################################################################################

#### Main program ##############################################################

OBJ=util.o transformations.o geometry.o singularities.o vertex.o element.o panel.o tripanel.o quadpanel.o vortex.o vortex_ring.o horseshoe_vortex.o
HSHOE=test_hshoe
DOUBLET=test_doublet
DOUBLET2=test_doublet2
TRIDOUBLET=test_tri_doublet
ROTDOUBLET=test_doublet_rotation
ROTDOUBLET2=test_doublet_rotation2
POINTDOUBLET=test_point_doublet
ONPANDOUBLET=test_doublet_onpanel
ONPANDOUBLET2=test_doublet_onpanel2
SOURCE=test_source
ROTSOURCE=test_source_rotation
ONPANSOURCE=test_source_onpanel
MIRROR=test_mirror
POTENTIAL=test_quad_potential
POTENTIAL2=test_tri_potential
SRCDIR=../../src
INCLUDE=-I../../include -I/usr/include/eigen3
#INCLUDE=-I../../include -I/data/dprosser/locally_installed/include/eigen3 -I/data/dprosser/locally_installed/include
#LDFLAGS=-L/data/dprosser/locally_installed/lib64

################################################################################

#### Preprocessor variables ####################################################

ifeq ($(CXXFLAGS), $(DEBUGFLAGS))
  PREPROC=-DDEBUG
else
  PREPROC=-UDEBUG
endif

################################################################################

all: $(HSHOE) $(DOUBLET) $(DOUBLET2) $(TRIDOUBLET) $(ROTDOUBLET) $(ROTDOUBLET2) $(POINTDOUBLET) $(ONPANDOUBLET) $(ONPANDOUBLET2) $(SOURCE) $(ROTSOURCE) $(ONPANSOURCE) $(MIRROR) $(POTENTIAL) $(POTENTIAL2)

$(HSHOE): $(OBJ) test_hshoe.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(PREPROC) -o $(HSHOE) $(OBJ) test_hshoe.o

$(DOUBLET): $(OBJ) test_doublet.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(PREPROC) -o $(DOUBLET) $(OBJ) test_doublet.o

$(DOUBLET2): $(OBJ) test_doublet2.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(PREPROC) -o $(DOUBLET2) $(OBJ) test_doublet2.o

$(TRIDOUBLET): $(OBJ) test_tri_doublet.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(PREPROC) -o $(TRIDOUBLET) $(OBJ) test_tri_doublet.o

$(ROTDOUBLET): $(OBJ) test_doublet_rotation.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(PREPROC) -o $(ROTDOUBLET) $(OBJ) test_doublet_rotation.o

$(ROTDOUBLET2): $(OBJ) test_doublet_rotation2.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(PREPROC) -o $(ROTDOUBLET2) $(OBJ) test_doublet_rotation2.o

$(POINTDOUBLET): $(OBJ) test_point_doublet.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(PREPROC) -o $(POINTDOUBLET) $(OBJ) test_point_doublet.o

$(ONPANDOUBLET): $(OBJ) test_doublet_onpanel.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(PREPROC) -o $(ONPANDOUBLET) $(OBJ) test_doublet_onpanel.o

$(ONPANDOUBLET2): $(OBJ) test_doublet_onpanel2.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(PREPROC) -o $(ONPANDOUBLET2) $(OBJ) test_doublet_onpanel2.o

$(SOURCE): $(OBJ) test_source.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(PREPROC) -o $(SOURCE) $(OBJ) test_source.o

$(ROTSOURCE): $(OBJ) test_source_rotation.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(PREPROC) -o $(ROTSOURCE) $(OBJ) test_source_rotation.o

$(ONPANSOURCE): $(OBJ) test_source_onpanel.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(PREPROC) -o $(ONPANSOURCE) $(OBJ) test_source_onpanel.o

$(MIRROR): $(OBJ) test_mirror.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(PREPROC) -o $(MIRROR) $(OBJ) test_mirror.o

$(POTENTIAL): $(OBJ) test_quad_potential.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(PREPROC) -o $(POTENTIAL) $(OBJ) test_quad_potential.o

$(POTENTIAL2): $(OBJ) test_tri_potential.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(PREPROC) -o $(POTENTIAL2) $(OBJ) test_tri_potential.o

clean: 
	rm -f *.o

util.o: $(SRCDIR)/util.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/util.cpp

clo_parser.o: $(SRCDIR)/clo_parser.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/clo_parser.cpp

algorithms.o: $(SRCDIR)/algorithms.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/algorithms.cpp

transformations.o: $(SRCDIR)/transformations.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/transformations.cpp

geometry.o: $(SRCDIR)/geometry.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/geometry.cpp

singularities.o: $(SRCDIR)/singularities.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/singularities.cpp

vertex.o: $(SRCDIR)/vertex.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/vertex.cpp

element.o: $(SRCDIR)/element.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/element.cpp

panel.o: $(SRCDIR)/panel.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/panel.cpp

tripanel.o: $(SRCDIR)/tripanel.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/tripanel.cpp

quadpanel.o: $(SRCDIR)/quadpanel.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/quadpanel.cpp

sectional_object.o: $(SRCDIR)/sectional_object.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/sectional_object.cpp

airfoil.o: $(SRCDIR)/airfoil.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/airfoil.cpp

section.o: $(SRCDIR)/section.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/section.cpp

vortex.o: $(SRCDIR)/vortex.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/vortex.cpp

horseshoe_vortex.o: $(SRCDIR)/horseshoe_vortex.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/horseshoe_vortex.cpp

vortex_ring.o: $(SRCDIR)/vortex_ring.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/vortex_ring.cpp

wake.o: $(SRCDIR)/wake.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/wake.cpp

wing.o: $(SRCDIR)/wing.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/wing.cpp

aircraft.o: $(SRCDIR)/aircraft.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/aircraft.cpp

test_hshoe.o: test_hshoe.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) test_hshoe.cpp

test_doublet.o: test_doublet.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) test_doublet.cpp

test_doublet2.o: test_doublet2.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) test_doublet2.cpp

test_tri_doublet.o: test_tri_doublet.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) test_tri_doublet.cpp

test_doublet_rotation.o: test_doublet_rotation.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) test_doublet_rotation.cpp

test_doublet_rotation2.o: test_doublet_rotation2.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) test_doublet_rotation2.cpp

test_point_doublet.o: test_point_doublet.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) test_point_doublet.cpp

test_doublet_onpanel.o: test_doublet_onpanel.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) test_doublet_onpanel.cpp

test_doublet_onpanel2.o: test_doublet_onpanel2.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) test_doublet_onpanel2.cpp

test_source.o: test_source.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) test_source.cpp

test_source_rotation.o: test_source_rotation.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) test_source_rotation.cpp

test_source_onpanel.o: test_source_onpanel.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) test_source_onpanel.cpp

test_mirror.o: test_mirror.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) test_mirror.cpp

test_quad_potential.o: test_quad_potential.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) test_quad_potential.cpp

test_tri_potential.o: test_tri_potential.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) test_tri_potential.cpp
