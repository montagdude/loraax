#### Basic compiler flags ######################################################

CXX=g++
DEBUGFLAGS=-g -Wall
GPPFLAGS=-O2 -fopenmp
#CXXFLAGS=$(DEBUGFLAGS)
CXXFLAGS=$(GPPFLAGS)

################################################################################

#### Main program ##############################################################

OBJ=util.o clo_parser.o algorithms.o settings.o transformations.o geometry.o singularities.o vertex.o element.o panel.o tripanel.o quadpanel.o sectional_object.o airfoil.o section.o wake_strip.o wake.o wing.o aircraft.o loraax.o
EXE=loraax
SRCDIR=../src
#INCLUDE=-I../include -I/usr/include/eigen3
INCLUDE=-I../include -I/data/dprosser/locally_installed/include/eigen3 -I/data/dprosser/locally_installed/include
LDFLAGS=-L/data/dprosser/locally_installed/lib64
LIBS=-lxfoil -ltinyxml2

################################################################################

#### Preprocessor variables ####################################################

ifeq ($(CXXFLAGS), $(DEBUGFLAGS))
  PREPROC=-DDEBUG
else
  PREPROC=-UDEBUG
endif

################################################################################

all: $(EXE)

$(EXE): $(OBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(PREPROC) -o $(EXE) $(OBJ) $(LIBS)

clean: 
	rm -f *.o

allclean: clean
	rm -f $(EXE)

util.o: $(SRCDIR)/util.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/util.cpp

clo_parser.o: $(SRCDIR)/clo_parser.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/clo_parser.cpp

algorithms.o: $(SRCDIR)/algorithms.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/algorithms.cpp

settings.o: $(SRCDIR)/settings.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/settings.cpp

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

wake_strip.o: $(SRCDIR)/wake_strip.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/wake_strip.cpp

wake.o: $(SRCDIR)/wake.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/wake.cpp

wing.o: $(SRCDIR)/wing.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/wing.cpp

aircraft.o: $(SRCDIR)/aircraft.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/aircraft.cpp

loraax.o: $(SRCDIR)/loraax.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/loraax.cpp
