#### Basic compiler flags ######################################################

CXX=g++
DEBUGFLAGS=-g -Wall
GPPFLAGS=-O3 -flto -funroll-loops -fopenmp
CXXFLAGS=$(DEBUGFLAGS)

################################################################################

#### Main program ##############################################################

OBJ=util.o clo_parser.o algorithms.o settings.o transformations.o geometry.o singularities.o vertex.o face.o triface.o quadface.o airfoil.o section.o panelcode.o
EXE=panelcode
SRCDIR=../src
INCLUDE=-I../include -I/usr/include/eigen3
#INCLUDE=-I../include -I/data/dprosser/locally_installed/include/eigen3 -I/data/dprosser/locally_installed/include
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

face.o: $(SRCDIR)/face.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/face.cpp

triface.o: $(SRCDIR)/triface.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/triface.cpp

quadface.o: $(SRCDIR)/quadface.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/quadface.cpp

airfoil.o: $(SRCDIR)/airfoil.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/airfoil.cpp

section.o: $(SRCDIR)/section.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/section.cpp

panelcode.o: $(SRCDIR)/panelcode.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/panelcode.cpp
