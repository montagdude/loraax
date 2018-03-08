#### Basic compiler flags ######################################################

CXX=g++
DEBUGFLAGS=-g -Wall
GPPFLAGS=-O3 -flto -funroll-loops -fopenmp
CXXFLAGS=$(DEBUGFLAGS)

################################################################################

#### Main program ##############################################################

OBJ=util.o settings.o transformations.o geometry.o singularities.o face.o triface.o quadface.o airfoil.o panelcode.o
EXE=panelcode
SRCDIR=../src
INCLUDE=-I../include -I/usr/include/eigen3

################################################################################

#### Preprocessor variables ####################################################

ifeq ($(CXXFLAGS), $(DEBUGFLAGS))
  PREPROC=-UNDEBUG
else
  PREPROC=-DNDEBUG
endif

################################################################################

all: $(EXE)

$(EXE): $(OBJ)
	$(CXX) $(CXXFLAGS) $(PREPROC) -o $(EXE) $(OBJ)

clean: 
	rm -f *.o

allclean: clean
	rm -f $(EXE)

util.o: $(SRCDIR)/util.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/util.cpp

settings.o: $(SRCDIR)/settings.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/settings.cpp

transformations.o: $(SRCDIR)/transformations.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/transformations.cpp

geometry.o: $(SRCDIR)/geometry.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/geometry.cpp

singularities.o: $(SRCDIR)/singularities.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/singularities.cpp

face.o: $(SRCDIR)/face.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/face.cpp

triface.o: $(SRCDIR)/triface.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/triface.cpp

quadface.o: $(SRCDIR)/quadface.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/quadface.cpp

airfoil.o: $(SRCDIR)/airfoil.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/airfoil.cpp

panelcode.o: $(SRCDIR)/panelcode.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/panelcode.cpp
