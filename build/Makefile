#### Basic compiler flags ######################################################

CXX=g++
DEBUGFLAGS=-g -Wall
GPPFLAGS=-O3 -flto -funroll-loops -fopenmp
CXXFLAGS=$(DEBUGFLAGS)

################################################################################

#### CPP Libs ##################################################################

LIBINSTALLDIR=/home/dan/cpplibs/lib
LIBINCLUDEDIR=/home/dan/cpplibs/include

################################################################################

#### Main program ##############################################################

OBJ=util.o settings.o transformations.o geometry.o singularities.o node.o edge.o face.o triface.o quadface.o body.o wakeline.o wake.o wake_relaxation.o file_io.o panelcode.o
EXE=panelcode
SRCDIR=../src
INCLUDE=-I../include -I$(LIBINCLUDEDIR) -I/usr/include/eigen3
LDFLAGS=-L$(LIBINSTALLDIR)
LIBS=-lnamelist

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
	$(CXX) $(CXXFLAGS) $(PREPROC) -o $(EXE) $(OBJ) $(LDFLAGS) $(LIBS) 

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

node.o: $(SRCDIR)/node.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/node.cpp

edge.o: $(SRCDIR)/edge.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/edge.cpp

face.o: $(SRCDIR)/face.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/face.cpp

triface.o: $(SRCDIR)/triface.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/triface.cpp

quadface.o: $(SRCDIR)/quadface.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/quadface.cpp

body.o: $(SRCDIR)/body.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/body.cpp

wakeline.o: $(SRCDIR)/wakeline.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/wakeline.cpp

wake.o: $(SRCDIR)/wake.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/wake.cpp

wake_relaxation.o: $(SRCDIR)/wake_relaxation.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/wake_relaxation.cpp

file_io.o: $(SRCDIR)/file_io.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/file_io.cpp

panelcode.o: $(SRCDIR)/panelcode.cpp
	$(CXX) -c $(CXXFLAGS) $(PREPROC) $(INCLUDE) $(SRCDIR)/panelcode.cpp
