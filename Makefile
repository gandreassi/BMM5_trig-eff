CXX=g++
CXXFLAGS=-O3 -fPIC -Wall

ROOT_FLAGS=$(shell $(ROOTSYS)/bin/root-config --cflags)
ROOT_INCLUDE=$(shell $(ROOTSYS)/bin/root-config  --incdir)
FITTER_INCLUDE=./src/

BUILDDIR=.
OBJ_DIR=.

INCLUDE=-I $(ROOT_INCLUDE) -I $(FITTER_INCLUDE)

ROOT_LIBS    = $(shell $(ROOTSYS)/bin/root-config --libs) -lTreePlayer -lMinuit -lXMLIO -lMLP -lRIO -lTMVA 
GLIBS := $(shell root-config --glibs)
FLAGS = $(CXXFLAGS)
FLAGS += ${ROOT_FLAGS}
FLAGS_ROOFIT = $(FLAGS)
FLAGS_ROOFIT += -lRooFitCore -lRooFit -lRooStats -lFoam


MACRO=fitter
EXE=trigeff

all: $(MACRO)

$(MACRO): main.cpp
	@echo "---> Making s_plot..."
	$(CXX) $(CXXFLAGS) $(GLIBS) $(FLAGS_ROOFIT) $(ROOT_LIBS) $(INCLUDE) $^ -o $(EXE).exe


clean:
	rm -f $(BUILDDIR)/$(EXE).exe
