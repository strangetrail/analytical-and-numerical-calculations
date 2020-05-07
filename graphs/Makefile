
myPATH=./fffc/
include ./fffc/mflibfffc.mk

currSHELL:=/usr/bin/

CXXFLAGS_DEFAULT:=-I/usr/include/freetype2 -I. -I./fffc -g3 -gdwarf-4 -ggdb -Wall -Wpedantic -Wextra
LDLIBS:=-lm -lGL -lX11 -lfreetype -lSDL2 -lGLEW -lGLU

ifneq ($(CXXFLAGS),)
	CXXFLAGS :=$(CXXFLAGS) $(CXXFLAGS_DEFAULT)
else
	CXXFLAGS :=$(CXXFLAGS_DEFAULT)
endif

.PHONY: all clean

all: alllibfffc graph
	
graph: graph.test2.o common/shader_utils.o
	$(currSHELL)g++ $(CXXFLAGS) -o graph.elf graph.test2.o shader_utils.o $(myPATH)fffc.o $(myPATH)libCZTcstm.o $(myPATH)libFFTcstm.o $(myPATH)libUtils.o $(LDLIBS) $(LDFLAGS)

graph.test2.o: graph.test2.cpp
	$(currSHELL)g++ -c $(CXXFLAGS) graph.test2.cpp

common/shader_utils.o: common/shader_utils.cpp
	$(currSHELL)g++ -c $(CXXFLAGS) common/shader_utils.cpp

clean: cleanlibfffc
	$(currSHELL)rm -ifv common/shader_utils.o graph.test2.o shader_utils.o graph.elf *output*txt *Output*txt *.aux *.log *.ps *.dvi *.toc *.los *.pdf *.tex *.pyc *.png *.bak *.so *sympy.c 1
	$(currSHELL)rm -frv build
