myPATH=./fffc/
include ./fffc/mflibfffc.mk
#LDLIBS=-lm -lglut -lGLEW -lGL -lGLU -lXmu -lX11 -L/usr/lib/libGLEW.so
LDLIBS=-lm -lglut -lGLEW -lGL -lGLU -lX11 -L/usr/lib/libGLEW.so -lfreetype
CXXFLAGS=-I/usr/include/freetype2
#
#
#CSTMFLAGS=-g3 -ggdb3
#
#
currSHELL=/usr/bin/
all: alllibfffc graph
graph: graph.o common/shader_utils.o
	$(currSHELL)g++ $(FLAGS) -o graph.elf graph.o shader_utils.o $(myPATH)fffc.o $(myPATH)libCZTcstm.o $(myPATH)libFFTcstm.o $(myPATH)libUtils.o $(LDLIBS)
graph.o: graph.cpp
	$(currSHELL)g++ -c -g3 $(FLAGS) $(CXXFLAGS) -pedantic -Wall -Wextra -ggdb3 graph.cpp
common/shader_utils.o: common/shader_utils.cpp
	$(currSHELL)g++ -c -g3 $(FLAGS) -pedantic -Wall -Wextra -ggdb3 common/shader_utils.cpp
.PHONY: all clean
clean: cleanlibfffc
	$(currSHELL)rm -ifv common/shader_utils.o graph.o shader_utils.o graph.elf *output*txt *Output*txt *.aux *.log *.ps *.dvi *.toc *.los *.pdf *.tex *.pyc *.png *.bak *.so *sympy.c 1
	$(currSHELL)rm -frv build
