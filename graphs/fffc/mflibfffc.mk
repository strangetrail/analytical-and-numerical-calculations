#LDLIBS=-lm
#ifndef $(myPATH)
#	myPATH = 
#endif
###DIRECT = $(myPATH)direct
###INVERSE = $(myPATH)inverse
#FLAGS = -Iinclude -std=c++11
#
#
#MKFLAGS_DEBUG:=-gdwarf-4 -ggdb3 -g3
#MKFLAGS_EXE_HARDENING:=-fpie -Wl,-pie
#
#
##FLAGS = -std=c99
#
#
MKCXXFLAGS_DEFAULT:=-pedantic -Wall -Wextra
MKCXXFLAGS:=-std=c++14 $(CXXFLAGS) $(CXXFLAGS_DEFAULT)
#
#
##CFLAGS = -pedantic -Wall -Wextra
#CFLAGS = -march=native -ggdb3
#DEBUGFLAGS   = -O0 -D _DEBUG
#RELEASEFLAGS = -O2 -D NDEBUG -combine -fwhole-program
#TARGET  = example.so
#SOURCES = $(shell echo src/*.c)
#HEADERS = $(shell echo include/*.h)
#OBJECTS = $(SOURCES:.c=.o)
#PREFIX = $(DESTDIR)/usr/local
#BINDIR = $(PREFIX)/bin
#FLAGS        = # -std=gnu99 -Iinclude
#CFLAGS       = -fPIC -g #-pedantic -Wall -Wextra -ggdb3
#LDFLAGS      = -shared
#all: $(TARGET)
#$(TARGET): $(OBJECTS)
# $(CC) $(FLAGS) $(CFLAGS) $(DEBUGFLAGS) -o $(TARGET) $(OBJECTS)
###alllibfffc: $(DIRECT)FFT.elf $(INVERSE)FFT.elf $(myPATH)libfffc.so
alllibfffc: $(myPATH)fffc.o $(myPATH)libCZTcstm.o $(myPATH)libFFTcstm.o $(myPATH)libUtils.o
###$(DIRECT)FFT.elf: $(DIRECT)FFT.o
###	g++ -o $(DIRECT)FFT.elf $(DIRECT)FFT.o
###$(DIRECT)FFT.o: $(DIRECT)FFT.cpp
###	g++ -c -g $(DIRECT)FFT.cpp -o $(DIRECT)FFT.o
###$(INVERSE)FFT.elf: $(INVERSE)FFT.o
###	g++ -o $(INVERSE)FFT.elf $(INVERSE)FFT.o
###$(INVERSE)FFT.o: $(INVERSE)FFT.cpp
###	g++ -c -g $(INVERSE)FFT.cpp -o $(INVERSE)FFT.o
##$(myPATH)libfffc.so: $(myPATH)fffc.o $(myPATH)libCZTcstm.o $(myPATH)libFFTcstm.o
##	g++ $(FLAGS) -shared -o $(myPATH)libfffc.so $(myPATH)fffc.o $(myPATH)libCZTcstm.o $(myPATH)libFFTcstm.o
$(myPATH)fffc.o: $(myPATH)fffc.cpp
	g++ -c $(MKCXXFLAGS) $(myPATH)fffc.cpp -o $(myPATH)fffc.o
$(myPATH)libCZTcstm.o: $(myPATH)libCZTcstm.cpp
	g++ -c $(MKCXXFLAGS) $(myPATH)libCZTcstm.cpp -o $(myPATH)libCZTcstm.o
$(myPATH)libFFTcstm.o: $(myPATH)libFFTcstm.cpp
	g++ -c $(MKCXXFLAGS) $(myPATH)libFFTcstm.cpp -o $(myPATH)libFFTcstm.o
$(myPATH)libUtils.o: $(myPATH)libUtils.cpp
	g++ -c $(MKCXXFLAGS) $(myPATH)libUtils.cpp -o $(myPATH)libUtils.o
cleanlibfffc:
	rm -vf $(DIRECT)FFT.elf $(DIRECT)FFT.o $(INVERSE)FFT.elf $(INVERSE)FFT.o $(myPATH)fffc.elf $(myPATH)fffc.o $(myPATH)*output*.* $(myPATH)check.m $(myPATH)res1.txt $(myPATH)res2.txt $(myPATH)libfffc.so $(myPATH)fffc.o $(myPATH)libCZTcstm.o $(myPATH)libFFTcstm.o $(myPATH)libUtils.o
