PATH2SRCS := src
PATH2INCS := inc

RM := rm -ifv
CC := gcc
CXX := g++

INCLUDES := -Iinclude -Iinc

DEBUGFLAGS := -gdwarf-4 -ggdb3

FLAGS := -pedantic -Wall -Wextra  # -Wno-write-strings
CPPFLAGS := -O0 $(INCLUDES)  # -O3 -Ofast -fno-omit-frame-pointer -fno-align-functions -fno-align-loops  # -S -O2
CFLAGS := -std=c99 -fpic
CXXFLAGS := -std=c++11 -fPIC

LDLIBS := -lgsl -lgslcblas -lm

all: quadpack_test

quadpack_test: quadpack_test.o
	$(CC) -shared $(FLAGS) -o $@ $+ $(LDLIBS)
quadpack_test.o: src/quadpack_test.c
	$(CC) -c -g3 $(FLAGS) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

.PHONY: all clean

clean:
	$(RM) *.aux *.log *.ps *.dvi *.toc *.los *.pdf *.tex *.png *.bak *.LISP *.UNLISP *.fasl 1 *.elf *.so *.o
