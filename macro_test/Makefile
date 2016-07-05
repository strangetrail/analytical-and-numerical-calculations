LDLIBS=-lm
FLAGS = -std=c++11 -gdwarf-4 -O0
currSHELL=/usr/bin/
all: test
test: test.o
	$(currSHELL)g++ $(FLAGS) -o test.elf test.o $(LDLIBS)
test.o: test.cpp
	$(currSHELL)g++ -c -g3 $(FLAGS) -pedantic -Wall -Wextra -ggdb3 test.cpp
.PHONY: all clean
clean:
	$(currSHELL)rm -ifv test.o test.elf
