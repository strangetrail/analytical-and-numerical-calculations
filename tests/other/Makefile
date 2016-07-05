currSHELL=/usr/bin/
FLAGS = -gdwarf-4
all: com
com: com.o
	$(currSHELL)g++ $(FLAGS) -o com.elf com.o
com.o: com.cpp
	$(currSHELL)g++ -c -g3 $(FLAGS) -pedantic -Wall -Wextra -ggdb3 com.cpp
.PHONY: all clean
clean:
	$(currSHELL)rm -ifv com.o com.elf
