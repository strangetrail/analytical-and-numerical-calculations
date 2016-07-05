myPATH=./
myLibPATH=../fffc/
LDLIBS=-lm
currSHELL=/usr/bin/
FLAGS = -std=c++11 -gdwarf-4 -O0
all: parallelepiped
parallelepiped: parallelepiped.o $(myLibPATH)libUtils.o
	$(currSHELL)g++ $(FLAGS) -o parallelepiped.elf parallelepiped.o $(myLibPATH)libUtils.o $(LDLIBS)
parallelepiped.o: parallelepiped.cpp
	$(currSHELL)g++ -c -g3 $(FLAGS) -pedantic -Wall -Wextra -ggdb3 parallelepiped.cpp
$(myLibPATH)libUtils.o: $(myLibPATH)libUtils.cpp
	g++ -c -g3 $(FLAGS) $(CFLAGS) $(myLibPATH)libUtils.cpp -o $(myLibPATH)libUtils.o
.PHONY: all clean
clean:
	$(currSHELL)rm -ifv parallelepiped.o parallelepiped.elf
	$(currSHELL)rm -ifv $(myLibPATH)libUtils.o
