FLAGS = -std=c++11 -gdwarf-4 -O0
CFLAGS = -pedantic -Wall -Wextra -Iinclude -fPIC -ggdb3
alllibutils: $(myPATH)libutils.so
$(myPATH)libutils.so: $(myPATH)utils.o
	g++ $(FLAGS) -shared -o $(myPATH)libutils.so $(myPATH)utils.o
$(myPATH)utils.o: $(myPATH)libUtils.cpp
	g++ -c -g3 $(FLAGS) $(CFLAGS) $(myPATH)libUtils.cpp -o $(myPATH)utils.o
cleanlibutils:
	rm -vf $(myPATH)libutils.so $(myPATH)utils.o
