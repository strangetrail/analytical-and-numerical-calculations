CXX=g++

CXXFLAGS= -D_FORTIFY_SOURCE=2 \
          -D_GLIBCXX_ASSERTIONS \
          -fasynchronous-unwind-tables \
          -fexceptions \
          -fpie -Wl,-pie\
          -fplugin=annobin \
          -fstack-clash-protection \
          -fstack-protector-strong \
          -fcf-protection \
          -O3 \
          -pipe \
          -Wall -Werror=format-security \
          -specs=/usr/lib/rpm/redhat/redhat-hardened-cc1 \
          -specs=/usr/lib/rpm/redhat/redhat-annobin-cc1 \
          -Wl,-z,defs \
          -Wl,-z,now \
          -Wl,-z,relro \
          -Wextra \
          -Wpedantic
_CXXFLAGS= -O0 \
          -Wall \
          -Wextra \
          -Wpedantic \
          -g3 -ggdb -ggdb3 -gdwarf-4

LDFLAGS= -fpie -Wl,-pie\
         -specs=/usr/lib/rpm/redhat/redhat-hardened-ld \
         -Wl,-z,defs \
         -Wl,-z,now \
         -Wl,-z,relro
_LDFLAGS=

LDLIBS= -lm \
        -lGL \
        -lX11 \
        -lfreetype \
        -lSDL2 \

OBJECTS= graph.o fffc.o libCZTcstm.o libFFTcstm.o libUtils.o shader_utils.o

.PHONY: all clean
.PRECIOUS: %.o

all: graph.elf

%.elf: $(OBJECTS)
	$(CXX) -std=c++17 $(LDFLAGS) -o $@ $(OBJECTS) $(LDLIBS)

%.o: %.cpp
	$(CXX) -std=c++17 $(CXXFLAGS) -I. -o $@ -c $<

clean:
	rm -f *.{elf,o}
