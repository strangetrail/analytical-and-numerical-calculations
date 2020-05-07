include ../../common/mfcommon.mk


FLAGS    += -fPIC $(WARNINGFLAGS)
CPPFLAGS += -O3 -Ofast $(INCLUDES)


LDLIBS :=-lgsl -lgslcblas -lm


ifeq ($(codesrcsuffix),)
	codesrcsuffix :=__export
endif


all: integrand_export.hpp libquadqags2d.so

integrand_export.hpp: $(PATH2INCS)/$(codesrcname)$(codesrcsuffix).hpp
	$(LN) -f ./$(codesrcname)$(codesrcsuffix).hpp $(PATH2INCS)/$@

libquadqags2d.so: quadqags2d_internals.o           \
                  quadqags2d.o                     \
                  $(codesrcname)$(codesrcsuffix).o  
	$(CXX) -shared $(FLAGS) $(CXXFLAGS) -o $@ $+ $(LDLIBS)

quadqags2d.o: $(PATH2SRCS)/quadqags2d.cpp
	$(CXX) -c -g3 $(FLAGS) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

quadqags2d_internals.o: $(PATH2SRCS)/quadqags2d_internals.cpp
	$(CXX) -c -g3 $(FLAGS) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

$(codesrcname)$(codesrcsuffix).o: $(PATH2SRCS)/$(codesrcname)$(codesrcsuffix).cpp
	$(CXX) -c -g3 $(FLAGS) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

clean: clean_common clean_quadqags2d

clean_quadqags2d:
	$(RM) $(PATH2INCS)/*$(codesrcsuffix).hpp $(PATH2SRCS)/*$(codesrcsuffix).cpp
