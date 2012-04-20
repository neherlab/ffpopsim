#############################################################################
#
# Licence:	
# Author:	Richard Neher, Boris Shraiman, Fabio Zanini
# Date:		2012/04/19
#
# Description:
# ------------
# Makefile for the PopGenLib library.
#
##==========================================================================
CXX=g++
libraries =  -lgsl -lgslcblas -lHandyTools -lPopGenLib -lboost_system -lboost_filesystem

# FIXME: include handytools into popgenlib (it's basically a few functions)
rlibpath = ../../

LDFLAGS = -L$(rlibpath)PopGenLib/src/ -L$(rlibpath)HandyTools/src/ -L/usr/ -O2
CXXFLAGS = -I$(rlibpath)PopGenLib/src/ -I$(rlibpath)HandyTools/src/ -Wall -O2 -c -fPIC -g3
##==========================================================================
TESTS_LOWD = tests_lowd
TESTS_HIGHD = tests_highd

TESTS = $(TESTS_LOWD) $(TESTS_HIGHD)
##==========================================================================
$(TESTS_HIGHD): $(TESTS_HIGHD:%=%.o) 
	${CXX} ${LDFLAGS} -o $@ $^ ${libraries} 

$(TESTS_HIGHD:%=%.o): $(TESTS_HIGHD:%=%.cpp) 
	${CXX} ${CXXFLAGS} -o $@ $<

clean:
	rm -rf *.o $(TESTS)

#############################################################################
