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
# Compiler and flags
CXX = g++
CXXFLAGS= -O2 -g3 -fPIC

##==========================================================================
PROGRAM := PopGenLib

HEADER_GENERIC = popgen.h
HEADER_LOWD = popgen_lowd.h
HEADER_HIGHD = popgen_highd.h

SOURCE_GENERIC = sample
SOURCE_LOWD = hypercube haploid_gt_dis
SOURCE_HIGHD = hypercube_function haploid_clone

##==========================================================================
# Rules for generating the final executable.
$(PROGRAM): $(SOURCE_LOWD:%=%.o) $(SOURCE_HIGHD:%=%.o)
	ar rcs libPopGenLib.a $^

##==========================================================================
# Rules for generating object files.
$(SOURCE_LOWD:%=%.o): $(SOURCE_GENERIC:%=%.cpp) $(SOURCE_LOWD:%=%.cpp) $(HEADER_LOWD) $(HEADER_GENERIC)
	$(CXX) ${CXXFLAGS} -o $@ -c $(SOURCE_GENERIC:%=%.cpp) $(SOURCE_LOWD:%=%.cpp)

$(SOURCE_HIGHD:%=%.o): $(SOURCE_GENERIC:%=%.cpp) $(SOURCE_HIGHD:%=%.cpp) $(HEADER_HIGHD) $(HEADER_GENERIC)
	$(CXX) ${CXXFLAGS} -o $@ -c $(SOURCE_GENERIC:%=%.cpp) $(SOURCE_HIGHD:%=%.cpp)

#
#$(SOURCE_GENERIC:%=%.o) $(SOURCE_LOWD:%=%.o) $(SOURCE_HIGHD:%=%.o): $(SOURCE_GENERIC:%=%.cpp) $(SOURCE_LOWD:%=%.cpp) $(SOURCE_HIGHD:%=%.cpp) $(HEADER_LOWD) $(HEADER_HIGHD) $(HEADER_GENERIC)
#	$(CXX) ${CXXFLAGS} -c $^
#
##==========================================================================
# Rules for cleaning.
clean:
	rm -rf *.a *.o

#############################################################################
