#############################################################################
#
# Licence:	
# Author:	Richard Neher, Fabio Zanini
# Date:		2012/05/15
#
# Description:
# ------------
# Makefile for the FFPopSim library.
#
# The first section deals with platform-specific options and paths.
# In particular, there are two variables that can differ across the audience:
#
# 1. PYTHON_PREFIX: the root of your Python installation. Two common choices
#                   are listed, but you should check your own path. One way
#                   of doing this is to run:
#
#                   which python
#
#                   in a shell, and copy the first chunk of the path.
#
# 2. PYTHON_LD_FLAGS_PLATFORM: the stategy used to build shared libraries.
#                   There are only two options here, depending on whether you
#                   are on an Apple computer or not.
#
# This file tries to guess reasonable values for both variables, by including
# Makefile_guesses. However, if you know their value on your system, please
# set it explicitely below.
#
# There are four main recipes:
# - src: library compilation and (static) linking
# - doc: documentation
# - tests: test cases compilation and linking against the library
# - python: python bindings
#
##==========================================================================
# PLATFORM-DEPENDENT OPTIONS
#
# Please edit the following lines to your needs. If you do not have doxygen
# or python, just leave that line untouched or delete it.
##==========================================================================
CC := gcc
CXX := g++
SWIG := swig
DOXY := doxygen
PYTHON := python2.7

# Try to guess the platform and whether Doxygen and/or Python are installed
include Makefile_guesses

# Apple users would usually use (and possibly edit) the first line
# Linux users the second line
#PYTHON_PREFIX := /Library/Frameworks/EPD64.framework/Versions/Current
#PYTHON_PREFIX := /usr

# Apple users would usually use (and possibly edit) the first line
# Linux users the second line
#PYTHON_LD_FLAGS_PLATFORM := -dynamiclib -flat_namespace -undefined suppress
#PYTHON_LD_FLAGS_PLATFORM := -shared

# PYTHON_PREFIX is actually only used in the following two variables.
# You can also directly edit these if you know what you are doing.
PYTHON_INCLUDES = $(PYTHON_PREFIX)/include/$(PYTHON)
NUMPY_INCLUDES = $(PYTHON_PREFIX)/lib/$(PYTHON)/site-packages/numpy/core/include

############################################################################
# !! DO NOT EDIT BELOW THIS LINE !!
############################################################################
##==========================================================================
# OVERVIEW
##==========================================================================
SRCDIR = src
DOCDIR = doc
TESTSDIR = tests
PYBDIR = $(SRCDIR)/python
PKGDIR = pkg
PFLDIR = profile

.PHONY : all clean src tests doc python profile clean-src clean-doc clean-tests clean-python clean-profile
all: src tests $(doc) $(python)
clean: clean-src clean-doc clean-tests clean-python clean-profile

# Profile flag to enable profiling with gprof.
# (Un)Comment the next line to switch off (on) profiling.
PROFILEFLAGS := -pg

##==========================================================================
# SOURCE
##==========================================================================
SRC_CXXFLAGS= -O2 -fPIC $(PROFILEFLAGS)

LIBRARY := libFFPopSim.a

HEADER_GENERIC = ffpopsim_generic.h
SOURCE_GENERIC = sample.cpp
OBJECT_GENERIC = $(SOURCE_GENERIC:%.cpp=%.o)

HEADER_LOWD = $(HEADER_GENERIC) ffpopsim_lowd.h
SOURCE_LOWD = hypercube_lowd.cpp haploid_lowd.cpp
OBJECT_LOWD = $(SOURCE_LOWD:%.cpp=%.o)

HEADER_HIGHD = $(HEADER_GENERIC) ffpopsim_highd.h
SOURCE_HIGHD = hypercube_highd.cpp haploid_highd.cpp
OBJECT_HIGHD = $(SOURCE_HIGHD:%.cpp=%.o)

HEADER_HIV = hivpopulation.h
SOURCE_HIV = $(HEADER_HIV:%.h=%.cpp)
OBJECT_HIV = $(SOURCE_HIV:%.cpp=%.o)
SOURCE_HIVGENE = hivgene.cpp
OBJECT_HIVGENE = $(SOURCE_HIVGENE:%.cpp=%.o)

OBJECTS = $(OBJECT_GENERIC) $(OBJECT_LOWD) $(OBJECT_HIGHD) $(OBJECT_HIV) $(OBJECT_HIVGENE)

# Recipes
src: $(SRCDIR)/$(LIBRARY)

$(SRCDIR)/$(LIBRARY): $(OBJECTS:%=$(SRCDIR)/%)
	ar rcs $@ $^
	mkdir -p $(PKGDIR)/lib
	cp $@ $(PKGDIR)/lib/
	mkdir -p $(PKGDIR)/include
	cp $(HEADER_GENERIC:%=$(SRCDIR)/%) $(PKGDIR)/include/
	cp $(HEADER_LOWD:%=$(SRCDIR)/%) $(PKGDIR)/include/
	cp $(HEADER_HIGHD:%=$(SRCDIR)/%) $(PKGDIR)/include/
	cp $(HEADER_HIV:%=$(SRCDIR)/%) $(PKGDIR)/include/

$(OBJECT_GENERIC:%=$(SRCDIR)/%): $(SOURCE_GENERIC:%=$(SRCDIR)/%)
	$(CXX) $(SRC_CXXFLAGS) -c -o $@ $(@:.o=.cpp)

$(OBJECT_LOWD:%=$(SRCDIR)/%): $(SOURCE_LOWD:%=$(SRCDIR)/%) $(HEADER_LOWD:%=$(SRCDIR)/%)
	$(CXX) $(SRC_CXXFLAGS) -c -o $@ $(@:.o=.cpp)

$(OBJECT_HIGHD:%=$(SRCDIR)/%): $(SOURCE_HIGHD:%=$(SRCDIR)/%) $(HEADER_HIGHD:%=$(SRCDIR)/%)
	$(CXX) $(SRC_CXXFLAGS) -c -o $@ $(@:.o=.cpp)

$(OBJECT_HIV:%=$(SRCDIR)/%): $(SOURCE_HIV:%=$(SRCDIR)/%) $(HEADER_HIV:%=$(SRCDIR)/%)
	$(CXX) $(SRC_CXXFLAGS) -c -o $@ $(@:.o=.cpp)

$(OBJECT_HIVGENE:%=$(SRCDIR)/%): $(SOURCE_HIVGENE:%=$(SRCDIR)/%) $(HEADER_HIV:%=$(SRCDIR)/%)
	$(CXX) $(SRC_CXXFLAGS) -c -o $@ $(@:.o=.cpp)

clean-src:
	cd $(SRCDIR); rm -rf $(LIBRARY) *.o *.h.gch
	cd $(PKGDIR); rm -rf lib include

##==========================================================================
# DOCUMENTATION
##==========================================================================
DOXYFILE   = $(DOCDIR)/Doxyfile

# Recipes
doc:
	$(DOXY) $(DOXYFILE)

clean_doc:
	rm -rf $(DOCDIR)/latex $(DOCDIR)/html

##==========================================================================
# TESTS
##==========================================================================
TESTS_CXXFLAGS = -I$(SRCDIR) -Wall -O2 -c -fPIC
TESTS_LDFLAGS = -O2 $(PROFILEFLAGS)
TEST_LIBDIRS = -L$(CURDIR)/$(SRCDIR)
TESTS_LIBS = -lFFPopSim -lgsl -lgslcblas

TESTS_LOWD = lowd
TESTS_HIGHD = highd

TESTS_SOURCE_LOWD = $(TESTS_LOWD:%=%.cpp)
TESTS_SOURCE_HIGHD = $(TESTS_HIGHD:%=%.cpp)

TESTS_HEADER_LOWD = $(TESTS_LOWD:%=%.h)
TESTS_HEADER_HIGHD = $(TESTS_HIGHD:%=%.h)

TESTS_OBJECT_LOWD = $(TESTS_LOWD:%=%.o)
TESTS_OBJECT_HIGHD = $(TESTS_HIGHD:%=%.o)

# Recipes
tests: $(SRCDIR)/$(LIBRARY) $(TESTS_LOWD:%=$(TESTSDIR)/%) $(TESTS_HIGHD:%=$(TESTSDIR)/%)

$(TESTS_LOWD:%=$(TESTSDIR)/%): $(TESTS_OBJECT_LOWD:%=$(TESTSDIR)/%) $(SRCDIR)/$(LIBRARY)
	$(CXX) $(TESTS_LDFLAGS) $^ $(TEST_LIBDIRS) $(TESTS_LIBS) -o $@

$(TESTS_OBJECT_LOWD:%=$(TESTSDIR)/%): $(TESTS_SOURCE_LOWD:%=$(TESTSDIR)/%) $(TESTS_HEADER_LOWD:%=$(TESTSDIR)/%)
	$(CXX) $(TESTS_CXXFLAGS) -c $(@:.o=.cpp) -o $@

$(TESTS_HIGHD:%=$(TESTSDIR)/%): $(TESTS_OBJECT_HIGHD:%=$(TESTSDIR)/%) $(OBJECT_HIV:%=$(SRCDIR)/%) $(SRCDIR)/$(LIBRARY)
	$(CXX) $(TESTS_LDFLAGS) $^ $(TEST_LIBDIRS) $(TESTS_LIBS) -o $@

$(TESTS_OBJECT_HIGHD:%=$(TESTSDIR)/%): $(TESTS_SOURCE_HIGHD:%=$(TESTSDIR)/%) $(TESTS_HEADER_HIGHD:%=$(TESTSDIR)/%)
	$(CXX) $(TESTS_CXXFLAGS) -c $(@:.o=.cpp) -o $@

clean-tests:
	cd $(TESTSDIR); rm -rf *.o $(TESTS_LOWD) $(TESTS_HIGHD)

##==========================================================================
# PROFILE
##==========================================================================
PROFILE_CXXFLAGS = -I$(SRCDIR) -Wall -O2 -c -fPIC $(PROFILEFLAGS)
PROFILE_LDFLAGS = -O2 $(PROFILEFLAGS)
PROFILE_LIBDIRS = -L$(CURDIR)/$(SRCDIR)
PROFILE_LIBS = -lFFPopSim -lgsl -lgslcblas

PROFILE = profile
PROFILE_SOURCE = $(PROFILE:%=%.cpp)
PROFILE_OBJECT = $(PROFILE:%=%.o)

# Recipes
profile: $(SRCDIR)/$(LIBRARY) $(PROFILE:%=$(PFLDIR)/%)

$(PROFILE:%=$(PFLDIR)/%): $(PROFILE_OBJECT:%=$(PFLDIR)/%) $(SRCDIR)/$(LIBRARY)
	$(CXX) $(PROFILE_LDFLAGS) $^ $(PROFILE_LIBDIRS) $(PROFILE_LIBS) -o $@

$(PROFILE_OBJECT:%=$(PFLDIR)/%): $(PROFILE_SOURCE:%=$(PFLDIR)/%)
	$(CXX) $(PROFILE_CXXFLAGS) -c $(@:.o=.cpp) -o $@

##==========================================================================
# PYTHON BINDINGS
##==========================================================================
SWIGFLAGS = -c++ -python -O -castmode -keyword

SWIG_INTERFACE = FFPopSim.i
SWIG_WRAP = $(SWIG_INTERFACE:%.i=%_wrap.cpp)
SWIG_WRAP_OBJECT = $(SWIG_WRAP:%.cpp=%.o)
SWIG_OBJECT = $(SWIG_INTERFACE:%.i=_%.so)
SWIG_PYMODULE = $(SWIG_INTERFACE:%.i=%.py)
SWIG_PYCMODULE = $(SWIG_INTERFACE:%.i=%.pyc)
SWIG_SUPPORT_1 = ffpopsim_highd.i
SWIG_SUPPORT_2 = hivpopulation.i
SWIG_SUPPORT_3 = ffpopsim_lowd.i
SWIG_SUPPORT_4 = ffpopsim_generic.i

PYTHON_CFLAGS = -O2 -fPIC -I$(SRCDIR) -I$(PYTHON_INCLUDES) -I$(NUMPY_INCLUDES)
PYTHON_LDFLAGS= -O2 -fPIC $(PYTHON_LD_FLAGS_PLATFORM)
PYTHON_LIBDIRS = -L$(CURDIR)/$(SRCDIR)
PYTHON_LIBS = -lFFPopSim -lgsl -lgslcblas

# Recipes
python: $(SWIG_PYMODULE:%=$(PYBDIR)/%) $(SWIG_PYMCODULE:%=$(PYBDIR)/%) $(SWIG_OBJECT:%=$(PYBDIR)/%) 

$(SWIG_PYMODULE:%=$(PYBDIR)/%) $(SWIG_PYMCODULE:%=$(PYBDIR)/%) $(SWIG_OBJECT:%=$(PYBDIR)/%): $(SWIG_WRAP:%=$(PYBDIR)/%) $(SWIG_SOURCE:%=$(PYBDIR)/%) $(DISTUTILS_SETUP:%=$(PYBDIR)/%) $(SRCDIR)/$(LIBRARY)
	$(CC) $(PYTHON_CFLAGS) -c $(SWIG_WRAP:%=$(PYBDIR)/%) -o $(SWIG_WRAP_OBJECT:%=$(PYBDIR)/%)
	$(CXX) $(PYTHON_LDFLAGS) $(SWIG_WRAP_OBJECT:%=$(PYBDIR)/%) $(PYTHON_LIBDIRS) $(PYTHON_LIBS) -o $(SWIG_OBJECT:%=$(PYBDIR)/%)
	mkdir -p $(PKGDIR)/python
	cp -f $(SWIG_PYMODULE:%=$(PYBDIR)/%) $(PKGDIR)/python/
	cp -f $(SWIG_OBJECT:%=$(PYBDIR)/%) $(PKGDIR)/python/

$(SWIG_WRAP:%=$(PYBDIR)/%): $(SWIG_HEADER_HIV:%=$(PYBDIR)/%) $(SWIG_INTERFACE:%=$(PYBDIR)/%) $(SRCDIR)/$(LIBRARY) $(SWIG_SUPPORT_1:%=$(PYBDIR)/%) $(SWIG_SUPPORT_2:%=$(PYBDIR)/%) $(SWIG_SUPPORT_3:%=$(PYBDIR)/%) $(SWIG_SUPPORT_4:%=$(PYBDIR)/%)
	$(SWIG) $(SWIGFLAGS) -o $@ $(SWIG_INTERFACE:%=$(PYBDIR)/%)

clean-python:
	cd $(PYBDIR); rm -rf $(SWIG_WRAP) $(SWIG_OBJECT) $(SWIG_WRAP_OBJECT) $(SWIG_PYMODULE) $(SWIG_PYCMODULE)
	cd $(TESTSDIR); rm -rf $(SWIG_OBJECT) $(SWIG_PYMODULE) $(SWIG_PYCMODULE)
	cd $(PKGDIR)/python; rm -rf $(SWIG_OBJECT) $(SWIG_PYMODULE) $(SWIG_PYCMODULE)

#############################################################################
