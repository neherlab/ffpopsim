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
# The first section deals with platform-specific programs, options, and paths.
# If you are only compiling the C++ library and do not need the Python
# bindings, please comment the lines that define the PYTHON, SWIG, and SPHINX
# variables.
#
# The second section of this Makefile, below the clause within !!, is where
# the make recipes are listed and specified. Modify that part of the file
# only if you know what you are doing and want to change some technical detail
# of the dependecy chain.
#
# There are five main recipes:
# - src: library compilation and (static) linking
# - doc: documentation
# - tests: test cases compilation and linking against the library
# - python: python bindings
# - python-doc: python documentation
#
# The compiled library, the C++ include files, the Python bindings, and the
# documentation are put into the pkg folder after building.
#
#############################################################################

##==========================================================================
# PLATFORM-DEPENDENT OPTIONS
#
# Please edit the following lines to your needs. If you do not have doxygen
# or python, just leave that line untouched or delete it.
##==========================================================================
CC := gcc
CXX := g++
DOXY := doxygen
# Comment the following lines to avoid building the Python 2.7 bindings
PYTHON := python2.7
SWIG := swig
SPHINX := sphinx-build2

# Lower this number if you prefer to use mildly optimized code only.
OPTIMIZATION_LEVEL := 1

############################################################################
# !! DO NOT EDIT BELOW THIS LINE !!
############################################################################
##==========================================================================
# OVERVIEW
##==========================================================================
SRCDIR = src
DOCDIR = doc
TESTSDIR = tests
PYBDIR := $(SRCDIR)/python
PYDOCDIR := $(DOCDIR)/python
PKGDIR = pkg
PFLDIR = profile
DISTUTILS_SETUP := setup.py

# Can we compile Python bindings?
ifdef PYTHON
    python := python
endif

# List all explicit recipes
.PHONY : all src tests doc python python-doc profile swig clean clean-all clean-src clean-doc clean-tests clean-python clean-python-doc clean-profile clean-swig
all: src tests $(python)
clean: clean-src clean-tests clean-python clean-profile
clean-all: clean clean-doc clean-python-doc clean-swig

# Profile flag to enable profiling with gprof.
# (Un)Comment the next line to switch off (on) profiling.
#PROFILEFLAGS := -pg

##==========================================================================
# SOURCE
##==========================================================================
SRC_CXXFLAGS= -O$(OPTIMIZATION_LEVEL) -fPIC $(PROFILEFLAGS)

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
DOXYFILE   = $(DOCDIR)/cpp/Doxyfile

# Recipes
doc:
	$(DOXY) $(DOXYFILE)
	mkdir -p $(PKGDIR)/doc/cpp
	cp -rf $(DOCDIR)/cpp/html $(PKGDIR)/doc/cpp/

clean-doc:
	cd $(DOCDIR)/cpp; rm -rf latex html
	cd $(PKGDIR)/doc; rm -rf cpp

##==========================================================================
# TESTS
##==========================================================================
TESTS_CXXFLAGS = -I$(SRCDIR) -Wall -O$(OPTIMIZATION_LEVEL) -c -fPIC
TESTS_LDFLAGS = -O$(OPTIMIZATION_LEVEL) $(PROFILEFLAGS)
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
PROFILE_CXXFLAGS = -I$(SRCDIR) -Wall -O$(OPTIMIZATION_LEVEL) -c -fPIC $(PROFILEFLAGS)
PROFILE_LDFLAGS = -O$(OPTIMIZATION_LEVEL) $(PROFILEFLAGS)
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
# PYTHON BINDINGS AND DOCUMENTATION
##==========================================================================
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

# Recipes
python: $(SWIG_PYMODULE:%=$(PYBDIR)/%) $(SWIG_PYMCODULE:%=$(PYBDIR)/%) $(SWIG_OBJECT:%=$(PYBDIR)/%) $(DISTUTILS_SETUP)

$(SWIG_OBJECT:%=$(PYBDIR)/%): $(SWIG_WRAP:%=$(PYBDIR)/%) $(SWIG_PYMODULE:%=$(PYBDIR)/%) $(SOURCE_LOWD:%=$(SRCDIR)/%) $(SOURCE_HIGHD:%=$(SRCDIR)/%) $(SOURCE_HIV:%=$(SRCDIR)/%) $(SOURCE_HIVGENE:%=$(SRCDIR)/%) $(SOURCE_GENERIC:%=$(SRCDIR)/%)
	CFLAGS='-O$(OPTIMIZATION_LEVEL)' $(PYTHON) setup.py build_ext --inplace
	rm -rf build
	cp -f $(SWIG_PYMODULE:%=$(PYBDIR)/%) $(PKGDIR)/python/
	cp -f $(SWIG_OBJECT:%=$(PYBDIR)/%) $(PKGDIR)/python/

clean-python:
	cd $(PYBDIR); rm -rf $(SWIG_OBJECT) $(SWIG_WRAP_OBJECT) $(SWIG_PYCMODULE)
	cd $(TESTSDIR); rm -rf $(SWIG_OBJECT) $(SWIG_PYMODULE) $(SWIG_PYCMODULE)
	cd $(PKGDIR)/python; rm -rf $(SWIG_OBJECT) $(SWIG_PYMODULE) $(SWIG_PYCMODULE)

python-doc:
	cd $(PYDOCDIR); $(MAKE) SPHINXBUILD=$(SPHINX) html
	mkdir -p $(PKGDIR)/doc/python
	cp -rf $(PYDOCDIR)/build/html $(PKGDIR)/doc/python/

clean-python-doc:
	cd $(PYDOCDIR); rm -rf build
	cd $(PKGDIR)/doc; rm -rf python

##==========================================================================
# SWIG (USED FOR PYTHON BINDINGS)
##==========================================================================
SWIGFLAGS = -c++ -python -O -castmode -keyword

swig: $(SWIG_WRAP:%=$(PYBDIR)/%) $(SWIG_PYMODULE:%=$(PYBDIR)/%)

$(SWIG_WRAP:%=$(PYBDIR)/%) $(SWIG_PYMODULE:%=$(PYBDIR)/%): $(SWIG_HEADER_HIV:%=$(PYBDIR)/%) $(SWIG_INTERFACE:%=$(PYBDIR)/%) $(SWIG_SUPPORT_1:%=$(PYBDIR)/%) $(SWIG_SUPPORT_2:%=$(PYBDIR)/%) $(SWIG_SUPPORT_3:%=$(PYBDIR)/%) $(SWIG_SUPPORT_4:%=$(PYBDIR)/%)
	$(SWIG) $(SWIGFLAGS) -o $(SWIG_WRAP:%=$(PYBDIR)/%) $(SWIG_INTERFACE:%=$(PYBDIR)/%)

clean-swig:
	cd $(PYBDIR); rm -rf $(SWIG_WRAP) $(SWIG_PYMODULE)
#############################################################################
