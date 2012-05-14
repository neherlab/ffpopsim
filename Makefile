#############################################################################
#
# Licence:	
# Author:	Richard Neher, Fabio Zanini
# Date:		2012/04/19
#
# Description:
# ------------
# Makefile for the PopGenLib library.
#
# There are four sections:
# - src: library compilation and (static) linking
# - doc: documentation
# - tests: test cases compilation and linking against the library
# - python: python bindings
#
##==========================================================================
SRCDIR = src
DOCDIR = doc
TESTSDIR = tests
PYBDIR = $(SRCDIR)/python

# System folders
LINKER_FOLDER = /usr
NUMPY_INCLUDES = /usr/lib/python2.7/site-packages/numpy/core/include

.PHONY : all clean src doc tests python clean-src clean-doc clean-tests clean-python
all: src doc tests python
clean: clean-src clean-doc clean-tests clean-python

##==========================================================================
# SOURCE
##==========================================================================
CXX = g++
CXXFLAGS= -O2 -g3 -fPIC

LIBRARY := libPopGenLib.a

HEADER_GENERIC = popgen.h
SOURCE_GENERIC = sample.cpp
OBJECT_GENERIC = $(SOURCE_GENERIC:%.cpp=%.o)

HEADER_LOWD = $(HEADER_GENERIC) popgen_lowd.h
SOURCE_LOWD = hypercube.cpp haploid_gt_dis.cpp
OBJECT_LOWD = $(SOURCE_LOWD:%.cpp=%.o)

HEADER_HIGHD = $(HEADER_GENERIC) popgen_highd.h
SOURCE_HIGHD = hypercube_function.cpp haploid_clone.cpp
OBJECT_HIGHD = $(SOURCE_HIGHD:%.cpp=%.o)

HEADER_HIV = hivpopulation.h
SOURCE_HIV = $(HEADER_HIV:%.h=%.cpp)
OBJECT_HIV = $(SOURCE_HIV:%.cpp=%.o)

OBJECTS = $(OBJECT_GENERIC) $(OBJECT_LOWD) $(OBJECT_HIGHD) $(OBJECT_HIV)

# Recipes
src: $(SRCDIR)/$(LIBRARY)

$(SRCDIR)/$(LIBRARY): $(OBJECTS:%=$(SRCDIR)/%)
	ar rcs $@ $^

$(OBJECT_GENERIC:%=$(SRCDIR)/%): $(SOURCE_GENERIC:%=$(SRCDIR)/%)
	$(CXX) $(CXXFLAGS) -c -o $@ $(@:.o=.cpp)

$(OBJECT_LOWD:%=$(SRCDIR)/%): $(SOURCE_LOWD:%=$(SRCDIR)/%) $(HEADER_LOWD:%=$(SRCDIR)/%)
	$(CXX) $(CXXFLAGS) -c -o $@ $(@:.o=.cpp)

$(OBJECT_HIGHD:%=$(SRCDIR)/%): $(SOURCE_HIGHD:%=$(SRCDIR)/%) $(HEADER_HIGHD:%=$(SRCDIR)/%)
	$(CXX) $(CXXFLAGS) -c -o $@ $(@:.o=.cpp)

$(OBJECT_HIV:%=$(SRCDIR)/%): $(SOURCE_HIV:%=$(SRCDIR)/%) $(HEADER_HIV:%=$(SRCDIR)/%)
	$(CXX) $(CXXFLAGS) -c -o $@ $(@:.o=.cpp)

clean-src:
	cd $(SRCDIR); rm -rf $(LIBRARY) *.o *.h.gch

##==========================================================================
# DOCUMENTATION
##==========================================================================
DOXYFILE   = $(DOCDIR)/Doxyfile
DOXY       = doxygen

# Recipes
doc:
	$(DOXY) $(DOXYFILE)

clean-doc:
	rm -rf $(DOCDIR)/latex $(DOCDIR)/html

##==========================================================================
# TESTS
##==========================================================================
TESTS_LIBRARIES = -lPopGenLib -lgsl -lgslcblas

TESTS_LDFLAGS = -L$(SRCDIR) -L$(LINKER_FOLDER) -O2
TESTS_CXXFLAGS = -I$(SRCDIR) -Wall -O2 -c -fPIC -g3

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
	$(CXX) $(TESTS_LDFLAGS) -o $@ $^ $(TESTS_LIBRARIES)

$(TESTS_OBJECT_LOWD:%=$(TESTSDIR)/%): $(TESTS_SOURCE_LOWD:%=$(TESTSDIR)/%) $(TESTS_HEADER_LOWD:%=$(TESTSDIR)/%)
	$(CXX) $(TESTS_CXXFLAGS) -c -o $@ $(@:.o=.cpp)

$(TESTS_HIGHD:%=$(TESTSDIR)/%): $(TESTS_OBJECT_HIGHD:%=$(TESTSDIR)/%) $(OBJECT_HIV:%=$(SRCDIR)/%) $(SRCDIR)/$(LIBRARY)
	$(CXX) $(TESTS_LDFLAGS) -o $@ $^ $(TESTS_LIBRARIES)

$(TESTS_OBJECT_HIGHD:%=$(TESTSDIR)/%): $(TESTS_SOURCE_HIGHD:%=$(TESTSDIR)/%) $(TESTS_HEADER_HIGHD:%=$(TESTSDIR)/%)
	$(CXX) $(TESTS_CXXFLAGS) -c -o $@ $(@:.o=.cpp)

clean-tests:
	cd $(TESTSDIR); rm -rf *.o $(TESTS_LOWD) $(TESTS_HIGHD)

##==========================================================================
# PYTHON BINDINGS
##==========================================================================
SWIG = swig
SWIGFLAGS = -c++ -python -O -castmode -keyword

SWIG_INTERFACE = PopGenLib.i
SWIG_WRAP = $(SWIG_INTERFACE:%.i=%_wrap.cpp)
SWIG_OBJECT = $(SWIG_INTERFACE:%.i=_%.so)
SWIG_PYMODULE = $(SWIG_INTERFACE:%.i=%.py)
SWIG_PYCMODULE = $(SWIG_INTERFACE:%.i=%.pyc)
SWIG_SUPPORT_1 = popgen_highd.i
SWIG_SUPPORT_2 = hivpopulation.i
SWIG_SUPPORT_3 = popgen_lowd.i
SWIG_SUPPORT_4 = popgen.i

## The following syntax is really black magic... look into the distutils source code and cry
PYTHON_SETUP = setup.py
PYTHON_LIBRARIES = -l'PopGenLib gsl gslcblas'
PYTHON_INCLUDES = -I$(SRCDIR):$(NUMPY_INCLUDES)
PYTHON_LIBDIRS = -L$(CURDIR)/$(SRCDIR)
PYTHON_FLAGS = --inplace $(PYTHON_INCLUDES) $(PYTHON_LIBRARIES) $(PYTHON_LIBDIRS)

# Recipes
python: $(SWIG_OBJECT:%=$(PYBDIR)/%)

$(SWIG_PYMODULE:%=$(PYBDIR)/%) $(SWIG_PYMCODULE:%=$(PYBDIR)/%) $(SWIG_OBJECT:%=$(PYBDIR)/%): $(SWIG_WRAP:%=$(PYBDIR)/%) $(SWIG_SOURCE:%=$(PYBDIR)/%) $(DISTUTILS_SETUP:%=$(PYBDIR)/%) $(SRCDIR)/$(LIBRARY)
	cd $(PYBDIR); python2 $(PYTHON_SETUP) build_ext $(PYTHON_FLAGS)

$(SWIG_WRAP:%=$(PYBDIR)/%): $(SWIG_HEADER_HIV:%=$(PYBDIR)/%) $(SWIG_INTERFACE:%=$(PYBDIR)/%) $(SRCDIR)/$(LIBRARY) $(SWIG_SUPPORT_1:%=$(PYBDIR)/%) $(SWIG_SUPPORT_2:%=$(PYBDIR)/%) $(SWIG_SUPPORT_3:%=$(PYBDIR)/%) $(SWIG_SUPPORT_4:%=$(PYBDIR)/%)
	$(SWIG) $(SWIGFLAGS) -o $@ $(SWIG_INTERFACE:%=$(PYBDIR)/%)

clean-python:
	cd $(PYBDIR); rm -rf $(SWIG_WRAP) $(SWIG_OBJECT) $(SWIG_PYMODULE) $(SWIG_PYCMODULE)

#############################################################################
