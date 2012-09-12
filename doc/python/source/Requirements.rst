.. _Requirements:

Requirements
============

FFpopSim is developed to work on 32 or 64 bit machines running Linux or Mac OSX.
The basic library is written in C++ and can be used and extended independently
of the Python bindings.

.. note:: The Enthought Python Distribution (EPD_) is a widely used and
          well-maintained Python environment that provides all necessary
          Python packages for running FFPopSim (other than GSL_ and BOOST_).
          A basic EPD version is available for free at the following website:
          http://www.enthought.com/products/epd_free.php 

Runtime Requirements
--------------------

- Python_ 2.7+ (no Python 3 support yet) are the only supported Python versions.
  Older Python versions will *never* be supported, but Python 3 might become so
  in the future.

- NumPy_ is used extensively in this library. We suggest to import numpy
  explicitely before using the library (but it will work in any case).

- matplotlib_ is used in the plot functions. As long as you do not call those,
  you can live without it. However, we suggest to import it explicitely before
  using the library.

- GSL_ shared libraries, from the GNU Scientific library

Building Requirements
---------------------
In order to build the Python bindings to FFPopSim, you need the following programs:

   - a C++ compiler, e.g. GCC_
   - Python_ 2.7+ (no Python 3 support yet), including header files
   - NumPy_, including header files and shared libraries
   - GSL_, the GNU Scientific library
   - BOOST_, the C++ extension library
   - an implementation of Make, e.g. `GNU Make`_
   - distutils_, a library for installing Python packages

In addition, if you modify the sources and want to regenerate the Python bindings, you
will need the following programs:

   - SWIG_, the Simplified Wrapper and Interface Generator

Finally, if you want to rebuild the documentation, you will need the following programs:

   - Sphinx_, the Python documentation generator, for Python 2.x

The building process has been tested on Python 2.7, Numpy 1.6, gcc 4.7, gsl 1.15, boost
1.50. The regeneration part has been tested on SWIG 2.0. The documentation has been
created with Sphinx 1.1.

.. _GCC: http://gcc.gnu.org/
.. _GSL: http://www.gnu.org/software/gsl/
.. _BOOST: http://www.boost.org/
.. _SWIG: http://www.swig.org/
.. _Python: http://www.python.org/
.. _NumPy: http://numpy.scipy.org/
.. _matplotlib: http://matplotlib.sourceforge.net/
.. _Sphinx: http://sphinx.pocoo.org/
.. _GNU Make: http://www.gnu.org/software/make/
.. _distutils: http://docs.python.org/library/distutils.html
.. _EPD: http://www.enthought.com/products/epd.php
