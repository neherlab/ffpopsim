.. _Requirements:

Requirements
============

FFpopSim is developed to work on 32 or 64 bit machines running Linux or Mac OSX.
The basic library is written in C++ and can be used and extended independently
of the Python bindings.

The Python bindings are also distributed as ready-to-use binary files (see
:doc:`Install`); you can still build the library yourself though, if you prefer
to do so.

- On Linux, FFPopSim is expected to be compatible with all distributions, provided
  they are up to date (glibc 2.14 recommended).

- On Mac OSX, only Intel CPUs are expected to work, and only **Mac OSX 10.6 or
  later**. If you have an earlier Mac computer, you can still try to build the
  library yourself (see :doc:`Install`), but that might fail.


Runtime Requirements
--------------------

- Python_ 2.7 (no Python 3 support yet): older Python versions will *never* be
  supported, but Python 3 might become so in the future. If you have only Python
  2.6 or earlier, consider using the EPD_ Python distribution or updating your
  system.

- NumPy_ 1.6: if your Python distribution has only NumPy 1.5 or earlier,
  consider using the EPD_ Python distribution, building FFPopSim from source, or
  updating your system. It is recommended to import numpy explicitely before
  using the library, as shown in the examples.

- matplotlib_ is used in the plot functions. As long as you do not call those
  functions, you can live without it. However, it is recommended to import it
  explicitely before using the library, as shown in the examples.

.. note:: The Enthought Python Distribution (EPD_) is a widely used and
          well-maintained Python environment that provides all necessary
          Python packages for running FFPopSim, including a recent NumPy
          version (but no GSL_ and BOOST_). A basic EPD version is available
          for free at the following website:

          http://www.enthought.com/products/epd_free.php


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
