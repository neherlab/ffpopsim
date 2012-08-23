.. _Requirements:

Requirements
============
FFpopSim is developed to work on 32 or 64 bit machines running Linux or Mac OSX. The basic library is written in C++ and can be used and extended independently of the Python bindings. This page will focus on the Python interface.

Runtime Requirements
--------------------
- Python_ 2.7+ (no Python 3 support yet) are the only supported Python versions. Older Python versions will *never* be supported, but Python 3 might become so in the future.
- NumPy_ is used extensively in this library. We suggest to import numpy explicitely before using the library (but it will work in any case).
- matplotlib_ is used in the plot functions. As long as you do not call those, you can live without it. However, we suggest to import it explicitely before using the library.

.. _Python: http://www.python.org/
.. _NumPy: http://numpy.scipy.org/
.. _matplotlib: http://matplotlib.sourceforge.net/

Building Requirements
---------------------
In order to build the Python bindings to FFPopSim, you need the following programs:

- Python_ 2.7+ (no Python 3 support yet), including header files
- NumPy_, and matplotlib_, including header files and shared libraries
- SWIG_
- a C++ compiler

The building process has been tested on SWIG 2.0, Numpy 1.6, matplotlib 1.1, and g++ 4.7.

We are trying to reduce the requirements, lifting e.g. the need for SWIG, and making the building less manual. This work is in progress.

.. _SWIG: http://www.swig.org/
.. _Python: http://www.python.org/
.. _NumPy: http://numpy.scipy.org/
.. _matplotlib: http://matplotlib.sourceforge.net/

