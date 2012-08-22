.. FFPopSim documentation master file, created by
   sphinx-quickstart2 on Wed Aug 22 09:37:23 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to FFPopSim's documentation!
====================================

.. toctree::
   :maxdepth: 1

   haploid_lowd

.. automodule:: FFPopSim

Requirements
------------
- NumPy_ is used extensively in this library, and is *strongly* recommended for all users. We suggest to import numpy explicitely before using the library (but it will work in any case).
- matplotlib_ is used in the plot functions. As long as you do not call those, you can live without it. However, we suggest to import it explicitely before using the library.

.. _NumPy: http://numpy.scipy.org/
.. _matplotlib: http://matplotlib.sourceforge.net/

**Note**: the Python interface does not offer the full functionality of the underlying library, nor as much
flexibility. If you need to perform a peculiar type of simulations that is not accessible by this Python
interface, consider subclassing in C++ (or Python, but it might be slow).


Building Requirements
---------------------
In order to build the Python bindings to FFPopSim, you need the following programs:

- Python_ 2.7+ (no Python 3 support yet), including header files
- NumPy_, and matplotlib_, including header files and shared libraries
- SWIG_
- a C++ compiler

The building process has been tested on SWIG 2.0, Numpy 1.6, matplotlib 1.1, and g++ 4.7.

.. _SWIG: http://www.swig.org/
.. _Python: http://www.python.org/
.. _NumPy: http://numpy.scipy.org/
.. _matplotlib: http://matplotlib.sourceforge.net/

For the impatient ones...
-------------------------
An effective way to discover all available methods is to import FFPopSim from
an interactive shell (e.g. iPython), create a population, and use TAB autocompletion:

.. sourcecode:: ipython

    In [1]: import FFPopSim as h
    In [2]: pop = h.haploid_lowd(5, 2000)
    In [3]: pop.      <--- TAB


Reference Manual
----------------
FFPopSim exposes most classes and functions to Python, but not all. Here is a list of the currently supported ones:

- :ref:`haploid_lowd <haploid_lowd>`;

.. - :ref:`haploid_highd <haploid_highd>`;
   - :ref:`hivpopulation <hivpopulation>`;


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

