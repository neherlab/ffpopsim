.. _First steps with FFPopSim:

First steps with FFPopSim
=========================
FFPopSim is supposed to be easy to use. This page is meant to help new users to familiarize themselves with the library by means of examples. For a complete reference of classes and functions, please see the :ref:`Contents` page.

For the impatient ones...
-------------------------
An effective way to discover all available methods is to import FFPopSim from
the iPython_ interactive shell, create a population, and use TAB autocompletion:

.. sourcecode:: ipython

    In [1]: import FFPopSim as h
    In [2]: pop = h.haploid_lowd(5)   #create a population with 5 loci
    In [3]: pop.      <--- TAB

Importing FFPopSim
------------------
FFPopSim is a single Python module. As such, you can import it with the python ``import`` statement,
provided the files ``FFPopSim.py`` and ``_FFPopSim.so`` are in a folder in your ``PYTHONPATH``.
If you wish to perform a system-wide installation of FFPopSim, call the make recipe ``python-install``
*as a superuser*:

.. code-block:: bash

   $ sudo make python-install

.. note:: if this sounds new to you, just put those two files into your current directory, from which
          you plan to call the Python interpreter. ``import`` statements first look in the current
          folder for modules.

We recommend to import Numpy_ and matplotlib_ together with FFPopSim. In short, all your scripts should
begin with the following piece of code::

   import numpy as np
   import matplotlib.pyplot as plt
   import FFPopSim

The selective import of parts of FFPopSim using ``from FFPopSim import <xxx>`` is discouraged and its results
are untested.

Examples
--------
See the main page for :ref:`examples <examples_main>`.


.. _iPython: http://ipython.org/
.. _NumPy: http://numpy.scipy.org/
.. _matplotlib: http://matplotlib.sourceforge.net/
