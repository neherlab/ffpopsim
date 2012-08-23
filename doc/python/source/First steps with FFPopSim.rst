.. _First steps with FFPopSim:

First steps with FFPopSim
=========================
FFPopSim is supposed to be easy to use. This page is meant to help new users to familiarize themselves with the library by means of examples. For a complete reference of classes and functions, please see the :ref:`Contents` page.

For the impatient ones...
-------------------------
An effective way to discover all available methods is to import FFPopSim from
an interactive shell (e.g. iPython), create a population, and use TAB autocompletion:

.. sourcecode:: ipython

    In [1]: import FFPopSim as h
    In [2]: pop = h.haploid_lowd(5, 2000)
    In [3]: pop.      <--- TAB

Examples
--------
The source code for all examples can be found in the ``examples`` folder.

.. toctree::
   :maxdepth: 1

   Examples lowd
   Examples highd
