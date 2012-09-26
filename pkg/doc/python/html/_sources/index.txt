.. FFPopSim documentation master file, created by
   sphinx-quickstart2 on Wed Aug 22 09:37:23 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to FFPopSim
====================================
.. automodule:: FFPopSim

Documentation
-----------------
.. toctree::
   :maxdepth: 1

   Requirements
   Install
   First steps with FFPopSim
   Contents

In addition, the underlying C++ library is documented `here
<../../cpp/html/index.html>`_. Note that some objects and functions have
slightly different names in C++ and Python.

.. _examples_main:

Examples
---------------------
Usage examples of FFPopSim can be found at the following pages.

The descriptions focus on FFPopSim and tend to ignore aesthetic aspects of the
scripts such as figures, labels, *et similia*. This also means that glueing
together the various code chunks found at those pages will not produce exactly
the same figures; neither is this necessary at all, because the full source
code for all examples (and more) can be found in the `examples
<../../../../examples>`_ folder.

.. note:: examples are ordered by increasing complexity.

Low-dimensional examples
^^^^^^^^^^^^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 1

   examples/lowd/Decay_of_linkage_disequilibrium
   examples/lowd/Mutation-selection
   examples/lowd/Time_complexity_and_scaling
   examples/lowd/Valley_crossing
   examples/lowd/Fitness_wave


High-dimensional examples
^^^^^^^^^^^^^^^^^^^^^^^^^
.. toctree::
   :maxdepth: 1

   examples/highd/Genetic_drift
   examples/highd/Genetic_drift_versus_genetic_draft
   examples/highd/Mutation-selection
   examples/highd/Condensation_of_genotypes_driven_by_epistasis


Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

