.. _Contents:

.. currentmodule:: FFPopSim

Contents
========
The functionality of :mod:`FFPopSim` is described in detail in the following sections.

Population Class Overviews
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Simulations in :mod:`FFPopSim` are centered around population classes. These methods of these classes are divided in categories and described at the following pages:

.. toctree::
   :maxdepth: 1

   haploid_lowd_categories
   haploid_highd_categories
   hivpopulation_categories

The genealogies at individual loci can be tracked within :class:`haploid_highd`, using :class:`multi_locus_genealogy`. Internally, this class itself contains instances of :class:`rooted_tree`. See :class:`haploid_highd`, the examples, or the reference pages below for more information on genealogy tracking.

Reference Pages
^^^^^^^^^^^^^^^
The classes used in :mod:`FFPopSim` are listed in the following pages for reference purposes:

.. toctree::
   :maxdepth: 1

   clone
   genotype_value_pair
   haploid_lowd
   haploid_highd
   hivgene
   hivpopulation
   index_value_pair
   multi_locus_genealogy
   polymorphism
   rooted_tree
   stat
   tree_key
   tree_step
   tree_node
   tree_edge

In addition, the following helper functions are defined:

.. toctree::
   :maxdepth: 1

   binarify
   integerify
