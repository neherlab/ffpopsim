.. _Examples highd:

Examples for the high-dimensional package
=========================================

Condensation of genotypes driven by epistasis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Epistatis and rugged landscapes favour more clonal populations compared to smooth and monotonic landscapes. This process is observed in this example.

.. warning:: the careful reader will notice that the argument parsing part of the example is skipped in this guide.

First, we load the FFPopSim module and the number crunchung and plotting facilities::

   import numpy as np
   import matplotlib.pyplot as plt
   import FFPopSim as ffpop

Second, we set up the population::

   L=64
   pop=ffpop.haploid_highd(L)
   pop.outcrossing_rate=0
   pop.set_random_epistasis(0.05)
   pop.recombination_model = ffpop.FREE_RECOMBINATION
   pop.set_allele_frequencies(np.ones(L)*0.5, 10000)

Third, we let the population evolve and collect statistics on fitness, clone size, and participation ratio along the way::

   pfit = pop.get_fitness_statistics()
   popstat = []
   for gen in xrange(1, 200):
       #append current statistics to the list
       pfit = pop.get_fitness_statistics()
       popstat.append([gen,pfit.mean, pfit.variance, pop.participation_ratio, pop.number_of_clones])
       
       #evolve for dt generations and clean up
       pop.evolve()
       pop.unique_clones()
       pop.calc_stat()
   
   popstat=np.array(popstat)

Fourth, we plot some interesting observables:

.. image:: figures/examples/condensation_participation_ratio.png
.. image:: figures/examples/condensation_number_of_clones.png
