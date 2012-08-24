.. _Examples lowd:

Examples for the low-dimensional package
========================================

Decay of linkage
^^^^^^^^^^^^^^^^
In recombining populations, genetic linkage decays with time and distance. This script shows the decay curves of a simple population. 

First, we load the FFPopSim module and the number crunchung and plotting facilities::

   import numpy as np
   import matplotlib.pyplot as plt
   import FFPopSim as h

.. note:: you might need to add the directory of the ``FFPopSim.py`` file to your Python ``PATH`` for this to work. If you decide to copy files around instead, please make sure that ``FFPopSim.py`` and ``_FFPopSim.so`` always be in the same folder.

Second, we set up the population::

   # specify parameters
   N = 500000                          # Population size
   L = 4                               # number of loci
   mu = 0.0                            # no new mutations
   r = 0.01                            # recombination rate
   
   ### set up
   pop = h.haploid_lowd(L)             # produce an instance of haploid_lowd with L loci
   pop.carrying_capacity = N           # set the steady-state population size
   pop.set_recombination_rates(r)      # assign the recombination rates
   pop.set_mutation_rates(mu)          # assign the mutation rate
   
   # initialize the population with N/2 individuals with genotypes 0, that is ----
   # and N/2 with the opposite genotype, that is ++++
   pop.set_genotypes([0, 2**L-1],[N/2, N/2])

Third, we let the population evolve in little steps and we track the linkage disequilibrium via the ``get_LD`` function::

   max_gen = 50
   LD_trajectories = [[pop.generation,pop.get_LD(0,1), pop.get_LD(0,2), pop.get_LD(0,3)]]
   for ii in range(max_gen):
       pop.evolve(5)               #N/10 generations between successive samples
       LD_trajectories.append([pop.generation, pop.get_LD(0,1), pop.get_LD(0,2), pop.get_LD(0,3)])
   LD_trajectories=np.array(LD_trajectories)

Fourth, we plot the resulting linkage disequilibrium curves::

   cols = ['r', 'b', 'g', 'm', 'c']
   for ii in range(LD_trajectories.shape[1]-1):
       plt.plot(LD_trajectories[:,0], LD_trajectories[:,ii+1], color=cols[ii], label=r'$D_{1'+str(ii+1)+'}$')
       plt.plot(LD_trajectories[:,0], np.exp(-LD_trajectories[:,0]* r * ii), ls='--', color=cols[ii])
   
   plt.legend()
   plt.xlabel('Time [generations]')
   plt.ylabel('LD $D_{ij}$')
   
   plt.ion()
   plt.show()

The typical plot we obtain is the following:

.. image:: figures/examples/decay_of_LD.png







