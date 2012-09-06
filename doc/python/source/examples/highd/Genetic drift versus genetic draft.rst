Genetic drift versus genetic draft
==================================
The next examples explores the interplay between genetic drift and
genetic draft, i.e., the effect of linked selection on the
trajectories of neutral alleles. The basic script is the same as
above, only that we now set a fitness landscape and change the effects
of some mutations from being deleterious to beneficial during the
simulation. This generates selective sweeps.

After importing the relevant modules, we build the population::

   L=256                                           # simulate 256 loci
   
   pop = h.haploid_highd(L)                        # produce an instance of haploid_lowd with L loci
   pop.carrying_capacity = 50000                   # set the average population size to 50000
   pop.outcrossing_rate = 1                        # make the species obligate outcrossing
   pop.crossover_rate = 0.02/pop.L                 # set the crossover rate of the segment to 2 centimorgans
   pop.mutation_rate = 0.1/pop.carrying_capacity   # per locus mutation rate equal to 0.1/N

    
In addition, we set the selection coefficients to 0 for most loci, but
make every 10th locus strongly deleterious::

    m=10
    selection_coefficients = 0.0*np.ones(pop.L)                 # most loci are neutral
    selection_coefficients[::m] = -0.1                          # every m-th locus is strongly deleterious
    pop.set_trait_additive(selection_coefficients,0)	        # trait 0 is by default fitness
    
Neutral loci are set the frequency 1/2, while the deleterious ones to
frequency 0. We initialize the population with those allele frequencies,
in linkage equilibrium::

    initial_allele_frequencies = 0.5*np.ones(pop.L)
    initial_allele_frequencies[::m] = 0
    pop.set_allele_frequencies(initial_allele_frequencies, pop.carrying_capacity)

Next, we start evolving and track the allele frequencies as we go
along. Every 200 generations, we pick a random locus from the
deleterious ones and make it beneficial.

::

    #evolve for 2000 generations and track the allele frequencies
    maxgen = 2000
    allele_frequencies = [pop.get_allele_frequencies()]
    tp = [pop.generation]
    while pop.generation<maxgen:
        pop.evolve(10)                                                  # procede 10 generations
    	if (pop.generation%200==0):                                     # every 200 generations, make one of the deleterious mutations beneficial
            print "generation:", pop.generation, 'out of', maxgen
            selection_coefficients[m*np.random.randint(0,25)] = 0.01
            pop.set_trait_additive(selection_coefficients)              # update fitness function

    allele_frequencies.append(pop.get_allele_frequencies())             # save the allele frequencies
    tp.append(pop.generation)                                           # and the associated generation

We now plot the frequency trajectories of all selected loci. Those
that become beneficial in the process have risen quickly to high
frequencies. When they sweep, they influence the trajectories of
linked neutral loci, of which also a few trajectories are shown.

::

    for locus in xrange(0,pop.L,m):         #plot the allele frequency trajectories of the selected mutations
        plt.plot(tp, allele_frequencies[:,locus], c=cm.cool(locus),lw=2)
    
    for locus in xrange(5,pop.L,50):        #plot a few neutral trajectories
        plt.plot(tp, allele_frequencies[:,locus], c=cm.cool(locus), lw=2)
    
.. image:: ../../figures/examples/drift_and_draft.png
