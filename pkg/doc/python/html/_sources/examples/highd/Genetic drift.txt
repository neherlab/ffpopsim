Genetic drift
=============
As a very basic example, lets simulate how allele frequencies change
due to genetic drift. If we want to track a large number of loci, we
can use haploid_hd. First, we import the necessary tools, including
FFPopSim::

    import numpy as np
    from matplotlib import pyplot as plt
    import FFPopSim as h

Next, we specify the parameters of the population and set up the
population.

::

   L=256                                           # simulate 256 loci
   
   pop = h.haploid_highd(L)                        # produce an instance of haploid_lowd with L loci
   pop.carrying_capacity = 5000                    # set the average population size to 50000
   pop.outcrossing_rate = 1                        # make the species obligate outcrossing
   pop.crossover_rate = 0.02/pop.L                 # set the crossover rate of the segment to 2 centimorgans
   pop.mutation_rate = 0.1/pop.carrying_capacity   # per locus mutation rate equal to 0.1/N
   initial_allele_frequencies = 0.5*np.ones(pop.L) # define some initial allele frequencies
   
   # initialize the population in LD with the specified allele frequencies
   pop.set_allele_frequencies(initial_allele_frequencies, pop.carrying_capacity)

We now have a population of size 5000 where each of the 256 loci is
present at frequency 1/2. Next, we want to evolve the population and
track the frequencies of the alleles.

::

    #evolve for 2000 generations and track the allele frequencies
    maxgen = 2000
    #make lists to store allele frequencies and time points
    allele_frequencies = [pop.get_allele_frequencies()] 
    tp = [pop.generation]
    
    while pop.generation<maxgen:
        pop.evolve(10)                                          #procede 10 generations
        allele_frequencies.append(pop.get_allele_frequencies()) #save the allele frequencies
        tp.append(pop.generation)                               #and the associated generation

The array *allele_frequencies* now contains the frequencies of 256
loci every 10 generations. We can plot a subset of these as follows::

    for locus in xrange(5,pop.L,50):        #plot a few neutral trajectories
        plt.plot(tp, allele_frequencies[:,locus], c=cm.cool(locus), lw=2)

.. image:: ../../figures/examples/genetic_drift.png
