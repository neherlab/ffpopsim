.. _Examples highd:

Examples for the high-dimensional package
=========================================
.. contents:: `Table of contents`
   :depth: 2

Genetic Drift
^^^^^^^^^^^^^
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

   #simulate 256 loci
   L=256
   
   ### set up
   pop = h.haploid_highd(L)                        #produce an instance of haploid_lowd with L loci
   pop.carrying_capacity = 5000                    #set the average population size to 50000
   pop.outcrossing_rate = 1                        #make the species obligate outcrossing
   pop.crossover_rate = 0.02/pop.L                 #set the crossover rate of the segment to 2 centimorgans
   pop.mutation_rate = 0.1/pop.carrying_capacity   #per locus mutation rate equal to 0.1/N
   initial_allele_frequencies = 0.5*np.ones(pop.L) #define some initial allele frequencies
   
   #initialize the population in LD with the specified allele frequencies
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

    #convert to an array to enable slicing
    allele_frequencies=np.array(allele_frequencies)

The array *allele_frequencies* now contains the frequencies of 256
loci every 10 generations. We can plot a subset of these as follows::

    plt.figure()
    for locus in xrange(5,pop.L,50):        #plot a few neutral trajectories
        plt.plot(tp, allele_frequencies[:,locus], c=cm.cool(locus), lw=2)
	
    plt.title('Genetic Drift')
    plt.xlabel('Time [generations]')
    plt.ylabel('Allele frequencies')

.. image:: figures/examples/genetic_drift.png


Genetic drift versus genetic draft
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The next examples explores the interplay between genetic drift and
genetic draft, i.e., the effect of linked selection on the
trajectories of neutral alleles. The basic script is the same as
above, only that we now set a fitness landscape and change the effects
of some mutations from being deleterious to beneficial during the
simulation. This generates selective sweeps. As before, we have::

    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import FFPopSim as h
    
    L=256
    
    ### set up
    pop = h.haploid_highd(L)                        #produce an instance of haploid_lowd with L loci
    pop.carrying_capacity = 50000                   #set the average population size to 50000
    pop.outcrossing_rate = 1                        #make the species obligate outcrossing
    pop.crossover_rate = 0.02/pop.L                 #set the crossover rate of the segment to 2 centimorgans
    pop.mutation_rate = 0.1/pop.carrying_capacity   #per locus mutation rate equal to 0.1/N

    
In addition, we set the selection coefficients to 0 for most loci, but
make every 10th locus strongly deleterious::

    m=10
    selection_coefficients = 0.0*np.ones(pop.L)     #most loci are neutral
    selection_coefficients[::m] = -0.1              #every m-th locus is strongly deleterious
    pop.set_trait_additive(selection_coefficients,0)	#trait 0 is by default fitness
    
Neutral loci are set the frequency 1/2, while the deleterious ones to
frequency 0::

    initial_allele_frequencies = 0.5*np.ones(pop.L) #define some initial allele frequencies
    initial_allele_frequencies[::m] = 0

    #initialize the population in linkage equilibrium with the specified allele frequencies
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
        pop.evolve(10)                              #procede 10 generations
    	if (pop.generation%200==0):                 #every 200 generations, make one of the deleterious mutations beneficial
            print "generation:", pop.generation, 'out of', maxgen
            selection_coefficients[m*np.random.randint(0,25)] = 0.01
            pop.set_trait_additive(selection_coefficients)      #update fitness function

    allele_frequencies.append(pop.get_allele_frequencies()) #save the allele frequencies
    tp.append(pop.generation)                               #and the associated generation

    #convert to an array to enable slicing
    allele_frequencies=np.array(allele_frequencies)

We now plot the frequency trajectories of all selected loci. Those
that become beneficial in the process have risen quickly to high
frequencies. When they sweep, they influence the trajectories of
linked neutral loci, of which also a few trajectories are shown.

::

    plt.figure()
    for locus in xrange(0,pop.L,m):         #plot the allele frequency trajectories of the selected mutations
        plt.plot(tp, allele_frequencies[:,locus], c=cm.cool(locus),lw=2)
    
    for locus in xrange(5,pop.L,50):        #plot a few neutral trajectories
        plt.plot(tp, allele_frequencies[:,locus], c=cm.cool(locus), lw=2)
    
    plt.title('Drift and Draft')
    plt.xlabel('Time [generations]')
    plt.ylabel('Allele frequencies')

.. image:: figures/examples/drift_and_draft.png


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
