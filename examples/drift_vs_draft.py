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

m=10
selection_coefficients = 0.0*np.ones(pop.L)     #most loci are neutral
selection_coefficients[::m] = -0.1              #every m-th locus is strongly deleterious

initial_allele_frequencies = 0.5*np.ones(pop.L) #define some initial allele frequencies
initial_allele_frequencies[::m] = 0

#initialize the population in LD with the specified allele frequencies
pop.set_allele_frequencies(initial_allele_frequencies, pop.carrying_capacity)
pop.set_trait_additive(selection_coefficients)


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

plt.figure()
for locus in xrange(0,pop.L,m):         #plot the allele frequency trajectories of the selected mutations
    plt.plot(tp, allele_frequencies[:,locus], c=cm.cool(locus),lw=2)

for locus in xrange(5,pop.L,50):        #plot a few neutral trajectories
    plt.plot(tp, allele_frequencies[:,locus], c=cm.cool(locus), lw=2)

plt.title('Drift and Draft')
plt.xlabel('Time [generations]')
plt.ylabel('Allele frequencies')
