# vim: fdm=indent
'''
author:     Fabio Zanini
date:       28/05/12
content:    Immune escape example for FFPopGen
'''
# Import module
import sys
sys.path.append('../pkg/python')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import FFPopSim as h



# Globals
L = 4                           # Number of loci
N = 1e10                        # Population size
mu = 1e-5                       # Mutation rate
r = [1e-4] * (L-1)              # Recombination rates
f = [0.3, 0.2, 0.1, 0.05]       # Fitness (main/additive effects)




# Script
if __name__ == '__main__':

    # Create population
    pop = h.haploid_lowd(L)

    # Set genotypes, recombination and mutation rates, fitness landscape
    pop.set_genotypes([0],[N])
    pop.set_recombination_rates(r)
    pop.set_mutation_rate(mu)
    pop.set_fitness_additive(f)

    # Evolve the population
    times = []
    genotype_frequencies = []
    while pop.get_genotype_frequency(0b1111) < 0.99 and pop.generation<1e7:
        pop.evolve()
        times.append(pop.generation)
        genotype_frequencies.append(pop.get_genotype_frequencies())
    genotype_frequencies = np.asarray(genotype_frequencies)

    # Check what genotypes reach high frequencies
    ind = []
    for i in xrange(1<<L):
        if (genotype_frequencies[:,i] > 0.01).any():
            ind.append(i)
    n_plots = len(ind)

    # Plot
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    colors = [cm.jet(int(255.0 * i / (n_plots))) for i in xrange(n_plots)]
    lstyles = [[1,0],
               [2, 4, 7, 4],
               [7, 4],
               [2, 2],
               [5, 2],
               [2, 4, 14, 4],
               [6, 6]]
    for i, ii in enumerate(ind):
        line = ax.plot(times, genotype_frequencies[:,ii], c=colors[i],
                       lw=2.5,
                       label='{:04d}'.format(int(bin(ii)[2:])))
        if i < len(lstyles):
            line[0].set_dashes(lstyles[i])
    ax.set_xlabel('Generations')
    ax.set_ylabel('Genotype frequency')
    ax.set_title('Genotype frequencies in HIV immune escape', fontsize=14)
    ax.legend(loc=1, bbox_to_anchor = (1.0, 0.8))

    # Show the plot
    plt.ion()
    fig.show()
