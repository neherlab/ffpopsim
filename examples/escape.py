# vim: fdm=indent
'''
author:     Fabio Zanini, Richard Neher
date:       28/05/12
content:    Example of immune escape in HIV using a few effective loci and
            haploid_lowd

Note: this example shows some alternative Python syntaxes from most other
example scripts, such as the following:

    - if __name__ == '__main__'
    - 1<<L
    - fig.add_subplot(1, 1, 1)
    - '{:04d}'.format(...)
    - line[0].set_dashes
    - bbox_to_anchor = (...)

They are mainly cosmetic changes for easy scripts such as this one, but might
become relevant for more complex programs. Please refer to the Python
documentation for more information on those topics.
'''
# Import modules (setting the path should not be necessary when the module is
# installed in the PYTHONPATH)
import sys
sys.path.insert(0, '../pkg/python')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import FFPopSim as h


# specify parameters
L = 4                           # number of loci
N = 1e10                        # population size
mu = 1e-5                       # mutation rate
r = 1e-4                        # recombination rate
f = [0.3, 0.2, 0.1, 0.05]       # fitness (additive effects)



# Script
if __name__ == '__main__':

    pop = h.haploid_lowd(L)     # produce an instance of haploid_lowd with L loci
    pop.set_genotypes([0],[N])  # start with a wildtype population
    
    pop.set_recombination_rates(r)  # set recombination rate
    pop.set_mutation_rates(mu)      # set mutation rate
    pop.set_fitness_additive(f)     # set fitness landscape

    # evolve while tracking genotype frequencies
    times = []
    genotype_frequencies = []
    while (pop.get_genotype_frequency(0b1111) < 0.99) and (pop.generation < 1e7):
        pop.evolve()

        # get genotype frequencies and time
        genotype_frequencies.append(pop.get_genotype_frequencies())
        times.append(pop.generation)

    genotype_frequencies = np.asarray(genotype_frequencies)

    # check what genotypes reach high frequencies
    ind = []
    for i in xrange(1<<L):
        if (genotype_frequencies[:,i] > 0.01).any():
            ind.append(i)

    # plot
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    n_plots = len(ind)
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

    plt.ion()
    fig.show()
