# vim: fdm=indent
'''
author:     Fabio Zanini
date:       07/09/12
content:    Example of tracking the fitness distribution wave along the simulation.
            This example also shows subclussing of FFPopSim classes.
'''
# Import module
import sys
sys.path.append('../pkg/python')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import FFPopSim as h



# Subclass haploid_lowd
class haploid_lowd_track(h.haploid_lowd):
    '''This class tracks the fitness distribution'''

    def __init__(self, *args, **kwargs):
        super(haploid_lowd_track, self).__init__(*args, **kwargs)

        # Prepare the tracking list
        self.fitness_wave = []

    def evolve(self, *args, **kwargs):
        ret = super(haploid_lowd_track, self).evolve(*args, **kwargs)
        self.fitness_wave.append((self.generation,
                                  self.get_fitness_histogram(density=True,
                                                             n_sample=50000,
                                                             bins=30)))
        return ret

    def plot_tracks(self):
        '''Plot the tracked histograms'''
        l = len(self.fitness_wave)
        colors = cm.jet([int(255.0 * i / l) for i in xrange(l)])
        for i, (g, dist) in enumerate(self.fitness_wave):
            x = dist[1]
            x = 0.5 * (x[:-1] + x[1:])
            y = dist[0] + 1e-5
            plt.plot(x, y, c=colors[i], lw=2, label='gen '+str(g))
        plt.xlabel('Fitness')
        plt.legend(loc=1)
        plt.yscale('log')
        plt.ylim(1e-1)
        plt.title('Fitness distribution over time')
        plt.ion()
        plt.show()




# specify parameters
N = 500000                          # Population size
L = 12                               # number of loci
mu = 1e-5                           # mutation rate
r = 0.01                            # recombination rate


# script
if __name__ == '__main__':

    ### set up
    pop = haploid_lowd_track(L)                             # produce an instance of haploid_lowd with L loci
    pop.carrying_capacity = N                               # set the steady-state population size
    pop.set_recombination_rates(r, h.SINGLE_CROSSOVER)      # assign the recombination rates
    pop.set_mutation_rates(mu)                              # assign the mutation rate
    
    # initialize the population with N/2 individuals with genotypes 0, that is ----
    # and N/2 with the opposite genotype, that is ++++
    pop.set_wildtype(N)

    # set a fitness landscape
    pop.set_fitness_additive(1e-3 * np.linspace(1, 4, L))

    # Track 5 distributions, with 200 generations between each other
    for i in xrange(5):
        pop.evolve(200)

    # Plot the fitness wave
    pop.plot_tracks()

