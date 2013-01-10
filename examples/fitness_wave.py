# vim: fdm=indent
'''
author:     Fabio Zanini, Richard Neher
date:       07/09/12
content:    Example of tracking the fitness distribution wave along the
            simulation.
            
Note: This example also shows subclassing of FFPopSim classes. In this case, the
syntax construct if __name__ == '__main__' is actually useful to import the
subclass in other scripts.
'''
# Import modules (setting the path should not be necessary when the module is
# installed in the PYTHONPATH)
import sys
sys.path.insert(0, '../pkg/python')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import FFPopSim as h


# subclass haploid_lowd
class haploid_lowd_track(h.haploid_lowd):
    '''This class tracks the fitness distribution'''

    def __init__(self, *args, **kwargs):
        '''Initialize the class instance'''
        super(haploid_lowd_track, self).__init__(*args, **kwargs)

        # prepare the tracking list
        self.fitness_wave = []


    def evolve(self, *args, **kwargs):
        '''Evolve the population and keep track of the fitness histogram'''
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

        # plot each histogram
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
L = 12                              # number of loci
mu = 1e-5                           # mutation rate
r = 0.01                            # recombination rate


# script
if __name__ == '__main__':
    print "This script illustrates subclassing of FFPopSim."
    print "In addition to FFPopSim, this class tracks the fitness \ndistribution and allows plotting of its history.\n"


    # set up population
    pop = haploid_lowd_track(L)                             # produce an instance of haploid_lowd with L loci
    pop.carrying_capacity = N                               # set the steady-state population size
    pop.set_recombination_rates(r, h.SINGLE_CROSSOVER)      # assign the recombination rates
    pop.set_mutation_rates(mu)                              # assign the mutation rate
    
    pop.set_wildtype(N)                                     # initialize a wildtype population

    pop.set_fitness_additive(1e-3 * np.linspace(1, 4, L))   # set additive fitness landscape
    pop.status()
    
    # track fitness distribution withfive points and 200 generations between each other
    for i in xrange(5):
        pop.evolve(200)

    # plot the fitness wave
    pop.plot_tracks()
