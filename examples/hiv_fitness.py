# vim: fdm=indent
'''
author:     Fabio Zanini, Richard Neher
date:       25/04/12
content:    Example of the hivpopulation subclass. This example is similar to
            fitness_wave.py, but exploits the specialized subclass hivpopulation
            for simulating HIV populations.
'''
# Import modules (setting the path should not be necessary when the module is
# installed in the PYTHONPATH)
import sys
sys.path.insert(0, '../pkg/python')

import numpy as np
import matplotlib.pyplot as plt
import FFPopSim as h


# specify parameters
N = 10000                       # population size
adaptive_fraction = 0.01        # fraction of beneficial/adaptive sites
effect_size_adap = 0.03         # mean selection coefficient of beneficial alleles


# script
if __name__ == '__main__':

    pop = h.hivpopulation(N)    # create population with default parameters
    
    # set random replication/fitness landscape
    pop.set_replication_landscape(adaptive_fraction=adaptive_fraction,
                                  effect_size_adaptive=effect_size_adap)
    
    # prepare figure
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    colors = ['b', 'r', 'g', 'cyan']
    
    # evolve and plot histograms
    x0 = pop.get_fitness_statistics().mean
    for i in xrange(4):
        pop.evolve(250)
        h = pop.get_fitness_histogram()
        x = h[1][:-1] - x0
        y = h[0]
        w = (x[1:] - x[:-1]).mean()
        ax.bar(x, y,width=w,
               color=colors[i%len(colors)], alpha=0.8,
               label=str(pop.generation))

    ax.set_xlabel('Fitness (relative to founder)')
    ax.set_title('Fitness distribution of the population')
    ax.legend(loc=2)

    plt.ion()
    fig.show()
