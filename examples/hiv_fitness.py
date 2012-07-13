# vim: fdm=indent
'''
author:     Fabio Zanini
date:       25/04/12
content:    Test script for the python bindings
'''

# Import module
import sys
sys.path.append('../pkg/python')

import numpy as np
import matplotlib.pyplot as plt
import FFPopSim as h



# Globals
N = 10000                   # Population size
adaptive_fraction = 0.01    # Fraction of beneficial/adaptive sites
effect_size_adaptive = 0.03 # Mean selection coefficient of beneficial alleles



# Script
if __name__ == '__main__':

    # Create population
    pop = h.hivpopulation(N)
    
    # Set replication/fitness landscape
    pop.set_replication_landscape(adaptive_fraction=adaptive_fraction,
                                  effect_size_adaptive=effect_size_adaptive)
    
    # Prepare figure
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    colors = ['b', 'r', 'g', 'cyan']
    
    # Evolve and plot histograms
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

    # Show the plot
    plt.ion()
    fig.show()
