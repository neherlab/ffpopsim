#!/usr/bin/env python
"""
Simple proof of concept for FFPopSim library
Simulates evolution of a haploid population with a few loci
"""
import sys
sys.path.insert(0, 'pkg/python')

import FFPopSim
import numpy as np

print("=" * 60)
print("FFPopSim Proof of Concept")
print("=" * 60)

# Create a haploid population with 4 loci
L = 4  # number of loci
N = 1000  # population size

print(f"\nCreating population with {L} loci and {N} individuals...")
pop = FFPopSim.haploid_lowd(L)

# Set fitness landscape to favor genotypes with more 1s
print("\nSetting fitness landscape...")
pop.set_fitness_additive([0.05, 0.08, 0.05, 0.08])  # additive benefits for derived alleles

# Initialize with a mixed population
# Add wildtype (0000) and a few mutant genotypes
print("Initializing population with genetic diversity...")
genotypes = [0b0000, 0b0001, 0b0010, 0b1111]
counts = [400, 200, 200, 200]
pop.set_genotypes(genotypes, counts)

print(f"Initial population size: {pop.N}")
print(f"Number of loci: {pop.L}")

# Evolve for a few generations
generations = 100
print(f"\nEvolving for {generations} generations...")

for gen in range(generations):
    pop.evolve()
    if gen % 20 == 0:
        stats = pop.get_fitness_statistics()
        print(f"  Generation {gen:3d}: mean fitness = {stats.mean:.4f}, variance = {stats.variance:.4f}, pop size = {pop.N}")

# Final statistics
stats = pop.get_fitness_statistics()
print(f"\nFinal generation {generations}:")
print(f"  Mean fitness: {stats.mean:.4f}")
print(f"  Fitness variance: {stats.variance:.4f}")
print(f"  Population size: {pop.N}")

# Get allele frequencies
print("\nFinal allele frequencies:")
freqs = pop.get_allele_frequencies()
for i, freq in enumerate(freqs):
    print(f"  Locus {i}: {freq:.3f}")

print("\n" + "=" * 60)
print("Success! FFPopSim is working correctly.")
print("=" * 60)
