/**
 * @file mainpage.cpp
 * @brief Main page for the documentation.
 * @author Richard Neher, Boris Shraiman, Fabio Zanini
 * @version 
 * @date 2012-04-19
 *
 * @mainpage PopGenLib Documentation
 *
 * PopGenLib is a C++ library used to simulate evolution of genome populations.
 *
 * Currently, two kinds of populations are supported, depending on their genome size.
 * 1. For very short genomes, the full space of possible genomes is generated and simulations can be performed with any population size.
 * 2. For longer genomes, the dimension of the genotype space increases too rapidly, and only the observed genotypes are monitored.
 * If a new genotype is created either by mutation or by recombination, a larger amount of resources is required.
 * Population sizes up to \f$10^5\f$ or \f$10^6\f$ for a genome length of \f$10^4\f$ can be modeled this way.
 */

