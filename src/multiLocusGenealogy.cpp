/*
 * multiLocusGenealogy.cpp
 *
 *  Created on: Oct 14, 2012
 *      Author: richard
 */

#include "ffpopsim_highd.h"

multiLocusGenealogy::multiLocusGenealogy() {
	// TODO Auto-generated constructor stub

}

multiLocusGenealogy::~multiLocusGenealogy() {
	// TODO Auto-generated destructor stub
}

void multiLocusGenealogy::add_generation(double baseline){
	for (unsigned int locusIndex = 0; locusIndex<trees.size(); locusIndex++){
		trees[locusIndex].add_generation(newGenerations[locusIndex], baseline);
	}
}

void multiLocusGenealogy::track_locus(int newLocus){
	loci.push_back(newLocus);
	rootedTree tempTree;
	trees.push_back(tempTree);
	vector <node_t> tempGeneration;
	newGenerations.push_back(tempGeneration);
}

/**
 * @brief: allocates memory for a sufficient number of slots in the genealogy
 *
 * @params: the total number of slots required (should be population.size() in most cases)
 */
int multiLocusGenealogy::extend_storage(int n) {
	node_t temp_ancestor;
	temp_ancestor.own_key.age=0; temp_ancestor.own_key.index=-1;
	temp_ancestor.fitness = 0;
	temp_ancestor.number_of_offspring = 0;
	temp_ancestor.clone_size = 0;
	temp_ancestor.crossover[0] = 0;
	temp_ancestor.crossover[1] = RT_VERYLARGE;
	for (unsigned int locus=0; locus<loci.size(); locus++){
		temp_ancestor.parent_node = trees[locus].get_MRCA();
		newGenerations[locus].resize(n, temp_ancestor);
	}
	return 0;
}
