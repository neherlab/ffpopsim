/*
 * multi_locus_genealogy.cpp
 *
 *  Created on: Oct 14, 2012
 *      Author: richard
 */

#include "ffpopsim_highd.h"

multi_locus_genealogy::multi_locus_genealogy() {
	// TODO Auto-generated constructor stub

}

multi_locus_genealogy::~multi_locus_genealogy() {
	// TODO Auto-generated destructor stub
}

void multi_locus_genealogy::add_generation(double baseline){
	for (unsigned int locusIndex = 0; locusIndex<trees.size(); locusIndex++){
		trees[locusIndex].add_generation(newGenerations[locusIndex], baseline);
	}
}

void multi_locus_genealogy::track_locus(int new_locus){
	loci.push_back(new_locus);
	rooted_tree temp_tree;
	trees.push_back(temp_tree);
	vector <node_t> temp_generation;
	newGenerations.push_back(temp_generation);
}

/**
 * @brief: allocates memory for a sufficient number of slots in the genealogy
 *
 * @params: the total number of slots required (should be population.size() in most cases)
 */
int multi_locus_genealogy::extend_storage(int n) {
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
