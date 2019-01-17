/**
 * @file test_genealogies.cpp
 * @brief Tests for the genealogies in the high-dimensional simulation library.
 * @author Richard Neher, Boris Shraiman, Fabio Zanini
 * @version
 * @date 2012-04-20
 */

/* Include directives */
#include "ffpopsim_highd.h"
#define HIGHD_BADARG -1354341
#define NOTHING 1e-10

/* Be verbose? */
#define HIGHD_VERBOSE 1


/* Test evolution */
int genealogy() {
	int L = 2000;
	int N = 200;

	haploid_highd pop(L);

	pop.set_mutation_rate(0.0001);
	pop.outcrossing_rate = 1e-2;
	pop.crossover_rate = 0;
	pop.recombination_model = CHROMOSOMES;
	vector <int> gen_loci;
	gen_loci.push_back(0);
	gen_loci.push_back(L-1);
	pop.track_locus_genealogy(gen_loci);

	pop.set_wildtype(N/2);		// start with a population of the right size
	boost::dynamic_bitset<> nonWT;
	nonWT.resize(L);
	for (int i=0; i<L;i++){
		nonWT[i]=1;
	}
	// pop.add_genotype(nonWT, N/2);

	vector <int> loci;
	for(int i=0; i< L; i++) {
		loci.assign(1, i);
		pop.add_fitness_coefficient(1e-5, loci);
		loci.clear();
	}

	// pop.evolve(5);
	// for (int i=0; i<pop.get_number_of_clones(); i++){
	// 	cout <<"gt "<<i<<" "<<pop.get_genotype_string(i)<<"\n";
	// }

	stat_t fitstat;
	int err=0;
	gsl_histogram *SFS = gsl_histogram_alloc(20);
	gsl_histogram_set_ranges_uniform(SFS,0,1);
	for (int i=0; i<10000; i++) {
		pop.evolve();
		// pop.calc_stat();
		// fitstat = pop.get_fitness_statistics();
		// cerr  <<"generations: "<<i<<" af: "<<pop.get_allele_frequency(5)<<'\t'<<pop.get_allele_frequency(50)<<'\t'<<fitstat.mean<<'\t'<<fitstat.variance<<'\n';
	}
	// cerr <<"check trees: "<<pop.genealogy.trees.size()<<endl;
	// for (unsigned int genlocus=0; genlocus<gen_loci.size(); genlocus++){
	// 	err =pop.genealogy.trees[genlocus].check_tree_integrity();
	// }
	vector <int> clones;
	vector <tree_key_t> clone_keys;
	rooted_tree subtree;
	tree_key_t temp;
	int subtree_size = 10;
	for (unsigned int genlocus=0; genlocus<gen_loci.size(); genlocus++){
		cerr <<"PRINT ENTIRE TREE"<<endl;
		cerr <<pop.genealogy.trees[genlocus].print_newick()<<endl;
		clones.clear();
		clone_keys.clear();
		pop.random_clones(subtree_size,&clones);
		for (vector <int>::iterator clone=clones.begin(); clone!=clones.end(); clone++){
			temp.age=pop.get_generation()-1;
			temp.index = *clone;
			clone_keys.push_back(temp);
			cerr <<"gt "<<*clone<<" "<<pop.get_genotype_string(*clone)<<"\n";
		}
		pop.genealogy.trees[genlocus].construct_subtree(clone_keys, subtree);
		cerr <<"PRINT SUB TREE WITH "<<subtree_size<<" LEAFS"<<endl;
		cerr<<subtree.print_newick()<<endl;
		cerr <<'\n';
	}
	return err;
}

/* MAIN */
int main(int argc, char **argv){

	int status= 0;
	if (argc > 1) {
		cout<<"Usage: "<<argv[0]<<endl;
		status = 1;
	} else {
		status += genealogy();
	}
	cout<<"Number of errors: "<<status<<endl;
	return status;
}


