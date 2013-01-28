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

// test neutral evolution
int kingman_coalescent() {
	int L = 1000;
	int N = 100;
	int nbins = 20;
	haploid_highd pop(L);

	pop.set_mutation_rate(1e-3);
	pop.outcrossing_rate = 1e-1;
	pop.crossover_rate = 1e-2;
	pop.recombination_model = CROSSOVERS;
	vector <int> gen_loci;
	gen_loci.push_back(100);
	gen_loci.push_back(500);
	gen_loci.push_back(900);
	pop.track_locus_genealogy(gen_loci);

	pop.set_wildtype(N);		// start with a population of the right size

	stat_t fitstat;
	gsl_histogram *SFS = gsl_histogram_alloc(nbins);
	gsl_histogram_set_ranges_uniform(SFS,0,1);
	int nrun=500,T[4];
	T[0]=0;T[1]=0; T[2]=0; T[3]=0;
	pop.evolve(4*N);
	for (int n=0; n<nrun; n++){
		cout <<n<<" out of "<<nrun<<endl;
		pop.evolve(N);

		vector <int> clones;
		vector <tree_key_t> clone_keys;
		rooted_tree subtree;
		tree_key_t temp;
		for (unsigned int genlocus=0; genlocus<gen_loci.size(); genlocus++){
			T[0]+=pop.get_generation() - pop.genealogy.trees[genlocus].get_MRCA().age-1;
			pop.genealogy.trees[genlocus].SFS(SFS);

			for (int cs = 2; cs<5; cs++){
				clones.clear();
				clone_keys.clear();
				pop.random_clones(cs,&clones);
				for (vector <int>::iterator clone=clones.begin(); clone!=clones.end(); clone++){
					temp.age=pop.get_generation()-1;
					temp.index = *clone;
					clone_keys.push_back(temp);
				}
				pop.genealogy.trees[genlocus].construct_subtree(clone_keys, subtree);
				T[cs-1]+=pop.get_generation() - subtree.get_MRCA().age-1;
			}
		}
	}

	double norm = 1.0*N*nrun*gen_loci.size();

	cout <<"SITE FREQUENCY SPECTRUM"<<endl;
	for (int i=0; i<nbins; i++){
		cout <<(i+0.5)/nbins<<'\t'<<gsl_histogram_get(SFS,i)/norm*nbins<<'\t'<<2*nbins/(i+0.5)<<endl;
	}

	cout <<"COALESCENCE TIMES"<<endl;
	cout <<T[0]/norm<<'\t'<<T[1]/norm<<'\t'<<T[2]/norm<<'\t'<<T[3]/norm<<endl;
	return 0;
}


// test neutral evolution
int large_populations() {
	int L = 1000;
	int N = 10000;
	int nbins = 20;
	haploid_highd pop(L,0,1,true);

	pop.outcrossing_rate = 1e-1;
	pop.crossover_rate = 1e-3;
	pop.recombination_model = CROSSOVERS;
	vector <int> gen_loci;
	gen_loci.push_back(100);
	gen_loci.push_back(500);
	gen_loci.push_back(900);
	pop.track_locus_genealogy(gen_loci);

	vector <int>loci;
	for(int i=0; i< L; i++) {
		loci.assign(1, i);
		pop.add_fitness_coefficient(1e-2, loci);
		loci.clear();
	}
	pop.set_wildtype(N);		// start with a population of the right size
	stat_t fitstat;
	gsl_histogram *SFS = gsl_histogram_alloc(nbins);
	double bins[nbins+1];
	bins[0]=0;bins[nbins]=1.0;
	for (int b=1; b<nbins; b++){bins[b] = 1.0/(N*exp(-2*b*log(N)/nbins)+1);}
	gsl_histogram_set_ranges(SFS,bins,nbins+1);
	int nrun=100,T[4];
	T[0]=0;T[1]=0; T[2]=0; T[3]=0;
	pop.evolve(500);
	for (int n=0; n<nrun; n++){
		cout <<n<<" out of "<<nrun<<endl;
		pop.evolve(10);

		vector <int> clones;
		vector <tree_key_t> clone_keys;
		rooted_tree subtree;
		tree_key_t temp;
		for (unsigned int genlocus=0; genlocus<gen_loci.size(); genlocus++){
			T[0]+=pop.get_generation() - pop.genealogy.trees[genlocus].get_MRCA().age-1;
			pop.genealogy.trees[genlocus].SFS(SFS);

			for (int cs = 2; cs<5; cs++){
				clones.clear();
				clone_keys.clear();
				pop.random_clones(cs,&clones);
				for (vector <int>::iterator clone=clones.begin(); clone!=clones.end(); clone++){
					temp.age=pop.get_generation()-1;
					temp.index = *clone;
					clone_keys.push_back(temp);
				}
				pop.genealogy.trees[genlocus].construct_subtree(clone_keys, subtree);
				T[cs-1]+=pop.get_generation() - subtree.get_MRCA().age-1;
			}
		}
	}

	double norm = 1.0*nrun*gen_loci.size();

	cout <<"SITE FREQUENCY SPECTRUM"<<endl;
	double lower,upper;
	for (int i=0; i<nbins; i++){
		gsl_histogram_get_range(SFS,i,&lower,&upper);
		cout <<0.5*(upper+lower)<<'\t'<<gsl_histogram_get(SFS,i)/norm/(upper-lower)<<'\t'<<2*nbins/(i+0.5)<<endl;
	}

	cout <<"COALESCENCE TIMES"<<endl;
	cout <<T[0]/norm<<'\t'<<T[1]/norm<<'\t'<<T[2]/norm<<'\t'<<T[3]/norm<<endl;
	return 0;
}


/* Test evolution */
int genealogy() {
	int L = 1000;
	int N = 50;

	haploid_highd pop(L);

	pop.set_mutation_rate(1e-3);
	pop.outcrossing_rate = 1e-2;
	pop.crossover_rate = 1e-2;
	pop.recombination_model = CROSSOVERS;
	vector <int> gen_loci;
	gen_loci.push_back(100);
	gen_loci.push_back(500);
	gen_loci.push_back(900);
	pop.track_locus_genealogy(gen_loci);

	pop.set_wildtype(N);		// start with a population of the right size
	vector <int> loci;
	for(int i=0; i< L; i++) {
		loci.assign(1, i);
		pop.add_fitness_coefficient(1e-2, loci);
		loci.clear();
	}


	stat_t fitstat;
	int err;
	gsl_histogram *SFS = gsl_histogram_alloc(20);
	gsl_histogram_set_ranges_uniform(SFS,0,1);
	for (int i=0; i< 500; i++) {
		pop.evolve();
		pop.calc_stat();
		fitstat = pop.get_fitness_statistics();
		cerr  <<"generations: "<<i<<" af: "<<pop.get_allele_frequency(5)<<'\t'<<pop.get_allele_frequency(50)<<'\t'<<fitstat.mean<<'\t'<<fitstat.variance<<'\n';
	}
	cerr <<"check trees: "<<pop.genealogy.trees.size()<<endl;
	for (unsigned int genlocus=0; genlocus<gen_loci.size(); genlocus++){
		err =pop.genealogy.trees[genlocus].check_tree_integrity();
	}
	vector <int> clones;
	vector <tree_key_t> clone_keys;
	rooted_tree subtree;
	tree_key_t temp;
	for (unsigned int genlocus=0; genlocus<gen_loci.size(); genlocus++){
		cerr <<"PRINT ENTIRE TREE"<<endl;
		cerr <<pop.genealogy.trees[genlocus].print_newick()<<endl;
		for (int cs = 2; cs<5; cs++){
			clones.clear();
			clone_keys.clear();
			pop.random_clones(cs,&clones);
			for (vector <int>::iterator clone=clones.begin(); clone!=clones.end(); clone++){
				temp.age=pop.get_generation()-1;
				temp.index = *clone;
				clone_keys.push_back(temp);
			}
			pop.genealogy.trees[genlocus].construct_subtree(clone_keys, subtree);
			cerr <<"PRINT SUB TREE WITH "<<cs<<" LEAFS"<<endl;
			cerr<<subtree.print_newick()<<endl;
		}
		cerr <<'\n';

	}
	return err;
}

/* Test evolution */
int genealogy_infinite_sites() {
	int L = 1000;
	int N = 50;
	int err=0;
	haploid_highd pop(L,0,1,true);

	pop.outcrossing_rate = 1e-2;
	pop.crossover_rate = 1e-2;
	pop.recombination_model = CROSSOVERS;
	vector <int> gen_loci;
	gen_loci.push_back(100);
	gen_loci.push_back(500);
	gen_loci.push_back(900);
	err = pop.track_locus_genealogy(gen_loci);
	if (err) {
		cerr <<"ERROR setting up genealogies"<<endl;
		return -1;
	}
	cerr <<"number of trees: "<<pop.genealogy.trees.size()<<'\t'<<gen_loci.size()<<endl;

	pop.set_wildtype(N);		// start with a population of the right size
	vector <int> loci;
	for(int i=0; i< L; i++) {
		loci.assign(1, i);
		pop.add_fitness_coefficient(1e-2, loci);
		loci.clear();
	}


	stat_t fitstat;
	gsl_histogram *SFS = gsl_histogram_alloc(20);
	gsl_histogram_set_ranges_uniform(SFS,0,1);
	for (int i=0; i< 200; i++) {
		pop.evolve();
		pop.calc_stat();
		fitstat = pop.get_fitness_statistics();
		cerr <<"generations: "<<i<<" af: "<<pop.get_allele_frequency(5)<<'\t'<<pop.get_allele_frequency(50)<<'\t'<<pop.get_derived_allele_frequency(5)<<'\t'<<pop.get_derived_allele_frequency(50)<<'\t'<<fitstat.mean<<'\t'<<fitstat.variance<<'\n';
	}
	cerr <<"check trees: "<<pop.genealogy.trees.size()<<endl;
	for (unsigned int genlocus=0; genlocus<gen_loci.size(); genlocus++){
		cerr <<"at locus "<<gen_loci[genlocus]<<endl;
		err =pop.genealogy.trees[genlocus].check_tree_integrity();
	}
	vector <int> clones;
	vector <tree_key_t> clone_keys;
	rooted_tree subtree;
	tree_key_t temp;
	for (unsigned int genlocus=0; genlocus<gen_loci.size(); genlocus++){
		cerr <<"PRINT ENTIRE TREE"<<endl;
		cerr <<pop.genealogy.trees[genlocus].print_newick()<<endl;
		for (int cs = 2; cs<5; cs++){
			clones.clear();
			clone_keys.clear();
			pop.random_clones(cs,&clones);
			for (vector <int>::iterator clone=clones.begin(); clone!=clones.end(); clone++){
				temp.age=pop.get_generation()-1;
				temp.index = *clone;
				clone_keys.push_back(temp);
			}
			pop.genealogy.trees[genlocus].construct_subtree(clone_keys, subtree);
			cerr <<"PRINT SUB TREE WITH "<<cs<<" LEAFS"<<endl;
			cerr<<subtree.print_newick()<<endl;
		}
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
		status += genealogy_infinite_sites();
		//status += kingman_coalescent();
		//status += large_populations();
	}
	cout<<"Number of errors: "<<status<<endl;
	return status;
}


