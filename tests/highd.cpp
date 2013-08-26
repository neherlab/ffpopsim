/**
 * @file highd.cpp
 * @brief Tests for the high-dimensional simulation library.
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

/* Test generic library access */
int library_access() {
	index_value_pair_t ivp;
	ivp.index = 4;
	ivp.val = 4.5;
	if(HIGHD_VERBOSE) {
		cerr<<"index_value_pair_t: index = "<<ivp.index;
		cerr<<", val = "<<ivp.val<<endl;
	}

	coeff_single_locus_t csl(ivp.val, ivp.index);
	if(HIGHD_VERBOSE) {
		cerr<<"coeff_single_locus_t: locus = "<<csl.locus;
		cerr<<", value = "<<csl.value<<endl;
	}
	return 0;
}


/* Test sample initialization */
int sample_initialize() {
	int N = 5;

	sample sam;
	sam.set_up(N);

	for(int i=0; i< N; i++)
		sam.values[i] = i + 3;
	if(HIGHD_VERBOSE)
		cerr<<"Mean of sample: "<<sam.calc_mean()<<endl;
	return 0;
}

/* Test hypercube initialization */
int hc_initialize() {
	hypercube_highd hc;
	int status = hc.set_up(4,3);
	if(HIGHD_VERBOSE){
		cerr<<"First constructor: Memory allocation = "<<status<<endl;
		cerr<<"First constructor: Dimension: "<<hc.get_dim()<<endl;
	}
	hypercube_highd hc2(4,3);
	if(HIGHD_VERBOSE){
		cerr<<"Second constructor: Dimension: "<<hc2.get_dim()<<endl;
	}
	return 0;
}

/* Test setting coefficients on hypercube */
int hc_setting() {
	int L = 1000;

	hypercube_highd hc(L,3);

	double values[] = {0.1, 0.043, 1.3, -4.3};
	int myints[] = {16,2,77,29};
	vector<int> loci;
	for(int i=0; i<4; i++) {
		loci.assign(1, myints[i]);
		hc.add_coefficient(values[i], loci);
		loci.clear();
	}	
	
	if(HIGHD_VERBOSE){
		cerr<<"Coefficient indices and values: ";
		for(int i=0; i<4; i++) {
			cerr<<"("<<myints[i]<<","<<hc.get_additive_coefficient(myints[i])<<") ";
		}
		cerr<<endl;
	}
	return 0;
}

/* Test population initialization */
int pop_initialize() {
	int L = 1000;

	haploid_highd pop(L, 3, 1);
	if(HIGHD_VERBOSE)
		cerr<<"L = "<<pop.get_number_of_loci()<<endl;	
	return 0;	
}

/* Test evolution */
int pop_evolve() {
	int L = 1000;
	int N = 100;

	haploid_highd pop(L);

	pop.set_mutation_rate(1e-4);
	pop.outcrossing_rate = 1;
	pop.crossover_rate = 1e-2;
	pop.recombination_model = CROSSOVERS;
	pop.set_wildtype(N);		// start with a population of the right size

	vector <int> loci;
	for(int i=0; i< L; i++) {
		loci.assign(1, i);
		pop.add_fitness_coefficient(0.01, loci);
		loci.clear();
	}
	stat_t fitstat;
	cerr <<"    "<<"af5\t"<<"af50\t"<<"Fitness mean\t"<<"Fitness var"<<endl;
	cerr <<"---------------------------------------------"<<endl;
	for (int i=0; i< 10; i++) {
		pop.evolve();
		pop.calc_stat();
		fitstat = pop.get_fitness_statistics();
		cerr <<i<<": "<<pop.get_allele_frequency(5)<<'\t'<<pop.get_allele_frequency(50)<<'\t'<<fitstat.mean<<'\t'<<fitstat.variance<<'\n';
	}
	pop.calc_stat();

	stat_t fitness = pop.get_fitness_statistics();

	if(HIGHD_VERBOSE) {
		double af = 0;
		double tmp;
		for(int i=0; i< L; i++) {
			tmp = pop.get_allele_frequency(i);
			af = MAX(af, tmp);
		}

		cerr<<"Generation: "<<pop.get_generation()<<endl;
		cerr<<"Max allele freq: "<<af<<endl;
		cerr<<"Fitness mean and variance: "<<fitness.mean<<", "<<fitness.variance<<endl;
	}

	return 0;	
}


/* Test evolution */


/* Test random sampling */
int pop_sampling() {
	int L = 100;
	int N = 100;

	haploid_highd pop(L);

	pop.set_mutation_rate(1e-3);
	pop.outcrossing_rate = 1e-2;
	pop.crossover_rate = 1e-2;
	pop.recombination_model = CROSSOVERS;
	pop.set_wildtype(N);		// start with a population of the right size

	pop.evolve(10);

	vector <int> sample;
	pop.random_clones(10, &sample);

	if(HIGHD_VERBOSE) {
		cout<<"Number of clones: "<<pop.get_number_of_clones()<<endl;
		cout<<"Random individuals:";
		for(unsigned int i=0; i < sample.size(); i++)
			cout<<" "<<sample[i];
		cout<<endl;
	}


	return 0;
}

///* Test Hamming distance functions */
//int pop_Hamming() {
//	int L = 100;
//	int N = 100;
//
//	haploid_highd pop(L);
//
//	pop.mutation_rate = 1e-3;
//	pop.outcrossing_rate = 1e-2;
//	pop.crossover_rate = 1e-2;
//	pop.recombination_model = CROSSOVERS;
//	pop.set_wildtype(N);		// start with a population of the right size
//
//	pop.evolve(10);
//
//	vector <int> sample;
//	pop.random_clones(10, &sample);
//
//	// create chunks
//	vector <unsigned int *> chunks;
//	unsigned int * tmp;
//	tmp = new unsigned int [2]; tmp[0] = 0; tmp[1] = 10; chunks.push_back(tmp);
//	tmp = new unsigned int [2]; tmp[0] = 20; tmp[1] = 30; chunks.push_back(tmp);
//	tmp = new unsigned int [2]; tmp[0] = 50; tmp[1] = 80; chunks.push_back(tmp);
//
//	if(HIGHD_VERBOSE) {
//		cerr<<"Chunks: "<<endl;
//		for(size_t i=0; i < chunks.size(); i++)
//			cerr<<chunks[i][0]<<" "<<chunks[i][1]<<endl;
//	}
//
//	// measure distances
//	vector <int> distance1;
//	vector <int> distance2;
//	vector <int> distance3;
//	for(size_t i1=0; i1<sample.size(); i1++)
//		for(size_t i2=0; i2<i1; i2++) {
//			distance1.push_back(pop.distance_Hamming(((*(pop.current_pop))[sample[i1]]).genotype, ((*(pop.current_pop))[sample[i2]]).genotype));
//			distance2.push_back(pop.distance_Hamming(sample[i1], sample[i2]));
//			distance3.push_back(pop.distance_Hamming(sample[i1], sample[i2], &chunks, 2));
//		}
//
//
//	if(HIGHD_VERBOSE) {
//		cerr<<"Hamming distances:"<<endl;
//		for(size_t i=0; i < distance1.size(); i++)	cerr<<" "<<distance1[i]; cerr<<endl;
//		for(size_t i=0; i < distance2.size(); i++)	cerr<<" "<<distance2[i]; cerr<<endl;
//		for(size_t i=0; i < distance3.size(); i++)	cerr<<" "<<distance3[i]; cerr<<endl;
//	}
//
//	// free memory for the chunks
//	size_t n_chunks = chunks.size();
//	for(size_t i=0; i < n_chunks; i++) {
//		delete chunks.back();
//		chunks.pop_back();
//	}
//
//	return 0;
//}
//
///* Test divergence and diversity statistics */
//int pop_divdiv() {
//	int L = 100;
//	int N = 100;
//
//	haploid_highd pop(L);
//
//	pop.mutation_rate = 1e-3;
//	pop.outcrossing_rate = 1e-2;
//	pop.crossover_rate = 1e-2;
//	pop.recombination_model = CROSSOVERS;
//	pop.set_wildtype(N);		// start with a population of the right size
//
//	pop.evolve(10);
//
//	stat_t fitness = pop.get_fitness_statistics();
//	stat_t divergence = pop.get_divergence_statistics();
//	stat_t diversity = pop.get_diversity_statistics();
//
//	if(HIGHD_VERBOSE) {
//		cerr<<"Fitness: mean = "<<fitness.mean<<", variance = "<<fitness.variance<<endl;
//		cerr<<"Divergence: mean = "<<divergence.mean<<", variance = "<<divergence.variance<<endl;
//		cerr<<"Diversity: mean = "<<diversity.mean<<", variance = "<<diversity.variance<<endl;
//	}
//
//	return 0;
//}
//
//
///* Test divergence and diversity statistics */
//int pop_histograms() {
//	int L = 100;
//	int N = 100;
//	int err = 0;
//
//	haploid_highd pop(L);
//
//	pop.mutation_rate = 1e-3;
//	pop.outcrossing_rate = 1e-2;
//	pop.crossover_rate = 1e-2;
//	pop.recombination_model = CROSSOVERS;
//	pop.set_wildtype(N);		// start with a population of the right size
//
//	vector <int> loci;
//	for(int i=0; i< L/2; i++) {
//		loci.push_back(i);
//		pop.add_fitness_coefficient(0.003, loci);
//		loci.pop_back();
//	}
//
////	pop.evolve(10);
//	pop.calc_stat();
//
//	if(HIGHD_VERBOSE) {
//		cerr<<"Number of clones: "<<pop.get_number_of_clones()<<endl;
//		for(size_t i = 0; i < (unsigned int)(pop.get_number_of_clones()); i++)
//			cerr<<pop.get_fitness(i)<<" ";
//		cerr<<endl;
//	}
//
//	unsigned int bins = 10;
//	double leftbin, rightbin;
//	gsl_histogram *fhist;
//	gsl_histogram *dvhist;
//	gsl_histogram *dshist;
//	err = pop.get_fitness_histogram(&fhist);
//	if(err == 0) {
//		bins = gsl_histogram_bins(fhist);
//		if(HIGHD_VERBOSE) {
//			cerr<<"Fitness histogram:"<<endl;
//			for(size_t i=0; i<bins; i++) {gsl_histogram_get_range(fhist,i,&leftbin,&rightbin); cerr<<leftbin<<" ";}	cerr<<endl;
//			for(size_t i=0; i<bins; i++) cerr<<gsl_histogram_get(fhist,i)<<" "; cerr<<endl;
//		}
//		gsl_histogram_free(fhist);
//	}
//
//	err = pop.get_divergence_histogram(&dvhist);
//	if(err == 0) {
//		bins = gsl_histogram_bins(dvhist);
//		if(HIGHD_VERBOSE) {
//			cerr<<"Divergence histogram:"<<endl;
//			for(size_t i=0; i<bins; i++) {gsl_histogram_get_range(dvhist,i,&leftbin,&rightbin); cerr<<leftbin<<" ";}	cerr<<endl;
//			for(size_t i=0; i<bins; i++) cerr<<gsl_histogram_get(dvhist,i)<<" "; cerr<<endl;
//		}
//		gsl_histogram_free(dvhist);
//	}
//
//	err = pop.get_diversity_histogram(&dshist);
//	if(err == 0) {
//		bins = gsl_histogram_bins(dshist);
//		if(HIGHD_VERBOSE) {
//			cerr<<"Diversity histogram:"<<endl;
//			for(size_t i=0; i<bins; i++) {gsl_histogram_get_range(dshist,i,&leftbin,&rightbin); cerr<<leftbin<<" ";}	cerr<<endl;
//			for(size_t i=0; i<bins; i++) cerr<<gsl_histogram_get(dshist,i)<<" "; cerr<<endl;
//		}
//		gsl_histogram_free(dshist);
//	}
//
//	return err;
//}
//
//

/* MAIN */
int main(int argc, char **argv){

	int status= 0;
	if (argc > 1) {
		cout<<"Usage: "<<argv[0]<<endl;
		status = 1;
	} else {
//		status = library_access();
//		status += sample_initialize();
//		status += hc_initialize();
//		status += hc_setting();
//		status += pop_initialize();
		status += pop_evolve();
//		status += pop_sampling();
//		status += pop_Hamming();
//		status += pop_divdiv();
//		status += pop_histograms();

	}
	cout<<"Number of errors: "<<status<<endl;
	return status;
}


