/**
 * @file highd.cpp
 * @brief Tests for the high-dimensional simulation library.
 * @author Richard Neher, Boris Shraiman, Fabio Zanini
 * @version 
 * @date 2012-04-20
 */

/* Include directives */
#include "hivpopulation.h"
#define HIVPOP_BADARG -1354341
#define NOTHING 1e-10

/* Be verbose? */
#define HIV_VERBOSE 1

/* Test end-user subclass initialization */
int hiv_initialize() {

	int N = 1000;
	hivpopulation pop(N);

	cout<<"Env start: "<<pop.env.start<<", env end: "<<pop.env.end<<endl;


	if(HIV_VERBOSE)
		cerr<<"L = "<<pop.get_number_of_loci()<<", N = "<<pop.carrying_capacity<<endl;

	if(HIV_VERBOSE) {
		cerr<<"Mutation rate: "<<pop.get_mutation_rate()<<", ";
		cerr<<"coinfection rate: "<<pop.outcrossing_rate<<endl;
	}
	return 0;
}

/* Test end-user subclass evolution and output */
int hiv_evolve() {
	int N = 1000;
	ifstream model("hiv_model.dat", ifstream::in);

	hivpopulation pop(N, 0, 2e-5, 1e-3, 1e-3);

	if(HIV_VERBOSE) cerr<<"Reading model coefficients...";
	pop.read_replication_coefficients(model);
	if(HIV_VERBOSE) cerr<<"read!"<<endl;

	// err checks for extinction
	int err = pop.evolve(10);
	if(err==0) {
		vector <int> sample;
		pop.random_clones(10, &sample);

		if(HIV_VERBOSE) {
			cerr<<"Number of clones: "<<pop.get_number_of_clones()<<endl;
			cerr<<"Random individuals:";
			for(unsigned int i=0; i < sample.size(); i++)
				cout<<" "<<sample[i];
			cout<<endl;
		}
	}
	return bool(err);
}


/* Test evolution of different populations on the same landscape, to test random effects */
int hiv_multiple_evolution() {
	int N = 1000, err=0;
	ifstream model("hiv_model.dat", ifstream::in);

	for(size_t i=0; i < 5; i++) {
		if(HIV_VERBOSE) cerr<<"Population n. "<<i<<"...";
		hivpopulation pop(N, 0, 2e-5, 1e-3, 1e-3);

		if(HIV_VERBOSE) cerr<<"Reading model coefficients...";
		pop.read_replication_coefficients(model);
		if(HIV_VERBOSE) cerr<<"read!"<<endl;

		for(size_t j = 0; j < 5; j++)
			err += pop.evolve(50);

		if(HIV_VERBOSE) {
			int sites[] = {5, 5000, 9000};
			cerr<<"             \t";
			for(size_t j = 0; j < 3; j++) {
				cerr<<sites[j];
				if (j!=2)
					cerr<<"\t";
			}
			cerr<<endl<<"-----------------------------------------------";
			cerr<<endl<<"allele freqs:\t";
			for(size_t j = 0; j < 3; j++) {
				cerr<<pop.get_allele_frequency(sites[j]);
				if (j!=2)
					cerr<<"\t";
			}
			cerr<<endl;
		}
		if(HIV_VERBOSE) cerr<<"done"<<endl;
	}

	return err;
}


/* Test genes of HIV */
int hiv_genes() {
	int N = 1000, err=0;
	ifstream model("hiv_model.dat", ifstream::in);

	hivpopulation pop(N, 0, 2e-5, 1e-3, 1e-3);

	pop.read_replication_coefficients(model);

	if(HIV_VERBOSE) {
		cerr<<"pol:\t"<<pop.pol.start<<"\t"<<pop.pol.end<<endl;
		cerr<<"env:\t"<<pop.env.start<<"\t"<<pop.env.end<<endl;
		cerr<<"nef:\t"<<pop.nef.start<<"\t"<<pop.nef.end<<endl;
		cerr<<"tat:\t"<<pop.tat.start<<"\t"<<pop.tat.end<<"\t";
		cerr<<pop.tat.second_start<<"\t"<<pop.tat.second_end<<endl;
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
		status += hiv_initialize();
		status += hiv_evolve();
		status += hiv_multiple_evolution();
		status += hiv_genes();
	}
	cout<<"Number of errors: "<<status<<endl;
	return status;
}


