/**
 * @file highd.cpp
 * @brief Profile functions for the simulation library.
 * @author Richard Neher, Fabio Zanini
 * @version 
 * @date 2012-04-20
 */
/* Include directives */
#include <boost/progress.hpp>
#include "hivpopulation.h"

/* Defines */
#define HIGHD_BADARG -1354341
#define NOTHING 1e-10

/* Be verbose? */
#define PROFILE_VERBOSE 1

/* Declarations */
int hiv_profile();


/* MAIN */
int main(int argc, char **argv){
	boost::progress_timer timer;

	int status= 0;
	if (argc > 1) {
		cout<<"Usage: "<<argv[0]<<endl;
		status = 1;
	} else {
                for(int i=0; i < 5; i++)
		status += hiv_profile();
	}
	cout<<"Number of errors: "<<status<<endl;
	return status;
}


int hiv_profile() {
	int N = 1000, err=0;
	ifstream model;

	for(size_t i=0; i < 1; i++) {
		if(PROFILE_VERBOSE) cerr<<"Population n. "<<i<<"...";
		model.open("../tests/hiv_model.dat", ifstream::in);
		hivpopulation pop(N, 1, 2e-5, 1e-3, 1e-3);
	
		if(PROFILE_VERBOSE) cerr<<"Reading model coefficients...";
		pop.read_replication_coefficients(model);
		if(PROFILE_VERBOSE) cerr<<"read!"<<endl;
	
		for(size_t j = 0; j < 5; j++) {
			pop.get_max_fitness();
			err += pop.evolve(100);
		}
		if(PROFILE_VERBOSE) cerr<<"done"<<endl;

		cout<<"Fitnesses: ";
		for(size_t k=0; k!= 5; k++)	{
			cout<<pop.get_fitness(k)<<"  ";		
		}
		cout<<endl;
		model.close();
	}

	return err;
}
