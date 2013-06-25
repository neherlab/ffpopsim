/**
 * @file hivpopulation.cpp
 * @brief Implementation of an HIV population with drug treatment.
 * @author Richard Neher, Fabio Zanini
 * @version 
 * @date 2012-04-23
 *
 * Copyright (c) 2012-2013, Richard Neher,  Fabio Zanini
 * All rights reserved.
 *
 * This file is part of FFPopSim.
 *
 * FFPopSim is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * FFPopSim is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with FFPopSim. If not, see <http://www.gnu.org/licenses/>.
 */
#include "hivpopulation.h"

/**
 * @brief Construct a HIV population with certain parameters
 *
 * @param N_in number of viral particles 
 * @param rng_seed seed for the random number generator. If this is 0, time(NULL)+getpid() is used.
 * @param mutrate mutation rate in events / generation / site
 * @param coinfection_rate probability of coinfection of the same cell by two viral particles in events / generation
 * @param crossover_rate probability of template switching during coinfection in events / site
 *
 * *Note*: the genome length is 10000 (see HIVGENOME).
 * *Note*: exceptions are propagated from the base class constructor (haploid_highd).
 */
hivpopulation::hivpopulation(int N_in, int rng_seed, double mutation_rate_in, double coinfection_rate_in, double crossover_rate_in) : haploid_highd(HIVGENOME, rng_seed, 2),
	            gag(GAG_START, GAG_END),
	            pol(POL_START, POL_END),
	            env(ENV_START, ENV_END),
	            nef(NEF_START, NEF_END),
	            vif(VIF_START, VIF_END),
	            vpu(VPU_START, VPU_END),
	            vpr(VPR_START, VPR_END),
	            tat(TAT1_START, TAT1_END, TAT2_START, TAT2_END),
	            rev(REV1_START, REV1_END, REV2_START, REV2_END),
		    treatment(0) {

	outcrossing_rate = coinfection_rate_in;
	mutation_rate = mutation_rate_in;
	crossover_rate = crossover_rate_in;
	recombination_model = CROSSOVERS;

	//by default, create a population of size carrying capacity at 00...0 (this is cheap, a single clone)
        if(N_in > 0)
        	set_wildtype(N_in);
}

/**
 * @brief Destructor.
 *
 * Only calls the method of the base class (which manages its own memory).
 */
hivpopulation::~hivpopulation() {
}


/**
 * @brief Set the fitness from replication and resistance for a single clone
 *
 * @param tempgt clone whose fitness is calculated
 */
void hivpopulation::calc_individual_fitness_from_traits(clone_t *tempgt) {
	tempgt->fitness = tempgt->trait[0] + treatment * tempgt->trait[1];
}


int hivpopulation::read_replication_coefficients(istream &model){
	if (HIVPOP_VERBOSE){
		cerr<<"hivpopulation::read_selection_coefficients(): read coefficients ";
	}
	if (model.bad()){
		cerr<<"hivpopulation::read_selection_coefficients(): BAD MODEL STREAM!"<<endl;
		return HIVPOP_BADARG;
	}
	double val;
	vector <int> loci;
	vector<string> strs;
	string line;

	// reset the hypercube
	trait[0].reset();
	
	// read the stream
	while(!model.eof()){
		strs.clear();
		getline(model, line);
		boost::split(strs, line, boost::is_any_of("\t "));
		if (strs.size()>1){
			for (unsigned int entry=0; entry<strs.size()-1; entry++){
				loci.push_back(atoi(strs[entry].c_str()));
			}
			val=atof(strs.back().c_str());
			add_trait_coefficient(val, loci, 0);
			if (HIVPOP_VERBOSE) cerr<<loci[0]<<" "<<val<<"  "<<loci.size()<<endl;
			loci.clear();
		}
	}

	// update the replication and fitness of all clones
	update_traits();
	update_fitness();

	if (HIVPOP_VERBOSE) cerr<<"...done"<<endl;
	return 0;
}


int hivpopulation::read_resistance_coefficients(istream &model){
	if (HIVPOP_VERBOSE){
		cerr<<"hivpopulation::read_resistance_coefficients(): read coefficients ";
	}
	if (model.bad()){
		cerr<<"hivpopulation::read_resistance_coefficients(): BAD MODEL STREAM!"<<endl;
		return HIVPOP_BADARG;
	}
	double val, wt_resistance=0;
	vector <int> loci;
	vector<string> strs;
	string line;

	// reset the hypercube
	trait[1].reset();
	
	// read the stream
	while(!model.eof()){
		strs.clear();
		loci.clear();
		getline(model, line);
		boost::split(strs, line, boost::is_any_of("\t "));
		//cout <<"a "<<line<<"  "<<strs.size();
		if (strs.size()>1){
			for (unsigned int entry=0; entry<strs.size()-1; entry++){
				loci.push_back(atoi(strs[entry].c_str()));
				//cout<<loci.back()<<" "<<strs[entry].c_str()<<"  ";
			}
			val=atof(strs.back().c_str());
			add_trait_coefficient(val, loci,1);
			wt_resistance+=val*pow(-1.0,(double)loci.size());
			//cout <<loci.size()<<"  "<<wt_resistance<<endl;
		}
		//cout<<loci[0]<<" "<<val<<"  "<<loci.size()<<endl;
	}
	trait[1].hypercube_mean=-wt_resistance;

	// update the replication and fitness of all clones
	update_traits();
	update_fitness();

	if (HIVPOP_VERBOSE){
		cerr<<"...done"<<endl;
	}
	return 0;
}


int hivpopulation::write_genotypes(ostream &out, int sample_size, string gt_label, int start, int length){
	if (HIVPOP_VERBOSE) cerr<<"hivpopulation::write_genotypes()...";
	if (HIVPOP_VERBOSE) cerr<<"start = "<<start<<"...";
	if (HIVPOP_VERBOSE) cerr<<"length = "<<length<<"...";

	if (out.bad()){
		cerr<<"hivpopulation::write_genotypes(): BAD OUTPUT FILE!"<<endl;
		return HIVPOP_BADARG;
	}else{
		int gti;
		string temp;
		if (length <= 0)
			length = number_of_loci - start;

		produce_random_sample(sample_size);
		if (sample_size>get_population_size()){
			cerr<<"hivpopulation::write_genotypes(): requested sample size exceeds population size"<<endl;
			return HIVPOP_BADARG;
		}else{
			for (int s=0; s<sample_size; s++){
				gti=random_clone();
				out <<">GT-"<<gt_label<<"_"<<gti<<'\n';
				for (int i =start; i<start+length; i++ ){
					if (population[gti].genotype[i]) out <<'1';
					else out <<'0';
				}
				out<<'\n';
			}
		}
		if (HIVPOP_VERBOSE) cerr<<"...done."<<endl;
		return 0;
	}
}
