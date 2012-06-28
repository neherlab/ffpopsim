/**
 * @file hivpopulation.cpp
 * @brief Implementation of an HIV population with drug treatment.
 * @author Richard Neher, Fabio Zanini
 * @version 
 * @date 2012-04-23
 *
 * Copyright (c) 2012, Richard Neher, Boris Shraiman, Fabio Zanini
 * All rights reserved.

 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *    * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.

 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
hivpopulation::hivpopulation(int N_in, int rng_seed, double mutation_rate_in, double coinfection_rate_in, double crossover_rate_in) : haploid_highd(HIVGENOME, rng_seed, 2), env(ENV_START, ENV_END), treatment(0) {

	outcrossing_rate = coinfection_rate_in;
	mutation_rate = mutation_rate_in;
	crossover_rate = crossover_rate_in;
	recombination_model = CROSSOVERS;

	//by default, create a population of size carrying capacity at 00...0 (this is cheap, a single clone)
	set_wildtype(N_in);
}

/**
 * @brief Destructor.
 *
 * Only calls the method of the base class (which manages its own memory).
 */
hivpopulation::~hivpopulation() {
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
					if ((*current_pop)[gti].genotype[i]) out <<'1';
					else out <<'0';
				}
				out<<'\n';
			}
		}
		if (HIVPOP_VERBOSE) cerr<<"...done."<<endl;
		return 0;
	}
}
