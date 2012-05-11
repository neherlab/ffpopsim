/**
 * @file hivpopulation.cpp
 * @brief Implementation of an HIV population with drug treatment.
 * @author Richard Neher, Boris Shraiman, Fabio Zanini
 * @version 
 * @date 2012-04-23
 */
#include "popgen_highd.h"
#include "hivpopulation.h"

/**
 * @brief Default constructor.
 *
 * Only calls the method of the base class.
 */
hivpopulation::hivpopulation() {
}

/**
 * @brief Destructor.
 *
 * Only calls the method of the base class (which manages its own memory).
 */
hivpopulation::~hivpopulation() {
}


/**
 * @brief Construct a HIV population with certain parameters
 *
 * @param N_in number of viral particles 
 * @param rng_seed seed for the random number generator. If this is 0, time(NULL)+getpid() is used.
 * @param mutrate mutation rate in events / generation / site
 * @param coinfection_rate probability of coinfection of the same cell by two viral particles in events / generation
 * @param crossover_rate probability of template switching during coinfection in events / site
 *
 * @returns zero if successful, error codes otherwise
 *
 * Note: the genome length is 10000 (see HIVGENOME).
 */
int hivpopulation::set_up(int N_in, int rng_seed, double mutation_rate_in, double coinfection_rate_in, double crossover_rate_in){
	int err=set_up(N_in, HIVGENOME, rng_seed, 2); // we have exactly 2 traits
	outcrossing_rate = coinfection_rate_in;
	mutation_rate = mutation_rate_in;
	crossover_rate = crossover_rate_in;
	recombination_model = CROSSOVERS;
	treatment = 0;
	init_genotypes();
	return err;
}

int hivpopulation::read_selection_coefficients(istream &model){
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
	while(!model.eof()){
		strs.clear();
		getline(model, line);
		boost::split(strs, line, boost::is_any_of("\t "));
		if (strs.size()>1){
			for (unsigned int entry=0; entry<strs.size()-1; entry++){
				loci.push_back(atoi(strs[entry].c_str()));
			}
			val=atof(strs.back().c_str());
			add_fitness_coefficient(val, loci);
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
			add_fitness_coefficient(val, loci,1);
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
