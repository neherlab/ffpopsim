/**
 * @file hivpython.cpp
 * @brief Implementation of an HIV population with drug treatment.
 * @author Richard Neher, Boris Shraiman, Fabio Zanini
 * @version 
 * @date 2012-04-23
 */

#include "popgen.h"
#include "popgen_highd.h"
#include "hivpopulation.h"
#include "hivpython.h"
#include <iostream>
#include <fstream>

hivpython::hivpython() {
}


hivpython::~hivpython() {
}

int hivpython::evolve(int gen) {
	int err=haploid_clone::evolve(gen);
	if(err==0)
		haploid_clone::calc_stat();
	return err;
}

int hivpython::read_selection_coefficients(char *model){
	ifstream modelstream(model);
	return hivpopulation::read_selection_coefficients(modelstream);
}

int hivpython::read_resistance_coefficients(char *model){
	ifstream modelstream(model);
	return hivpopulation::read_resistance_coefficients(modelstream);
}

void hivpython::get_allele_frequencies(double *af) {
	for(size_t i=0; i < number_of_loci; i++)
		af[i] = haploid_clone::get_allele_frequency(i);
}

void hivpython::random_clone(unsigned short *seq) {
	boost::dynamic_bitset<> newgt = (*current_pop)[haploid_clone::random_clone()].genotype;
	for(size_t i=0; i < number_of_loci; i++)
		seq[i] = newgt[i];
}


//int hivpython::read_selection_coefficients(istream &model){
//	if (HIVPOP_VERBOSE){
//		cerr<<"hivpython::read_selection_coefficients(): read coefficients ";
//	}
//	if (model.bad()){
//		cerr<<"hivpython::read_selection_coefficients(): BAD MODEL STREAM!"<<endl;
//		return HIVPOP_BADARG;
//	}
//	double val;
//	vector <int> loci;
//	vector<string> strs;
//	string line;
//	while(!model.eof()){
//		strs.clear();
//		getline(model, line);
//		boost::split(strs, line, boost::is_any_of("\t "));
//		if (strs.size()>1){
//			for (unsigned int entry=0; entry<strs.size()-1; entry++){
//				loci.push_back(atoi(strs[entry].c_str()));
//			}
//			val=atof(strs.back().c_str());
//			add_fitness_coefficient(val, loci);
//			if (HIVPOP_VERBOSE) cerr<<loci[0]<<" "<<val<<"  "<<loci.size()<<endl;
//			loci.clear();
//		}
//	}
//	if (HIVPOP_VERBOSE) cerr<<"...done"<<endl;
//	return 0;
//}
//
//
//int hivpython::read_resistance_coefficients(istream &model){
//	if (HIVPOP_VERBOSE){
//		cerr<<"hivpython::read_resistance_coefficients(): read coefficients ";
//	}
//	if (model.bad()){
//		cerr<<"hivpython::read_resistance_coefficients(): BAD MODEL STREAM!"<<endl;
//		return HIVPOP_BADARG;
//	}
//	double val, wt_resistance=0;
//	vector <int> loci;
//	vector<string> strs;
//	string line;
//	while(!model.eof()){
//		strs.clear();
//		loci.clear();
//		getline(model, line);
//		boost::split(strs, line, boost::is_any_of("\t "));
//		//cout <<"a "<<line<<"  "<<strs.size();
//		if (strs.size()>1){
//			for (unsigned int entry=0; entry<strs.size()-1; entry++){
//				loci.push_back(atoi(strs[entry].c_str()));
//				//cout<<loci.back()<<" "<<strs[entry].c_str()<<"  ";
//			}
//			val=atof(strs.back().c_str());
//			add_fitness_coefficient(val, loci,1);
//			wt_resistance+=val*pow(-1.0,(double)loci.size());
//			//cout <<loci.size()<<"  "<<wt_resistance<<endl;
//		}
//		//cout<<loci[0]<<" "<<val<<"  "<<loci.size()<<endl;
//	}
//	trait[1].hypercube_mean=-wt_resistance;
//	if (HIVPOP_VERBOSE){
//		cerr<<"...done"<<endl;
//	}
//	return 0;
//}
//
//
//int hivpython::write_genotypes(ostream &out, int sample_size, string gt_label, int start, int length){
//	if (out.bad()){
//		cerr<<"hivpython::write_genotypes(): BAD OUTPUT FILE!"<<endl;
//		return HIVPOP_BADARG;
//	}else{
//		int gti;
//		int string_length;
//		string temp;
//		if (length>0)
//			string_length = length;
//		else
//			string_length = number_of_loci - start;
//
//		produce_random_sample(sample_size);
//		if (sample_size>get_pop_size()){
//			cerr<<"hivpython::write_genotypes(): requested sample size exceeds population size"<<endl;
//			return HIVPOP_BADARG;
//		}else{
//			for (int s=0; s<sample_size; s++){
//				gti=random_clone();
//				out <<">GT-"<<gt_label<<"_"<<gti<<'\n';
//				for (int i =start; i<start+length; i++ ){
//					if ((*current_pop)[gti].genotype[i]) out <<'1';
//					else out <<'0';
//				}
//				out<<'\n';
//				//out <<get_genotype_string(gti).substr(start,string_length)<<'\n';
//			}
//		}
//		return 0;
//	}
//}
//
