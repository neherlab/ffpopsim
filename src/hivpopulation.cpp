/**
 * @file hivpopulation.cpp
 * @brief Implementation of an HIV population with drug treatment.
 * @author Richard Neher, Boris Shraiman, Fabio Zanini
 * @version 
 * @date 2012-04-23
 */

#include "popgen.h"
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
int hivpopulation::set_up(int N_in, int rng_seed, double mutrate, double coinfection_rate, double crossover_rate){
	int err=set_up(N_in, HIVGENOME, rng_seed, 2); // we have exactly 2 traits
	outcrossing_probability = coinfection_rate;
	mutation_rate = mutrate;
	crossover_rate = crossover_rate;
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

/*
 * Calculate mean and variance of the divergence from the [00...0] bitset
 */
stat_t hivpopulation::get_divergence()
{
	int n_sample = 1000;
	
	// Cumulative partition the population according to clones
	vector <unsigned long> partition_cum = partition_cumulative();

	// Calculate divergence of the sample
	boost::dynamic_bitset<> genotype;
	stat_t divergence;
	double tmp;
	for (int i = 0; i < n_sample; i++, tmp = 0) {
		// Calculate divergence
		genotype = (*current_pop)[random_clone()].genotype;
		for (int j = 0; j < genotype.size(); tmp += genotype[j++]);
		divergence.mean += tmp;
		divergence.variance += tmp * tmp;
	}
	divergence.mean /= n_sample;
	divergence.variance /= n_sample;
	divergence.variance -= divergence.mean * divergence.mean;
	return divergence;
}

/*
 * Calculate diversity in the current population, i.e. Hamming distance between all pairs of sequences, with mean and variance.
 */
stat_t hivpopulation::get_diversity()
{
	int n_sample = 1000;
	
	// Calcolate random distances
	stat_t diversity;
	unsigned long c1, c2, tmp;
	for (int i = 0; i < n_sample; i++) {
		c1 = random_clone();
		c2 = random_clone();
		// Calculate distance if they belong to different clones
		if (c1 != c2 ) {
			tmp = distance_Hamming((*current_pop)[c1].genotype,(*current_pop)[c2].genotype);
			diversity.mean += tmp;
			diversity.variance += tmp * tmp;
		}
	}
	diversity.mean /= n_sample;
	diversity.variance /= n_sample;
	diversity.variance -= diversity.mean * diversity.mean;
	return diversity;
}

/*
 * Get histogram of divergence, to ingestivate more in detail than just mean and variance
 */
int hivpopulation::get_divergence_histogram(gsl_histogram **hist, int bins, int position, int begin, int end)
{
	/* Note: position is in {1,2,3} */
	int n_sample = 1000;
	int i,j;
	int Ltmp = number_of_loci;

	// Set begin, end, and jump in comparisons
	if (end < 0) end = Ltmp - 1;			// If end is not specified, compare till the end
	int jump = 1;
	if (position) {
		begin += ((Ltmp - Ltmp%3 + (position-1) - begin) % 3);	// Start at the first compatible position (avoid negative numbers)
		end -= ((3 + end - (position-1)) % 3);	// Stop at the last compatible position
		jump = 3;
	}
	if (end - begin < 0) return -2;			// Check boundaries

	// Calculate divergence of the sample
	vector <unsigned long> divergence(n_sample, 0);
	unsigned long tmp;
	boost::dynamic_bitset<> genotype;
	for (i = 0; i < n_sample; i++, tmp = 0) {
		// Calculate divergence
		genotype = (*current_pop)[random_clone()].genotype;
		for (j = begin, tmp = 0; j <= end; j+=jump)
			tmp += genotype[j];
		divergence[i] = tmp;
	}

	// Prepare the histogram
	unsigned long dmax = *max_element(divergence.begin(),divergence.end());
	unsigned long dmin = *min_element(divergence.begin(),divergence.end());

	/* Aliasing: see get_diversity_histogram */
	int width, binsnew;
	if (dmin == dmax)
		width = 1;
	else {
		width = (dmax - dmin) / (bins-1);
		width += ((dmax - dmin)%(bins-1))?1:0;
	}
	binsnew = ((dmax - dmin) / width) + 1;
	if (binsnew > bins) {
		cout<<"binsnew: "<<binsnew<<"\tbins: "<<bins<<endl;
		return -1;
	}

	// Fill and scale histogram
	*hist = gsl_histogram_alloc(binsnew); 
	gsl_histogram_set_ranges_uniform(*hist, dmin - 0.5 * width, dmax + 0.5 * width);
	for (i = 0; i < n_sample; i++)
		gsl_histogram_increment(*hist, divergence[i]);
	gsl_histogram_scale(*hist, 1/(double)n_sample);
	
	return 0;
}

/*
 * Get histogram of diversity, to ingestivate more in detail than just mean and variance
 */
int hivpopulation::get_diversity_histogram(gsl_histogram **hist, int bins, int position, int begin, int end)
{
	int n_sample = 1000;
	int i;

	// Calculate diversity of the sample
	vector <unsigned long> diversity(n_sample, 0);
	unsigned long c1, c2;
	clone_t *tempgt1, *tempgt2;
	for(i=0; i < n_sample; i++) {
		c1 = random_clone();
		c2 = random_clone();
		if (c1 != c2)
			diversity[i] = distance_Hamming((*current_pop)[c1].genotype,(*current_pop)[c2].genotype, position, begin, end);
	}

	// Prepare the histogram
	unsigned long dmax = *max_element(diversity.begin(),diversity.end());
	unsigned long dmin = *min_element(diversity.begin(),diversity.end());

	/* Aliasing: change the number of bins and the width as needed.
	 * This induces a (negative) bias in the last bin, but we cannot do anything for it;
	 * moreover, for memory reasons we must output a vector of the same size, so we set
	 * the rest equal to the top leftborder, and the counts at zero (note that binsnew <= bins) */
	int width, binsnew;
	if (dmin == dmax)
		width = 1;
	else {
		width = (dmax - dmin) / (bins -1);
		width += ((dmax - dmin)%(bins-1))?1:0;
	}
	binsnew = ((dmax - dmin) / width) + 1;
	if (binsnew > bins) {
		cout<<"binsnew: "<<binsnew<<"\tbins: "<<bins<<"\tdelta: "<<(dmax - dmin)<<endl;
		return -1;
	}
	
	// Fill and scale histogram
	*hist = gsl_histogram_alloc(binsnew); 
	gsl_histogram_set_ranges_uniform(*hist, dmin - 0.5 * width, dmax + 0.5 * width);
	for (i = 0; i < n_sample; i++)
		gsl_histogram_increment(*hist, diversity[i]);
	gsl_histogram_scale(*hist, 1/(double)n_sample);
	
	return 0;
}


/**
 * @brief calculate histogram of fitness from traits.
 *
 * There is a small problem here, namely that the fitness distribution has a horrible tail which
 * messes up the calculation of the bin width. Thus we first calculate the fitness average and
 * variance, and then set the bin width so that the deleterious tail is within 2 sigma or so.
 */
int hivpopulation::get_fitness_histogram(gsl_histogram **hist, int bins) {
	int n_sample = 1000;
	int i;

	// Calculate fitness of the sample
	double fitnesses[n_sample];
	for(i=0;i<n_sample;i++)
		fitnesses[i] = (*current_pop)[random_clone()].fitness;

	// Set the bins according to average and variance in fitness in the population
	calc_fitness_stat();
	double fitmean = fitness_stat.mean;
	double fitstd = sqrt(fitness_stat.variance);
	double histtail = 2 * fitstd;

	// Prepare histogram (here antialiasing is not a problem)
	double fmax = *max_element(fitnesses,fitnesses+n_sample);
	double fmin = fitmean - histtail;
	bins = min(n_sample / 30, bins);
	double width = (fmax - fmin) / (bins - 1);
	*hist = gsl_histogram_alloc(bins); 
	gsl_histogram_set_ranges_uniform(*hist, fmin - 0.5 * width, fmax + 0.5 * width);

	// Fill and scale histogram
	for (i = 0; i < n_sample; i++)
		gsl_histogram_increment(*hist, fitnesses[i]);
	gsl_histogram_scale(*hist, 1/(double)n_sample);

	return 0;
}


/**
 * @brief calculate the hamming distance between two sequences (not normalized)
 */
unsigned int hivpopulation::distance_Hamming(boost::dynamic_bitset<> genotype, boost::dynamic_bitset<> genotype1, int position, int begin, int end)
{
	unsigned int d = 0;
	int Ltmp = number_of_loci;
	// Set begin, end, and jump in comparisons
	if (end < 0) end = Ltmp - 1;			// If end is not specified, compare till the end
	int jump = 1;
	if (position) {
		begin += ((Ltmp - Ltmp%3 + (position-1) - begin) % 3);	// Start at the first compatible position (avoid negative numbers)
		end -= ((3 + end - (position-1)) % 3);	// Stop at the last compatible position
		jump = 3;
	}
	if (end - begin < 0) return -2;			// Check boundaries
	
	for (int i = begin; i < end; i+=jump)
		d += (genotype[i] != genotype1[i]);
	return d;
}

/*
 * Calculate the cumulative partition of sequences in the clones
 */
vector <unsigned long> hivpopulation::partition_cumulative()
{
	vector <unsigned long> partition_cum;
	unsigned long tmp;
	partition_cum.push_back((*current_pop)[0].clone_size);		
	for (int i = 1; i < get_number_of_clones(); i++) {
		partition_cum.push_back((*current_pop)[i].clone_size + partition_cum[i-1]);		
	}
	return partition_cum;
}


int hivpopulation::write_genotypes(ostream &out, int sample_size, string gt_label, int start, int length){
	if (out.bad()){
		cerr<<"hivpopulation::write_genotypes(): BAD OUTPUT FILE!"<<endl;
		return HIVPOP_BADARG;
	}else{
		int gti;
		int string_length;
		string temp;
		if (length>0)
			string_length = length;
		else
			string_length = number_of_loci - start;

		produce_random_sample(sample_size);
		if (sample_size>get_pop_size()){
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
				//out <<get_genotype_string(gti).substr(start,string_length)<<'\n';
			}
		}
		return 0;
	}
}

