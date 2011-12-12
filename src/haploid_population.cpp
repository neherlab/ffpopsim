

/*
 * haploid_population.cpp
 *
 *  Created on: Aug 28, 2008
 *      Author: neher
 */

#include "popgen.h"
int comp_labels_hp(const void *p1, const void *p2)
{
	label_gt_pair *gt1=(label_gt_pair*)p1;
	label_gt_pair *gt2=(label_gt_pair*)p2;
	if (gt1->label>gt2->label) return 1;
	else if  (gt1->label<gt2->label) return -1;
	else return 0;
}


haploid_population::haploid_population() {
	mem=false;
	clone_size_mem=false;
	cumulants_mem=false;
	circular=false;
}

haploid_population::~haploid_population() {
	if (mem) free_mem();
}

int haploid_population::set_up(int N_in,int L_in,  int rng_seed)
{
	if (N_in<1 or L_in <1)
	{
		cerr <<"haploid_population::set_up(): Bad Arguments!\n";
		return HP_BADARG;
	}
	number_of_loci=L_in;
	number_of_individuals=2*N_in;
	target_pop_size=N_in;
	number_of_words=int(ceil(double(number_of_loci)/WORD_LENGTH));  //number of words used to store genotypes

	if (rng_seed==0){
		seed=time(NULL)+getpid();
	}else{seed=rng_seed;}

	mem=false;
	return allocate_mem();
}

int haploid_population::allocate_mem()
{
	if (mem)
	{
		cerr <<"haploid_population::allocate_mem(): memory already allocated, freeing and reallocating ...!\n";
		free_mem();
	}

	if (HP_VERBOSE) cerr <<"haploid_population::allocate_mem(): genotypes: "<<number_of_individuals<<"  loci: "<<number_of_loci<<endl;

	number_of_bits=new int [number_of_words];
	for(int i=0; i<number_of_words; i++)	number_of_bits[i]=WORD_LENGTH;
	number_of_bits[number_of_words-1]=number_of_loci-(number_of_words-1)*WORD_LENGTH;

	//allocate all the memory
	if (HP_VERBOSE) cerr <<"genotypes*...";
	genotypes= new int [number_of_individuals*number_of_words];
	new_genotypes= new int [number_of_individuals*number_of_words];
	if (HP_VERBOSE) cerr <<"genotypes...";

	if (HP_VERBOSE) cerr <<"gametes...";
	sex_gametes=new int [number_of_individuals];
	asex_gametes=new int [number_of_individuals];
	genome = new int [number_of_loci+1];
	crossovers= new int [number_of_loci];
	for (int i=0; i<number_of_loci; i++) genome[i]=i;

	if (HP_VERBOSE) cerr <<"fitness...";
	fitness_values = new sample;
	new_fitness_values = new sample;
	fitness_values->set_up(number_of_individuals);
	new_fitness_values->set_up(number_of_individuals);
	if (HP_VERBOSE) cerr <<"outcrossing_rates...";
	outcrossing_rates=new sample;
	new_outcrossing_rates=new sample;
	outcrossing_rates->set_up(number_of_individuals);
	new_outcrossing_rates->set_up(number_of_individuals);
	outcrossing_rates->set_range(0,1);
	new_outcrossing_rates->set_range(0,1);
	if (HP_VERBOSE) cerr <<"genotype labels...";
	gt_labels=new sample;
	new_gt_labels=new sample;
	gt_labels->set_up(number_of_individuals);
	new_gt_labels->set_up(number_of_individuals);



	if (HP_VERBOSE) cerr <<"allele frequencies...";
	allele_frequencies =new double [number_of_loci];
	gamete_allele_frequencies =new double [number_of_loci];

	//array of numbers need to choose individuals from.
	numbers= new int [number_of_individuals];
	for (int i = 0; i < number_of_individuals; ++i) {
		numbers[i]=i;
	}

	//Random number generator
	evo_generator=gsl_rng_alloc(RNG);
	gsl_rng_set(evo_generator, seed);
	label_generator=gsl_rng_alloc(RNG);
	fitness.set_up(number_of_loci, gsl_rng_uniform_int(evo_generator, 10000000));
	//outcrossing_rates.set_up(number_of_loci, number_of_words, gsl_rng_uniform_int(rng_seed, 100000000));

	mem=true;
	if (HP_VERBOSE) cerr <<"done.\n";
	generation=0;
	return 0;
}

int haploid_population::allocate_clone_size_distribution(){
	if ((mem==true) and (clone_size_mem==false)){
		clone_size_distribution = new int [2*number_of_individuals];
		if (clone_size_distribution!=NULL) clone_size_mem=true;
		return 0;
	}
	else return HP_MEMERR;
}

int haploid_population::set_evolve_outcrossing_rates(){
	evolve_outcrossing_rates=true;
	outcrossing_rate.set_up(number_of_loci, gsl_rng_uniform_int(evo_generator, 100000000));
}

int haploid_population::free_mem()
{
	if (!mem)
	{
		cerr <<"haploid_population::free_mem(): No memory allocated!\n";
		return HP_BADARG;
	}
	delete [] genotypes;
	delete [] new_genotypes;
	delete [] sex_gametes;
	delete [] asex_gametes;
	delete [] allele_frequencies;
	delete [] gamete_allele_frequencies;
	mem=false;
	return 0;
}

// create genotypes at random from allele frequencies saved in nu[k] k=0..number_of_loci-1
int haploid_population::init_genotypes(double* nu, int n_o_genotypes)
{
	int i,w,l, locus;
	if (!mem)
	{
		cerr <<"haploid_population::init_genotypes(): Allocate and set up first!\n";
		return HP_MEMERR;
	}
	if (n_o_genotypes<0 or n_o_genotypes>=number_of_individuals)
	{
		cerr <<"haploid_population::init_genotypes(): number of genotypes has to be positive and smaller than max! got: "<<n_o_genotypes<<"!\n";
		return HP_BADARG;
	}
	if (HP_VERBOSE) cerr <<"haploid_population::init_genotypes(double* nu, int n_o_genotypes) ....";
	generation=0;
	if (n_o_genotypes==0) pop_size=target_pop_size;
	else pop_size=n_o_genotypes;
	fitness_values->number_of_values=pop_size;
	outcrossing_rates->number_of_values=pop_size;
	gt_labels->number_of_values=pop_size;
	for (i=0; i<pop_size; i++)
	{
		locus=0;	//loop over the number of words, use the cumulative locus index for nu[i]
		for (w=0; w<number_of_words; w++)
		{
			genotypes[index(i, w)]=0;		//put together the genotype: + state with prob nu[i]
			for (l=0; l<number_of_bits[w]; l++)
			{
				if (gsl_rng_uniform(evo_generator)<nu[locus])		genotypes[index(i,w)]+=(1<<l);
				locus++;
			}
		}
	}
	for (i=pop_size; i<number_of_individuals; i++)
	{
		for (w=0; w<number_of_words; w++) genotypes[index(i,w)]=NO_GENOTYPE;
	}
	//calculate its fitness and recombination rates
	calc_fit();
	calc_rec();
	calc_gt_labels();
	calc_stat();
	if (HP_VERBOSE) cerr <<"done.\n";
	//set generations counter to zero and calculate the statistics for the intial configuration
	generation=0;
	return 0;
}


//init with all allele frequencies equal to 0.5 -- n copies of each inital draw-default is 1
int haploid_population::init_genotypes_diverse(int n_o_genotypes, int no_copies)
{
	if (!mem)
	{
		cerr <<"haploid_population::init_genotypes_diverse(): Allocate and set up first!\n";
		return HP_MEMERR;
	}
	if (n_o_genotypes<0 or n_o_genotypes>number_of_individuals)
	{
		cerr <<"haploid_population::init_genotypes_diverse(): number of genotypes has to be positive and smaller than max! got: "<<n_o_genotypes<<"!\n";
		return HP_BADARG;
	}
	pop_size=0;
	if (HP_VERBOSE) cerr <<"haploid_population::init_genotypes_diverse(int n_o_genotypes, int no_copies) with population size "<<n_o_genotypes<<"...";
	generation=0;
	int i,j,w;
	for (i=0; i<n_o_genotypes;)
	{
		for(j=0;j<no_copies; j++);
		{
			for (w=0; w<number_of_words; w++)
			{		//random integer in 0:(1<<number_of_bits[k])-1 for each word
				genotypes[index(i,w)]=gsl_rng_uniform_int(evo_generator,(1<<number_of_bits[w]));
			}
			i++; pop_size++;
		}
	}
	fitness_values->number_of_values=pop_size;
	outcrossing_rates->number_of_values=pop_size;
	gt_labels->number_of_values=pop_size;
	for (i=pop_size; i<number_of_individuals; i++)
	{
		for (w=0; w<number_of_words; w++) genotypes[index(i,w)]=NO_GENOTYPE;
	}
	//calculate its fitness and recombination rates
	calc_fit();
	calc_gt_labels();
	calc_rec();
	calc_stat();	//make sure everything is calculated
	if (HP_VERBOSE) cerr <<"done."<<endl;
	return 0;
}


//init with genotypes with 0
int haploid_population::init_genotypes(int n_o_genotypes)
{
	if (!mem)
	{
		cerr <<"haploid_population::init_genotypes(): Allocate and set up first!\n";
		return HP_MEMERR;
	}
	if (n_o_genotypes<0 or n_o_genotypes>number_of_individuals)
	{
		cerr <<"haploid_population::init_genotypes(): number of genotypes has to be positive and smaller than max! got: "<<n_o_genotypes<<"!\n";
		return HP_BADARG;
	}
	pop_size=0;
	if (HP_VERBOSE) cerr <<"haploid_population::init_genotypes(int n_o_genotypes, int no_copies) with population size "<<n_o_genotypes<<"...";
	generation=0;
	int i,w;
	for (i=0; i<n_o_genotypes;i++)
	{
		for (w=0; w<number_of_words; w++)
		{		//random integer in 0:(1<<number_of_bits[k])-1 for each word
			genotypes[index(i,w)]=0;
		}
		pop_size++;
	}
	fitness_values->number_of_values=pop_size;
	gt_labels->number_of_values=pop_size;
	outcrossing_rates->number_of_values=pop_size;
	for (i=pop_size; i<number_of_individuals; i++)
	{
		for (w=0; w<number_of_words; w++) genotypes[index(i,w)]=NO_GENOTYPE;
	}
	//calculate its fitness and recombination rates
	calc_fit();
	calc_rec();
	calc_gt_labels();
	calc_stat();	//make sure everything is calculated
	if (HP_VERBOSE) cerr <<"done."<<endl;
	return 0;
}

/*
 * Calculate allele frequencies and mean and variance in fitness
 */
void haploid_population::calc_stat()
{
	if (HP_VERBOSE) cerr<<"haploid_population::calc_stat()...";
	fitness_values->calc_variance();
	if (evolve_outcrossing_rates) outcrossing_rates->calc_variance();
	int i,w,k,locus;
	for (locus=0; locus<number_of_loci; locus++) allele_frequencies[locus]=0;
	for (i=0; i<pop_size; i++)
	{
		locus=0; //use locus as a cumulative index across words
		for (w=0; w<number_of_words; w++)
		{
			for (k=0; k<number_of_bits[w]; k++)
			{
				if (genotypes[index(i,w)]&(1<<k)) allele_frequencies[locus]++;
				locus++;
			}
		}
	}
	//convert counts into frequencies
	for (locus=0; locus<number_of_loci; locus++)
	{
		allele_frequencies[locus]/=pop_size;
	}
	if (HP_VERBOSE) cerr<<"done.\n";
}

double haploid_population::get_multi_point_frequency(vector <int> loci)
{
	if (HP_VERBOSE) cerr<<"haploid_population::calc_multi_point_frequency()...";
	int *gt_mask=new int [number_of_words];
	int w,k,i;
	unsigned int  locus;
	double frequency=0;
	for (w=0; w<number_of_words; w++) gt_mask[w]=0;
	for (locus=0; locus<loci.size(); locus++)
	{
		w=locus_word(loci[locus]);
		k=locus_bit(loci[locus]);
		gt_mask[w]+=(1<<k);
		//cerr <<w<<" "<<k<<" "<<gt_mask[w]<<" ";
	}
	int temp;
	for (i=0; i<pop_size; i++)
	{
		temp=0;
		for (w=0; w<number_of_words; w++)
		{
			if ((genotypes[index(i,w)]&gt_mask[w])==gt_mask[w]) temp++;
		}
		if (temp==number_of_words) {
			/*for (w=0; w<number_of_words; w++)
			{
				cerr <<"word: "<<w<<" ";
				for (k=0; k<number_of_bits[w]; k++) {
					if ((genotypes[index(i,w)]&(1<<k))>0) cerr <<'1';
					else cerr <<'0';
				}
				cerr <<"  ";
				for (k=0; k<number_of_bits[w]; k++) {
					if ((gt_mask[w]&(1<<k))>0) cerr <<'1';
					else cerr <<'0';
				}
				cerr <<"  ";
				for (k=0; k<number_of_bits[w]; k++) {
					if (((genotypes[index(i,w)]&gt_mask[w])&(1<<k))>0) cerr <<'1';
					else cerr <<'0';
				}
				cerr <<endl;
			}*/
			frequency+=1.0;
		}
	}
	if (HP_VERBOSE) cerr<<"done.\n";
	delete [] gt_mask;
	//cerr <<" "<< frequency/pop_size<<endl;
	return frequency/pop_size;
}


int haploid_population::evolve(){
	int err=0;
	err+=select_gametes();
	err+=new_generation();
	return err;
}

/* TODO: disentangle crossovers and outcrossing processes
 *
 */
int haploid_population::select_gametes()
{
	double cpot=chemical_potential();
	//draw gametes according to parental fitness
	int total_offsprings=0;
	int os, o, i;
	n_asex_gametes=0;
	n_sex_gametes=0;
	for (i=0; i<pop_size; i++)
	{
		//poisson distributed random numbers -- mean exp(f)/bar{exp(f)}) (since death rate is one, the growth rate is (f-bar{f})
		os=gsl_ran_poisson(evo_generator, exp(fitness_values->values[i]-cpot));
		if (total_offsprings+os<number_of_individuals)
		{
			//gametes is a list of parents, each parent shows up os times
			for (o=total_offsprings; o<total_offsprings+os; o++)
			{
				if (recombination_model==FREE_RECOMBINATION)
				{
					if (gsl_rng_uniform(evo_generator)<outcrossing_rates->values[i]) {sex_gametes[n_sex_gametes]=i; n_sex_gametes++;}
					else {asex_gametes[n_asex_gametes]=i; n_asex_gametes++;}
				}
				else
				{
					sex_gametes[n_sex_gametes]=i; n_sex_gametes++;
				}
			}
			total_offsprings+=os;
		}
		else
		{
			cerr <<"haploid_population::evolve: maximal number of offspring reached: "<<total_offsprings<<" N: "<< number_of_individuals<<"\n";
			cerr <<os<<"  "<<endl;
			return HP_MEMERR;
		}
	}
	return 0;
}

int haploid_population::new_generation()
{
	//construct new generation
	int new_pop_size=0, parent1, parent2, n_o_c=1, i,w, stddev;
	//asexual offspring
	for (i = 0; i < n_asex_gametes; ++i) {
		parent1=asex_gametes[i];
		for (w = 0; w < number_of_words; ++w) {
			new_genotypes[index(i,w)]=genotypes[index(parent1,w)];
		}
		new_fitness_values->values[i]=fitness_values->values[parent1];
		new_outcrossing_rates->values[i]=outcrossing_rates->values[parent1];
		new_gt_labels->values[i]= gt_labels->values[parent1];
	}
	new_pop_size+=n_asex_gametes;
	stddev=sqrt(fitness_values->variance);
	if (n_sex_gametes>0)
	{
		//sexual offspring -- shuffle the set of gametes to ensure random mating
		gsl_ran_shuffle(evo_generator, sex_gametes, n_sex_gametes, sizeof(int));
		//if the number of sex_gametes is odd, the last one is unlucky
		for (i = 0; i < n_sex_gametes-1; i+=2) {
			parent1=sex_gametes[i];
			parent2=sex_gametes[i+1];
			//The recombination function store two new genotypes in new_genotypes[new_pop_size] and [new_pop_size+1]
			if (recombination_model==FREE_RECOMBINATION) { //no linkage
				recombine_free(parent1, parent2, new_pop_size);
				n_o_c=1;
			} else {		//linear chromosomes. returns number of crossovers.
			    n_o_c=recombine_crossover(parent1, parent2, new_pop_size);
			}
			//if number of crossovers is greater zero, things have to be recalculated
			if (n_o_c>0)
			{
				new_fitness_values->values[new_pop_size]=fitness.get_func_words(&new_genotypes[index(new_pop_size,0)], number_of_words);
				new_fitness_values->values[new_pop_size+1]=fitness.get_func_words(&new_genotypes[index(new_pop_size+1,0)], number_of_words);
				new_gt_labels->values[new_pop_size]=calc_label(&new_genotypes[index(new_pop_size,0)]);
				new_gt_labels->values[new_pop_size+1]=calc_label(&new_genotypes[index(new_pop_size+1,0)]);
				if (evolve_outcrossing_rates)
				{
					new_outcrossing_rates->values[new_pop_size]=outcrossing_rate.get_func_words(&new_genotypes[index(new_pop_size,0)], number_of_words);
					new_outcrossing_rates->values[new_pop_size+1]=outcrossing_rate.get_func_words(&new_genotypes[index(new_pop_size+1,0)], number_of_words);
				}
				else
				{
					new_outcrossing_rates->values[new_pop_size]=fixed_outcrossing_rate;
					new_outcrossing_rates->values[new_pop_size+1]=fixed_outcrossing_rate;
				}
			}else { //without crossovers, values can be copied
				new_fitness_values->values[new_pop_size]=fitness_values->values[parent1];
				new_outcrossing_rates->values[new_pop_size]=outcrossing_rates->values[parent1];
				new_gt_labels->values[new_pop_size]= gt_labels->values[parent1];
				new_fitness_values->values[new_pop_size+1]=fitness_values->values[parent2];
				new_outcrossing_rates->values[new_pop_size+1]=outcrossing_rates->values[parent2];
				new_gt_labels->values[new_pop_size+1]= gt_labels->values[parent2];
			}
			new_pop_size+=2;
		}
	}

	for (i = new_pop_size; i < pop_size; ++i) {
		for (w = 0; w < number_of_words; ++w) {
			new_genotypes[index(i,w)]=NO_GENOTYPE;
		}
	}
	//make the new generation the old one
	int *temp_gt;
	temp_gt=genotypes;
	genotypes=new_genotypes;
	new_genotypes=temp_gt;

	sample *temp_values;
	temp_values=fitness_values;
	fitness_values=new_fitness_values;
	new_fitness_values=temp_values;

	temp_values=outcrossing_rates;
	outcrossing_rates=new_outcrossing_rates;
	new_outcrossing_rates=temp_values;

	temp_values=gt_labels;
	gt_labels=new_gt_labels;
	new_gt_labels=temp_values;


	pop_size=new_pop_size;
	fitness_values->number_of_values=pop_size;
	outcrossing_rates->number_of_values=pop_size;
	gt_labels->number_of_values=pop_size;

	generation++;
	return 0;
}

//recombine two genotypes parent1 and parent2 to two new stored in new_genotypes
//at positions ng and ng+1
int haploid_population::recombine_free(int parent1, int parent2, int ng)
{
	int w, temp;
	for (w=0; w<number_of_words; w++)
	{
		temp=gsl_rng_uniform_int(evo_generator, 1<<number_of_bits[w]);
		new_genotypes[index(ng,w)]=(genotypes[index(parent1,w)])&(temp);		//pull out the mothers bits
		new_genotypes[index(ng,w)]+=(genotypes[index(parent2,w)])&(~temp);		//pull out the fathers bits (by negation of temp)
		new_genotypes[index(ng+1,w)]=(genotypes[index(parent1,w)])&(~temp);		//pull out the mothers bits the complement of the mother
		new_genotypes[index(ng+1,w)]+=(genotypes[index(parent2,w)])&(temp);		//complement of the father
	}
	return 1;
}


//recombine two genotypes gam1 and gam2 to a new one stored in new genotypes at position ng
int haploid_population::recombine_crossover(int parent1, int parent2, int ng)
{
	int w, temp;
	int strand=1, c=0, old_co;
	double crossover_rate=0.5*(outcrossing_rates->values[parent1]+outcrossing_rates->values[parent2]);
	//double crossover_rate=0.5*(outcrossing_rates->values[parent1]+outcrossing_rates->values[parent2]);
	int n_o_c=gsl_ran_poisson(evo_generator,number_of_loci*crossover_rate);
	if (circular)
	{
		n_o_c*=2;
		n_o_c=(n_o_c<number_of_loci)?n_o_c:number_of_loci;
		gsl_ran_choose(evo_generator,crossovers,n_o_c,genome,number_of_loci,sizeof(int));
	}
	else
	{
		n_o_c=(n_o_c<number_of_loci)?n_o_c:(number_of_loci-1);
		gsl_ran_choose(evo_generator,crossovers,n_o_c,(genome+1),number_of_loci-1,sizeof(int));
	}
	for (w=0; w<number_of_words; w++)
	{
		new_genotypes[index(ng,w)]=0;
		new_genotypes[index(ng+1,w)]=0;
		old_co=0;
		while((crossovers[c]< (w+1)*WORD_LENGTH)&&(c<n_o_c))	//loop over cross overs within chromosome
		{
			temp=(1<<(crossovers[c]-w*WORD_LENGTH))-(1<<old_co);
			if (strand){
				new_genotypes[index(ng,w)]+=(genotypes[index(parent1,w)])&(temp);		//pull out the mothers bits
				new_genotypes[index(ng+1,w)]+=(genotypes[index(parent2,w)])&(temp);		//pull out the mothers bits
			} else {
				new_genotypes[index(ng,w)]+=(genotypes[index(parent2,w)])&(temp);		//pull out the fathers bits
				new_genotypes[index(ng+1,w)]+=(genotypes[index(parent1,w)])&(temp);		//pull out the fathers bits
			}
			old_co=crossovers[c]-w*WORD_LENGTH;
			c++;
			strand=(strand+1)%2;
		}
		temp=(1<<(WORD_LENGTH))-(1<<old_co);			//assign the remainder of the chromosome to new genotype (complete chromosome if no crossover)
		if (strand)	{
			new_genotypes[index(ng,w)]+=(genotypes[index(parent1,w)])&(temp);		//pull out the mothers bits
			new_genotypes[index(ng+1,w)]+=(genotypes[index(parent2,w)])&(temp);		//pull out the mothers bits
		} else {
			new_genotypes[index(ng,w)]+=(genotypes[index(parent2,w)])&(temp);				//pull out the fathers bits
			new_genotypes[index(ng+1,w)]+=(genotypes[index(parent1,w)])&(temp);				//pull out the fathers bits
		}
	}
	return n_o_c;
}


//mutate all loci
void haploid_population::mutate()
{
	if (HP_VERBOSE)
		cerr <<"haploid_population::mutate() .....";
	int i, actual_n_o_mutations,locus=0;
	double mean_n_o_mutations=mutation_rate*pop_size;
	int *mutations=new int [int(ceil(mean_n_o_mutations))*5+100];
	for (locus = 0; locus<number_of_words; locus++) {
		actual_n_o_mutations=gsl_ran_poisson(evo_generator, mean_n_o_mutations);
		gsl_ran_choose(evo_generator, mutations, actual_n_o_mutations, numbers,pop_size, sizeof(int));
		for (i = 0; i < actual_n_o_mutations; ++i) {
			flip_single_locus(mutations[i], locus);
		}
	}
	delete mutations;
	if (HP_VERBOSE)
		cerr <<"done";
}

//flip a spin at locus in random individual.
int haploid_population::flip_single_locus(int locus)
{
	int individual = gsl_rng_uniform_int(evo_generator, pop_size);
	flip_single_locus(individual, locus);
	return individual;
}


//flip a spin at locus w,k (word and bit) in individual, assign new fitness and recombination rate.
void haploid_population::flip_single_locus(int individual, int locus)
{
	int k,w;
	k=locus_bit(locus);
	w=locus_word(locus);
	//check the state of the locus and add or subtract the relevant power of two
	if (genotypes[index(individual,w)] & (1<<k)) {genotypes[index(individual,w)]-=(1<<k);}
	else genotypes[index(individual,w)]+=(1<<k);

	calc_individual_fitness(individual);

	if (HP_VERBOSE) cerr <<"subpop::flip_single_spin(): mutated individual "<<individual<<" in word "<<w<<" bit "<<k<<endl;
}

void haploid_population::calc_individual_fitness(int individual)
{
	//calculate the new fitness value of the mutant
	fitness_values->values[individual]=fitness.get_func_words(&genotypes[index(individual,0)], number_of_words);
	gt_labels->values[individual]=calc_label(&genotypes[index(individual,0)]);
	//calculate the recombination rate
	if (evolve_outcrossing_rates) outcrossing_rates->values[individual]=outcrossing_rate.get_func_words(&genotypes[index(individual,0)], number_of_words);
	else outcrossing_rates->values[individual]=fixed_outcrossing_rate;
}

int haploid_population::bottleneck(int size_of_bottleneck)
{
	shuffle_genotypes();
	return remove_last_genotypes(pop_size-size_of_bottleneck);
}

int haploid_population::remove_last_genotypes(int m)
{
	int i,w;
	if (m<0 or m>pop_size)
	{
		cerr <<"haploid_population::remove_last_genotypes(): Bad argument: "<<m<<"\n";
		return HP_BADARG;
	}
	for (i = pop_size-m; i < pop_size; ++i) {
		for (w = 0; w < number_of_words; ++w) {
			genotypes[index(i,w)]=NO_GENOTYPE;
		}
	}
	pop_size-=m;
	fitness_values->number_of_values-=m;
	outcrossing_rates->number_of_values-=m;

	return 0;
}

string haploid_population::get_genotype_string(int i){
	stringstream gt;
	gt.str("");
	for (int w = 0; w < number_of_words; w++) {
		for (int k=0; k<number_of_bits[w]; k++){
			gt <<((genotypes[index(i,w)]&(1<<k))>0);
		}
	}
	return gt.str();
}

int haploid_population::add_genotypes(int *gt, int n)
{
	int err=0;
	for (int i = 0; i < n; ++i) {
		err=add_genotype(gt);
	}
	return err;
}

int haploid_population::add_genotype(int *gt)
{
	if (pop_size>=number_of_individuals)
	{
		cerr<<"haploid_population::add_genotype(): exhausted memory!\n";
		return HP_MEMERR;
	}
	for (int w = 0; w < number_of_words; ++w) {
		genotypes[index(pop_size,w)]=gt[w];
	}
	fitness_values->values[pop_size]=fitness.get_func_words(&genotypes[index(pop_size,0)], number_of_words);
	if (evolve_outcrossing_rates)
	{
		outcrossing_rates->values[pop_size]=outcrossing_rate.get_func_words(&genotypes[index(pop_size,0)], number_of_words);
	}
	else
	{
		outcrossing_rates->values[pop_size]=fixed_outcrossing_rate;
	}

	fitness_values->number_of_values++;
	outcrossing_rates->number_of_values++;
	pop_size++;
	return 0;
}

double haploid_population::chemical_potential()
{
	fitness_values->calc_mean();
	return fitness_values->mean+0.6931*(double(pop_size)/target_pop_size-1);
}

void haploid_population::calc_fit()
{
	if (HP_VERBOSE)
		cerr <<"haploid_population::calc_fit() .....";
	for (int i=0; i<pop_size; i++)
	{
		fitness_values->values[i]=fitness.get_func_words(&genotypes[index(i,0)], number_of_words);
		new_fitness_values->values[i]=fitness_values->values[i];
	}
	fitness_values->calc_variance();
	new_fitness_values->calc_variance();
	if (HP_VERBOSE)
		cerr <<"done"<<endl;
}

double haploid_population::calc_label(int* i){
	int w=0;
	double temp_label=0;
	int temp_seed=*i;
	gsl_rng_set(label_generator, temp_seed);
	for (w = 0; w < number_of_words-1; w++) {
		temp_label+=gsl_rng_uniform(label_generator);
		temp_seed+=*(i+w+1)+gsl_rng_uniform_int(label_generator, 1<<20);
		gsl_rng_set(label_generator,temp_seed);
	}
	temp_label+=gsl_rng_uniform(label_generator);
	return temp_label;
}

void haploid_population::calc_gt_labels()
{
	if (HP_VERBOSE)
		cerr <<"haploid_population::calc_gt_labels() .....";
	for (int i=0; i<pop_size; i++)
	{
		gt_labels->values[i]=calc_label(&genotypes[index(i,0)]);
		new_gt_labels->values[i]=gt_labels->values[i];
	}
	if (HP_VERBOSE)
		cerr <<"done"<<endl;
}

void haploid_population::calc_rec()
{
	if (HP_VERBOSE)
		cerr <<"haploid_population::calc_rec() .....";
	//calculate its fitness and recombination rates
	for (int i=0; i<pop_size; i++){
		if (evolve_outcrossing_rates) outcrossing_rates->values[i]=outcrossing_rate.get_func_words(&genotypes[index(i,0)], number_of_words);
		else outcrossing_rates->values[i]=fixed_outcrossing_rate;
		new_outcrossing_rates->values[i]=outcrossing_rates->values[i];
	}
	outcrossing_rates->calc_variance();
	new_outcrossing_rates->calc_variance();
	if (HP_VERBOSE)
		cerr <<"done"<<endl;
}

double haploid_population::get_max_fitness()
{
	double mf=fitness_values->values[0];
	for (int i = 0; i < pop_size; ++i) {
		if (mf<fitness_values->values[i]) mf=fitness_values->values[i];
	}
	return mf;
}


//function that shuffles the entire genotype distribution, this is useful for example before
//or after a migration step.
void haploid_population::shuffle_genotypes()
{
	int i, w;
	for (i = 0; i < pop_size; ++i) {
		asex_gametes[i]=i;
	}
	gsl_ran_shuffle(evo_generator, asex_gametes, pop_size, (sizeof(int)));
	for (i = 0; i < pop_size; ++i) {
		for (w = 0; w < number_of_words; ++w) {
			new_genotypes[index(i,w)]=genotypes[index(asex_gametes[i],w)];
		}
		new_fitness_values->values[i]=fitness_values->values[asex_gametes[i]];
		outcrossing_rates->values[i]=outcrossing_rates->values[asex_gametes[i]];
	}
	//make the new genotypes the old one
	int *temp_gt;
	temp_gt=genotypes;
	genotypes=new_genotypes;
	new_genotypes=temp_gt;
	sample *temp_values;
	temp_values=fitness_values;
	fitness_values=new_fitness_values;
	new_fitness_values=temp_values;
	temp_values=outcrossing_rates;
	outcrossing_rates=new_outcrossing_rates;
	new_outcrossing_rates=temp_values;
	fitness_values->number_of_values=pop_size;
	outcrossing_rates->number_of_values=pop_size;
}

//print allele frequencies to stream
int haploid_population::print_allele_frequencies(ostream &out)
{
	if (out.bad())
	{
		cerr <<"haploid_population::print_allele_frequencies: bad stream\n";
		return HP_BADARG;
	}
	calc_stat();
	out <<setw(10)<<generation;
	for (int l=0; l<number_of_loci; l++) out <<setw(15)<<allele_frequencies[l];
	out <<endl;
	return 0;
}

/*
 * Use genotype labels to calculate the clone size distribution.
 * -> sort labels Nlog(N), walk through labels and count how many consecutive labels are the same
 * each of these group corresponds to one clone size
 */
int haploid_population::calc_clone_size_distribution()
{
	if (clone_size_mem==false)
	{
		cerr <<"haploid_population::calc_clone_size_distribution(): allocate memory for genotype distribution first!\n";
		return HP_MEMERR;
	}
	if (HP_VERBOSE) cerr <<"haploid_population::calc_clone_size_distributions(): ...";
	label_gt_pair *sorted_labels=new label_gt_pair [pop_size];
	for (int i = 0; i < pop_size; ++i) {
		sorted_labels[i].gt=i;
		sorted_labels[i].label=gt_labels->values[i];

	}
	if (HP_VERBOSE) cerr <<"sorting...";
	qsort(sorted_labels, pop_size, sizeof(label_gt_pair), comp_labels_hp);
	if (HP_VERBOSE) cerr <<"done...";
	double label=sorted_labels[0].label;
	int clone_size=1, gt=sorted_labels[0].gt;
	number_of_clones=0;
	int temp_popsize=0;
	if (HP_VERBOSE) cerr <<"counting...";
	for (int i = 1; i < pop_size; i++) {
		if (label+HP_NOTHING<sorted_labels[i].label)
		{	//new clone
			clone_size_distribution[2*number_of_clones]=clone_size;
			clone_size_distribution[2*number_of_clones+1]=gt;
			number_of_clones++;
			temp_popsize+=clone_size;
			clone_size=1;
			gt=sorted_labels[i].gt;
			label=sorted_labels[i].label;
		}else{
			clone_size++;
		}
	}
	clone_size_distribution[2*number_of_clones]=clone_size;
	clone_size_distribution[2*number_of_clones+1]=gt;
	number_of_clones++;
	temp_popsize+=clone_size;
	if (temp_popsize!=pop_size) cerr <<"haploid_population::calc_clone_size_distribution(): Lost genotypes!"<<endl;
	if (HP_VERBOSE) cerr <<"done.\n";
	delete [] sorted_labels;
	return 0;
}


//function that reads the output from ms to initialize the population
int haploid_population::read_ms_sample(istream &gts, int initial_locus, int multiplicity){
	if (gts.bad()){
		cerr<<"haploid_population::read_ms_sample(): bad stream!\n";
		return HP_BADARG;
	}
	char *line= new char [2*number_of_loci+5000];
	bool gt=false;
	string header;
	int count=0;
	int segsites, site, locus,k,w;
	segsites=0;
	int *newgt=new int [number_of_words];
	if (mem){
		init_genotypes(0);
		pop_size=0;
		while(gts.eof()==false){
			gts.get(line, 2*number_of_loci+5000);
			gts.get();
			while (gts.peek()=='\n'){gts.get();};

			//cout <<count<<"  "<<gt<<" "<<line<<endl;
			count++;
			//if (count>10) break;
			if (gt and line[0]!='\0'){
				for (w=0; w<number_of_words; w++) newgt[w]=0;
				w=0;
				for (site=0; site<segsites; site++){
					if (line[site]=='1'){
						locus=initial_locus+site*(number_of_loci-initial_locus)/segsites;
						k=locus_bit(locus);
						w=locus_word(locus);
						newgt[w]+=(1<<k);
					}
				}
				add_genotypes(newgt, multiplicity);
			}else{
				header.assign(line);
				if (header.compare(0,2, "//")==0) {
					gts.get(line, 2*number_of_loci);
					gts.get();
					while (gts.peek()=='\n'){gts.get();};
					cerr <<count<<"  "<<gt<<" "<<line<<endl;
					header.assign(line);
					segsites=atoi(header.substr(9,header.size()-9).c_str());
					gts.get(line, 2*number_of_loci);
					gts.get();
					while (gts.peek()=='\n'){gts.get();};
					cerr <<count<<"  "<<gt<<" "<<line<<endl;
					gt=true;
				}
			}
		}
	}
	delete [] line;
	return 0;
}



