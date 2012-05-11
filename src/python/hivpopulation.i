/* renames and ignores */
%ignore set_up(int N_in, int L,  int rng_seed=0, int number_of_traits=1); /* This seems to work for mysterious reasons */


/**** HIVPOPULATION ****/
%extend hivpopulation {
/* evolve */
%ignore evolve;
%ignore _evolve;
%rename (evolve) _evolve;
int _evolve(int gen) {
        int err=$self->evolve(gen);
        if(err==0)
                $self->calc_stat();
        return err;
}

/* get allele frequencies */
void get_allele_frequencies(double ARGOUT_ARRAY1[HIVGENOME]) {
        for(size_t i=0; i < HIVGENOME; i++)
                ARGOUT_ARRAY1[i] = $self->get_allele_frequency(i);
}

/* get fitnesses of all clones */
void _get_fitnesses(int DIM1, double* ARGOUT_ARRAY1) {
        for(size_t i=0; i < DIM1; i++)
                ARGOUT_ARRAY1[i] = $self->get_fitness(i);
}
%pythoncode {
    def get_fitnesses(self):
        return self._get_fitnesses(self.get_number_of_clones())
}

/* get genotypes */
void _get_genotype(unsigned int i, short ARGOUT_ARRAY1[HIVGENOME]) {
        boost::dynamic_bitset<> newgt = (*($self->current_pop))[i].genotype;
        for(size_t i=0; i < HIVGENOME; i++)
                ARGOUT_ARRAY1[i] = newgt[i];
}
%pythoncode {
    def get_genotypes(self, ind=None):
        import numpy as np
        if ind is None:
            ind = np.arange(self.get_number_of_clones())
        else:
            ind = np.array(ind, ndmin=1)
        n = len(ind)
        genotypes = np.zeros((n, self.get_number_of_loci()), bool)
        for i, indi in enumerate(ind):
            genotypes[i] = self._get_genotype(indi)
        if n == 1:
            return genotypes[0]
        else:
            return genotypes
}

/* get random clones/genotypes */
%pythoncode {
    def random_genomes(self, n):
        import numpy as np
        genotypes = np.zeros((n, self.get_number_of_loci()), bool)
        for i in xrange(genotypes.shape[0]):
            genotypes[i] = self._get_genotype(self.random_clone())
        return genotypes
}
void random_clones(int DIM1, unsigned int * ARGOUT_ARRAY1) {
        vector <int> sample = vector <int>(0);
        int err = $self->random_clones(DIM1, &sample);
        if(!err)
                for(size_t i=0; i < DIM1; i++)
                        ARGOUT_ARRAY1[i] = sample[i];
}

/* divergence/diversity/fitness distributions and plot */
%pythoncode {
    def get_fitness_histogram(self, bins=10, n_sample=1000, **kwargs):
        '''Calculate the fitness histogram.'''
        import numpy as np
        fit = [self.get_fitness(self.random_clone()) for i in xrange(n_sample)]
        h = np.histogram(fit, bins=bins, **kwargs)
        return h
    
    
    def plot_fitness_histogram(self, axis=None, **kwargs):
        '''Plot a distribution of fitness in the population.'''
        import matplotlib.pyplot as plt
        fit = self.get_fitnesses();
    
        if axis is None:
            fig = plt.figure()
            axis = fig.add_subplot(111)
            axis.set_title('Fitness histogram')
            axis.set_xlabel('Fitness')
    
        axis.hist(fit, **kwargs)
    
    
    def get_divergence_histogram(self, bins=10, chunks=None, every=1, n_sample=1000, **kwargs):
        '''Get the divergence histogram restricted to those chunks of the genome.'''
        import numpy as np
    
        # Check chunks
        if chunks is not None:
            chunks = np.asarray(chunks)
            if (np.rank(chunks) != 2) or (chunks.shape[1] != 2):
                raise ValueError('Please input an N x 2 matrix with the chunks initial and (final+1) positions')
    
        # Get the random genotypes
        genotypes = self.random_genomes(n_sample)
    
        # Restrict to the chunks
        if chunks is not None:
            ind = np.zeros(genotypes.shape[1], bool)
            for chunk in chunks:
                inde = np.arange(chunk[1] - chunk[0])
                inde = inde[(inde % every) == 0] + chunk[0]
                ind[inde] = True
            genotypes = genotypes[:,ind]
    
        # Calculate divergence
        div = genotypes.sum(axis=1)
    
        # Calculate histogram
        h = np.histogram(div, bins=bins, **kwargs)
        return h 
    
    
    def plot_divergence_histogram(self, axis=None, **kwargs):
        import matplotlib.pyplot as plt
        import numpy as np
        n_sample = 1000
        genotypes = self.random_genomes(n_sample)
        div = genotypes.sum(axis=1)
     
        if axis is None:
            fig = plt.figure()
            axis = fig.add_subplot(111)
            axis.set_title('Divergence histogram')
            axis.set_xlabel('Divergence')
    
        if 'bins' not in kwargs:
            kwargs['bins'] = np.arange(10) * max(1, (div.max() + 1 - div.min()) / 10) + div.min()
    
        axis.hist(div, **kwargs)
    
    
    def get_diversity_histogram(self, bins=10, chunks=None, every=1, n_sample=1000, **kwargs):
        '''Get the diversity histogram restricted to those chunks of the genome.'''
        import numpy as np
    
        # Check chunks
        if chunks is not None:
            chunks = np.asarray(chunks)
            if (np.rank(chunks) != 2) or (chunks.shape[1] != 2):
                raise ValueError('Please input an N x 2 matrix with the chunks initial and (final+1) positions')
    
        # Get the random genotypes
        genotypes = self.random_genomes(2 * n_sample)
    
        # Restrict to the chunks
        if chunks is not None:
            ind = np.zeros(genotypes.shape[1], bool)
            for chunk in chunks:
                inde = np.arange(chunk[1] - chunk[0])
                inde = inde[(inde % every) == 0] + chunk[0]
                ind[inde] = True
            genotypes = genotypes[:,ind]
    
        # Calculate diversity
        genotypes1 = genotypes[:genotypes.shape[0] / 2]
        genotypes2 = genotypes[-genotypes1.shape[0]:]
        div = (genotypes1 != genotypes2).sum(axis=1)
    
        # Calculate histogram
        h = np.histogram(div, bins=bins, **kwargs)
        return h 


    def plot_diversity_histogram(self, axis=None, **kwargs):
        import matplotlib.pyplot as plt
        import numpy as np
        n_sample = 1000
        genotypes1 = self.random_genomes(n_sample)
        genotypes2 = self.random_genomes(n_sample)
        div = (genotypes1 != genotypes2).sum(axis=1)
    
        if axis is None:
            fig = plt.figure()
            axis = fig.add_subplot(111)
            axis.set_title('Diversity histogram')
            axis.set_xlabel('Diversity')
    
        if 'bins' not in kwargs:
            kwargs['bins'] = np.arange(10) * max(1, (div.max() + 1 - div.min()) / 10) + div.min()
    
        axis.hist(div, **kwargs)
}

/* read selection/resistance coefficients */
%ignore read_selection_coefficients;
%rename (read_selection_coefficients) _read_selection_coefficients;
int _read_selection_coefficients(char *model){
        ifstream modelstream(model);
        return $self->read_selection_coefficients(modelstream);
}

%ignore read_resistance_coefficients;
%rename (read_resistance_coefficients) _read_resistance_coefficients;
int _read_resistance_coefficients(char *model){
        ifstream modelstream(model);
        return $self->read_resistance_coefficients(modelstream);
}


/* write genotypes */
%ignore write_genotypes;
%rename (write_genotypes) _write_genotypes;
int _write_genotypes(char * filename, int sample_size, char * gt_label=NULL, int start=0, int length=0) {
        if(HIVPOP_VERBOSE >= 1) cerr<<"hivpython::write_genotypes(char * filename, int sample_size, char * gt_label, int start, int length)...";

        ofstream genotype_file(filename);        
        string gt_labelstring;
        if(gt_label == NULL)
                gt_labelstring = "";
        else
                gt_labelstring = string(gt_label);

        return $self->write_genotypes(genotype_file, sample_size, gt_labelstring, start, length);
}


%pythoncode {
    def write_genotypes_compressed(self, filename, sample_size, gt_label='', start=0, length=0):
        '''Write genotypes into a compressed archive.'''
        import numpy as np 
        if length <= 0:
            length = self.get_number_of_loci() - start
        d = {}
        for i in xrange(sample_size):
            rcl = self.random_clone()
            d['>'+str(i)+'_GT-'+gt_label+'_'+str(rcl)] = self._get_genotype(rcl)[start:start+length]
        np.savez_compressed(filename, **d)    
}

/* fitness landscape */
void calc_fitness_from_traits(int l) {
        $self->calc_fitness_from_traits(&((*($self->current_pop))[l]));
}

} /* extend hivpopulation */


