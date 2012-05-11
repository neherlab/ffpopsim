/* renames and ignores */
%ignore coeff_t;
%ignore coeff_single_locus_t;
%ignore hypercube_function;

/**** CLONE_T ****/
%extend clone_t {
%rename (_trait) trait;
%rename (_genotype) genotype;

/* number of traits */
int get_number_of_traits() {
        return ($self->trait).size();
}

/* get single trait */
void _get_trait(int DIM1, double* ARGOUT_ARRAY1) {
        for(size_t i=0; i<($self->trait).size(); i++)
                ARGOUT_ARRAY1[i] = ($self->trait)[i];
}
%pythoncode {
    @property
    def trait(self):
        return self._get_trait(self.get_number_of_traits())
}

/* genotype */
void _get_genotype(int DIM1, short* ARGOUT_ARRAY1) {
        for(size_t i=0; i < ($self->genotype).size(); i++) ARGOUT_ARRAY1[i] = ($self->genotype)[i];
}

int _get_genotype_length() {return ($self->genotype).size();}
%pythoncode {
    @property
    def genotype(self):
        import numpy as np
        return np.array(self._get_genotype(self._get_genotype_length()), bool)
}
}

/**** HAPLOID_CLONE ****/
%apply (double* IN_ARRAY1, int DIM1) {(double *nu, int n_o_genotypes=0)};
%extend haploid_clone {
/* Hamming distance (full Python reimplementation) */
%ignore distance_Hamming;
%pythoncode {
    def distance_Hamming(self, clone_gt1, clone_gt2, chunks=None, every=1):
        import numpy as np
        if np.isscalar(clone_gt1):
            genotypes = self.get_genotypes((clone_gt1, clone_gt2))
            clone_gt1 = genotypes[0]
            clone_gt2 = genotypes[1]

        if chunks is not None:
            ind = np.zeros(clones.shape[1], bool)
            for chunk in chunks:
                inde = np.arange(chunk[1] - chunk[0])
                inde = inde[(inde % every) == 0] + chunk[0]
                ind[inde] = True
            clone_gt1 = clone_gt1[ind]
            clone_gt2 = clone_gt2[ind]
        return (clone_gt1 != clone_gt2).sum()
}
}
