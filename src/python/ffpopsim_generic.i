/* renames and ignores */
%ignore SAMPLE_ERROR;
%ignore sample;


/**** INDEX_VALUE_PAIR_T ****/
%feature("autodoc", "Pair of an index and a value") index_value_pair_t;

%rename(index_value_pair) index_value_pair_t;
%extend index_value_pair_t {
const char* __str__() {
        static char buffer[255];
        sprintf(buffer,"index: %u, val: %.2e", $self->index, $self->val);
        return &buffer[0];
}

const char* __repr__() {
        static char buffer[255];
        sprintf(buffer,"(%u, %.2e)", $self->index, $self->val);
        return &buffer[0];
}

%feature("autodoc", "Index") index;
%feature("autodoc", "Value") val;

} /* extend index_value_pair_t */


/**** GENOTYPE_VALUE_PAIR_T ****/
%feature("autodoc", "Pair of a genotype and a value") genotype_value_pair_t;

%rename(genotype_value_pair) genotype_value_pair_t;
%extend genotype_value_pair_t {

/* string representations */
const char* __str__() {
        static char buffer[255];
        unsigned long L = ($self->genotype).size();
        if(L > 2)
                sprintf(buffer,"genotype: %d...%d, val: %.2e", (int)($self->genotype)[0], (int)($self->genotype)[L-1], $self->val);
        else if(L == 2)
                sprintf(buffer,"genotype: %d%d, val: %.2e", (int)($self->genotype)[0], (int)($self->genotype)[1], $self->val);
        else if(L == 1)
                sprintf(buffer,"genotype: %d, val: %.2e", (int)($self->genotype)[0], $self->val);
        else
                sprintf(buffer,"genotype: (empty), val: %.2e", $self->val);
        return &buffer[0];
}

const char* __repr__() {
        static char buffer[255];
        unsigned long L = ($self->genotype).size();
        if(L > 2)
                sprintf(buffer,"([%d, ..., %d], %.2e)", (int)($self->genotype)[0], (int)($self->genotype)[L-1], $self->val);
        else if(L == 2)
                sprintf(buffer,"([%d, %d], %.2e)", (int)($self->genotype)[0], (int)($self->genotype)[1], $self->val);
        else if(L == 1)
                sprintf(buffer,"([%d], %.2e)", (int)($self->genotype)[0], $self->val);
        else
                sprintf(buffer,"([], %.2e)", $self->val);
        return &buffer[0];
}

/* constructor */
%typemap(in) boost::dynamic_bitset<> genotype_in (boost::dynamic_bitset<> temp) {
        /* Ensure input is a Python sequence */
        PyObject *tmplist = PySequence_Fast($input, "I expected a sequence");
        unsigned long L = PySequence_Length(tmplist);

        /* Create boost::dynamic_bitset from Python list */
        temp.resize(L);
        long tmplong;
        for(size_t i=0; i < L; i++) {
                tmplong = PyInt_AsLong(PySequence_Fast_GET_ITEM(tmplist, i));
                if(tmplong < 0) {
                        PyErr_SetString(PyExc_ValueError, "Expecting an array of bool.");
                        SWIG_fail;
                }
                temp[i] = (bool)tmplong; 
        }      
        $1 = temp;
}


/* genotype */
%ignore genotype;
int _get_genotype_length() {return ($self->genotype).size();}
void _get_genotype(int DIM1, short* ARGOUT_ARRAY1) {
        for(size_t i=0; i < ($self->genotype).size(); i++) ARGOUT_ARRAY1[i] = ($self->genotype)[i];
}

void _set_genotype(boost::dynamic_bitset<> genotype_in) {$self->genotype = genotype_in;}

%pythoncode {
@property
def genotype(self):
    '''Genotype'''
    import numpy as np
    return np.array(self._get_genotype(self._get_genotype_length()), bool)


@genotype.setter
def genotype(self, genotype):
    self._set_genotype(genotype)
}

%feature("autodoc", "Value") val;

} /* extend genotype_value_pair_t */

/**** STAT_T ****/
%feature("autodoc", "Mean and variance of a statistical distribution") stat_t;

%rename(stat) stat_t;
%extend stat_t {
const char* __str__() {
        static char buffer[255];
        sprintf(buffer,"mean: %e, variance: %e", $self->mean, $self->variance);
        return &buffer[0];
}

const char* __repr__() {
        static char buffer[255];
        sprintf(buffer,"(%e, %e)", $self->mean, $self->variance);
        return &buffer[0];
}

%feature("autodoc", "Mean") mean;
%feature("autodoc", "Variance") variance;

} /* extend stat_t */
