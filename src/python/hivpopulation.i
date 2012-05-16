/**** HIVPOPULATION ****/
%define DOCSTRING_HIVPOPULATION
"Class for HIV population genetics (genome size = 10000).

This class is the main object for simulating the evolution of HIV.
The class offers a number of functions, but an example will explain the basic
idea:

#####################################
#   EXAMPLE SCRIPT                  #
#####################################
import numpy as np
import matplotlib.pyplot as plt
import PopGenLib as h

c = h.hivpopulation(2000)
c.init_genotypes() 
c.evolve(10)
c.plot_divergence_histogram()
plt.show()
#####################################

An effective way to discover all available methods is to import PopGenLib from
an interactive shell (e.g. iPython), create a population as above, and use TAB
autocompletion:

In [1]: import PopGenLib as h
In [2]: c = h.haploid_clone(5000, 2000)
In [3]: c.      <--- TAB

In addition to the haploid_clone class, this class offers functions for reading
fitness and drug resistance landscapes from a text file, and to save genomes as
plain text or in compressed numerical Python format.
"
%enddef
%feature("autodoc", DOCSTRING_HIVPOPULATION) hivpopulation;

%extend hivpopulation {


/* we do not need this and it conflicts with base class constructor calling */
%ignore set_up;

/* constructor */
%exception hivpopulation {
        try {
                $action
        } catch (int err) {
                PyErr_SetString(PyExc_ValueError,"Construction impossible. Please check input args.");
                SWIG_fail;
        }
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

} /* extend hivpopulation */


