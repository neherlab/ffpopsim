/**
 * Copyright (c) 2012-2013, Richard Neher, Fabio Zanini
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

/* ignore some classes */
%ignore coeff_t;
%ignore coeff_single_locus_t;
%ignore hypercube_highd;
%ignore step_t;
%ignore node_t;

/*****************************************************************************/
/* CLONE_T                                                                   */
/*****************************************************************************/
%feature("autodoc", "Clone of isogenetic individuals") clone_t;
%rename (clone) clone_t;
%extend clone_t {
/* string representations */
%feature("autodoc", "x.__str__() <==> str(x)") __str__;
%feature("autodoc", "x.__repr__() <==> repr(x)") __repr__;
const char* __str__() {
        static char buffer[255];
        sprintf(buffer,"clone: %u traits, genome size = %u",
                       (unsigned int)($self->trait).size(),
                       (unsigned int)($self->genotype).size());
        return &buffer[0];
}
const char* __repr__() {
        static char buffer[255];
        sprintf(buffer,"<clone>");
        return &buffer[0];
}

/* read/write attributes */
%feature("autodoc", "Genotype of the clone") genotype;
%feature("autodoc", "Fitness of the clone") fitness;
%feature("autodoc", "Number of individuals of the clone") clone_size;

/* traits */
%feature("autodoc", "Number of traits (read-only)") number_of_traits;
const int number_of_traits;

%ignore trait;
void _get_trait(int DIM1, double* ARGOUT_ARRAY1) {
        for(size_t i=0; i < (size_t)DIM1; i++)
                ARGOUT_ARRAY1[i] = ($self->trait)[i];
}
%pythoncode
%{
@property
def trait(self):
    '''Traits vector of the clone'''
    return self._get_trait(self.number_of_traits)
%}

} /* extend clone_t */
%{
const int clone_t_number_of_traits_get(clone_t *c) {
        return (const int) c->trait.size();
}
%} /* attributes of clone_t */

/*****************************************************************************/
/* POLY_T                                                                    */
/*****************************************************************************/
%feature("autodoc", "Polymorphism history") poly_t;
%rename (polymorphism) poly_t;
%extend poly_t {
/* string representations */
%feature("autodoc", "x.__str__() <==> str(x)") __str__;
%feature("autodoc", "x.__repr__() <==> repr(x)") __repr__;
const char* __str__() {
        static char buffer[255];
        sprintf(buffer,"polymorphism: birth = %d, sweep time = %d, effect = %f, fitness = %f, fitness variance = %f",
                       (int)($self->birth),
                       (int)($self->sweep_time),
                       (double)($self->effect),
                       (double)($self->fitness),
                       (double)($self->fitness_variance));
        return &buffer[0];
}
const char* __repr__() {
        static char buffer[255];
        sprintf(buffer,"polymorphism(b=%d, age=%d, e=%f, f=%f, fvar=%f)",
                       (int)($self->birth),
                       (int)($self->sweep_time),
                       (double)($self->effect),
                       (double)($self->fitness),
                       (double)($self->fitness_variance));
        return &buffer[0];
}

/* read/write attributes */
%feature("autodoc", "Birth generation") birth;
%feature("autodoc", "Sweep time [in generations]") sweep_time;
%feature("autodoc", "Fitness effect of the mutation") effect;
%feature("autodoc", "Relative fitness of the clone at birth") fitness;
%feature("autodoc", "Fitness variance of the population at birth") fitness_variance;

} /* extend polymorphism */

/*****************************************************************************/
/* TREE_KEY_T                                                                */
/*****************************************************************************/
%feature("autodoc", "Key for a phylogenetic tree, with index and age.") tree_key_t;
%rename (tree_key) tree_key_t;
%extend tree_key_t {

/* string representations */
%feature("autodoc", "x.__str__() <==> str(x)") __str__;
%feature("autodoc", "x.__repr__() <==> repr(x)") __repr__;
const char* __str__() {
        static char buffer[255];
        sprintf(buffer,"tree key: index = %d, age = %d", (int)($self->index), (int)($self->age));
        return &buffer[0];
}
const char* __repr__() {
        static char buffer[255];
        sprintf(buffer,"tree_key(%d, %d)", (int)($self->index), (int)($self->age));
        return &buffer[0];
}

/* constructor */
%feature("autodoc",
"Initialize new tree_key.

Parameters:
   - index: index of the key
   - age: age of the key
") tree_key_t;

/* __hash__ needs to be overloaded if __eq__ has been so to be used as a dict key */
const long __hash__() {
        return (long)(($self->index) * RT_VERYLARGE + ($self->age));
}

/* read/write attributes */
%feature("autodoc", "Index of the key") index;
%feature("autodoc", "Age [in generations]") age;

} /* extend tree_key_t */

/*****************************************************************************/
/* STEP_T                                                                    */
/*****************************************************************************/
%feature("autodoc", "Step in a phylogenetic tree search") step_t;
%rename (tree_step) step_t;
%extend step_t {

/* string representations */
%feature("autodoc", "x.__str__() <==> str(x)") __str__;
%feature("autodoc", "x.__repr__() <==> repr(x)") __repr__;
const char* __str__() {
        static char buffer[255];
        sprintf(buffer,"tree_step: pos = %d, step = %d", (int)($self->pos), (int)($self->step));
        return &buffer[0];
}
const char* __repr__() {
        static char buffer[255];
        sprintf(buffer,"tree_step(%d, %d)", (int)($self->pos), (int)($self->step));
        return &buffer[0];
}

/* constructor */
%feature("autodoc",
"Initialize new step.

Parameters:
   - pos: position
   - step: length of step
") step_t;

/* __hash__ needs to be overloaded if __eq__ has been so to be used as a dict key */
const long __hash__() {
        return (long)(($self->pos) * RT_VERYLARGE + ($self->step));
}

/* read/write attributes */
%feature("autodoc", "Position") pos;
%feature("autodoc", "Step [in generations]") step;

} /* extend step_t */

/*****************************************************************************/
/* NODE_T                                                                    */
/*****************************************************************************/
%rename (tree_node) node_t;
%feature("autodoc", "Node of a phylogenetic tree") node_t;
%extend node_t {

/* string representations */
%feature("autodoc", "x.__str__() <==> str(x)") __str__;
%feature("autodoc", "x.__repr__() <==> repr(x)") __repr__;
const char* __str__() {
        static char buffer[255];
        sprintf(buffer,"tree_node");
        return &buffer[0];
}
const char* __repr__() {
        static char buffer[255];
        sprintf(buffer,"<tree_node>");
        return &buffer[0];
}

/* cloak child_edges with a Pythonic flavour */
%rename (_child_edges) child_edges;
%pythoncode
%{
@property
def child_edges(self):
    '''Child edges of the node'''
    return list(self._child_edges)


@child_edges.setter
def child_edges(self, es):
    self._child_edges = list_tree_key(es)
%}


/* cloak crossover */
%ignore crossover;
int _get_crossover_chunk(int i) {return ($self->crossover)[i];}
void _set_crossover_chunk(int value, int i) {($self->crossover)[i] = value;}
%pythoncode
%{
@property
def crossover(self):
    '''Crossover of node'''
    return [self._get_crossover_chunk(i) for i in xrange(2)]

@crossover.setter
def crossover(self, value):
    if len(value) != 2:
        raise ValueError('Crossover is a pair of integers.')
    [self._set_crossover_chunk(value[i], i) for i in xrange(2)]
%}

/* cloak weight_distribution with a Pythonic flavour */
%rename (_weight_distribution) weight_distribution;
%pythoncode
%{
@property
def weight_distribution(self):
    '''Distribution of weights of this node'''
    return list(self._weight_distribution)

@weight_distribution.setter
def weight_distribution(self, distr):
    self._weight_distribution = vector_tree_step(distr)
%}

/* read/write attributes */
%feature("autodoc", "Parent tree key") parent_node;
%feature("autodoc", "Own tree key") own_key;
%feature("autodoc", "Number of offspring") number_of_offspring;
%feature("autodoc", "Fitness  of the clone represented by the node") fitness;
%feature("autodoc", "Size of the clone represented by the node") clone_size;

} /* extend node_t */

/*****************************************************************************/
/* EDGE_T                                                                    */
/*****************************************************************************/
%rename (tree_edge) edge_t;
%feature("autodoc", "Edge of a phylogenetic tree") edge_t;
%extend edge_t {

/* string representations */
%feature("autodoc", "x.__str__() <==> str(x)") __str__;
%feature("autodoc", "x.__repr__() <==> repr(x)") __repr__;
const char* __str__() {
        static char buffer[255];
        sprintf(buffer,"tree_edge");
        return &buffer[0];
}
const char* __repr__() {
        static char buffer[255];
        sprintf(buffer,"<tree_edge>");
        return &buffer[0];
}

/* cloak segment */
%ignore segment;
int _get_segment_chunk(int i) {return ($self->segment)[i];}
%pythoncode
%{
@property
def segment(self):
    '''Segment of edge'''
    return [self._get_segment_chunk(i) for i in xrange(2)]
%}

/* read/write attributes */
%feature("autodoc", "Parent tree key") parent_node;
%feature("autodoc", "Own tree key") own_key;
%feature("autodoc", "Number of offspring") number_of_offspring;
%feature("autodoc", "Edge length [in generations]") length;

} /* extend edge_t */

/*****************************************************************************/
/* ROOTED_TREE                                                               */
/*****************************************************************************/
%feature("autodoc",
"Rooted phylogenetic tree.

This class is used to represent the phylogenetic tree of a single locus.
It is possible to print the tree in Newick format, to get the subtree
spanned by some of the leaves, and to look at the tree nodes and edges.
") rooted_tree;
%extend rooted_tree {

/* string representations */
%feature("autodoc", "x.__str__() <==> str(x)") __str__;
%feature("autodoc", "x.__repr__() <==> repr(x)") __repr__;
const char* __str__() {
        static char buffer[255];
        sprintf(buffer,"genealogy tree with %u nodes", (unsigned int)($self->nodes).size());
        return &buffer[0];
}
const char* __repr__() {
        static char buffer[255];
        sprintf(buffer,"<rooted_tree(%u nodes)>", (unsigned int)($self->nodes).size());
        return &buffer[0];
}



/* ignore weird stuff */
%ignore SFS;

/* redundant */
%ignore get_MRCA;

/* TODO: allow modifications of the tree (?) */
%ignore reset;
%ignore add_generation;
%ignore add_terminal_node;
%ignore erase_edge_node;
%ignore bridge_edge_node;
%ignore update_leaf_to_root;
%ignore update_tree;
%ignore erase_child;
%ignore delete_extra_children;
%ignore delete_one_child_nodes;
%ignore check_node;
%ignore check_tree_integrity;
%ignore clear_tree;

%feature("autodoc",
"Measure the length of the external branches.

Returns:
   - length: the sum of the lengths of the external branches.
") external_branch_length;

%feature("autodoc",
"Measure the length of the branches.

Returns:
   - length: the sum of the lengths of all branches.
") total_branch_length;

%feature("autodoc",
"Recalculate the weight of some internal nodes.

Parameters:
   - subtree_root: the tree_key of the node whose hanging subtree is recalculated

Returns:
   - error code: zero if successful
") calc_weight_distribution;

/* ancestors at age */
%ignore ancestors_at_age;
vector <tree_key_t> _ancestors_at_age(int age, tree_key_t subtree_root) {
        vector <tree_key_t> ancestors;
        $self->ancestors_at_age(age, subtree_root, ancestors);
        return ancestors;
}
%pythoncode
%{
def ancestors_at_age(self, age, subtree):
    '''Find nodes in subtree younger than a certain age
    
    Parameters:
       - age: critical age to check
       - subtree: subtree to look for nodes in
    
    Returns:
       - ancestors: the ancestors at that age
    '''
    return list(self._ancestors_at_age(age, subtree))
%}

/* construct subtree */
%ignore construct_subtree;
%feature("autodoc",
"Create a subtree from a list of leaves.

Parameters:
   - leaves: the leaves used to contruct the subtree

Returns:
   - subtree: the subtree spanned by those leaves.

.. note:: leaves can be a Python list or a numpy array of tree_key, or a vector_tree_key.
") create_subtree_from_keys;
rooted_tree create_subtree_from_keys(vector <tree_key_t> leaves) {
        rooted_tree other;
        $self->construct_subtree(leaves, other);
        return other;
}

%feature("autodoc",
"Print the tree in Newick format.

Returns:
   - tree: string of the tree in Newick format.

.. note:: You can pipe the output of this function to a cStingIO.StringIO
          for further manipulations.
") print_newick;

%feature("autodoc",
"Print a subtree in Newick format.

Parameters:
   - subtree_root: tree_key of the root of the subtree to print

Returns:
   - subtree: string of the subtree in Newick format.

.. note:: You can pipe the output of this function to a cStingIO.StringIO
          for further manipulations.
") subtree_newick;

%feature("autodoc",
"Read from Newick string.

Returns:
   - zero if successful, nonzero otherwise.
") read_newick;


/* cloak edges with a Pythonic flavour */
%rename (_edges) edges;
%pythoncode
%{
@property
def edges(self):
    '''Edges of the tree'''
    return dict(self._edges)


@edges.setter
def edges(self, es):
    self._edges = map_key_edge(es)
%}

/* cloak nodes with a Pythonic flavour */
%rename (_nodes) nodes;
%pythoncode
%{
@property
def nodes(self):
    '''Nodes of the tree'''
    return dict(self._nodes)


@nodes.setter
def nodes(self, ns):
    self._nodes = map_key_node(ns)
%}

/* cloak leafs with a Pythonic flavour */
%rename (_leafs) leafs;
%pythoncode
%{
@property
def leafs(self):
    '''Leaves of the tree'''
    return list(self._leafs)


@leafs.setter
def leafs(self, leaves):
    self._leafs = vector_tree_key(leaves)
%}

%pythoncode
%{
def to_Biopython_tree(self):
    '''Convert the tree into Biopython format
    
    Returns:
       - tree: Biopython.Phylo phylogenetic tree representation of self
    '''
    from cStringIO import StringIO
    from Bio import Phylo
     
    treedata = self.print_newick()
    handle = StringIO(treedata)
    tree = Phylo.read(handle, "newick")
    return tree
%}

} /* extend rooted_tree */

/*****************************************************************************/
/* MULTI_LOCUS_GENEALOGY                                                     */
/*****************************************************************************/
%feature("autodoc", "Genealogy for multiple loci") multi_locus_genealogy;
%extend multi_locus_genealogy {

/* string representations */
%feature("autodoc", "x.__str__() <==> str(x)") __str__;
%feature("autodoc", "x.__repr__() <==> repr(x)") __repr__;
const char* __str__() {
        static char buffer[255];
        sprintf(buffer,"multi_locus_genealogy for %u loci", (unsigned int)($self->loci).size());
        return &buffer[0];
}

const char* __repr__() {
        static char buffer[255];
        sprintf(buffer,"<multi_locus_genealogy(%u)>", (unsigned int)($self->loci).size());
        return &buffer[0];
}

/* hide weird stuff */
%ignore add_generation;
%ignore extend_storage;

/* document functions */
%feature("autodoc", "Default constructor") multi_locus_genealogy;
%feature("autodoc",
"Start tracking a new locus.

Parameters:
   - new_locus: locus to be tracked

.. note:: the locus gets appended to the 'loci' array.
") track_locus;
%feature("autodoc", "Reset (empty) the genealogy.") reset;
%feature("autodoc", "Reset (empty) the genealogy but keep the loci indices.") reset_but_loci;

/* loci */
%ignore loci;
int _get_number_of_loci() {
        return ($self->loci).size();
}
void _get_loci(int DIM1, int* ARGOUT_ARRAY1) {
        for(size_t i = 0; i < ($self->loci).size(); i++)
                ARGOUT_ARRAY1[i] = ($self->loci)[i];
}
%pythoncode
%{
@property
def loci(self):
    '''The loci that are being tracked'''
    return self._get_loci(self._get_number_of_loci())
%}

/* trees */
%ignore trees;
%feature("autodoc",
"Get the genealogy tree for a certain locus.

Parameters:
   - locus: site whose tree is being returned

.. note:: if you want to know what loci are being tracked, look into the 'loci'
          attribute.
") get_tree;
%exception get_tree {
        try {
                $action
        } catch (int err) {
                PyErr_SetString(PyExc_ValueError,"Locus not found among the tracked ones.");
                SWIG_fail;
        }
}
rooted_tree get_tree(int locus) {
        vector<int>::iterator index;
        index = std::find(($self->loci).begin(), ($self->loci).end(), locus);
        if(index == ($self->loci).end()) {
                throw (int)RT_LOCUSNOTFOUND;
        } else
                return ($self->trees)[(int)(index - ($self->loci).begin())];
}

%exception _set_tree {
        try {
                $action
        } catch (int err) {
                PyErr_SetString(PyExc_ValueError,"Locus not found among the tracked ones.");
                SWIG_fail;
        }
}
void _set_tree(int locus, rooted_tree &tree) {
        vector<int>::iterator index;
        index = std::find(($self->loci).begin(), ($self->loci).end(), locus);
        if(index == ($self->loci).end()) {
                throw (int)RT_LOCUSNOTFOUND;
        } else
                ($self->trees)[(int)(index - ($self->loci).begin())] = tree;
}

/* we need newGenerations for the tree reconstruction */
%ignore newGenerations;
%exception _get_newGeneration {
        try {
                $action
        } catch (int err) {
                PyErr_SetString(PyExc_ValueError,"Locus not found among the tracked ones.");
                SWIG_fail;
        }
}
vector <node_t> _get_newGeneration(int locus) {
        int i = 0;
        for(vector< vector<node_t> >::iterator it=$self->newGenerations.begin(); it != $self->newGenerations.end(); it++, i++) {
                if ($self->loci[i] == locus) {
                        return $self->newGenerations[i];
                }
        }
        throw (int)RT_LOCUSNOTFOUND;
}

%exception _set_newGeneration {
        try {
                $action
        } catch (int err) {
                PyErr_SetString(PyExc_ValueError,"Locus not found among the tracked ones.");
                SWIG_fail;
        }
}
void _set_newGeneration(int locus, vector<node_t> newGeneration) {
        int i = 0;
        for(vector< vector<node_t> >::iterator it=$self->newGenerations.begin(); it != $self->newGenerations.end(); it++, i++) {
                if (($self->loci)[i] == locus) {
                        ($self->newGenerations)[i] = newGeneration;
                        return;
                }
        }
        throw (int)RT_LOCUSNOTFOUND;
}


} /* extend multi_locus_genealogy */

/*****************************************************************************/
/* HAPLOID_HIGHD                                                             */
/*****************************************************************************/
/* Note: C++ includes empty clones, but Python does not. The filter is
       performed by the _nonempty_clones list.                               */
/*****************************************************************************/
%define DOCSTRING_HAPLOID_HIGHD
"Class for high-dimensional population genetics (genomes larger than ~20 loci).

This class is the main object for simulating the evolution of populations with
many loci (more than ~20). The class offers a number of functions, but an
example will explain the basic idea::

   ######################################
   #  EXAMPLE SCRIPT FOR HAPLOID_HIGHD  #
   ######################################
   import numpy as np
   import matplotlib.pyplot as plt
   import FFPopSim as h
   c = h.haploid_highd(300)       # 300 loci
   pop.set_wildtype(1000)         # start with 1000 wildtype individuals
   pop.mutation_rate = 1e-4       # mutation rate per site per generation
   pop.outcrossing_rate = 1e-1    # probability of sexual reproduction per gen
   pop.crossover_rate = 1e-2      # probability of crossover per site per gen
   pop.evolve(100)                # evolve for 100 generations
   c.plot_divergence_histogram()
   plt.show()
   ######################################

Populations can have a number of phenotypic traits that contribute to the fitness
of each individual. The function that calculates fitness from the phenotype
identifies fitness with the first trait only by default. The user is, however,
free to subclass haploid_highd in C++ (as it is done in hivpopulation) and
implement their own phenotype -> fitness function.

In addition, the trait landscapes describe the genotype -> phenotype maps.
These can be set directly from Python (since the genotypic space has a finite
number of elements).

**Note**: fitness is not a phenotypic trait directly, but rather a function of *all*
phenotypic traits together. 
"
%enddef
%feature("autodoc", DOCSTRING_HAPLOID_HIGHD) haploid_highd;
%extend haploid_highd {
/* constructor */
%feature("autodoc",
"Construct a high-dimensional population with certain parameters.

Parameters:
   - L: number of loci
   - rng_seed: seed for the random generator. If zero (default) pick a random number
   - number_of_traits: number of phenotypic traits, defaults to one
   - all_polymorphic: option to use an infinite-sites model tracking ancestral alleles
                      (only available with a single phenotypic trait and zero mutation rate)
") haploid_highd;
%exception haploid_highd {
        try {
                $action
        } catch (int err) {
                PyErr_SetString(PyExc_ValueError,
                                "Construction impossible. Please check input args.");
                SWIG_fail;
        }
}

/* string representations */
%feature("autodoc", "x.__str__() <==> str(x)") __str__;
%feature("autodoc", "x.__repr__() <==> repr(x)") __repr__;
const char* __str__() {
        static char buffer[255];
        sprintf(buffer,"haploid_highd: L = %d, N = %d", $self->L(), $self->N());
        return &buffer[0];
}
const char* __repr__() {
        static char buffer[255];
        sprintf(buffer,"<haploid_highd(%d, %d)>", $self->L(), $self->N());
        return &buffer[0];
}

/* ignore hypercubes */
%ignore trait;

/* ignore stream stuff we never need in Python */
%ignore print_allele_frequencies;
%ignore get_genotype_string;
%ignore read_ms_sample;
%ignore read_ms_sample_sparse;

/* ignore weird functions using pointers */
%ignore get_pair_frequencies(vector < vector <int> > *loci);

/* read/write attributes */
%feature("autodoc", "is the genome circular?") circular;
%feature("autodoc", "current carrying capacity of the environment") carrying_capacity;
%feature("autodoc", "outcrossing rate (probability of sexual reproduction per generation)") outcrossing_rate;
%feature("autodoc", "crossover rate (probability of crossover per site per generation)") crossover_rate;
%feature("autodoc",
"Model of recombination to use

Available values:
   - FFPopSim.FREE_RECOMBINATION: free reassortment of all loci between parents
   - FFPopSim.CROSSOVERS: linear chromosome with crossover probability per locus
") recombination_model;
%feature("autodoc",
"Growth rate

This value is used to determine how fast a population converges to the
carrying capacity.

This parameter must be set strictly larger than 1 (very slow growth) and not
too big to avoid population explosion. The default is 2, which means that a
freely expanding population (N << carrying capacity) approximately doubles in
size every generation.

Note that when the population is shrinking, in order to avoid extinction, the
population decreases by ten times or so only. If you want a hard bottleneck,
use the bottleneck function.
") growth_rate;

/* mutation rate */
%rename(_get_mutation_rate) get_mutation_rate;
%rename(_set_mutation_rate) set_mutation_rate;
%pythoncode
%{
@property
def mutation_rate(self):
   '''mutation rate (per site per generation)'''
   return self._get_mutation_rate()

@mutation_rate.setter
def mutation_rate(self, m):
    if self.all_polymorphic:
        raise ValueError("You cannot set all_polymorphic and a nonzero mutation rate.")
    else:
        self._set_mutation_rate(m)
%}

/* do not expose the population, but rather only nonempty clones */
%ignore population;
%rename (_get_nonempty_clones) get_nonempty_clones;
%feature("autodoc", 
"Get a single clone

Parameters:
   - n: index of the clone

Returns:
   - clone: the n-th clone in the population
") get_clone;
clone_t get_clone(unsigned long n) {
        return $self->population.at(n);

}
%pythonprepend get_clone {
if len(args):
    args = list(args)
    if args[0] >= len(self._nonempty_clones):
        raise ValueError('The population has only '+str(len(self._nonempty_clones))+' clones.')
    args[0] = self._nonempty_clones[args[0]]
    args = tuple(args)
}

/* read only parameters */
%ignore get_number_of_loci;
%feature("autodoc", "Number of loci (read-only)") L;
%feature("autodoc", "Number of loci (read-only)") number_of_loci;
const int L;
const int number_of_loci;

%ignore get_population_size;
%feature("autodoc", "Population size (read-only)") N;
%feature("autodoc", "Population size (read-only)") population_size;
const int N;
const int population_size;

%ignore get_generation;
%ignore set_generation;
%feature("autodoc", "Current generation (read-only)") generation;
int generation;

%ignore get_number_of_clones;
%feature("autodoc", "Number of non-empty clones (read-only)") number_of_clones;
const int number_of_clones;

%ignore get_number_of_traits;
%feature("autodoc", "Number of traits (read-only)") number_of_traits;
const int number_of_traits;

%ignore get_max_fitness;
%feature("autodoc", "Maximal fitness in the population (read-only)") max_fitness;
const double max_fitness;

%ignore get_participation_ratio;
%feature("autodoc", "Participation ratio (read-only)") participation_ratio;
const double participation_ratio;

%ignore is_all_polymorphic;
%feature("autodoc", "All polymorphic? (read-only)") all_polymorphic;
const bool all_polymorphic;

%rename(_get_polymorphisms) get_polymorphisms;
%rename(_get_fixed_mutations) get_fixed_mutations;
%rename(_get_number_of_mutations) get_number_of_mutations;
%pythoncode
%{
@property
def polymorphisms(self):
    '''Polymorphisms from all_polymorphic (read-only)'''
    if not self.all_polymorphic:
        raise ValueError("all_polymorphic is not set.")
    return self._get_polymorphisms()


@property
def fixed_mutations(self):
    '''Fixed mutations from all_polymorphic (read-only)'''
    if not self.all_polymorphic:
        raise ValueError("all_polymorphic is not set.")
    return self._get_fixed_mutations()


@property
def number_of_mutations(self):
    '''Fixed mutations from all_polymorphic (read-only)'''
    if not self.all_polymorphic:
        raise ValueError("all_polymorphic is not set.")
    return self._get_number_of_mutations()
%}

/* trait weights */
%ignore set_trait_weights;
%pythonprepend _set_trait_weights {
if len(args) and (len(args[0]) != self.number_of_traits):
    raise ValueError('The weights must be a sequence of length equal to the number of traits.')
}
void _set_trait_weights(double* IN_ARRAY1, int DIM1) {
        /* call the C++ method */
        $self->set_trait_weights(IN_ARRAY1);
        $self->update_fitness();
}
%feature("autodoc",
"weight of each trait on fitness

.. note:: Fitness is updated automatically when the weights are changed.
") _get_trait_weights;
%pythonprepend _get_trait_weights {
args = tuple(list(args) + [self.number_of_traits])
}
void _get_trait_weights(double* ARGOUT_ARRAY1, int DIM1) {
        /* check trait number */
        if(DIM1 != $self->get_number_of_traits())
                throw HP_BADARG; 

        /* set the output array */
        for(size_t t=0; t < (size_t)DIM1; t++)
                ARGOUT_ARRAY1[t] = $self->get_trait_weight(t);
}
%pythoncode
%{
trait_weights = property(_get_trait_weights, _set_trait_weights)
%}

/* dump to file */
%pythoncode
%{
def dump(self, filename, format='bz2', include_genealogy=False):
    '''Dump a population to binary file, for later use.

    Parameters:
       - filename: the path to the file where to store the information
       - format: one of 'bz2' or 'plain'. Choose the former if you want compression.
       - include_genealogy: if True, the multi_locus_genealogy is stored as well (if present).

    .. note:: The population can be reloaded using the function FFPopSim.load_haploid_highd.
    '''

    try:
        import cPickle as pickle
    except:
        import pickle

    pop_dict = {}
    pop_dict['genotypes'] = self.get_genotypes()
    pop_dict['N'] = self.carrying_capacity
    pop_dict['L'] = self.L
    pop_dict['mu'] = self.mutation_rate
    pop_dict['crossover_rate'] = self.crossover_rate
    pop_dict['outcrossing_rate'] = self.outcrossing_rate
    pop_dict['circular'] = self.circular
    pop_dict['generation'] = self.generation
    pop_dict['clone_sizes'] = self.get_clone_sizes()
    pop_dict['recombination_model'] = self.recombination_model
    pop_dict['traits_additive'] = [self.get_trait_additive(i) for i in range(self.number_of_traits)]
    pop_dict['traits_epistasis'] = [self.get_trait_epistasis(i) for i in range(self.number_of_traits)]
    pop_dict['all_polymorphic']  = self.all_polymorphic
    pop_dict['ancestral'] = self.get_ancestral_states()
    pop_dict['trait_weights'] = self.trait_weights

    # Genealogy
    if include_genealogy and len(self.genealogy.loci):
        pop_dict['trees'] = {locus: self.genealogy.get_tree(locus).print_newick() for locus in self.genealogy.loci}
        pop_dict['_nonempty_clones'] = self._nonempty_clones

        # Save newGenerations as a non-SWIG object
        def serialize_leaf(leaf):
            serial = {}
            for key in ['clone_size', 'crossover', 'fitness', 'number_of_offspring']:
                serial[key] = getattr(leaf, key)
            serial['own_key'] = (leaf.own_key.index, leaf.own_key.age)
            serial['parent_node'] = (leaf.parent_node.index, leaf.parent_node.age)
            return serial
            
        newGenerations = []
        for locus in self.genealogy.loci:
            newGenerations.append(map(serialize_leaf, self.genealogy._get_newGeneration(locus)))
        pop_dict['_newGenerations'] = newGenerations

    
    with open(filename, 'wb') as f:
        dump = pickle.dumps(pop_dict, pickle.HIGHEST_PROTOCOL)

        # Try to compress if the user wishes so
        try:
            if format == 'bz2':
                import bz2
                dump = dump.encode('bz2')
        # Fallback on uncompressed
        except:
            import warnings
            warnings.warn('compression module ('+format+') not found. Defaulting to uncompressed file.')
            format = 'plain'

        # Dump to file
        f.write(dump)
%}

/* copy */
%pythoncode
%{
def copy(self, rng_seed=0):
    '''Copy population into new instance.
    
    Parameters:
       - rng_seed: random number to initialize the new population
    '''
    pop = haploid_highd(self.L, rng_seed=rng_seed, number_of_traits=self.number_of_traits)

    # Mutation and recombination
    pop.recombination_model =  self.recombination_model
    pop.outcrossing_rate = self.outcrossing_rate
    pop.crossover_rate = self.crossover_rate
    pop.mutation_rate = self.mutation_rate
    pop.circular = self.circular

    # Fitness
    for i in xrange(self.number_of_traits):
        pop.set_trait_additive(self.get_trait_additive(i), i)
        for coeff in self.get_trait_epistasis(i):
            pop.add_trait_coefficient(coeff[0], coeff[1], i)

    # Population parameters
    pop.carrying_capacity = self.carrying_capacity
    pop.set_genotypes(self.get_genotypes(), self.get_clone_sizes())    

    # Evolution
    pop.generation = self.generation
    
    return pop
%}

/* status function */
%pythoncode
%{
def status(self):
    '''Print a status list of the population parameters'''
    parameters = (('number of loci', 'L'),
                  ('circular', 'circular'),
                  ('number of traits', 'number_of_traits'),
                  ('population size', 'N'),
                  ('carrying capacity', 'carrying_capacity'),
                  ('generation', 'generation'),
                  ('outcrossing rate', 'outcrossing_rate'),
                  ('crossover rate', 'crossover_rate'),
                  ('recombination model', 'recombination_model'),
                  ('mutation rate', 'mutation_rate'),
                  ('participation ratio', 'participation_ratio'),
                  ('number of non-empty clones', 'number_of_clones'),
                 )
    lenmax = max(map(lambda x: len(x[0]), parameters))

    for (strin, name) in parameters:
        par = getattr(self, name)
        # Recombination model needs a conversion
        # (a very frequently used one, to be honest)
        if strin == 'recombination model':
            if par == 0:
                par = 'FREE_RECOMBINATION'
            else:
                par = 'CROSSOVERS'
        print ('{:<'+str(lenmax + 2)+'s}').format(strin)+'\t'+str(par)
%}

/* initialize wildtype */
%feature("autodoc",
"Initialize a population of wildtype individuals

Parameters:
   - N: the number of individuals

.. note:: the carrying capacity is set to the same value if still unset.
") set_wildtype;
%pythonappend set_wildtype {
self._nonempty_clones = _np.array(self._get_nonempty_clones())
return None
}

/* initalize frequencies */
%feature("autodoc",
"Initialize the population according to the given allele frequencies in linkage equilibrium.

Parameters:
   - frequencies: an array of length L with all allele frequencies
   - N: set the population size and, if still unset, the carrying
     capacity to this value
") set_allele_frequencies;
%pythonprepend set_allele_frequencies {
if len(args) and (len(args[0]) != self.L):
    raise ValueError('Please input an L dimensional list of allele frequencies.')
}
%exception set_allele_frequencies {
  $action
  if (result) {
     PyErr_SetString(PyExc_RuntimeError,"Error in the C++ function.");
     SWIG_fail;
  }
}
%pythonappend set_allele_frequencies {
self._nonempty_clones = _np.array(self._get_nonempty_clones())
return None
}

/* set genotypes */
%ignore set_genotypes(vector <genotype_value_pair_t> gt);
%feature("autodoc",
"Initialize population with fixed counts for specific genotypes.

Parameters:
   - genotypes: list of genotypes to set. Genotypes are lists of alleles,
     e.g. [[0,0,1,0], [0,1,1,1]] for genotypes 0010 and 0111   
   - counts: list of the number at which each of those genotypes it to be present

.. note:: the population size and, if unset, the carrying capacity will be set
          as the sum of the counts.

**Example**: if you want to initialize 200 individuals with genotype 001 and
             300 individuals with genotype 110, you can use
             ``set_genotypes([[0,0,1], [1,1,0]], [200, 300])``
") set_genotypes;
%pythonprepend set_genotypes {
if len(args) and (len(args) >= 2):
    genotypes = args[0]
    counts = args[1]
    genotypes = _np.array(genotypes, float, copy=False, ndmin=2)
    counts = _np.asarray(counts, float)
    if len(genotypes) != len(counts):
        raise ValueError('Genotypes and counts must have the same length')
    args = tuple([genotypes.ravel(), counts] + list(args[2:]))
}
%exception set_genotypes {
  $action
  if (result) {
     PyErr_SetString(PyExc_RuntimeError,"Error in the C++ function.");
     SWIG_fail;
  }
}
%pythonappend set_genotypes {
self._nonempty_clones = _np.array(self._get_nonempty_clones())
return None
}
%apply (int DIM1, double* IN_ARRAY1) {(int len1, double* genotypes), (int len2, double* counts)};
int set_genotypes(int len1, double* genotypes, int len2, double* counts) {
        /* We use a flattened array */
        len1 /= len2;
        vector<genotype_value_pair_t> gt;
        genotype_value_pair_t temp;
        for(size_t i = 0; i != (size_t)len2; i++) {
                temp.genotype = boost::dynamic_bitset<>(len1);
                for(size_t j=0; j < (size_t)len1; j++)
                        temp.genotype[j] = (bool)genotypes[i * len1 + j];
                temp.val = counts[i];
                gt.push_back(temp);
        }
        return $self->set_genotypes(gt);
}
%clear (int len1, double* genotypes);
%clear (int len2, double* counts);


/* set genotypes with ancestral state*/
%ignore set_genotypes_and_ancestral_state(vector <genotype_value_pair_t> gt, vector <int> anc_state);
%feature("autodoc",
"Initialize population with fixed counts for specific genotypes.

Parameters:
   - genotypes: list of genotypes to set. Genotypes are lists of alleles,
     e.g. [[0,0,1,0], [0,1,1,1]] for genotypes 0010 and 0111   
   - counts: list of the number at which each of those genotypes it to be present
   - ancestral state of the sample, a vector of 0 and 1
.. note:: the population size and, if unset, the carrying capacity will be set
          as the sum of the counts.

**Example**: if you want to initialize 200 individuals with genotype 001 and
             300 individuals with genotype 110, you can use
             ``set_genotypes([[0,0,1], [1,1,0]], [200, 300])``
") set_genotypes_and_ancestral_state;
%pythonprepend set_genotypes_and_ancestral_state {
if len(args) and (len(args) >= 3):
    genotypes = args[0]
    counts = args[1]
    anc_state = args[2]
    genotypes = _np.array(genotypes, float, copy=False, ndmin=2)
    counts = _np.asarray(counts, float)
    anc_state = _np.asarray(anc_state, float)
    if len(genotypes) != len(counts):
        raise ValueError('Genotypes and counts must have the same length')
    if (len(anc_state) != self.L):
        raise ValueError('Ancestral state vector must have length L')
    args = tuple([genotypes.ravel(), counts] + list(args[2:]))
}
%exception set_genotypes_and_ancestral_state {
  $action
  if (result) {
     PyErr_SetString(PyExc_RuntimeError,"Error in the C++ function.");
     SWIG_fail;
  }
}
%pythonappend set_genotypes_and_ancestral_state {
self._nonempty_clones = _np.array(self._get_nonempty_clones())
return None
}
%apply (int DIM1, double* IN_ARRAY1) {(int len1, double* genotypes), (int len2, double* counts), (int len3, double* anc_state)};
int set_genotypes_and_ancestral_state(int len1, double* genotypes, int len2, double* counts, int len3, double* anc_state) {
        /* We use a flattened array */
        len1 /= len2;
        vector<genotype_value_pair_t> gt;
        genotype_value_pair_t temp;
        for(size_t i = 0; i != (size_t)len2; i++) {
                temp.genotype = boost::dynamic_bitset<>(len1);
                for(size_t j=0; j < (size_t)len1; j++)
                        temp.genotype[j] = (bool)genotypes[i * len1 + j];
                temp.val = counts[i];
                gt.push_back(temp);
        }
                vector <int> ancestral_state($self->L(), 0);
                for (size_t locus=0; locus<len3; locus++){
                  ancestral_state[locus]=(anc_state[locus]<0.5)?0:1;
                }
        return $self->set_genotypes_and_ancestral_state(gt, ancestral_state);
}
%clear (int len1, double* genotypes);
%clear (int len2, double* counts);
%clear (int len3, double* anc_state);


/* evolve */
%feature("autodoc",
"Evolve for some generations.

Parameters:
   - gen: number of generations, defaults to one
") evolve;
%exception evolve {
  $action
  if (result) {
     PyErr_SetString(PyExc_RuntimeError,"Error in the C++ function.");
     SWIG_fail;
  }
}
%pythonappend evolve {
self.calc_stat()
self._nonempty_clones = _np.array(self._get_nonempty_clones())
return None
}

/* bottleneck */
%feature("autodoc",
"Make the population undergo a bottleneck

Parameters:
   - size_of_bottleneck: the number of individuals at the bottleneck
") bottleneck;
%pythonappend bottleneck {
self.calc_stat()
self._nonempty_clones = _np.array(self._get_nonempty_clones())
return None
}

/* flip single locus */
%feature("autodoc",
"
Take a random clone, flip the allele at the selected locus and create a new clone.

Parameters:
   - locus: locus to flip

Returns:
   - index: index of the new clone with the flipped locus
") flip_single_locus;

/* genealogy */
%feature("autodoc",
"
Track the genealogy of some loci.

Parameters:
   - loci: sites whose genealogy is being stored

Returns:
   - zero if successful
") track_locus_genealogy;
%exception track_locus_genealogy {
        $action
        if (result) {
                PyErr_SetString(PyExc_ValueError,"Track the genealogy before initializing the population.");
                SWIG_fail;
        }
}

/* implement multi_locus_genealogy as a read-only property */
%ignore genealogy;
%feature("autodoc",
"Genealogy of the tracked loci.

.. note:: This attribute is read-only.
") _get_genealogy;
multi_locus_genealogy _get_genealogy() {
        return $self->genealogy;
}
%pythoncode
%{
genealogy = property(_get_genealogy)
%}

/* statistics */
%feature("autodoc", "Calculate trait and fitness statistics for the population") calc_stat;

%feature("autodoc",
"Get the mean and variance of the divergence in the population.

Parameters:
   - n_sample: number of individuals to sample at random from the population

Returns:
   - stat: structure with mean and variance of divergence in the population
") get_divergence_statistics;

%feature("autodoc", 
"Get the mean and variance of the diversity in the population.

Parameters:
   - n_sample: number of individuals to sample at random from the population

Returns:
   - stat: structure with mean and variance of diversity in the population
") get_diversity_statistics;

%feature("autodoc",
"Get the mean and variance of a trait in the population.

Parameters:
   - t: number of the trait whose statistics are to be computed

Returns:
   - stat: structure with mean and variance of the trait in the population
") get_trait_statistics;

%feature("autodoc",
"Get the mean and variance of the fitness in the population.

Returns:
   - stat: structure with mean and variance of the fitness in the population
") get_fitness_statistics;

%feature("autodoc",
"Get the covariance of two traits in the population.

Parameters:
   - t1: first trait
   - t2: second trait

Returns:
   - cov: the covariance of the two traits
") get_trait_covariance;

/* get allele frequencies */
%feature("autodoc", "Get all allele frequencies") get_allele_frequencies;
%pythonprepend get_allele_frequencies {
args = tuple(list(args) + [self.L])
}
void get_allele_frequencies(double* ARGOUT_ARRAY1, int DIM1) {
        for(size_t i=0; i < (size_t)$self->get_number_of_loci(); i++)
                ARGOUT_ARRAY1[i] = $self->get_allele_frequency(i);
}

%feature("autodoc",
"Get the frequency of the + allele at the selected locus

Parameters:
   - locus: locus whose frequency of the + allele is to be returned

Returns:
   - frequency: allele frequency in the population
") get_allele_frequency;

/* get allele frequencies */
%feature("autodoc", "Get all derived allele frequencies") get_derived_allele_frequencies;
%pythonprepend get_derived_allele_frequencies {
args = tuple(list(args) + [self.L])
}
void get_derived_allele_frequencies(double* ARGOUT_ARRAY1, int DIM1) {
        if ($self->is_all_polymorphic()){
                for(size_t i=0; i < (size_t)$self->get_number_of_loci(); i++)
                        ARGOUT_ARRAY1[i] = $self->get_derived_allele_frequency(i);
    }
}


%feature("autodoc",
"Get the frequency of the derived allele at the selected locus

Parameters:
   - locus: locus whose frequency of the derived allele is to be returned

Returns:
   - frequency: allele frequency in the population
") get_derived_allele_frequency;

/* get ancestral state of all loci */
%feature("autodoc", "Get ancestral state of all loci") get_ancestral_states;
%pythonprepend get_ancestral_states {
args = tuple(list(args) + [self.L])
}
void get_ancestral_states(double* ARGOUT_ARRAY1, int DIM1) {
  for(size_t i=0; i < (size_t)$self->get_number_of_loci(); i++)
	ARGOUT_ARRAY1[i] = $self->get_ancestral_state(i);
}


%feature("autodoc",
"Get the joint frequency of two + alleles

Parameters:
    - locus1: first locus
    - locus2: second locus

Returns:
    - the joint frequency of the + alleles
") get_pair_frequency;

%feature("autodoc",
"Get :math:`\\chi_i` of an allele

Parameters:
    - locus: locus whose chi is to be computed

Returns:
    - the chi of that allele, :math:`\\chi_i := \\left<s_i\\right>`, where :math:`s_i \\in \{\\pm1\}`.
") get_chi;

%feature("autodoc",
"Get :math:`\\chi_i` of a derived allele

Parameters:
    - locus: locus whose chi is to be computed

Returns:
    - the chi of that derived allele, :math:`\\chi_i := \\left<s_i\\right>`, where :math:`s_i \\in \{\\pm1\}`.
") get_derived_chi;

%feature("autodoc",
"Get :math:`\\chi_{ij}`

Parameters:
    - locus1: first locus
    - locus2: second locus

Returns:
    - the linkage disequilibiurm between them, i.e. :math:`\\chi_{ij} := \\left<s_i s_j\\right> - \\chi_i \\cdot \\chi_j`.
") get_chi2;

%feature("autodoc",
"Get linkage disequilibrium

Parameters:
    - locus1: first locus
    - locus2: second locus

Returns:
    - the linkage disequilibiurm between them, i.e. :math:`LD := 1 / 4 \\left[\\left<s_i s_j\\right> - \\chi_i \\cdot \\chi_j\\right]`.
") get_LD;

%feature("autodoc",
"Get moment of two alleles in the -/+ basis

Parameters:
    - locus1: first locus
    - locus2: second locus

Returns:
    - the second moment, i.e. :math:`\\left<s_i s_j\\right>`, where :math:`s_i, s_j \in \{-1, 1\}`.
") get_moment;

/* add genotypes */
%feature("autodoc",
"Add new individuals to the population with certain genotypes

Parameters:
   - genotype: genotype to add to the population (Boolean list)
   - n: number of new individuals carrying that genotype
") add_genotype;
%pythonappend add_genotype {
    self._nonempty_clones = _np.array(self._get_nonempty_clones())
}

/* get single locus effects */
%feature("autodoc",
"Get an array with the additive coefficients of all loci of a trait. 

Parameters:
   - t: number of the trait

Returns:
   - coefficients: array of additive coefficients for the selected trait
") get_trait_additive;
%pythonprepend get_trait_additive {
if (len(args) > 1) and (args[1] >= self.number_of_traits):
    raise ValueError("There are only "+str(self.number_of_traits)+" traits.")
args = tuple([self.L] + list(args))
}
void get_trait_additive(double* ARGOUT_ARRAY1, int DIM1, int t=0) {
        /* Initialize to zero */
        for(size_t i=0; i < (size_t)DIM1; i++)
                ARGOUT_ARRAY1[i] = 0;

        /* Add any coefficient you found */
        hypercube_highd *trait = &(($self->trait)[t]);
        coeff_single_locus_t * coeff;
        for(size_t i=0; i < trait->coefficients_single_locus.size(); i++) {
                coeff = &(trait->coefficients_single_locus[i]);
                ARGOUT_ARRAY1[coeff->locus] += coeff->value;
        }
}

/* update functions we need in hivpopulation */
%rename (_update_traits) update_traits;
%rename (_update_fitness) update_fitness;

/* clear trait/fitness coefficients */
%feature("autodoc",
"Clear a trait landscape.

Parameters:
   - t: number of the trait to be cleared
") clear_trait;
%feature("autodoc", "Shortcut for clear_trait when there is only one trait") clear_fitness;
%exception clear_fitness {
        try {
                $action
        } catch (int err) {
                PyErr_SetString(PyExc_ValueError,"Fitness depends only on traits, not on the genome directly.");
                SWIG_fail;
        }
}
%feature("autodoc", "Clear all trait landscapes") clear_traits;

/* set single locus effects */
%feature("autodoc",
"Set the additive part of a trait

Parameters:
   - coefficients: array of coefficients for the trait (of length L). All previous additive coefficents are erased
   - t: number of the trait to set
") set_trait_additive;
%pythonprepend set_trait_additive {
if (len(args) > 1) and (args[1] >= self.number_of_traits):
    raise ValueError("There are only "+str(self.number_of_traits)+" traits.")
if len(args) and (len(args[0]) != self.L):
    raise ValueError("L coefficients expected.")
}
void set_trait_additive(int DIM1, double* IN_ARRAY1, int t=0) {
        /* reset trait landscape */
        $self->trait[t].reset_additive();
        
        /* set the new coefficients */
        vector <int> loci(1,0);
        for(size_t i = 0; i < (size_t)DIM1; i++) {
                if(abs(IN_ARRAY1[i]) > HP_NOTHING) {
                        loci[0] = i;
                        $self->add_trait_coefficient(IN_ARRAY1[i], loci, t);
                }
        }

        /* update the population */
        $self->update_traits();
        $self->update_fitness();
}

%feature("autodoc", "Shortcut for set_trait_additive when there is only one trait") set_fitness_additive;
%pythonprepend set_fitness_additive {
if len(args) and (len(args[0]) != self.L):
    raise ValueError("L coefficients expected.")
}
void set_fitness_additive(int DIM1, double *IN_ARRAY1) {
        /* reset trait landscape */
        $self->trait[0].reset_additive();
        
        /* set the new coefficients */
        vector <int> loci(1,0);
        for(size_t i = 0; i < (size_t)DIM1; i++) {
                if(abs(IN_ARRAY1[i]) > HP_NOTHING) {
                        loci[0] = i;
                        $self->add_trait_coefficient(IN_ARRAY1[i], loci, 0);
                }
        }

        /* update the population */
        $self->update_traits();
        $self->update_fitness();
}

%feature("autodoc",
"Add a coefficient to the trait landscape.
 
Parameters:
   - value: value of the coefficient
   - loci: array/list of loci indexed by the coefficient.
   - t: number of the trait to be changed

**Example**: to set a second-order epistatic term :math:`t_{ij} = 0.1`, use ``add_trait_coefficient(0.1, [i, j])``.

.. warning:: the -/+ basis is used throughout the library. If you are used to the 0/1 basis, keep in mind that the interaction series-expansion is different.
") add_trait_coefficient;

%feature("autodoc", "Shortcut for add_trait_coefficient when there is only one trait") add_fitness_coefficient;

/* get epistatic terms */
/* Note: the output is as compatible as possible with add_trait_coefficient,
         i.e. a tuple of (value, loci), where loci is a tuple itself. */
%feature("autodoc",
"Get the epistatic terms of the genotype/phenotype map of the chosen trait.

Parameters:
   - t: trait number

Returns:
   - coefficients: tuple of coefficients, with a value and a tuple of loci

.. note:: This function is designed to work well in conjunction with add_trait_coefficient.
") get_trait_epistasis;

/* random epistasis */
%feature("autodoc",
"Set a random epistatic term in the genotype-phenotype map. This is meant as an approximation to multi-locus epistasis to which many locus sets contribute. It assigns to each genotype a reprodrucible fitness component drawn from a Gaussian distribution.

Parameters:
   - epistasis_std: standard deviation of the random epistatic terms
   - t: trait number

.. note:: the epistatic terms will be Gaussian distributed around zero with the given standard deviation.
") set_random_trait_epistasis;

%feature("autodoc", "Shortcut for set_random_trait_epistasis when there is only one trait") set_random_epistasis;

/* fitness of clones */
%pythonprepend get_fitness {
if len(args) and (args[0] >= self.number_of_clones):
    raise ValueError('The population has only '+str(self.number_of_clones)+' clones.')
if len(args):
    args = list(args)
    args[0] = self._nonempty_clones[args[0]]
    args = tuple(args)
}
%feature("autodoc",
"Get the fitness of an individual

Parameters:
   - n: index of the clone whose fitness is to be computed

Returns:
   - fitness: fitness value of that clone
") get_fitness;

%pythoncode
%{
def get_fitnesses(self):
    '''Get the fitness of all clones.'''
    f = _np.zeros(self.number_of_clones)
    for i in xrange(self.number_of_clones):
        f[i] = self.get_fitness(i)
    return f
%}

/* traits of clones */
%pythonprepend get_trait {
if (len(args) > 1) and (args[1] >= self.number_of_traits):
    raise ValueError("There are only "+str(self.number_of_traits)+" traits.")
if len(args) and (args[0] >= self.number_of_clones):
    raise ValueError('The population has only '+str(self.number_of_clones)+' clones.')
if len(args):
    args = list(args)
    args[0] = self._nonempty_clones[args[0]]
    args = tuple(args)
}
%feature("autodoc",
"Get a trait of an individual

Parameters:
   - n: index of the clone whose trait is to be computed
   - t: trait to be computed

Returns:
   - trait: value of that trait for that clone
") get_trait;

%pythoncode
%{
def get_traits(self):
    '''Get all traits from all clones'''
    t = _np.zeros((self.number_of_clones, self.number_of_traits))
    for i in xrange(self.number_of_clones):
        for j in xrange(self.number_of_traits):
            t[i, j] = self.get_trait(i, j)
    return t
%}

/* get clone sizes */
%pythonprepend get_clone_size {
if len(args) and (args[0] >= self.number_of_clones):
    raise ValueError('The population has only '+str(self.number_of_clones)+' clones.')
if len(args):
    args = list(args)
    args[0] = self._nonempty_clones[args[0]]
    args = tuple(args)
}
%feature("autodoc", 
"Get the size of a clone

Parameters:
   - n: index of the clone

Returns:
   - size: size of the selected clone
") get_clone_size;

%pythoncode
%{
def get_clone_sizes(self):
    '''Get the size of all clones.'''
    s = _np.zeros(self.number_of_clones, int)
    for i in xrange(self.number_of_clones):
        s[i] = self.get_clone_size(i)
    return s
%}

/* get genotypes */
%pythonprepend get_genotype {
if len(args) and (args[0] >= self.number_of_clones):
    raise ValueError('The population has only '+str(self.number_of_clones)+' clones.')
if len(args):
    args = list(args)
    args[0] = self._nonempty_clones[args[0]]
    args = tuple(args)
}
boost::dynamic_bitset<> get_genotype(int n) {
        return $self->population[n].genotype;
}
%feature("autodoc",
"Get a genotype from the population

Parameters:
   - n: index of the clone whose genotype is to be returned

Returns:
   - genotype: Boolean array of the genotype
") get_genotype;

%pythoncode
%{
def get_genotypes(self):
    '''Get all genotypes of the population.

    Return:
       - genotypes: boolean 2D array with the genotypes

    .. note:: this function does not return the sizes of each clone.
    '''
    genotypes = _np.zeros((self.number_of_clones, self.number_of_loci), bool)
    for i in xrange(self.number_of_clones):
        genotypes[i] = self.get_genotype(i)
    return genotypes
%}

/* unique clones */
%feature("autodoc",
"Recompress the clone structure

During its evolution, identical clones might be generated by different routes at
different times. This function merges any such duplicates into unique clones with
the size equal to the sum of the sizes of the duplicates.
") unique_clones;

/* Hamming distance (full Python reimplementation) */
%ignore distance_Hamming;
%pythoncode
%{
def distance_Hamming(self, clone_gt1, clone_gt2, chunks=None, every=1):
    '''Calculate the Hamming distance between two genotypes

    Parameters:
       - clone_gt1: index of the clone corresponding to the first genotype
       - clone_gt2: index of the clone corresponding to the second genotype
       - chunks: list of pairs delimiting the genetic areas to include
       - every: do the comparison only on certain sites

    **Example**: to calculate the distance between the first two clones
    limited to third codon positions between locus 90 and 200, use:
    ``distance_Hamming(0, 1, chunks=[92, 200], every=3)``.
    '''
    if _np.isscalar(clone_gt1):
        genotypes = self.get_genotypes((clone_gt1, clone_gt2))
        clone_gt1 = genotypes[0]
        clone_gt2 = genotypes[1]

    if chunks is not None:
        ind = _np.zeros(clones.shape[1], bool)
        for chunk in chunks:
            inde = _np.arange(chunk[1] - chunk[0])
            inde = inde[(inde % every) == 0] + chunk[0]
            ind[inde] = True
        clone_gt1 = clone_gt1[ind]
        clone_gt2 = clone_gt2[ind]
    return (clone_gt1 != clone_gt2).sum()
%}

/* get random clones/genotypes */
%pythoncode
%{
def random_genomes(self, n):
    '''Get a sample of random genomes from the population

    Parameters:
       - n: number of random genomes to compute

    Returns:
       - gts: (n x L) bool matrix with the n genotypes
    '''

    L = self.number_of_loci
    genotypes = _np.zeros((n, L), bool)
    for i in xrange(genotypes.shape[0]):
        genotypes[i] = self.get_genotype(self.random_clone())
    return genotypes
%}

%feature("autodoc",
"Get a random clone

Returns:
   - clone: index of the random clone
") random_clone;
%pythonappend random_clone {
val = (self._nonempty_clones == val).nonzero()[0][0]
}
%ignore random_clones;
%pythoncode
%{
def random_clones(self, n):
    '''Get random clones
    
    Parameters:
       - n: number of random clones to return
    
    Returns:
       - clones: clone indices
    '''
    return _np.array([self.random_clone() for i in xrange(n)], int)
%}

/* divergence/diversity/fitness distributions and plot */
%ignore get_divergence_histogram;
%ignore get_diversity_histogram;
%ignore get_fitness_histogram;
%pythoncode
%{
def get_fitness_histogram(self, n_sample=1000, **kwargs):
    '''Calculate the fitness histogram of a population sample.

    Parameters:
       - n_sample: number of individuals to sample

    Returns:
       - h: numpy.histogram of fitness in the population
    '''

    fit = [self.get_fitness(self.random_clone()) for i in xrange(n_sample)]
    h = _np.histogram(fit, **kwargs)
    return h
    
    
def plot_fitness_histogram(self, axis=None, n_sample=1000, **kwargs):
    '''Plot a distribution of fitness of a population sample.

    Parameters:
       - axis: an axis to use. A new figure is created by default
       - n_sample: number of individuals to sample
       - kwargs: further optional keyword arguments to matplotlib.pyplot.hist

    Returns:
       - return value of axis.hist(...)
    '''

    import matplotlib.pyplot as plt
    fit = [self.get_fitness(self.random_clone()) for i in xrange(n_sample)]
    
    if axis is None:
        fig = plt.figure()
        axis = fig.add_subplot(111)
        axis.set_title('Fitness histogram')
        axis.set_xlabel('Fitness')
    return axis.hist(fit, **kwargs)
    
    
def get_divergence_histogram(self, bins=10, chunks=None, every=1, n_sample=1000, **kwargs):
    '''Get the divergence histogram restricted to those chunks of the genome.

    Parameters:
       - bins: number or array of bins to be used in the histogram (see also numpy.histogram)
       - chunks: restrict analysis to some chunk in the genome. It must be an n x 2 matrix with
                 the initial and (final+1) positions of the chunks
       - every: restrict analysis to every X positions. For instance, if every third site is neutral,
                this argument can be used to only look at those neutral sites
       - n_sample: number of individuals to sample
       - kwargs: further optional keyword arguments to numpy.histogram

    Returns:
       - h: numpy.histogram of divergence in the population
    '''

    # Check chunks
    if chunks is not None:
        chunks = _np.asarray(chunks)
        if (_np.rank(chunks) != 2) or (chunks.shape[1] != 2):
            raise ValueError('Please input an N x 2 matrix with the chunks initial and (final+1) positions')
    
    # Get the random genotypes
    genotypes = self.random_genomes(n_sample)
    
    # Restrict to the chunks
    if chunks is not None:
        ind = _np.zeros(genotypes.shape[1], bool)
        for chunk in chunks:
            inde = _np.arange(chunk[1] - chunk[0])
            inde = inde[(inde % every) == 0] + chunk[0]
            ind[inde] = True
        genotypes = genotypes[:,ind]
    
    # Calculate divergence
    div = genotypes.sum(axis=1)
    
    # Calculate histogram
    return _np.histogram(div, bins=bins, **kwargs)
    
    
def plot_divergence_histogram(self, axis=None, n_sample=1000, **kwargs):
    '''Plot the divergence histogram of a population sample.

    Parameters:
       - axis: an axis to use. A new figure is created by default
       - n_sample: number of individuals to sample
       - kwargs: further optional keyword arguments to matplotlib.pyplot.hist
    
    Returns:
       - return value of axis.hist(...)
    '''

    import matplotlib.pyplot as plt
    genotypes = self.random_genomes(n_sample)
    div = genotypes.sum(axis=1)
     
    if axis is None:
        fig = plt.figure()
        axis = fig.add_subplot(111)
        axis.set_title('Divergence histogram')
        axis.set_xlabel('Divergence')
    
    if 'bins' not in kwargs:
        kwargs['bins'] = _np.arange(10) * max(1, (div.max() + 1 - div.min()) / 10) + div.min()
    return axis.hist(div, **kwargs)
    
    
def get_diversity_histogram(self, bins=10, chunks=None, every=1, n_sample=1000, **kwargs):
    '''Get the diversity histogram restricted to those chunks of the genome.

    Parameters:
       - bins: number or array of bins to be used in the histogram (see also numpy.histogram)
       - chunks: restrict analysis to some chunk in the genome. It must be an n x 2 matrix with
                 the initial and (final+1) positions of the chunks
       - every: restrict analysis to every X positions. For instance, if every third site is neutral,
                this argument can be used to only look at those neutral sites
       - n_sample: number of individuals to sample
       - kwargs: further optional keyword arguments to numpy.histogram

    Returns:
       - h: numpy.histogram of diversity in the population
    '''

    # Check chunks
    if chunks is not None:
        chunks = _np.asarray(chunks)
        if (_np.rank(chunks) != 2) or (chunks.shape[1] != 2):
            raise ValueError('Please input an N x 2 matrix with the chunks initial and (final+1) positions')
    
    # Get the random genotypes
    genotypes = self.random_genomes(2 * n_sample)
    
    # Restrict to the chunks
    if chunks is not None:
        ind = _np.zeros(genotypes.shape[1], bool)
        for chunk in chunks:
            inde = _np.arange(chunk[1] - chunk[0])
            inde = inde[(inde % every) == 0] + chunk[0]
            ind[inde] = True
        genotypes = genotypes[:,ind]
    
    # Calculate diversity
    genotypes1 = genotypes[:genotypes.shape[0] / 2]
    genotypes2 = genotypes[-genotypes1.shape[0]:]
    div = (genotypes1 != genotypes2).sum(axis=1)
    
    # Calculate histogram
    return _np.histogram(div, bins=bins, **kwargs)


def plot_diversity_histogram(self, axis=None, n_sample=1000, **kwargs):
    '''Plot the diversity histogram of a population sample.

    Parameters:
       - axis: an axis to use. A new figure is created by default
       - n_sample: number of individuals to sample
       - kwargs: further optional keyword arguments to matplotlib.pyplot.hist
    
    Returns:
       - return value of axis.hist(...)
    '''

    import matplotlib.pyplot as plt
    genotypes1 = self.random_genomes(n_sample)
    genotypes2 = self.random_genomes(n_sample)
    div = (genotypes1 != genotypes2).sum(axis=1)
    
    if axis is None:
        fig = plt.figure()
        axis = fig.add_subplot(111)
        axis.set_title('Diversity histogram')
        axis.set_xlabel('Diversity')
    
    if 'bins' not in kwargs:
        kwargs['bins'] = _np.arange(10) * max(1, (div.max() + 1 - div.min()) / 10) + div.min()
    return axis.hist(div, **kwargs)
%}


/* add a tree to the mlg at the selected locus (used to load populations from files) */
%exception _set_tree_in_genealogy {
        try {
                $action
        } catch (int err) {
                PyErr_SetString(PyExc_ValueError,"Locus not found among the tracked ones.");
                SWIG_fail;
        }
}
void _set_tree_in_genealogy(int locus, rooted_tree tree) {
        multi_locus_genealogy* own_genealogy = &($self->genealogy);
        vector<int>::iterator index;
        index = std::find((own_genealogy->loci).begin(), (own_genealogy->loci).end(), locus);
        if(index == (own_genealogy->loci).end()) {
                throw (int)RT_LOCUSNOTFOUND;
        } else
                (own_genealogy->trees)[(int)(index - (own_genealogy->loci).begin())] = tree;
}
/* add a newGeneration to the mlg at the selected locus (used to load populations from files) */
%exception _set_newGeneration_in_genealogy {
        try {
                $action
        } catch (int err) {
                PyErr_SetString(PyExc_ValueError,"Locus not found among the tracked ones.");
                SWIG_fail;
        }
}
void _set_newGeneration_in_genealogy(int locus, vector<node_t> newGeneration) {
        multi_locus_genealogy* own_genealogy = &($self->genealogy);
        int i = 0;
        for(vector< vector<node_t> >::iterator it=own_genealogy->newGenerations.begin(); it != own_genealogy->newGenerations.end(); it++, i++) {
                if ((own_genealogy->loci)[i] == locus) {
                        (own_genealogy->newGenerations)[i] = newGeneration;
                        return;
                }
        }
        throw (int)RT_LOCUSNOTFOUND;
}

} /* extend haploid_highd */

%{
const int haploid_highd_L_get(haploid_highd *h) {
  return (const int) h->get_number_of_loci();
}
const int haploid_highd_number_of_loci_get(haploid_highd *h) {
  return (const int) h->get_number_of_loci();
}

const int haploid_highd_N_get(haploid_highd *h) {
  return (const int) h->get_population_size();
}
const int haploid_highd_population_size_get(haploid_highd *h) {
  return (const int) h->get_population_size();
}

int haploid_highd_generation_get(haploid_highd *h) {
  return (const int) h->get_generation();
}
void haploid_highd_generation_set(haploid_highd *h, int g) {
  h->set_generation(g);
}

const int haploid_highd_number_of_clones_get(haploid_highd *h) {
  return (const int) h->get_number_of_clones();
}

const int haploid_highd_number_of_traits_get(haploid_highd *h) {
  return (const int) h->get_number_of_traits();
}

const double haploid_highd_max_fitness_get(haploid_highd *h) {
  return (const double) h->get_max_fitness();
}

const double haploid_highd_participation_ratio_get(haploid_highd *h) {
  return (const double) h->get_participation_ratio();
}

const bool haploid_highd_all_polymorphic_get(haploid_highd *h) {
  return (const bool) h->is_all_polymorphic();
}
%}

/* load haploid_highd from file */
%pythoncode
%{
def load_haploid_highd(filename, gen_loci=[], include_genealogy=False):
    '''Load a population from a compressed pickle file

    Parameters:
       - filename: the path of the pickle file
       - gen_loci: start tracking these loci in the population
       - include_genealogy: load the old genealogy if present
    '''

    try:
        import cPickle as pickle
    except:
        import pickle

    # Try the compressed format first
    try:
        import bz2
        with bz2.BZ2File(filename, 'rb') as f:
            pop_dict = pickle.load(f)
    # Fallback on uncompressed
    except:
        with open(filename, 'rb') as f:
            pop_dict = pickle.load(f)
        

    pop = haploid_highd(pop_dict['L'],
                        all_polymorphic=pop_dict['all_polymorphic'],
                        number_of_traits=len(pop_dict['traits_additive']))
    pop.carrying_capacity = pop_dict['N']
    if pop.all_polymorphic == False:
        pop.mutation_rate = pop_dict['mu']
    pop.crossover_rate = pop_dict['crossover_rate']
    pop.outcrossing_rate = pop_dict['outcrossing_rate']
    pop.circular = pop_dict['circular']

    pop.recombination_model = pop_dict['recombination_model']
    for i in range(pop.number_of_traits):
        pop.set_trait_additive(pop_dict['traits_additive'][i], i)
        for (value, loci) in pop_dict['traits_epistasis'][i]:
            pop.add_trait_coefficient(value, loci, i)
    
    pop.trait_weights=pop_dict['trait_weights']

    # Load the genealogy and track new loci if wished
    if include_genealogy and 'trees' in pop_dict:
        old_loci = pop_dict['trees'].keys()
    else:
        old_loci = []
    all_loci = list(set(list(old_loci) + list(gen_loci)))

    if len(all_loci):
        pop.track_locus_genealogy(all_loci)

    pop.generation = pop_dict['generation']
            
    # Initialize the population
    # Note: if the tree is recovered from the past, we insert empty clones to
    # keep the labels of the leaves 
    if include_genealogy and 'trees' in pop_dict:
        import numpy as np
        _nonempty_clones = pop_dict['_nonempty_clones']
        _maxclone = _nonempty_clones.max()
        genotypes = np.zeros((_maxclone + 1, pop.L), bool)
        clone_sizes = np.zeros(_maxclone + 1, int)
        genotypes[_nonempty_clones] = pop_dict['genotypes']
        clone_sizes[_nonempty_clones] = pop_dict['clone_sizes']
        pop.set_genotypes_and_ancestral_state(genotypes, 
                                              clone_sizes, 
                                              pop_dict['ancestral'])
        for (locus, tree_s) in pop_dict['trees'].iteritems():
            tree = rooted_tree()
            tree.read_newick(tree_s)
            pop._set_tree_in_genealogy(locus, tree)

        # Make nodes out of non-SWIG objects for newGenerations
        def deserialize_leaf(serial):
            leaf = tree_node()
            for key in ['clone_size', 'crossover', 'fitness', 'number_of_offspring']:
                setattr(leaf, key, serial[key])
            leaf.own_key = tree_key(*serial['own_key'])
            leaf.parent_node = tree_key(*serial['parent_node'])
            return leaf

        for i, locus in enumerate(old_loci):
            pop._set_newGeneration_in_genealogy(locus, map(deserialize_leaf, pop_dict['_newGenerations'][i]))


    else:
        pop.set_genotypes_and_ancestral_state(pop_dict['genotypes'], 
                                              pop_dict['clone_sizes'], 
                                              pop_dict['ancestral'])
    

    return pop
%}
/*****************************************************************************/
