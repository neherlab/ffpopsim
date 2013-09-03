#ifndef MULTILOCUSGENEALOGY_H
#define MULTILOCUSGENEALOGY_H

#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <list>
#include <boost/dynamic_bitset.hpp>
//#include <gsl/gsl_histogram.h>

#define HCF_MEMERR -131545
#define HCF_BADARG -131546
#define HCF_VERBOSE 0
#define WORDLENGTH 28 	//length used to chop bitsets into words

#ifndef rooted_tree_H_
#define rooted_tree_H_
#define RT_VERBOSE 0
#define RT_VERYLARGE 10000000
#define RT_CHILDNOTFOUND -35343
#define RT_NODENOTFOUND -35765
#define RT_LOCUSNOTFOUND -35762
#define RT_FITNESS_MISSING -35722
#define RT_CROSSOVER_MISSING -35721
#define RT_SEGMENT_MISSING -35720
#define RT_ERROR_PARSING 1


using namespace std;

/*
 *	@brief a class that implements a rooted tree to store genealogies
 *
 *	Nodes and edges are stored as maps with a key that holds the age (rather the time) when the node lived
 *	and the index in the population at that time. The nodes themselves are sufficient to reconstruct the tree
 *	since they contain keys of parents and children
 */

struct tree_key_t
{

    int location;
    int age;
    int index;



    bool operator == (const tree_key_t &other)  {return (age == other.age) && (index == other.index) && (location == other.location);}
    bool operator != (const tree_key_t &other)  {return (age != other.age) || (index != other.index) || (location != other.location);}
    bool operator < (const tree_key_t &other) const
    {
                if(age < other.age) return true;
                else if (age > other.age) return false;
                else
                {
                    if (location < other.location) return true;
                    else if ( location > other.location) return false;
                    else {return (index < other.index); }
                }

    }
    bool operator > (const tree_key_t &other) const
    {
                if(age > other.age) return true;
                else if (age < other.age) return false;
                else
                {
                    if (location > other.location) return true;
                    else if ( location < other.location) return false;
                    return (index < other.index);
                }

    }
    tree_key_t(int index = 0, int age = 0, int location = 0) : index (index), age (age), location (location)  {};

};

std::ostream& operator<< ( std::ostream& os, const tree_key_t& key );

struct step_t {
    int pos;
    int step;
    bool operator<(const step_t &other) const {
                if(pos < other.pos) return true;
                else return false;
        }
    bool operator>(const step_t &other) const {
                if(pos > other.pos) return true;
                else return false;
        }
    bool operator==(const step_t &other) const {
                if(pos == other.pos) return true;
                else return false;
        }
        step_t(int pos=0, int step=0) : pos(pos), step(step) {};
};


struct node_t {
    tree_key_t parent_node;
    list < tree_key_t > child_edges;
    double fitness;
    tree_key_t own_key;
    vector <step_t> weight_distribution;
    vector < double > traits;
    boost::dynamic_bitset<> allele_freqs;
    int number_of_offspring;
    int clone_size;
    int crossover[2];
};


struct edge_t {
    tree_key_t parent_node;
    tree_key_t own_key;
    int segment[2];
    int length;
    int number_of_offspring;
};

struct poly_t {
    int birth;
    int sweep_time;
    double effect;
    double fitness;
    double fitness_variance;
    poly_t(int b=0, int age=0, double e=0, double f=0, double fvar=0) :
                birth(b), sweep_time(age), effect(e), fitness(f), fitness_variance(fvar) {};
};


class rooted_tree {
public:
    map < tree_key_t , edge_t > edges;
    map < tree_key_t , node_t > nodes;
    vector <tree_key_t> leafs;
    tree_key_t root;
    tree_key_t MRCA;

    rooted_tree();
    virtual ~rooted_tree();
    void reset();
    void add_generation(vector <node_t> &new_generation, double mean_fitness);
    int add_terminal_node(node_t &newNode);
    tree_key_t erase_edge_node(tree_key_t to_be_erased);
    tree_key_t bridge_edge_node(tree_key_t to_be_bridged);
    int external_branch_length();
    int total_branch_length();
    int ancestors_at_age(int age, tree_key_t subtree_root, vector <tree_key_t> &ancestors);
    int update_leaf_to_root(tree_key_t leaf);
    void update_tree();
    int calc_weight_distribution(tree_key_t subtree_root);
    void SFS(gsl_histogram *sfs);

    tree_key_t get_MRCA(){return MRCA;};
    int erase_child(map <tree_key_t,node_t>::iterator Pnode, tree_key_t to_be_erased);
    int delete_extra_children(tree_key_t subtree_root);
    int delete_one_child_nodes(tree_key_t subtree_root);
    bool check_node(tree_key_t node);
    int check_tree_integrity();
    void clear_tree();

        // print tree or subtrees
    string print_newick();
    string subtree_newick(tree_key_t root);
    string print_genotypes();
    string genotypes_newick(tree_key_t root);
    string print_weight_distribution(tree_key_t node_key);
	int read_newick(string newick_string);

private:
    int parse_label(std::string label, int *index, int *clone_size, int *branch_length);
    int parse_subtree(tree_key_t &parent_key, std::string &tree_s);
        // construct subtrees
    int construct_subtree(vector <tree_key_t> subtree_leafs, rooted_tree &other);


};

#endif /* rooted_tree_H_ */

/*
 * @brief short wrapper class that handles trees at different places in the genome
 *
 * the class contains a vector of rooted_tree instances that hold the genealogy
 *  in different places. In addition, there is a rooted_tree called subtree
 *  that is used on demand
 *
 *  Created on: Oct 14, 2012
 *      Author: richard
 */

class multi_locus_genealogy_parent {
public:
    vector <int> loci;				//vector of loci (positions on a genome) whose genealogy is to be tracked
    vector <rooted_tree> trees;                     //vector of rooted trees (one per locus)
    vector < vector < node_t > > newGenerations;	//used by the evolving class to store the new generation

    multi_locus_genealogy_parent();
    virtual ~multi_locus_genealogy_parent();
    void track_locus(int new_locus);
    void reset(){loci.clear(); trees.clear(); newGenerations.clear();}
    void reset_but_loci(){for(unsigned int i=0; i<loci.size(); i++){trees[i].reset();newGenerations[i].clear();}}
    virtual void add_generation(double mean_fitness){};
    virtual void add_generation(double baseline, vector < vector < node_t > >  new_Generations){};
    virtual int extend_storage(int n);
};

class multi_locus_genealogy : public multi_locus_genealogy_parent{

public:


   void add_generation(double baseline);

};


class multi_locus_genealogy_2 : public multi_locus_genealogy_parent{

public:
    multi_locus_genealogy_2();
    multi_locus_genealogy_2(vector < vector < node_t > >  new_Generations);
    void add_generation(double baseline, vector < vector < node_t > >  new_Generations);

};



#endif /* MULTILOCUSGENEALOGY_H_ */
