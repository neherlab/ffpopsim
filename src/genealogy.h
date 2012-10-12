/*
 * genealogy.h
 *
 *  Created on: Sep 21, 2012
 *      Author: richard
 */

#ifndef GENEALOGY_H_
#define GENEALOGY_H_
#define GEN_VERBOSE 0
#define GEN_VERYLARGE 10000000
#define GEN_CHILDNOTFOUND -35343

#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <list>
#include <gsl/gsl_histogram.h>


using namespace std;

struct key_t {
	int index;
	int age;
	bool operator==(const key_t &other)  {return (age == other.age) && (index == other.index);}
	bool operator!=(const key_t &other)  {return (age != other.age) || (index != other.index);}
	bool operator<(const key_t &other) const {
                if(age < other.age) return true;
                else if (age > other.age) return false;
                else { return (index<other.index); }
        }
	bool operator>(const key_t &other) const {
                if(age > other.age) return true;
                else if (age < other.age) return false;
                else { return (index>other.index); }
        }
};

struct node_t{
	key_t parent_node;
	list < key_t > child_edges;
	double fitness;
	key_t own_key;
	int number_of_offspring;
	int clone_size;
	int crossover[2];
};

struct edge_t{
	key_t parent_node;
	key_t own_key;
	int segment[2];
	int length;
	int number_of_offspring;
};

class genealogy{
	map < key_t , edge_t > edges;
	map < key_t , node_t > nodes;
	map < key_t , edge_t > subtree_edges;
	map < key_t , node_t > subtree_nodes;
	vector <key_t> leafs;
	key_t root;
	key_t MRCA;
	key_t subtree_MRCA;

public:
	genealogy();
	~genealogy();

	void reset();
	void add_generation(vector <node_t> &new_generation, double mean_fitness);
	key_t erase_edge_node(key_t to_be_erased);
	key_t bridge_edge_node(key_t to_be_bridged);
	int external_branch_length();
	int total_branch_length();
	void clear_tree();
	void update_leaf_to_root(key_t leaf_key);
	void update_tree();
	void SFS(gsl_histogram *sfs);
	key_t get_MRCA(){return MRCA;};
	bool check_node(key_t node);
	int erase_child(map <key_t,node_t>::iterator Pnode, key_t to_be_erased);
	int construct_subtree(vector <key_t> subtree_leafs);
	void delete_extra_children_in_subtree(key_t);
	string print_newick();
	string subtree_newick(key_t root);
};

#endif /* GENEALOGY_H_ */
