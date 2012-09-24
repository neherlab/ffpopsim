/*
 * genealogy.h
 *
 *  Created on: Sep 21, 2012
 *      Author: richard
 */

#ifndef GENEALOGY_H_
#define GENEALOGY_H_
#define GEN_VERBOSE 1
#include <map>
#include <vector>
#include <iostream>
#include <list>
#include <gsl/gsl_histogram.h>

using namespace std;

struct key_t {
	int index;
	int age;
	bool operator==(const key_t &other)  {return (age == other.age) && (index == other.index);}
	bool operator!=(const key_t &other)  {return (age != other.age) || (index != other.index);}
	bool operator<(const key_t &other)  {
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
	key_t parent_edge;
	list < key_t > child_edges;
	double fitness;
	key_t own_key;
	int number_of_offspring;
	int crossover[2];
};

struct edge_t{
	key_t parent_node;
	key_t child_node;
	key_t own_key;
	int segment[2];
	int length;
	int number_of_offspring;
};

class genealogy{
	map < key_t , edge_t > edges;
	map < key_t , node_t > nodes;
	vector <key_t> leafs;

public:
	genealogy();
	~genealogy();

	void add_generation(vector <node_t> &new_generation);
	key_t erase_edge_node(key_t to_be_erased);
	key_t bridge_edge_node(key_t to_be_bridged);
	int external_branch_length();
	int total_branch_length();
	void clear_tree();
};

#endif /* GENEALOGY_H_ */
