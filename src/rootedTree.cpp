/*
 * tree.cpp
 *
 *  Created on: Oct 14, 2012
 *      Author: richard
 *
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

#include "ffpopsim_highd.h"
/*
 * this overloads the ostream operator to output keys of nodes and edges
 */
std::ostream& operator<< ( std::ostream& os, const tree_key_t& key )
  {
    os <<"age: "<< key.age << " index "
      << key.index;
    return os;
  }


rooted_tree::rooted_tree() {
	reset();
}
rooted_tree::~rooted_tree() {
}

/*
 * @brief deletes everything and reinserts a root and an MRCA
 */
void rooted_tree::reset(){
	node_t root_node;
	node_t mrca_node;
 	edge_t to_root;

 	nodes.clear();
 	edges.clear();
 	leafs.clear();

 	//the root node will never be touched, the MRCA moves up with the tree
	root.age=-3;
	root.index=0;
	MRCA.age=-2;
	MRCA.index=0;

	//set up a trivial root node with the only child MRCA
	root_node.own_key = root;
	root_node.parent_node = root;
	root_node.clone_size=-1;
	root_node.number_of_offspring=-1;
	root_node.crossover[0]=0;
	root_node.crossover[1]=RT_VERYLARGE;
	root_node.child_edges.clear();
	root_node.child_edges.push_back(MRCA);

	//set up a trivial MRCA node to which any genotypes are to be attached upon initialization
	mrca_node.own_key = MRCA;
	mrca_node.parent_node = root;
	mrca_node.clone_size=-1;
	mrca_node.number_of_offspring=-1;
	mrca_node.crossover[0]=0;
	mrca_node.crossover[1]=RT_VERYLARGE;
	mrca_node.child_edges.clear();

	//produce an edge that links the MRCA to root.
	to_root.own_key = MRCA;
	to_root.length=mrca_node.own_key.age-root_node.own_key.age;
	to_root.segment[0]=0;
	to_root.segment[1]=RT_VERYLARGE;
	to_root.number_of_offspring=-1;
	to_root.parent_node = root;

	//as of now, the only leaf is the MRCA. root, MRCA and the edge are added to respective maps
	leafs.push_back(MRCA);
	nodes.insert(pair<tree_key_t,node_t>(root,root_node));
	nodes.insert(pair<tree_key_t,node_t>(MRCA,mrca_node));
	edges.insert(pair<tree_key_t,edge_t>(to_root.own_key,to_root));
}

/*
 * @brief takes a vector of nodes from the evolution class and knits it into the genealogy.
 * @params new_generation
 * @params double mean_fitness is the baseline with respect to which node fitness is measured.
 * TODO nothing is done with the mean_fitness argument
 */
void rooted_tree::add_generation(vector <node_t> &new_generation, double mean_fitness){
	if (RT_VERBOSE){
		cerr <<"rooted_tree::add_generation(). Number of leafs to add: "<<new_generation.size()<<endl;
	}

	vector<node_t>::iterator new_leaf = new_generation.begin();
	vector <tree_key_t> new_leafs;
	new_leafs.reserve(new_generation.size());
	//add new leafs to a temporary vector and to the tree
	for (; new_leaf!=new_generation.end(); new_leaf++){
		if (new_leaf->clone_size>0){	//restrict to those with positive clone sizes, all others are dummies
			add_terminal_node(*new_leaf);
			new_leafs.push_back(new_leaf->own_key);
		}
	}

	if (RT_VERBOSE){
		cerr <<"rooted_tree::add_generation(). added leafs... erase dead ends"<<endl;
		cerr <<"rooted_tree::add_generation(). rooted_tree size: "<<edges.size()<<" edges, "<<nodes.size()<<" nodes "<<endl;
	}

	//go over the previous leafs and delete those that leave no offspring, bridge those with one
	map <tree_key_t,edge_t>::iterator edge_pos = edges.end();
	map <tree_key_t,node_t>::iterator node_pos = nodes.end();
	tree_key_t parent_key;
	for(vector <tree_key_t>::iterator old_leaf_key=leafs.begin(); old_leaf_key!=leafs.end(); old_leaf_key++){
		node_pos = nodes.find(*old_leaf_key);
		if (node_pos!=nodes.end()){
			if (node_pos->second.child_edges.size()==0){			//no offspring -> delete
				parent_key = erase_edge_node(*old_leaf_key);
				while (nodes[parent_key].child_edges.size()==0){ 	//continue until offspring found
					parent_key = erase_edge_node(parent_key);
				}
				//bridge upstream nodes that have exactly one offspring until either multiple offsprings
				//are encountered or the parent is the root.
				while (nodes[parent_key].child_edges.size()==1 and parent_key!=root){
					parent_key = bridge_edge_node(parent_key);
				}
			}else if (node_pos->second.child_edges.size()==1){		//bridge nodes that have exactly one offspring
				parent_key = bridge_edge_node(*old_leaf_key);
				//bridge upstream nodes that have exactly one offspring until either multiple offsprings
				//are encountered or the parent is the root.
				while (nodes[parent_key].child_edges.size()==1 and root!=parent_key){
					parent_key = bridge_edge_node(parent_key);
				}
			}
		}else{
			cerr <<"rooted_tree::add_generation(). did not find old leaf "<<*old_leaf_key<<endl;
		}
	}

	if (RT_VERBOSE){
		cerr <<"rooted_tree::add_generation(). done "<<endl;
	}

	//copy the temporary leafs to the actual leafs
	leafs = new_leafs;
	if (RT_VERBOSE){
		cerr <<"rooted_tree::add_generation(). done "<<endl;
		cerr <<"rooted_tree::add_generation(). rooted_tree size: "<<edges.size()<<" edges, "<<nodes.size()<<" nodes "<<endl;
	}

	//make sure all nodes and edges are up-to-date in terms of offspring
	//number and other attributes (this is not necessary, could be called when needed)
	update_tree();
	return;
}

/*
 * @brief basic operation that attached a node to the tree
 * @params node_t reference to the newly added node
 */
int rooted_tree::add_terminal_node(node_t &new_node){
	edge_t new_edge;
	tree_key_t new_key;
	new_key = new_node.own_key;
	new_node.child_edges.clear();		//no kids
	//produce an edge that links the new leaf to its parent
	new_edge.own_key=new_key;			//reset the associated edge (has key of child node)
	new_edge.parent_node=new_node.parent_node;
	new_edge.number_of_offspring=1;
	new_edge.segment[0]=new_node.crossover[0];
	new_edge.segment[1]=new_node.crossover[1];
	new_edge.length=new_key.age-new_node.parent_node.age;
	nodes[new_node.parent_node].child_edges.push_back(new_key);			//add node as child of parent
	edges.insert(pair<tree_key_t,edge_t>(new_key, new_edge));			//insert node and edge
	nodes.insert(pair<tree_key_t,node_t>(new_key, new_node));
	return 0;
}


/*
 * @brief basic operation that deletes a node that leafs no offspring
 * @params tree_key_t the key of the one to be erased
 */
tree_key_t rooted_tree::erase_edge_node(tree_key_t to_be_erased){
	if (RT_VERBOSE) {
		cerr <<"rooted_tree::erase_edge_node(). ..."<<to_be_erased<<endl;
	}
	//iterators pointing to the elements that are to be erased, found by their key (argument)
	map <tree_key_t,node_t>::iterator Enode = nodes.find(to_be_erased);
	map <tree_key_t,edge_t>::iterator Eedge = edges.find(to_be_erased);

	if (Enode->second.child_edges.size()>0) {
		cerr <<"rooted_tree::erase_edge_node(): attempting to erase non-terminal node"<<endl;
	}

	//find the parents of the node that is going to be erased
	tree_key_t parent_key = Eedge->second.parent_node;
	map <tree_key_t,edge_t>::iterator Pedge = edges.find(parent_key);
	map <tree_key_t,node_t>::iterator Pnode = nodes.find(parent_key);

	//update the parents offspring numbers by subtracting what is lost when erasing the child
	Pnode->second.number_of_offspring-=Eedge->second.number_of_offspring;
	Pedge->second.number_of_offspring-=Eedge->second.number_of_offspring;

	//take out the child and check for good result
	if (erase_child(Pnode, to_be_erased)==RT_CHILDNOTFOUND) {
		cerr <<"rooted_tree::erase_edge_node(): child not found"<<endl;
	}
	//update the register of nodes and edges.
	nodes.erase(to_be_erased);
	edges.erase(to_be_erased);

	if (RT_VERBOSE) {
		cerr <<"rooted_tree::erase_edge_node(). done"<<endl;
	}

	//return the parent key, since the upstream function might want to check
	//whether the parent needs to be erased too
	return parent_key;
}


/*
 * @brief erases a child from the list of children of a node
 * @params <tree_key_t,node_t>::iterator Pnode parent node
 * @params tree_key_t to_be_erased
 */
int rooted_tree::erase_child(map <tree_key_t,node_t>::iterator Pnode, tree_key_t to_be_erased){
	for (list <tree_key_t>::iterator child = Pnode->second.child_edges.begin();child!=Pnode->second.child_edges.end(); child++){
		if (*child == to_be_erased){Pnode->second.child_edges.erase(child); return 0;}
	}
	return RT_CHILDNOTFOUND;
}

/*
 * @brief erases a node that has one offspring. wires that nodes child to that node parent
 * @params tree_key_t of node to be deleted
 */
tree_key_t rooted_tree::bridge_edge_node(tree_key_t to_be_bridged) {
	if (RT_VERBOSE){
		cerr <<"rooted_tree::bridge_edge_node(). ..."<<to_be_bridged<<endl;
	}

	//node and edge to be erased
	map <tree_key_t,node_t>::iterator Enode = nodes.find(to_be_bridged);
	map <tree_key_t,edge_t>::iterator Eedge = edges.find(to_be_bridged);

	//consistency checks. only nodes with exactly one child can be bridged.
	if (Enode->second.child_edges.size()!=1 or to_be_bridged==root){
		cerr <<"rooted_tree::bridge_edge_node(): attempting to bridge branched node or bridge root"<<endl;
	}

	// determine the parent key
	tree_key_t parent_key = Eedge->second.parent_node;

	//find the one and only child of the node to be bridged
	map <tree_key_t,edge_t>::iterator child_edge = edges.find(Enode->second.child_edges.front());
	map <tree_key_t,node_t>::iterator child_node = nodes.find(Enode->second.child_edges.front());
	//rewire: attach the child of the node to be erased to the parent of that node.
	child_edge->second.parent_node = parent_key;
	child_node->second.parent_node = parent_key;

	//update the part of the chromosome transmitted along this edge to the minimum of the two edges
	//cerr  <<child_edge->second.segment[0]<<" "<<Eedge->second.segment[0]<<"  "<<child_edge->second.segment[1]<<" "<<Eedge->second.segment[1]<<endl;
	child_edge->second.segment[0]=(child_edge->second.segment[0]<Eedge->second.segment[0])?(Eedge->second.segment[0]):(child_edge->second.segment[0]);
	child_edge->second.segment[1]=(child_edge->second.segment[1]>Eedge->second.segment[1])?(Eedge->second.segment[1]):(child_edge->second.segment[1]);
	child_edge->second.length+=Eedge->second.length;	//add length of edges
	//cerr  <<child_edge->second.segment[0]<<"  "<<child_edge->second.segment[1]<<endl;

	//find the parent node: add new child, delete old
	map <tree_key_t,node_t>::iterator Pnode = nodes.find(parent_key);
	Pnode->second.child_edges.push_back(child_edge->first);
	if (erase_child(Pnode, to_be_bridged)==RT_CHILDNOTFOUND){
		cerr <<"rooted_tree::bridge_edge_node(). child not found. index "<<to_be_bridged<<endl;
	}
	//erase bidged node from register
	nodes.erase(to_be_bridged);
	edges.erase(to_be_bridged);
	//if the MRCA was bridged, reset the MRCA key to the child key
	if (to_be_bridged == MRCA){MRCA = child_edge->first;}

	if (RT_VERBOSE){
		cerr <<"rooted_tree::bridge_edge_node(). done"<<endl;
	}

	//return the parent key, since the upstream function might want to check
	//whether the parent needs to be erased too
	return parent_key;
}

/*
 * @brief walk over all leafs and rebuild the number of ancestors leaf to root
 */
void rooted_tree::update_tree() {
	clear_tree();
	for (vector <tree_key_t>::iterator leaf=leafs.begin(); leaf!=leafs.end(); leaf++){
		update_leaf_to_root(*leaf);
	}
}

/*
 * @brief start at a leaf and add its population to all nodes in its lineage to the MRCA
 * @params tree_key_t to the leaf at which this is supposed to start
 * Note that before calling this for every leaf, it is assumed that the tree is reset.
 * This could also be done recursively starting from the root.
 */
int rooted_tree::update_leaf_to_root(tree_key_t leaf_key) {
	if (RT_VERBOSE){
		cerr <<"rooted_tree::update_leaf_to_root(). key:"<<leaf_key<<endl;
	}
	//find the leaf and check for its presence
	map <tree_key_t,node_t>::iterator leaf_node = nodes.find(leaf_key);
	map <tree_key_t,edge_t>::iterator leaf_edge = edges.find(leaf_key);
	if (leaf_node == nodes.end() or leaf_edge == edges.end()){
		cerr <<"rooted_tree::update_leaf_to_root(). leaf not found"<<endl;
		return RT_NODENOTFOUND;
	}

	//count the number of offspring along the lineage in the variable increment
	int increment = leaf_node->second.clone_size;
	leaf_edge->second.number_of_offspring = increment;
	map <tree_key_t,node_t>::iterator parent_node = nodes.find(leaf_edge->second.parent_node);
	map <tree_key_t,edge_t>::iterator parent_edge = edges.find(leaf_edge->second.parent_node);
	//walk up the tree and add this leafs contribution to every node and to the parent is the root.
	while (root != parent_node->first){
		parent_node->second.number_of_offspring+=increment;
		parent_edge->second.number_of_offspring+=increment;

		//TODO: this seems ugly
		leaf_node = parent_node;
		leaf_edge = parent_edge;
		//find parents.
		parent_node = nodes.find(leaf_edge->second.parent_node);
		parent_edge = edges.find(leaf_edge->second.parent_node);
		if (parent_node==nodes.end()){
			cerr <<"rooted_tree::update_leaf_to_root(). key:"<<leaf_key<<endl;
			cerr <<"rooted_tree::update_leaf_to_root(): key not found: "<<leaf_edge->second.parent_node<<" root: "<<root<<endl;
			break;
		}
	}
	if (RT_VERBOSE){
		cerr <<"rooted_tree::update_leaf_to_root(). done"<<endl;
		cerr <<"rooted_tree::update_leaf_to_root(): total of "<< nodes.find(MRCA)->second.number_of_offspring<<" offspring "<<increment<<endl;
	}
	return 0;
}

/*
 * @brief loop over all edges and calculate the site frequency spectrum of all derived neutral mutations
 * @params gsl_histogram* the histogram is accumulated and assumed to have bins between 0 and 1
 * NOTE that this assumes number_of_offspring counts on the tree are up-to-date
 */
void rooted_tree::SFS(gsl_histogram *sfs) {
	map <tree_key_t,edge_t>::iterator edge = edges.begin();
	int total_pop = nodes[MRCA].number_of_offspring;
	for (; edge!=edges.end(); edge++)
		gsl_histogram_accumulate(sfs, 1.0*edge->second.number_of_offspring/total_pop, edge->second.length);
}

/*
 * @brief calculate external branch length
 */
int rooted_tree::external_branch_length(){
	map <tree_key_t,edge_t>::iterator edge = edges.begin();
	int branchlength = 0;
	for (vector <tree_key_t>::iterator leaf=leafs.begin(); leaf!=leafs.end();leaf++)
		branchlength+=edges[*leaf].length;
	return branchlength;
}

/*
 * @brief calculate total branch length
 */
int rooted_tree::total_branch_length() {
	map <tree_key_t,edge_t>::iterator edge = edges.begin();
	int branchlength = 0;
	for (; edge!=edges.end(); edge++)
		branchlength+=edge->second.length;
	return branchlength;
}

/*
 * @brief delete all number of offspring values and set the leafs to clone size
 */
void rooted_tree::clear_tree() {
	if (RT_VERBOSE)
		cerr <<"rooted_tree::clear_tree()..."<<endl;

	//loop over all nodes and edges and reset.
	map <tree_key_t,node_t>::iterator node = nodes.begin();
	map <tree_key_t,edge_t>::iterator edge = edges.begin();

	for (; node!=nodes.end(); node++)
		node->second.number_of_offspring=0;
	for (; edge!=edges.end(); edge++)
		edge->second.number_of_offspring=0;

	//set the leaf number_of_offspring to the clone size (which is never erased as it is a node attribute)
	for (vector<tree_key_t>::iterator leaf=leafs.begin(); leaf!=leafs.end(); leaf++) {
		node = nodes.find(*leaf);
		if (node==nodes.end()) {
			cerr <<"rooted_tree::clear_tree(). key of leaf not found: "<<*leaf<<endl;
			break;
		}
		node->second.number_of_offspring=node->second.clone_size;
	}
	if (RT_VERBOSE) {
		cerr <<"rooted_tree::clear_tree(). done"<<endl;
	}
}

/*
 * @brief return tree in newick format as string
 */
string rooted_tree::print_newick() {
	return subtree_newick(MRCA)+";";
}

/*
 * @brief return newick string of a subtree.
 * @params tree_key_t root of the subtree that is to be returned.
 */
string rooted_tree::subtree_newick(tree_key_t root){
	stringstream tree_str;
	map <tree_key_t,node_t>::iterator root_node = nodes.find(root);
	map <tree_key_t,edge_t>::iterator edge = edges.find(root);
	if (root_node->second.child_edges.size()>0){
		list <tree_key_t>::iterator child = root_node->second.child_edges.begin();
		tree_str.str();
		tree_str <<"("<< subtree_newick(*child);
		child++;
		for (;child!=root_node->second.child_edges.end(); child++){
			tree_str<<","+subtree_newick(*child);
		}
		tree_str<<")";
	}
	tree_str<<root.index<<'_'<<root_node->second.clone_size<<":"<<edge->second.length;
	//tree_str<<root.index<<'_'<<root.age<<":"<<edge->second.length;
	return tree_str.str();
}

/*
 * @brief returns true of the tree_key_t node_key is part of the tree, false otherwise
 */
bool rooted_tree::check_node(tree_key_t node_key){
	map <tree_key_t,node_t>::iterator node = nodes.find(node_key);
	if (node == nodes.end()){
		return false;
	}else { return true;}
}

/*
 * @brief prune a tree such that only a subset of its nodes remain
 * @params vector <tree_key_t> subtree_leafs the leaf to retain
 * @params rooted_tree reference to tree that is to be pruned
 *
 * note that instead of pruning, this functions builds up a new tree from
 * scratch by walking through the super tree
 */
int rooted_tree::construct_subtree(vector <tree_key_t> subtree_leafs, rooted_tree &other){
	if (RT_VERBOSE){
		cerr <<"rooted_tree::construct_subtree()..."<<endl;
	}
	//the MRCA needs to go since we will rediscover it later. the tree does contain a root
	other.reset();
        other.nodes.erase(other.MRCA);
        other.edges.erase(other.MRCA);
        other.leafs.clear();

	//add all new leafs and make a set of leaves to be added next round. (set has unique elements)
	set <tree_key_t> new_nodes;
	new_nodes.clear();
	map <tree_key_t,node_t>::iterator node;
	map <tree_key_t,edge_t>::iterator edge;
	for (vector <tree_key_t>::iterator leaf=subtree_leafs.begin(); leaf!=subtree_leafs.end(); leaf++){
		if (check_node(*leaf)) {
			node = nodes.find(*leaf);
			other.nodes.insert(*node);
			other.edges.insert(*edges.find(*leaf));
			other.leafs.push_back(*leaf);
			new_nodes.insert(node->second.parent_node);
		}else{
			cerr <<"rooted_tree::construct_subtree(). leaf does not exist"<<endl;
			return RT_NODENOTFOUND;
		}

	}
	if (RT_VERBOSE){
		cerr <<"rooted_tree::construct_subtree(). added leafs"<<endl;
	}


	other.MRCA = other.root; //set MRCA to some definite value, doesn't matter
	//repeatedly add parents of all nodes until nothing else is added, or one node
	//is added and that one is the MRCA (set within the loop)
	while((new_nodes.size() > 1 or other.MRCA != *new_nodes.begin()) and new_nodes.size() > 0) {
		set <tree_key_t> temp = new_nodes;
		new_nodes.clear();
		for (set <tree_key_t>::iterator node_key=temp.begin(); node_key!=temp.end(); node_key++){
			if (check_node(*node_key)){
				node = nodes.find(*node_key);
				other.nodes.insert(*node);
				other.edges.insert(*edges.find(*node_key));
				if (node->second.parent_node!=root) {
					new_nodes.insert(node->second.parent_node);
				} else {	//of the parent is the root, that node is the MRCA
					other.MRCA = *node_key;
				}
			}
			else{
				cerr <<"rooted_tree::construct_subtree(). internal node did not exist: age "<<*node_key<<endl;
			}
		}
	}
	if (RT_VERBOSE){
		cerr <<"rooted_tree::construct_subtree(). added internal nodes"<<endl;
	}

	//wire the root and the MRCA together
	if (check_node(other.MRCA)){
		node = other.nodes.insert(*nodes.find(MRCA)).first;
		edge = other.edges.insert(*edges.find(MRCA)).first;
		edge->second.length = other.MRCA.age - other.root.age;
		node->second.parent_node = other.root;
		edge->second.parent_node = other.root;
		node = other.nodes.find(other.root);
		node->second.child_edges.clear();
		node->second.child_edges.push_back(other.MRCA);
	}

	//we have copied nodes from the superTree, whcih probably have more children than necessary
	//delete those
	other.delete_extra_children(other.MRCA);
	//bridge all nodes with a single child
	other.delete_one_child_nodes(other.MRCA);
	//
	other.update_tree();

	if (RT_VERBOSE){
		cerr <<"genealogy::construct_subtree(). done"<<endl;
	}

	return 0;
}

/*
 * @brief function that recursively deletes all children that no longer exist
 * @params takes as an argument the key of the root of the subtree that is to be handled
 */
int rooted_tree::delete_extra_children(tree_key_t subtree_root){
	if (RT_VERBOSE){
		cerr <<"rooted_tree::delete_extra_children(). node age "<<subtree_root<<endl;
	}

	map <tree_key_t,node_t>::iterator node = nodes.find(subtree_root);
	if (node == nodes.end()){
		cerr <<"rooted_tree::delete_extra_children(): subtree root not found! age: "<<subtree_root<<endl;
		return RT_NODENOTFOUND;
	}

	//loop over children.
	list <tree_key_t>::iterator child = node->second.child_edges.begin();
	while(child!=node->second.child_edges.end()){
		if (check_node(*child)){	//of child exists, apply function to child
			delete_extra_children(*child);
			child++;
		}else{						//otherwise, delete
			child = node->second.child_edges.erase(child);
		}
	}
	if (RT_VERBOSE){
		cerr <<"rooted_tree::delete_extra_children(). done"<<endl;
	}
	return 0;
}

/*
 * @brief function that recursively deletes all nodes that have only one child.
 * @params key of subtree root
 */
int rooted_tree::delete_one_child_nodes(tree_key_t subtree_root){
	if (RT_VERBOSE){
		cerr <<"rooted_tree::delete_one_child_nodes(). node age "<<subtree_root<<endl;
	}

	map <tree_key_t,node_t>::iterator node = nodes.find(subtree_root);
	if (node == nodes.end()){
		cerr <<"rooted_tree::delete_one_child_nodes(): subtree root not found! age: "<<subtree_root<<endl;
		return RT_NODENOTFOUND;
	}
	//make a copy of the children (they are modified by the recursive call)
	//loop over copy of children
	list <tree_key_t> tempchildren = node->second.child_edges;
	for (list <tree_key_t>::iterator child = tempchildren.begin(); child!=tempchildren.end(); child++){
		delete_one_child_nodes(*child);	//apply function to children
	}
	//if only one child, bridge this node
	if (node->second.child_edges.size()==1){
		bridge_edge_node(subtree_root);
	}

	if (RT_VERBOSE){
		cerr <<"rooted_tree::delete_one_child_nodes(). done"<<endl;
	}
	return 0;
}

/*
 * @brief function that checks the presence of all leafs, the congruent wiring between edges and nodes etc
 */
int rooted_tree::check_tree_integrity(){
	if (RT_VERBOSE){
		cerr <<"rooted_tree::check_tree_integrity()..."<<endl;
	}

	int err=0;
	map <tree_key_t,node_t>::iterator node;
	map <tree_key_t,edge_t>::iterator edge;
	//make sure all leafs do not have any children
	for (vector <tree_key_t>::iterator leaf=leafs.begin(); leaf!=leafs.end(); leaf++){
		if (check_node(*leaf)){
			node = nodes.find(*leaf);
			if (node->second.child_edges.size()==0){
				cerr <<"leaf "<<*leaf<<" found and has no children OK"<<endl;
			}else{
				err++;
				cerr <<"leaf "<<*leaf<<" found and has children "<<node->second.child_edges.size()<<" ERROR!"<<endl;
			}
		}else{
			err++;
			cerr <<"leaf "<<*leaf<<" not found ERROR!"<<endl;
		}
	}

	//check whether root has only one child and whether that is MRCA
	if (check_node(root)){
		node = nodes.find(root);
		if (node->second.child_edges.size()==1){
			cerr <<"root "<<root<<" found and has one child OK"<<endl;
			if (MRCA == (node->second.child_edges.front())){
				cerr <<"child is MRCA, OK: "<<MRCA<<endl;
			}else{
				err++;
				cerr <<"child not MRCA, ERROR"<<endl;
			}
		}else{
			err++;
			cerr <<"root "<<root<<" found and has "<<node->second.child_edges.size()<<" children. ERROR!"<<endl;
		}
	}else{
		err++;
		cerr <<"root "<<root<<" not found ERROR!"<<endl;
	}

	//check internal nodes
	unsigned int nnodes=0, nedges=0;
	for (node = nodes.begin(); node!=nodes.end(); node++){
		if (root != node->first){
			nedges+=node->second.child_edges.size();
			if (node->second.child_edges.size()==1){
				err++;
				cerr <<"node "<<node->first<<" is degenerate (only one child)! ERROR"<<endl;
			}
			nnodes++;
			edge=edges.find(node->first);
			if (edge!=edges.end()){
				if (edge->second.parent_node != node->second.parent_node){
					err++;
					cerr <<"edge and node "<<node->first<<" do not have the same parent! ERROR"<<endl;
				}
			}else{
				err++;
				cerr <<"edge "<<node->first<<" does not exist!"<<endl;
			}
		}else{nnodes++;nedges++;}
	}
	if ( nnodes!=nodes.size() ){
		err++;
		cerr <<"number of nodes encountered does not equal the size of nodes. ERROR"<<endl;
	}
	if ( nedges!=edges.size() ){
		err++;
		cerr <<"number of edges encountered does not equal the size of edges. "<<nedges<<" vs "<<edges.size()<<"  ERROR"<<endl;
	}
	if (err==0){
		cerr <<"Tree OK!"<<endl;
	}else{
		cerr <<"Tree messed up! "<< err<<" error(s) found!"<<endl;
	}
	return err;
}


/*
 * @brief recursive function that calculates what chunks of the genome
 *        are inherited by what number of individuals
 * @params key of the node whose descendants are to be investigated
 */
int rooted_tree::calc_weight_distribution(tree_key_t subtree_root){
	map <tree_key_t,node_t>::iterator node = nodes.find(subtree_root);
	if (RT_VERBOSE){
		cerr <<"rooted_tree::calc_weight_distribution(): of  "<<subtree_root<<endl;
	}
	if (node == nodes.end()){
		cerr <<"rooted_tree::calc_weight_distribution(): node "<<subtree_root<<" does not exist!"<<endl;
		return RT_NODENOTFOUND;
	}else{
		node->second.weight_distribution.clear();
		//if we are dealing with a terminal set the weight distribution to [0,1]*clone_size
		if (node->second.child_edges.size()==0){
			step_t temp_step;
			temp_step.pos = node->second.crossover[0];	//step up at the left end
			temp_step.step = node->second.clone_size;
			node->second.weight_distribution.push_back(temp_step);
			temp_step.pos = node->second.crossover[1];	//step down at the right end
			temp_step.step = -node->second.clone_size;
			node->second.weight_distribution.push_back(temp_step);
		}else{	//if an internal node, we loop of children
			map <tree_key_t,node_t>::iterator child_node;	//iterators used for finding
			map <tree_key_t,edge_t>::iterator child_edge;
			step_t temp_step,cumulative_step;
			//loop over children
			for (list <tree_key_t>::iterator child=node->second.child_edges.begin(); child!=node->second.child_edges.end();child++)
			{
				//recursively call yourself
				calc_weight_distribution(*child);
				//find the child
				child_node=nodes.find(*child);
				child_edge=edges.find(*child);
				cumulative_step.pos=0; cumulative_step.step=0;
				//intersect the childs weight with the edge leading up to it
				for (vector <step_t>::iterator child_step = child_node->second.weight_distribution.begin();
						child_step != child_node->second.weight_distribution.end();child_step++)
				{
					temp_step = *child_step;
					//move the position of the step up or down to the left or right edge of the
					//transmitted segment, depending on whether the step is +/-
					if (temp_step.step>0){
						if(temp_step.pos<child_edge->second.segment[0]){
							temp_step.pos = child_edge->second.segment[0];
						}
					}else if (temp_step.step<0){
						if(temp_step.pos>child_edge->second.segment[1]){
							temp_step.pos = child_edge->second.segment[1];
						}
					}
					//add step to previous step if at the same location.
					//otherwise, add previous step and remember new one
					//this works, since steps are sorted by position
					if (cumulative_step.pos!=temp_step.pos){
						node->second.weight_distribution.push_back(cumulative_step);
						cumulative_step = temp_step;
					}else{
						cumulative_step.step +=temp_step.step;
					}
				}
				node->second.weight_distribution.push_back(cumulative_step);
			}
			//sort the steps according to position. This makes the above algorithm possible
			sort(node->second.weight_distribution.begin(), node->second.weight_distribution.end());
		}
		if (RT_VERBOSE){
			cerr <<"rooted_tree::calc_weight_distribution(): done"<<endl;
		}
		return 0;
	}
}

/*
 * @brief return a string (table, in fact) representation of the weight distribution that can be plotted
 * @params key of node whose weight distribution is to be returned.
 */
string rooted_tree::print_weight_distribution(tree_key_t node_key){
	map <tree_key_t,node_t>::iterator node = nodes.find(node_key);
	stringstream WD_str;
	if (node == nodes.end()){
		cerr <<"rooted_tree::print_weight_distribution(): node "<<node_key<<" does not exist!"<<endl;
		return "bad node";
	}else{
		int cumulative = 0;
		for (vector <step_t>::iterator step=node->second.weight_distribution.begin();
				step!=node->second.weight_distribution.end();step++){
			WD_str<<step->pos<<'\t'<<0<<'\t'<<cumulative<<'\n';
			cumulative+=step->step;
			WD_str<<step->pos<<'\t'<<step->step<<'\t'<<cumulative<<'\n';
		}
	}
	return WD_str.str();
}

/*
 * @brief find nodes in a subtree younger than age
 * @params int age (the time at which the tree is to be slcied)
 * @params tree_key_t subtree_root the root of the subtree that is to be examined
 * @params reference to tree_key_t vector. vector into which the found ancestors are pushed
 */
int rooted_tree::ancestors_at_age(int age, tree_key_t subtree_root, vector <tree_key_t> &ancestors){
	if (RT_VERBOSE){
		cerr <<"rooted_tree::ancestors_at_age(): looking for ancestors in subtree with root "<<subtree_root<<endl;
	}
	int nanc=0;
	if (subtree_root.age>age){
		ancestors.push_back(subtree_root);
		return 1;
	}else{
		map <tree_key_t,node_t>::iterator node = nodes.find(subtree_root);
		for (list <tree_key_t>::iterator child = node->second.child_edges.begin(); child!=node->second.child_edges.end(); child++){
			nanc+=ancestors_at_age(age, *child, ancestors);
		}
	}
	return nanc;
}

/**
 * @brief Parse a label from the newick representation
 *
 * @param label the string with the label
 * @param index where to store the index
 * @param clone_size where to store the clone_size
 * @param branch_length where to store the branch_length
 *
 * @returns zero if successful, error codes otherwise
 */
int rooted_tree::parse_label(std::string label, int *index, int *clone_size, int *branch_length) {
	std::vector<std::string> strs;
	strs = boost::split(strs, label, boost::is_any_of(":_"));
	if(strs.size() < 3) {
		if (RT_VERBOSE) std::cerr<<"Error parsing "<<label<<"."<<std::endl;
		return RT_ERROR_PARSING;
	}
	std::stringstream(strs[0]) >> (*index);
	std::stringstream(strs[1]) >> (*clone_size);
	std::stringstream(strs[2]) >> (*branch_length);

	return 0;
}

/**
 * @brief Parse a Newick substring and add the nodes to the tree
 *
 * @param parent_key the parent key to hang the subtree from
 * @param tree_s the Newick substring to parse
 *
 * @returns zero if successful, error codes otherwise
 */
int rooted_tree::parse_subtree(tree_key_t &parent_key, std::string &tree_s) {

		std::string own_label;
		tree_key_t own_key;
		node_t own_node;
		int index, clone_size, branch_length;
		int status = 0;

		// Is this a leaf?
		size_t close_posn = tree_s.find_last_of(")");
		if (close_posn == std::string::npos) {
			// Leaf
			own_label = tree_s;
			if (RT_VERBOSE) std::cerr<<"Leaf";
		} else {
			// Internal node
			own_label = tree_s.substr(close_posn + 1, tree_s.length() - close_posn - 1);
			if (RT_VERBOSE) std::cerr<<"Internal node";
		}

		if (RT_VERBOSE) {
			std::cerr<<" of tree_key_t("<<parent_key.index<<", "<<parent_key.age<<"): ";
			std::cerr<<tree_s<<std::endl;
		}

		// Parse own label
		status = parse_label(own_label,
				     &index,
				     &clone_size,
				     &branch_length);
		if (status) {
			if (RT_VERBOSE) std::cerr<<"Error while parsing label: "<<own_label<<std::endl;
			return RT_ERROR_PARSING;	
		}
		own_key.index = index;
		own_key.age = parent_key.age + branch_length;

		// Make the current node
		own_node.parent_node = parent_key;
		own_node.own_key = own_key;
		own_node.number_of_offspring = 0;
		own_node.clone_size = clone_size;
		own_node.fitness = RT_FITNESS_MISSING;
		own_node.crossover[0] = RT_CROSSOVER_MISSING;
		own_node.crossover[1] = RT_CROSSOVER_MISSING;
		add_terminal_node(own_node);

		// Only if not a leaf, find out about subtrees
		if (close_posn == std::string::npos) {
			leafs.push_back(own_key);
			return status;
		}
		std::vector<std::string> subtrees;
        	std::string subtext;
		int plevel = 0;
		int prev = 1;
		char c;
		size_t posn = 1;
		for(std::string::iterator it=tree_s.begin() + 1; posn != close_posn; it++, posn++) {
			c = (*it);
			if (c == '(') {
				plevel++;
			} else if (c == ')') {
				plevel--;
			} else if ((c == ',') and (plevel == 0)) {
				subtext = tree_s.substr(prev, posn - prev);
				subtrees.push_back(subtext);
				prev = posn + 1;	
			}
		}
		subtext = tree_s.substr(prev, close_posn - prev);
		subtrees.push_back(subtext);
		
		// OK, now you have the subtrees with a lot of string copying around
		if (RT_VERBOSE) {
			std::cerr<<"Subtrees of tree_key_t("<<own_key.index<<", "<<own_key.age<<"): ";
			for(std::vector<std::string>::iterator it=subtrees.begin(); it != subtrees.end(); it++) {
				std::cerr<<(*it)<<"\t";
			}
			std::cerr<<std::endl;
		}

		// Add recursively all children
		for(std::vector<std::string>::iterator it=subtrees.begin(); it != subtrees.end(); it++) {
				status = parse_subtree(own_key, (*it));
			}

		//FIXME
		return status;
}

/**
 * @brief Read a newick strong and make the tree from there
 *
 * @param tree_s the newick string
 *
 * @returns zero if successful, error codes otherwise
 */
int rooted_tree::read_newick(string tree_s) {

	int status = 0;
        if(RT_VERBOSE) std::cout<<"Tree in:  "<<tree_s<<std::endl;

        // Check for open and closed parentheses
        unsigned long nopen = 0, nclose = 0;
	for(std::string::iterator i=tree_s.begin(); i != tree_s.end(); i++) {
		if ((*i) == '(')
			nopen++;
		else if ((*i) == ')')
			nclose++;
	}
	if (nopen != nclose) {
		if (RT_VERBOSE) std::cerr<<"Newick file corrupted (parentheses not matching)"<<std::endl;
		return 1;
	}

	// Parse
	tree_key_t tmp_key;
	node_t tmp_node;
	edge_t tmp_edge;
	std::string tmp_s;
	int index, clone_size, branch_length;
	std::string root_label;

	// is it a trivial tree?
	size_t close_posn = tree_s.find_last_of(")");
	if (close_posn == std::string::npos)
		root_label = tree_s.substr(0, tree_s.length() - 1);
	else
		root_label = tree_s.substr(close_posn + 1, tree_s.length() - close_posn - 2);

	// The first level of parsing is special because:
	// 1. there is a single node (the root)
	// 2. the node must be added to the tree is a special way
	// 3. the label ends with a semicolon and a newline (probably, FIXME)
	if (RT_VERBOSE) std::cerr<<"Root label: "<<root_label<<std::endl;
	status = parse_label(root_label,
			     &index,
			     &clone_size,
			     &branch_length);
	if (status) {
		if (RT_VERBOSE) std::cerr<<"Error while parsing root label"<<std::endl;
		return RT_ERROR_PARSING;	
	}
	// Make the MRCA node
	reset();
	leafs.clear();
	tmp_key.index = index;
	tmp_key.age = branch_length + root.age;
	tmp_node.parent_node = MRCA;
	tmp_node.own_key = tmp_key;
	tmp_node.number_of_offspring = 0; // updated later
	tmp_node.clone_size = clone_size;
	tmp_node.fitness = RT_FITNESS_MISSING;
	tmp_node.crossover[0] = RT_CROSSOVER_MISSING;
	tmp_node.crossover[1] = RT_CROSSOVER_MISSING;
	add_terminal_node(tmp_node);
	bridge_edge_node(MRCA);

	// is it a trivial tree?
	if (close_posn == std::string::npos) {
		leafs.push_back(tmp_key);
	} else {
		// Find out about subtrees
		std::vector<std::string> subtrees;
		std::vector<std::string>::iterator it;
        	std::string subtext;
		int plevel = 0;
		int prev = 1;
		char c;
		size_t posn = 1;
		for(std::string::iterator its=tree_s.begin() + 1; posn != close_posn; its++, posn++) {
			c = (*its);
			if (c == '(')
				plevel++;
			else if (c == ')')
				plevel--;
			else if ((c == ',') and (plevel == 0)) {
				subtext = tree_s.substr(prev, posn - prev);
				subtrees.push_back(subtext);
				prev = posn + 1;
			}
		}
		subtext = tree_s.substr(prev, close_posn - prev);
		subtrees.push_back(subtext);

		// OK, now you have the subtrees with a lot of string copying around
		if (RT_VERBOSE) {
			std::cerr<<"Subtrees:"<<std::endl;
			for(it=subtrees.begin(); it != subtrees.end(); it++)
				std::cerr<<(*it)<<std::endl;
		}

		// Add recursively all children
		for(it=subtrees.begin(); it != subtrees.end(); it++)
				status = parse_subtree(tmp_key, (*it));
	}

	// update number of offspring
	update_tree();

	// Now the tree is done. Let us print newick to test
	if (RT_VERBOSE) std::cout<<"Tree out: "<<print_newick()<<std::endl;

	return status;
}
