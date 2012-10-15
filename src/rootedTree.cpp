/*
 * tree.cpp
 *
 *  Created on: Oct 14, 2012
 *      Author: richard
 */

#include "ffpopsim_highd.h"

std::ostream& operator<< ( std::ostream& os, const key_t& key )
  {
    os <<"age: "<< key.age << " index "
      << key.index;
    return os;
  }


rootedTree::rootedTree() {
	reset();
}
rootedTree::~rootedTree() {
}


void rootedTree::reset(){
	node_t root_node;
	node_t mrca_node;
 	edge_t to_root;

 	nodes.clear();
 	edges.clear();
 	leafs.clear();

	root.age=-2;
	root.index=0;
	MRCA.age=-1;
	MRCA.index=0;

	root_node.own_key = root;
	root_node.parent_node = root;
	root_node.clone_size=-1;
	root_node.number_of_offspring=-1;
	root_node.crossover[0]=0;
	root_node.crossover[1]=RT_VERYLARGE;
	root_node.child_edges.clear();
	root_node.child_edges.push_back(MRCA);

	mrca_node.own_key = MRCA;
	mrca_node.parent_node = root;
	mrca_node.clone_size=-1;
	mrca_node.number_of_offspring=-1;
	mrca_node.crossover[0]=0;
	mrca_node.crossover[1]=RT_VERYLARGE;
	mrca_node.child_edges.clear();

	to_root.own_key = MRCA;
	to_root.length=mrca_node.own_key.age-root_node.own_key.age;
	to_root.segment[0]=0;
	to_root.segment[1]=RT_VERYLARGE;
	to_root.number_of_offspring=-1;
	to_root.parent_node = root;

	nodes.insert(pair<key_t,node_t>(root,root_node));
	nodes.insert(pair<key_t,node_t>(MRCA,mrca_node));
	edges.insert(pair<key_t,edge_t>(to_root.own_key,to_root));
}


void rootedTree::add_generation(vector <node_t> &new_generation, double mean_fitness){
	if (RT_VERBOSE){
		cerr <<"rootedTree::add_generation(). Number of leafs to add: "<<new_generation.size()<<endl;
	}

	vector<node_t>::iterator new_leaf = new_generation.begin();
	vector <key_t> new_leafs;
	new_leafs.reserve(new_generation.size());
	//add new leafs
	for (; new_leaf!=new_generation.end(); new_leaf++){
		if (new_leaf->clone_size>0){
			add_terminal_node(*new_leaf);
			new_leafs.push_back(new_leaf->own_key);
		}
	}

	if (RT_VERBOSE){
		cerr <<"rootedTree::add_generation(). added leafs... erase dead ends"<<endl;
		cerr <<"rootedTree::add_generation(). rootedTree size: "<<edges.size()<<" edges, "<<nodes.size()<<" nodes "<<endl;
	}

	map <key_t,edge_t>::iterator edge_pos = edges.end();
	map <key_t,node_t>::iterator node_pos = nodes.end();
	key_t parent_key;

	for(vector <key_t>::iterator old_leaf_key=leafs.begin(); old_leaf_key!=leafs.end(); old_leaf_key++){
		node_pos = nodes.find(*old_leaf_key);
		if (node_pos!=nodes.end()){
			if (node_pos->second.child_edges.size()==0){
				parent_key = erase_edge_node(*old_leaf_key);
				while (nodes[parent_key].child_edges.size()==0){
					parent_key = erase_edge_node(parent_key);
				}
				while (nodes[parent_key].child_edges.size()==1 and parent_key!=root){
					parent_key = bridge_edge_node(parent_key);
				}
			}else if (node_pos->second.child_edges.size()==1){
				parent_key = bridge_edge_node(*old_leaf_key);
				while (nodes[parent_key].child_edges.size()==1){
					parent_key = bridge_edge_node(parent_key);
				}
			}
		}else{
			cerr <<"rootedTree::add_generation(). did not find old leaf"<<old_leaf_key->age<<" "<<old_leaf_key->index<<endl;
		}
	}

	if (RT_VERBOSE){
		cerr <<"rootedTree::add_generation(). done "<<endl;
	}

	leafs = new_leafs;

	if (RT_VERBOSE){
		cerr <<"rootedTree::add_generation(). done "<<endl;
		cerr <<"rootedTree::add_generation(). rootedTree size: "<<edges.size()<<" edges, "<<nodes.size()<<" nodes "<<endl;
	}

	update_tree();
	return;
}

int rootedTree::add_terminal_node(node_t &newNode){
	edge_t newEdge;
	key_t newKey;
	newKey = newNode.own_key;
	newNode.child_edges.clear();
	newEdge.own_key=newKey;
	newEdge.parent_node=newNode.parent_node;
	newEdge.number_of_offspring=1;
	newEdge.segment[0]=newNode.crossover[0];
	newEdge.segment[1]=newNode.crossover[1];
	newEdge.length=newKey.age-newNode.parent_node.age;
	nodes[newNode.parent_node].child_edges.push_back(newKey);
	edges.insert(pair<key_t,edge_t>(newKey, newEdge));
	nodes.insert(pair<key_t,node_t>(newKey, newNode));
	return 0;
}

key_t rootedTree::erase_edge_node(key_t to_be_erased){
	if (RT_VERBOSE){
		cerr <<"rootedTree::erase_edge_node(). ..."<<to_be_erased.age<<" "<<to_be_erased.index<<endl;
	}

	map <key_t,node_t>::iterator Enode = nodes.find(to_be_erased);
	map <key_t,edge_t>::iterator Eedge = edges.find(to_be_erased);

	if (Enode->second.child_edges.size()>0){
		cerr <<"rootedTree::erase_edge_node(): attempting to erase non-terminal node"<<endl;
	}

	key_t parent_key = Eedge->second.parent_node;
	map <key_t,edge_t>::iterator Pedge = edges.find(parent_key);
	map <key_t,node_t>::iterator Pnode = nodes.find(parent_key);

	Pnode->second.number_of_offspring-=Eedge->second.number_of_offspring;
	Pedge->second.number_of_offspring-=Eedge->second.number_of_offspring;


	if (erase_child(Pnode, to_be_erased)==RT_CHILDNOTFOUND){
		cerr <<"rootedTree::erase_edge_node(): child not found"<<endl;
	}

	nodes.erase(to_be_erased);
	edges.erase(to_be_erased);

	if (RT_VERBOSE){
		cerr <<"rootedTree::erase_edge_node(). done"<<endl;
	}

	return parent_key;
}

int rootedTree::erase_child(map <key_t,node_t>::iterator Pnode, key_t to_be_erased){
	for (list <key_t>::iterator child = Pnode->second.child_edges.begin();child!=Pnode->second.child_edges.end(); child++){
		if (*child == to_be_erased){Pnode->second.child_edges.erase(child); return 0;}
	}
	return RT_CHILDNOTFOUND;
}

key_t rootedTree::bridge_edge_node(key_t to_be_bridged){
	if (RT_VERBOSE){
		cerr <<"rootedTree::bridge_edge_node(). ..."<<to_be_bridged.age<<" "<<to_be_bridged.index<<endl;
	}
	map <key_t,node_t>::iterator Enode = nodes.find(to_be_bridged);
	map <key_t,edge_t>::iterator Eedge = edges.find(to_be_bridged);

	if (Enode->second.child_edges.size()!=1 or to_be_bridged==root){
		cerr <<"rootedTree::bridge_edge_node(): attempting to bridge branched node or bridge root"<<endl;
	}

	key_t parent_key = Eedge->second.parent_node;
	map <key_t,edge_t>::iterator Pedge = edges.find(Enode->second.child_edges.front());
	map <key_t,node_t>::iterator ChildNode = nodes.find(Enode->second.child_edges.front());
	Pedge->second.parent_node = Eedge->second.parent_node;
	ChildNode->second.parent_node = Eedge->second.parent_node;
	Pedge->second.segment[0]=(Pedge->second.segment[0]<Eedge->second.segment[0])?(Eedge->second.segment[0]):(Pedge->second.segment[0]);
	Pedge->second.segment[1]=(Pedge->second.segment[1]>Eedge->second.segment[0])?(Eedge->second.segment[1]):(Pedge->second.segment[1]);
	Pedge->second.length+=Eedge->second.length;
	Pedge->second.parent_node = Eedge->second.parent_node;

	map <key_t,node_t>::iterator Pnode = nodes.find(Eedge->second.parent_node);
	Pnode->second.child_edges.push_back(Pedge->first);
	if (erase_child(Pnode, to_be_bridged)==RT_CHILDNOTFOUND){
		cerr <<"rootedTree::bridge_edge_node(). child not found. index "<<to_be_bridged.index<<" age "<<to_be_bridged.age<<endl;
	}
	nodes.erase(to_be_bridged);
	edges.erase(to_be_bridged);
	if (to_be_bridged == MRCA){MRCA = Pedge->first;}

	if (RT_VERBOSE){
		cerr <<"rootedTree::bridge_edge_node(). done"<<endl;
	}

	return parent_key;
}

void rootedTree::update_tree(){
	clear_tree();
	for (vector <key_t>::iterator leaf=leafs.begin(); leaf!=leafs.end(); leaf++){
		update_leaf_to_root(*leaf);
	}
}

int rootedTree::update_leaf_to_root(key_t leaf_key){
	if (RT_VERBOSE){
		cerr <<"rootedTree::update_leaf_to_root(). key:"<<leaf_key.index<<" "<<leaf_key.age<<endl;
	}
	map <key_t,node_t>::iterator leaf_node = nodes.find(leaf_key);
	map <key_t,edge_t>::iterator leaf_edge = edges.find(leaf_key);
	if (leaf_node == nodes.end() or leaf_edge == edges.end()){
		cerr <<"rootedTree::update_leaf_to_root(). leaf not found"<<endl;
		return RT_NODENOTFOUND;
	}

	int increment = leaf_node->second.clone_size;
	leaf_edge->second.number_of_offspring = increment;
	map <key_t,node_t>::iterator parent_node = nodes.find(leaf_edge->second.parent_node);
	map <key_t,edge_t>::iterator parent_edge = edges.find(leaf_edge->second.parent_node);
	while (root != parent_node->first){
		parent_node->second.number_of_offspring+=increment;
		parent_edge->second.number_of_offspring+=increment;

		leaf_node = parent_node;
		leaf_edge = parent_edge;
		parent_node = nodes.find(leaf_edge->second.parent_node);
		parent_edge = edges.find(leaf_edge->second.parent_node);
		if (parent_node==nodes.end()){
			cerr <<"rootedTree::update_leaf_to_root(). key:"<<leaf_key.index<<" "<<leaf_key.age<<endl;
			cerr <<"rootedTree::update_leaf_to_root(): key not found: "<<leaf_edge->second.parent_node.index<<" "<<leaf_edge->second.parent_node.age<<" root: "<<root.index<<" "<<root.age <<endl;
			break;
		}
	}
	if (RT_VERBOSE){
		cerr <<"rootedTree::update_leaf_to_root(). done"<<endl;
		cerr <<"rootedTree::update_leaf_to_root(): total of "<< nodes.find(MRCA)->second.number_of_offspring<<" offspring "<<increment<<endl;
	}
	return 0;
}

void rootedTree::SFS(gsl_histogram *sfs){
	map <key_t,edge_t>::iterator edge = edges.begin();
	int total_pop = nodes[MRCA].number_of_offspring;
	for (; edge!=edges.end(); edge++){
		gsl_histogram_accumulate(sfs, 1.0*edge->second.number_of_offspring/total_pop, edge->second.length);
	}
}

int rootedTree::external_branch_length(){
	map <key_t,edge_t>::iterator edge = edges.begin();
	int branchlength = 0;
	for (vector <key_t>::iterator leaf=leafs.begin(); leaf!=leafs.end();leaf++){
		branchlength+=edges[*leaf].length;
	}
	return branchlength;
}

int rootedTree::total_branch_length(){
	map <key_t,edge_t>::iterator edge = edges.begin();
	int branchlength = 0;
	for (; edge!=edges.end(); edge++){
		branchlength+=edge->second.length;
	}
	return branchlength;
}


void rootedTree::clear_tree(){
	if (RT_VERBOSE){
		cerr <<"rootedTree::clear_tree()..."<<endl;
	}
	map <key_t,node_t>::iterator node = nodes.begin();
	map <key_t,edge_t>::iterator edge = edges.begin();

	for (; node!=nodes.end(); node++){
		node->second.number_of_offspring=0;
	}
	for (; edge!=edges.end(); edge++){
		edge->second.number_of_offspring=0;
	}


	for (vector<key_t>::iterator leaf=leafs.begin(); leaf!=leafs.end(); leaf++){
		node = nodes.find(*leaf);
		if (node==nodes.end()){
			cerr <<"rootedTree::clear_tree(). key of leaf not found"<<endl;
			break;
		}
		node->second.number_of_offspring=1;
	}
	if (RT_VERBOSE){
		cerr <<"rootedTree::clear_tree(). done"<<endl;
	}
}

string rootedTree::print_newick(){
	return subtree_newick(MRCA)+";";
}


string rootedTree::subtree_newick(key_t root){
	stringstream tree_str;
	map <key_t,node_t>::iterator root_node = nodes.find(root);
	map <key_t,edge_t>::iterator edge = edges.find(root);
	if (root_node->second.child_edges.size()>0){
		list <key_t>::iterator child = root_node->second.child_edges.begin();
		tree_str.str();
		tree_str <<"("<< subtree_newick(*child);
		child++;
		for (;child!=root_node->second.child_edges.end(); child++){
			tree_str<<","+subtree_newick(*child);
		}
		tree_str<<")";
	}
	//tree_str<<root.index<<'_'<<root_node->second.clone_size<<":"<<edge->second.length;
	tree_str<<root.index<<'_'<<root.age<<":"<<edge->second.length;
	return tree_str.str();
}

bool rootedTree::check_node(key_t node_key){
	map <key_t,node_t>::iterator node = nodes.find(node_key);
	if (node == nodes.end()){
		return false;
	}else { return true;}
}

int rootedTree::construct_subtree(vector <key_t> subtree_leafs, rootedTree &superTree){
	if (RT_VERBOSE){
		cerr <<"rootedTree::construct_subtree()..."<<endl;
	}
	reset();nodes.erase(MRCA);edges.erase(MRCA);

	set <key_t> new_nodes;
	new_nodes.clear();
	map <key_t,node_t>::iterator node;
	map <key_t,edge_t>::iterator edge;

	for (vector <key_t>::iterator leaf=subtree_leafs.begin(); leaf!=subtree_leafs.end(); leaf++){
		if (superTree.check_node(*leaf)){
			node = superTree.nodes.find(*leaf);
			nodes.insert(*node);
			edges.insert(*superTree.edges.find(*leaf));
			leafs.push_back(*leaf);
			new_nodes.insert(node->second.parent_node);
		}else{
			cerr <<"rootedTree::construct_subtree(). leaf does not exist"<<endl;
			return RT_NODENOTFOUND;
		}

	}
	if (RT_VERBOSE){
		cerr <<"rootedTree::construct_subtree(). added leafs"<<endl;
	}

	MRCA=root;
	while((new_nodes.size()>1 or MRCA!= *new_nodes.begin()) and new_nodes.size()>0){
		set <key_t> temp = new_nodes;
		new_nodes.clear();
		for (set <key_t>::iterator node_key=temp.begin(); node_key!=temp.end(); node_key++){
			if (superTree.check_node(*node_key)){
				node = superTree.nodes.find(*node_key);
				nodes.insert(*node);
				edges.insert(*superTree.edges.find(*node_key));
				if (node->second.parent_node!=root) {new_nodes.insert(node->second.parent_node);}
				else {
					MRCA = *node_key;
					//cerr <<"ran into root "<<node_key->age<<" "<<node_key->index<<endl;
				}
			}
			else{
				cerr <<"rootedTree::construct_subtree(). internal node did not exist: age "<<node_key->age<<" index "<<node_key->index<<endl;
			}
		}
	}
	if (RT_VERBOSE){
		cerr <<"rootedTree::construct_subtree(). added internal nodes"<<endl;
	}

	//key_t lastnode = *new_nodes.begin();
	//MRCA = lastnode;
	if (superTree.check_node(MRCA)){
		node = nodes.insert(*superTree.nodes.find(MRCA)).first;
		edge = edges.insert(*superTree.edges.find(MRCA)).first;
		edge->second.length = MRCA.age - root.age;
		node->second.parent_node=root;
		edge->second.parent_node=root;
		node = nodes.find(root);
		node->second.child_edges.clear();
		node->second.child_edges.push_back(MRCA);
	}


	delete_extra_children(MRCA);
	delete_one_child_nodes(MRCA);

	update_tree();

	if (RT_VERBOSE){
		cerr <<"genealogy::construct_subtree(). done"<<endl;
	}

	return 0;
}


int rootedTree::delete_extra_children(key_t subtree_root){
	if (RT_VERBOSE){
		cerr <<"rootedTree::delete_extra_children(). node age "<<subtree_root.age<<" index "<<subtree_root.index<<endl;
	}

	map <key_t,node_t>::iterator node = nodes.find(subtree_root);
	if (node == nodes.end()){
		cerr <<"rootedTree::delete_extra_children(): subtree root not found! age: "<<subtree_root.age<<" index: "<<subtree_root.index<<endl;
		return RT_NODENOTFOUND;
	}

	list <key_t>::iterator child = node->second.child_edges.begin();
	while(child!=node->second.child_edges.end()){
		if (check_node(*child)){
			delete_extra_children(*child);
			child++;
		}else{
			child = node->second.child_edges.erase(child);
		}
	}
	if (RT_VERBOSE){
		cerr <<"rootedTree::delete_extra_children(). done"<<endl;
	}
	return 0;
}

int rootedTree::delete_one_child_nodes(key_t subtree_root){
	if (RT_VERBOSE){
		cerr <<"rootedTree::delete_one_child_nodes(). node age "<<subtree_root.age<<" index "<<subtree_root.index<<endl;
	}

	map <key_t,node_t>::iterator node = nodes.find(subtree_root);
	if (node == nodes.end()){
		cerr <<"rootedTree::delete_one_child_nodes(): subtree root not found! age: "<<subtree_root.age<<" index: "<<subtree_root.index<<endl;
		return RT_NODENOTFOUND;
	}
	list <key_t> tempchildren = node->second.child_edges;
	for (list <key_t>::iterator child = tempchildren.begin(); child!=tempchildren.end(); child++){
		delete_one_child_nodes(*child);
	}
	if (node->second.child_edges.size()==1){
		bridge_edge_node(subtree_root);
	}

	if (RT_VERBOSE){
		cerr <<"rootedTree::delete_one_child_nodes(). done"<<endl;
	}
	return 0;
}

int rootedTree::check_tree_integrity(){
	if (RT_VERBOSE){
		cerr <<"rootedTree::check_tree_integrity()..."<<endl;
	}

	int err=0;
	map <key_t,node_t>::iterator node;
	map <key_t,edge_t>::iterator edge;
	//make sure all leafs do not have any children
	for (vector <key_t>::iterator leaf=leafs.begin(); leaf!=leafs.end(); leaf++){
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
				cerr <<"child is MRCA, OK"<<endl;
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
		cerr <<"number of edges encountered does not equal the size of edges."<<nedges<<" vs "<<edges.size()<<"  ERROR"<<endl;
	}
	if (err==0){
		cerr <<"Tree OK!"<<endl;
	}else{
		cerr <<"Tree messed up! "<< err<<" error(s) found!"<<endl;
	}
	return err;
}
