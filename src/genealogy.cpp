#import "genealogy.h"


genealogy::genealogy(){
	reset();
}

void genealogy::reset(){
	node_t root_node;
	node_t mrca_node;
 	edge_t to_root;

	root.age=-2;
	root.index=0;
	MRCA.age=-1;
	MRCA.index=0;

	root_node.own_key = root;
	root_node.parent_node = root;
	root_node.clone_size=-1;
	root_node.number_of_offspring=-1;
	root_node.crossover[0]=0;
	root_node.crossover[1]=GEN_VERYLARGE;
	root_node.child_edges.clear();
	root_node.child_edges.push_back(MRCA);

	mrca_node.own_key = MRCA;
	mrca_node.parent_node = root;
	mrca_node.clone_size=-1;
	mrca_node.number_of_offspring=-1;
	mrca_node.crossover[0]=0;
	mrca_node.crossover[1]=GEN_VERYLARGE;
	mrca_node.child_edges.clear();

	to_root.own_key = MRCA;
	to_root.length=mrca_node.own_key.age-root_node.own_key.age;
	to_root.segment[0]=0;
	to_root.segment[1]=GEN_VERYLARGE;
	to_root.number_of_offspring=-1;
	to_root.parent_node = root;

	nodes.insert(pair<key_t,node_t>(root,root_node));
	nodes.insert(pair<key_t,node_t>(MRCA,mrca_node));
	edges.insert(pair<key_t,edge_t>(to_root.own_key,to_root));
}

genealogy::~genealogy(){}

void genealogy::add_generation(vector <node_t> &new_generation, double mean_fitness){
	if (GEN_VERBOSE){
		cerr <<"genealogy::add_generation(). Number of leafs to add: "<<new_generation.size()<<endl;
	}

	vector<node_t>::iterator new_leaf = new_generation.begin();
	edge_t new_edge;
	key_t new_key;
	vector <key_t> new_leafs;
	new_leafs.reserve(new_generation.size());
	//add new leafs
	for (; new_leaf!=new_generation.end(); new_leaf++){
		if (new_leaf->clone_size>0){
			new_key= new_leaf->own_key;
			new_leaf->child_edges.clear();
			new_edge.own_key=new_key;
			new_edge.parent_node=new_leaf->parent_node;
			new_edge.number_of_offspring=1;
			new_edge.segment[0]=new_leaf->crossover[0];
			new_edge.segment[1]=new_leaf->crossover[1];
			new_edge.length=1;
			new_leafs.push_back(new_key);
			nodes[new_leaf->parent_node].child_edges.push_back(new_key);
			edges.insert(pair<key_t,edge_t>(new_key, new_edge));
			nodes.insert(pair<key_t,node_t>(new_key, *new_leaf));
		}
	}

	if (GEN_VERBOSE){
		cerr <<"genealogy::add_generation(). added leafes... erase dead ends"<<endl;
		cerr <<"genealogy::add_generation(). genealogy size: "<<edges.size()<<" edges, "<<nodes.size()<<" nodes "<<endl;
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
			cerr <<"genealogy::add_generation(). did not find old leaf"<<old_leaf_key->age<<" "<<old_leaf_key->index<<endl;
		}
	}

	if (GEN_VERBOSE){
		cerr <<"genealogy::add_generation(). done "<<endl;
	}

	leafs.clear();
	leafs.reserve(new_leafs.size());
	for(vector <key_t>::iterator new_leaf_key=new_leafs.begin(); new_leaf_key!=new_leafs.end(); new_leaf_key++){
		leafs.push_back(*new_leaf_key);
	}

	if (GEN_VERBOSE){
		cerr <<"genealogy::add_generation(). done "<<endl;
		cerr <<"genealogy::add_generation(). genealogy size: "<<edges.size()<<" edges, "<<nodes.size()<<" nodes "<<endl;
	}

	clear_tree();
	update_tree();
	return;
}



key_t genealogy::erase_edge_node(key_t to_be_erased){
	if (GEN_VERBOSE){
		cerr <<"genealogy::erase_edge_node(). ..."<<to_be_erased.age<<" "<<to_be_erased.index<<endl;
	}

	map <key_t,node_t>::iterator Enode = nodes.find(to_be_erased);
	map <key_t,edge_t>::iterator Eedge = edges.find(to_be_erased);

	if (Enode->second.child_edges.size()>0){
		cerr <<"genealogy::erase_edge_node(): attempting to erase non-terminal node"<<endl;
	}

	key_t parent_key = Eedge->second.parent_node;
	map <key_t,edge_t>::iterator Pedge = edges.find(parent_key);
	map <key_t,node_t>::iterator Pnode = nodes.find(parent_key);

	Pnode->second.number_of_offspring-=Eedge->second.number_of_offspring;
	Pedge->second.number_of_offspring-=Eedge->second.number_of_offspring;


	if (erase_child(Pnode, to_be_erased)==GEN_CHILDNOTFOUND){
		cerr <<"genealogy::erase_edge_node(): child not found"<<endl;
	}

	nodes.erase(to_be_erased);
	edges.erase(to_be_erased);

	if (GEN_VERBOSE){
		cerr <<"genealogy::erase_edge_node(). done"<<endl;
	}

	return parent_key;
}

int genealogy::erase_child(map <key_t,node_t>::iterator Pnode, key_t to_be_erased){
	list <key_t>::iterator child = Pnode->second.child_edges.begin();
	for (;child!=Pnode->second.child_edges.end(); child++){
		if (*child == to_be_erased){Pnode->second.child_edges.erase(child); return 0;}
	}
	return GEN_CHILDNOTFOUND;
}

key_t genealogy::bridge_edge_node(key_t to_be_bridged){
	if (GEN_VERBOSE){
		cerr <<"genealogy::erase_bridged_node(). ..."<<to_be_bridged.age<<" "<<to_be_bridged.index<<endl;
	}
	map <key_t,node_t>::iterator Enode = nodes.find(to_be_bridged);
	map <key_t,edge_t>::iterator Eedge = edges.find(to_be_bridged);

	if (Enode->second.child_edges.size()!=1){
		cerr <<"genealogy::erase_bridged_node(): attempting to bridge branched node"<<endl;
	}

	key_t parent_key = Eedge->second.parent_node;
	map <key_t,edge_t>::iterator Pedge = edges.find(Enode->second.child_edges.front());
	Pedge->second.parent_node = Eedge->second.parent_node;
	Pedge->second.segment[0]=(Pedge->second.segment[0]<Eedge->second.segment[0])?(Eedge->second.segment[0]):(Pedge->second.segment[0]);
	Pedge->second.segment[1]=(Pedge->second.segment[1]>Eedge->second.segment[0])?(Eedge->second.segment[1]):(Pedge->second.segment[1]);
	Pedge->second.length+=Eedge->second.length;
	Pedge->second.parent_node = Eedge->second.parent_node;

	map <key_t,node_t>::iterator Pnode = nodes.find(Eedge->second.parent_node);
	Pnode->second.child_edges.push_back(Pedge->first);
	if (erase_child(Pnode, to_be_bridged)==GEN_CHILDNOTFOUND){
		cerr <<"genealogy::erase_bridged_node(). child not found"<<endl;
	}
	nodes.erase(to_be_bridged);
	edges.erase(to_be_bridged);
	if (to_be_bridged == MRCA){MRCA = Pedge->first;}

	if (GEN_VERBOSE){
		cerr <<"genealogy::bridge_edge_node(). done"<<endl;
	}

	return parent_key;
}

void genealogy::update_tree(){
	clear_tree();
	for (vector <key_t>::iterator leaf=leafs.begin(); leaf!=leafs.end(); leaf++){
		update_leaf_to_root(*leaf);
	}
}

void genealogy::update_leaf_to_root(key_t leaf_key){
	if (GEN_VERBOSE){
		cerr <<"genealogy::update_leaf_to_root(). key:"<<leaf_key.index<<" "<<leaf_key.age<<endl;
	}
	map <key_t,node_t>::iterator leaf_node = nodes.find(leaf_key);
	map <key_t,edge_t>::iterator leaf_edge = edges.find(leaf_key);
	int increment = leaf_node->second.number_of_offspring;
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
			cerr <<"genealogy::update_leaf_to_root(): key not found: "<<leaf_edge->second.parent_node.index<<" "<<leaf_edge->second.parent_node.age<<" root: "<<root.index<<" "<<root.age <<endl;
			break;
		}
	}
	if (GEN_VERBOSE){
		cerr <<"genealogy::update_leaf_to_root(). done"<<endl;
		cerr <<"genealogy::update_leaf_to_root(): total of "<< nodes.find(MRCA)->second.number_of_offspring<<" offspring "<<increment<<endl;
	}
}

void genealogy::SFS(gsl_histogram *sfs){
	map <key_t,edge_t>::iterator edge = edges.begin();
	int total_pop = nodes[root].number_of_offspring;
	for (; edge!=edges.end(); edge++){
		gsl_histogram_accumulate(sfs, 1.0*edge->second.number_of_offspring/total_pop, edge->second.length);
	}
}

int genealogy::external_branch_length(){
	map <key_t,edge_t>::iterator edge = edges.begin();
	int branchlength = 0;
	for (vector <key_t>::iterator leaf=leafs.begin(); leaf!=leafs.end();leaf++){
		branchlength+=edges[*leaf].length;
	}
	return branchlength;
}

int genealogy::total_branch_length(){
	map <key_t,edge_t>::iterator edge = edges.begin();
	int branchlength = 0;
	for (; edge!=edges.end(); edge++){
		branchlength+=edge->second.length;
	}
	return branchlength;
}


void genealogy::clear_tree(){
	if (GEN_VERBOSE){
		cerr <<"genealogy::clear_tree()..."<<endl;
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
			cerr <<"genealogy::clear_tree(). key note found"<<endl;
			break;
		}
		node->second.number_of_offspring=1;
	}
	if (GEN_VERBOSE){
		cerr <<"genealogy::clear_tree(). done"<<endl;
	}
}

int genealogy::construct_subtree(vector <key_t> subtree_leafs){
	vector <key_t>::iterator leaf=subtree_leafs.begin();
	set <key_t> new_nodes;
	map <key_t,node_t>::iterator node = nodes.find(*leaf);
	map <key_t,edge_t>::iterator edge = edges.find(*leaf);
	int oldest_node_age = node->second.age;
	subtree_nodes.insert(node);
	subtree_edges.insert(edge);
	leaf++;
	new_nodes.clear();
	for (; leaf!=subtree_leafs.end(); leaf++){
		node = nodes.find(*leaf);
		edge = edges.find(*leaf);
		if (node->second.age<oldest_node_age) oldest_node_age = node->second.age;
		subtree_nodes.insert(node);
		subtree_edges.insert(edge);
		new_nodes.insert(node->second.parent_node);
	}

	while(new_nodes.size()>1){
		set <key_t> temp = new_nodes;
		new_nodes.clear();
		for (set <key_t>::iterator node_key=temp.begin(); node_key!=temp.end(); node_key++){
			node = nodes.find(*node_key);
			edge = edges.find(*node_key);
			if (node->second.age<oldest_node_age){
				if (node->second.parent_node.age<oldest_node_age) oldest_node_age = node->second.parent_node.age;
				subtree_nodes.insert(node);
				subtree_edges.insert(edge);
				new_nodes.insert(node->second.parent_node);
			}
		}
	}
	subtree_MRCA = *new_nodes.begin();
	node = nodes.find(subtree_MRCA);
	edge = edges.find(subtree_MRCA);
	subtree_nodes.insert(node);

	delete_extra_children_in_subtree(subtree_MRCA);

	return 0;
}

string genealogy::print_newick(){
	return subtree_newick(MRCA)+";";
}

string genealogy::subtree_newick(key_t root){
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
	tree_str<<root.index<<'_'<<root_node->second.clone_size<<":"<<edge->second.length;
	return tree_str.str();
}

bool genealogy::check_node(key_t node_key){
	map <key_t,node_t>::iterator node = nodes.find(node_key);
	if (node == nodes.end()){
		return false;
	}else { return true;}
}

