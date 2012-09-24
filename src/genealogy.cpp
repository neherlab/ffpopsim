#import "genealogy.h"


genealogy::genealogy(){
	node_t root;
	edge_t to_root;
}

genealogy::~genealogy(){}

void genealogy::add_generation(vector <node_t> &new_generation){
	if (GEN_VERBOSE){
		cerr <<"genealogy::add_generation(). Number of leafs to add:"<<new_generation.size()<<endl;
	}
	vector<node_t>::iterator new_leaf = new_generation.begin();
	edge_t new_edge;
	key_t new_key;
	vector <key_t> new_leafs;
	new_leafs.reserve(new_generation.size());
	//add new leafs
	for (; new_leaf!=new_generation.end(); new_leaf++){
		if (new_leaf->number_of_offspring>0){
			new_key= new_leaf->own_key;
			new_edge.child_node=new_key;
			new_edge.parent_node=new_leaf->parent_edge;
			new_edge.number_of_offspring=1;
			new_edge.segment[0]=new_leaf->crossover[0];
			new_edge.segment[1]=new_leaf->crossover[1];
			new_edge.length=1;
			new_leafs.push_back(new_key);
			nodes[new_leaf->parent_edge].child_edges.push_back(new_key);
			edges.insert(<new_key, new_edge>);
			nodes.insert(new_key, *new_leaf);
		}
	}

	if (GEN_VERBOSE){
		cerr <<"genealogy::add_generation(). added leafes... erase dead ends"<<endl;
	}


	map <key_t,edge_t>::iterator edge_pos = edges.end();
	map <key_t,node_t>::iterator node_pos = nodes.end();
	key_t parent_key;
	for(vector <key_t>::iterator old_leaf_key=leafs.begin(); old_leaf_key!=leafs.end(); old_leaf_key++){
		node_pos = nodes.find(*old_leaf_key);
		if (node_pos->second.child_edges.size()==0){
			parent_key = erase_edge_node(*old_leaf);
			while (nodes[parent_key].child_edges.size()==0){
				parent_key = erase_edge_node(parent_key);
			}
			while (nodes[parent_key].child_edges.size()==1){
				parent_key = bridge_edge_node(parent_key);
			}
		}else if (node_pos->child_edges.size()==1){
			parent_key = bridge_edge_node(*old_leaf);
			while (nodes[parent_key].child_edges.size()==1){
				parent_key = bridge_edge_node(parent_key);
			}
		}
	}


	leafs.clear();
	leafs.reserve(new_leafs.size());
	for(vector <key_t>::iterator new_leaf_key=new_leafs.begin(); new_leaf_key!=new_leafs.end(); new_leaf_key++){
		leafs.push_back(*new_leaf_key);
	}

	if (GEN_VERBOSE){
		cerr <<"genealogy::add_generation(). done "<<endl;
	}

	return;
}


key_t genealogy::erase_edge_node(key_t to_be_erased){
	if (GEN_VERBOSE){
		cerr <<"genealogy::erase_edge_node(). ..."<<to_be_erased.age<<" "<<to_be_erased.index<<endl;
	}

	map <key_t,node_t>::iterator Enode = nodes.find(to_be_erased);
	map <key_t,edge_t>::iterator Eedge = edges.find(to_be_erased);

	if (Enode->second.child_edges.size()>1){
		cerr <<"genealogy::erase_edge_node(): attempting to erase non-terminal node"<<endl;
	}

	key_t parent_key = Eedge->second.parent_node;
	map <key_t,edge_t>::iterator Pedge = edges.find(parent_key);
	map <key_t,node_t>::iterator Pnode = nodes.find(parent_key);

	Pnode->second.number_of_offspring-=Eedge->second.number_of_offspring;
	Pedge->second.number_of_offspring-=Eedge->second.number_of_offspring;

	list <key_t>::iterator children = Pnode->second.child_edges.begin();
	for (;children!=Pnode->second.child_edges.end(); children++){
		if (*children == to_be_erased){Pnode->second.child_edges.pop(to_be_erased); break;}
	}
	nodes.erase(to_be_erased);
	edges.erase(to_be_erased);

	if (GEN_VERBOSE){
		cerr <<"genealogy::erase_edge_node(). done"<<endl;
	}

	return parent_key;
}


key_t genealogy::erase_bridged_node(key_t to_be_bridged){
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
	Pedge->second.parent_key = Eedge->second.parent_node;

	Pedge->second.segment[0]=(Pedge->second.segment[0]<Eedge->second.segment[0])?(Eedge->second.segment[0]):(Pedge->second.segment[0]);
	Pedge->second.segment[1]=(Pedge->second.segment[1]>Eedge->second.segment[0])?(Eedge->second.segment[1]):(Pedge->second.segment[1]);
	Pedge->length+=Eedge->length;
	nodes.erase(to_be_bridged);
	edges.erase(to_be_bridged);

	if (GEN_VERBOSE){
		cerr <<"genealogy::bridge_edge_node(). done"<<endl;
	}

	return parent_key;
}

void genealogy::update_leaf_to_root(key_t leaf_key){
	map <key_t,node_t>::iterator leaf_node = nodes.find(leaf_key);
	map <key_t,node_t>::iterator leaf_node = nodes.find(leaf_key);
	leaf_edge->number_of_offspring = leaf_node->number_of_offspring;
	map <key_t,node_t>::iterator parent_node = nodes.find(leaf_edge->parent);
	map <key_t,edge_t>::iterator parent_edge = nodes.find(leaf_edge->parent);
	while (parent_node!=root){
		parent_node->number_of_offspring+=leaf_edge->number_off_spring;
		parent_edge->number_of_offspring=parent_node->number_of_offspring;

		leaf_node = parent_node;
		edge_node = parent_edge;
		parent_node = nodes.find(leaf_edge->parent);
		parent_edge = nodes.find(leaf_edge->parent);
	}
}

void genealogy::SFS(gsl_histogram sfs){
	map <key_t,edge_t>::iterator edge = edges.begin();
	int total_pop = root.number_of_offspring;
	for (; edge!=edges.end(); edge++){
		gsl_histogram_accumulate(sfs, 1.0*edge->second.number_of_offspring/total_pop, edge->second.length);
	}
}

int genealogy::external_branch_length(){
	map <key_t,edge_t>::iterator edge = edges.begin();
	int branchlength = 0;
	for (vector <key_t>::iterator leaf=leafs.begin(); leaf!=leafs.end();leaf++){
		branchlength+=edges[*leaf]->length;
	}
	return branchlength;
}

int genealogy::total_branch_length(){
	map <key_t,edge_t>::iterator edge = edges.begin();
	int branchlength = 0;
	for (; edge!=edges.end(); edge++){
		branchlength+=edge->length;
	}
	return branchlength;
}


void genealogy::clear_tree(){
	map <key_t,node_t>::iterator node = nodes.begin();
	map <key_t,edge_t>::iterator edge = edges.begin();
	for (; node!=nodes.end(); node++){
		node->number_of_offspring=0;
	}
	for (; edge!=edge.end(); edge++){
		edge->number_of_offspring=0;
	}


	for (vector<node_t>::iterator leaf=leafs.begin(); leaf!=leafs.end(); leaf++){
		node = nodes.find(node_t, *leaf);
		node->second.number_of_offspring=1;
	}
}


